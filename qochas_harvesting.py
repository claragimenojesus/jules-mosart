#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import os
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.stats as stats
from shapely.geometry import mapping
import xesmf as xe

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])
SCENARIO_SH_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])
RC= str(os.environ['RC'])
SC = str(os.environ['SCENARIO'])

## jules output. change accordingly
#jules = xr.open_dataset(home+"/hydrodata/jules_runs/coupled_for_nbs/rcp85_ACCESS1-3_coupled_jules_oggm_00_99.nc",decode_coords="all").sel(time=slice("2000-01-01", "2018-12-31"))
jules = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_coupled_jules_oggm_00_99.nc"), decode_coords="all")

#mask >4000masl
mask=xr.open_dataset("/rds/general/user/cg2117/home/netcdf/mask_4000m.nc")['frac']

### Qochas parameter relations derived from medians from a Sierra Azul small dataset.
# We will produce linear spaces for storage volume but the range of storage volume would be delimited by the
# corresponding contribution area as this cannot exceed the cell size.
# Storage to contribution area is 13729.72/202500 = 0.06558. The cell size is an average of 1.56E7 m2.
# then the storage volume will be 1023030 M3, so 1E6 m3.
# [1] 13729.72 ## storage volume m3
# [1] 202500 ## contribution area m2 
# [1] 11652.64 ## qocha area m2
qochas_volumes=1E6
vol_area=11652.64/13729.72
vol_acc=202500/13729.72

# read unit hydrograph for shallow reinfiltration
uh = pd.read_csv(home+'/unit_hydro.csv').fillna(0)
uh["conv"]=uh["conv"]/uh["conv"].sum()

## hydraulic conductivity
satcon = xr.open_dataset("/rds/general/user/cg2117/home/netcdf/jules_soil_props_2015_rosetta3_ESA_rahu_modified_v2.nc")['satcon']

# cell area
cell_area = xr.open_dataset("/rds/general/user/cg2117/home/netcdf/gridcell_area_rahu_v2.nc")['area']

## masking cells where there exists flows over the entire time
masking=jules.surf_roff.sum(dim="time")
masking=masking.where(masking>0)

## start scenarios array
qochas_scenarios_flow= np.zeros((qochas_volumes.size,12))
available_storage=np.zeros_like(qochas_scenarios_flow)
et_month=np.zeros_like(qochas_scenarios_flow)


qochas_cap=st_max
qochas_cap=np.multiply(xr.full_like(mask,st_max).where(masking>0),mask)
qochas_area = qochas_cap*vol_area
qochas_acc = qochas_cap*vol_acc

### Zero arrays for various intermediate variables
St=xr.zeros_like(jules.surf_roff) #Storage at time step t
Qav = xr.zeros_like(satcon[0,...]) # storage volume + contributing runoff - losses (Et + drainage). or available storage
R = xr.zeros_like(jules.surf_roff) # Potential recharge at time step t
Qin = xr.zeros_like(satcon[0,...])
ET = xr.zeros_like(jules.surf_roff) 


### Qocha model
# cell_area=15526711 #grid_area[0,j,i].values for grid area
for t in range(jules.sizes["time"]-1):
    #### Calculate losses
    # Recharge limited by the hydraulic conductivity or the available water
    R[t,...]=np.minimum(86.4*satcon[0,...]*qochas_area,St[t,...])
    St[t,...]=np.add(St[t,...],-R[t,...])
    ## calculate ET
    ET[t,...]=np.minimum(np.abs(jules.fao_et0[t,...]) * qochas_area * 86.4,
                            St[t,...])
    St[t,...]=np.add(St[t,...],-ET)
    ## Calculate available storage
    Qav=np.add(qochas_cap,-St[t,...])
    # water that actually enters the qocha is limited by the storage capacity
    # Qin = (xr.where(Qav>qochas_cap, qochas_cap-St[t,...], qochas_acc*jules.surf_roff[t,...]*86.4))
    # Calculate inflow. This is considering overflow, i.e. cannot flow more than available storage
    Qin= np.minimum(qochas_acc*jules.surf_roff[t,...]*86.4,Qav)
    ## correcting flow
    jules.surf_roff[t,...]= (jules.surf_roff[t,...]-(np.divide(Qin.values/86.4,cell_area[0,...])))
    #update storage
    St[t+1,...] = np.add(Qin,St[t,...])
# add shallow subsurface delay because of re infiltration
shallow_asnp = R.values
shallow_asnp = np.append(shallow_asnp, np.zeros((365,75,102)),axis=0)
shallow_asnp1 = np.zeros_like(shallow_asnp)
for i in range(jules.sizes["lon"]):
    for j in range(jules.sizes["lat"]):
        for t in range(jules.sizes["time"]):
            if shallow_asnp[t,j,i] > 0:
                shallow_asnp1[t:t+uh["conv"].size,j,i]+= (shallow_asnp[t,j,i]*uh["conv"]).values
#replace flow
jules.sub_surf_roff.values += np.divide(shallow_asnp1[:jules.sizes["time"],...]/86.4,cell_area[0,...].values)
jules.runoff.values = jules.surf_roff.values + jules.sub_surf_roff.values

## accumulate surface
surf_roff=jules.surf_roff.sel(time=jules.surf_roff.time.dt.year.isin(np.arange(2080,2101)),drop=True).resample(time="ME").mean().groupby("time.month").mean()
surf_roff=np.multiply(surf_roff,cell_area[0,...]).sum(dim=["lat","lon"])
surf_roff=surf_roff/1000
   
# Subsurface partitioning
q_deep = xr.zeros_like(jules.sub_surf_roff)
# calculate deep percolation. Hydraulic conductivity is multiplied by 1000 (m/s to mm/s )
q_deep[...]=np.minimum(jules.sub_surf_roff,1000*hyd_cond0)
## correct shallow subsurface flow
jules.sub_surf_roff.values = np.add(jules.sub_surf_roff,-q_deep)
sub_surf_roff=jules.sub_surf_roff.sel(time=jules.sub_surf_roff.time.dt.year.isin(np.arange(2080,2101)),drop=True).resample(time="ME").mean().groupby("time.month").mean()
sub_surf_roff=np.multiply(sub_surf_roff,cell_area[0,...]).sum(dim=["lat","lon"])
sub_surf_roff=sub_surf_roff/1000
## JULES gridded results as timeseries
deep_ts = np.multiply(q_deep,cell_area[0,...]).sum(dim=["lat","lon"])
deep_ts = deep_ts/1000
deep_ts = deep_ts.values
## convolve timeseries
[deeper_mod,uh2]= unit_hydro(V,k0,n0,np.arange(deep_ts.size),deep_ts)
pd_deep=pd.DataFrame({"deeper":deeper_mod},index=q_deep.time.values)["2080-01-01":]
pd_deep=pd_deep.groupby(pd_deep.index.month).mean()
total_flow=surf_roff.values+sub_surf_roff.values+pd_deep.to_numpy().flatten()
qochas_scenarios_flow[counter,:]=total_flow
available_storage[counter,:]=St.sel(time=jules.surf_roff.time.dt.year.isin(np.arange(2080,2101)),drop=True).resample(time="ME").mean().groupby("time.month").mean().mean(dim=["lat","lon"])
et_month[counter,:]=ET.sel(time=jules.surf_roff.time.dt.year.isin(np.arange(2080,2101)),drop=True).resample(time="ME").mean().groupby("time.month").mean().mean(dim=["lat","lon"])

    counter+=1

    
np.savetxt("scaling_qochas_imp_85.csv",qochas_scenarios_flow,delimiter=',')
np.savetxt("scaling_qochas_sto_85.csv",available_storage,delimiter=',')
np.savetxt("scaling_qochas_ET_85.csv",et_month,delimiter=',')
