#import libraries
import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import os
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import mapping
import matplotlib.colors as colors
# import xesmf as xe

rcp45=xr.open_dataset("/home/clara/rcp45_ACCESS1-0_degradation_coupled_jules_oggm_00_99.nc",decode_coords="all") #4.5
rcp45=rcp45.sel(time=rcp45.time.dt.year.isin(np.arange(2078,2101)),drop=True)
# rcp85=xr.open_dataset("/home/clara/rcp85_ACCESS1-0_degradation_coupled_jules_oggm_00_99.nc",decode_coords="all") #8.5
# rcp85=rcp85.sel(time=rcp85.time.dt.year.isin(np.arange(2078,2101)),drop=True)
#mask
mask=xr.open_dataset("/home/clara/NbS_paper/mask_4000m.nc")['frac']

# read unit hydrograph
uh = pd.read_csv('/home/jec18/unit_hydro.csv').fillna(0)
uh["conv"]=uh["conv"]/uh["conv"].sum()
## hydraulic conductivity
hyd_cond = xr.open_dataset("/home/jec18/gis_layer/groundwater/hydrogeo_k.nc").fillna(0)
## I am copying the hydraulic conductivity values across time dimension. This allows us to do array operations in xarray without dimension problems
hyd_cond0 = xr.zeros_like(rcp45.sub_surf_roff)

for t in range(rcp45.sizes["time"]):
    hyd_cond0[t,...]=hyd_cond.Band1.values
catchment= gpd.read_file("/home/jec18/dem/km105/km105.shp")

## unit hydrograph
def unit_hydro(V,k,n,t,er):
    from scipy.special import gamma
    import numpy as np
    from scipy.sparse import spdiags

    #inputs
    # V = volume of unit hydrograph
    # k,n = parameters from Nash (1957)
    # t = unit hydrograph time
    # er = effective rainfall time-series
    
    # outputs
    # DGW = deeper groundwater flow
    # uh = unit hydrograph
    
    uh = np.array(V/(k*gamma(n))*np.multiply(np.exp(-t/k),np.power((t/k),n-1)))
    uh = uh/np.sum(uh)
    
    # from Ana Mijic as from Anthony script
    m = er.size
    r = t.size
    n = r+m-1
    
    pmat = spdiags(np.multiply(np.ones(r),np.row_stack(er)),np.arange(0,-er.size,-1),n,r)
    
    DGW = pmat*uh
    DGW = DGW[0:er.size]
    
    return [DGW, uh]

# k0 = 18.59720523 for KM105, and k0= 22.32687392 for PISAC
# n0 = 22.57835457 for KM105, and n0= 18.8277447 for PISAC
k0=18.59720523
n0=22.57835457
V=1

## masking cells where there exists flows over the entire time
masking=rcp45.surf_roff.sum(dim="time")
masking=masking.where(masking>0)
area_catchment=catchment.to_crs("EPSG:24893").area
area_catchment=area_catchment/1000/1000 #catchment area is in km2
area_catchment=area_catchment[0]
cell_area = xr.open_dataset("/home/clara/rahu_data/netcdf/gridcell_area_rahu_v2.nc")['area']
min_diversion=4E-6 ## minimum flow to start diversion
ratio_amuna = 3.5E-5/8 ## 1 amuna per 2 sq km carry 3.5 L/s/km2. Then for our cell area, the density of amuna is around 8 per cell.

## the diversion maximum was calculate as 3.5e-5 on top of the minimum flow to start to operate
## this was equal to 8 amunas per cell. Maximum is 100 amunas 
## 35 l/s/km2 was for 2 km2. it's equivalent was around 120000 m2 of greened area.
## if we would like to get all greened (maximum scenario), this would equate to the grid cell size / greened area =
amunas_linspace=np.linspace(0,100,num=101)
amunas_scenarios_flow= np.zeros((amunas_linspace.size,12))

## et area will vary linearly with the amount of water diverted
et_ratio = 1.2E6/75 #75l/s equated to 1.2 km2 of greened area 

counter=0
for div_ratio in amunas_linspace:
    jules=rcp45.copy(deep=True)
    area=xr.zeros_like(jules.surf_roff)
    ET0=xr.zeros_like(jules.surf_roff)
    ## calculating diverted water
    diversion=np.multiply(jules.surf_roff,mask)
    diversion= div_ratio/100*diversion.where(diversion.time.dt.month.isin([1,2,3,4,12])) # kg/m2/s
    area= np.multiply(diversion,cell_area[0,...])*et_ratio # this is the greened area
    ## Correct surface flow
    jules["surf_roff"]= jules.surf_roff-diversion.fillna(0)
    # Correct subsurface flow
    # Calculate ET
    ET0=np.divide(np.multiply(area,jules.fao_et0), #kg/s
                  cell_area[0,...]) #kg/m2/s
    # Calculate infiltration
    infilt=diversion.fillna(0)-ET0.fillna(0) ##infiltration is now in kg/m2/s
    ## unit hydrograph i,j wise for shallow subsurface flow
    shallow_asnp = infilt.values
    shallow_asnp = np.append(shallow_asnp,np.zeros((365,75,102)),axis=0)
    shallow_asnp1 = np.zeros_like(shallow_asnp)
    for i in range(jules.sizes["lon"]):
        for j in range(jules.sizes["lat"]):
            for t in range(jules.sizes["time"]):
                if shallow_asnp[t,j,i] > 0:
                    shallow_asnp1[t:t+uh["conv"].size,j,i]+= shallow_asnp[t,j,i]*uh["conv"].values
                    
    #replace flow
    jules.sub_surf_roff.values +=shallow_asnp1[:jules.sizes["time"],...] # no need to divide by area
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
    
    amunas_scenarios_flow[counter,:]=total_flow
    counter+=1
amunas_scenarios_flow  

np.savetxt("scaling_amunas_45.csv",amunas_scenarios_flow,delimiter=',')
