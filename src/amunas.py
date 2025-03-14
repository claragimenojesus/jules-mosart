#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
import xesmf as xe

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
#JULES_PATH = str(os.environ['JULES_PATH_SH'])
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])
SCENARIO_SH_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])
RC = str(os.environ['RC'])
SC = str(os.environ['SCENARIO'])

## import jules runs (note this is 2000 to 2099)
## scenario rcp45_ACCESS1-0
#jules = xr.open_dataset(home+"/hydrodata/jules_runs/coupled_for_nbs/rcp85_ACCESS1-3_coupled_jules_oggm_00_99.nc",decode_coords="all")
jules = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_coupled_jules_oggm_00_99.nc"), decode_coords="all")

# import high res dem
dem = rxr.open_rasterio(home+'/sa_dem_3s.tif').isel(band=0).rio.clip_box(minx=-74,miny=-16,maxx=-70,maxy=-12)


## regrid dem
threshold = 4000
regridder=xe.Regridder(dem,jules,"bilinear")
dr_out=regridder((dem>threshold).astype(float),keep_attrs=True)


### We will deviating water according to Ochoa Tocachi et al (2019)
### i.e. the diversion operation rules is: will diverge whenever flow is greater than 4 l/s/km2 but no more than 39.9 l/s/km2. 10^6 l/s/km2 = kg/m2/s. So it is 4E-5 and 4E-6
## where operator will limit the diversion for december to april
## If else statement will shut off diversion if flow is too low
## np.minimum function will limit the diversion to the minimum between the diversion (flow-minimum flow for supporting water use) and the channel capacity
operator=jules.surf_roff.where(jules.surf_roff.time.dt.month.isin([1,2,3,4,12]))
operator=operator.where(operator>4E-6)
diversion=np.minimum(operator-4E-6,4E-5)
diversion=diversion*dr_out

## load UH
uh = pd.read_csv(home+'/unit_hydro.csv').fillna(0)
uh["conv"]=uh["conv"]/uh["conv"].sum()


## Correct surface flow
a=jules.surf_roff-diversion.fillna(0)
jules.surf_roff.values=a.values

# Correct subsurface flow
# Calculate ET
operator_et=jules.fao_et0.where(diversion>0)
et=np.minimum(operator_et*120000*8/15526711,diversion) ## note that we need to adjust 15526711 to cell size
# Calculate infiltration
infilt=diversion.fillna(0)-et.fillna(0)
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
jules.sub_surf_roff.values += shallow_asnp1[:jules.sizes["time"],...]
jules.runoff.values = jules.surf_roff.values + jules.sub_surf_roff.values
jules.esoil_gb.values += et.values
# Save netcdf
jules.to_netcdf(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_"+"amunas_coupled_jules_oggm_00_99.nc"))
