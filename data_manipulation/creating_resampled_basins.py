#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
import glob
from tqdm import tqdm
import os

path="/home/clara/project/NbS/jules-output/rcp45/ACCESS1-0/rcp45_ACCESS1-0_degradation_coupled_jules_oggm_00_99.nc"
#path = "/mnt/scratch/clara/project/rahu/jules-output/rcp45/ACCESS1-0/rcp45_ACCESS1-0_amunas_coupled_jules_oggm_00_99.nc"
hyd_cond = xr.open_dataset("/home/jec18/gis_layer/groundwater/hydrogeo_k.nc").fillna(0)

hist = xr.open_dataset(path)[["surf_roff","sub_surf_roff","melt"]].sel(time=slice("2080-01-01","2099-12-31"),drop=True)

## I am copying the hydraulic conductivity values across time dimension. This allows us to do array operations in xarray without dimension problems
hyd_cond0 = xr.zeros_like(hist.sub_surf_roff)
for t in range(hist.sizes["time"]):
    hyd_cond0[t,...]=hyd_cond.Band1.values

# substracting deep drainage
q_deep = xr.zeros_like(hist.sub_surf_roff)
q_deep[...]=np.minimum(hist.sub_surf_roff,1000*hyd_cond0)
hist.sub_surf_roff.values = np.add(hist.sub_surf_roff,-q_deep)

# substracting glacier melt from surface runoff
hist["surf_roff"]-=hist["melt"]

hist_month= hist.resample(time="M").mean().groupby("time.month").mean()
hist_dry=hist_month.sel(month=slice(5,11)).mean(dim="month")*3600*24 #dry season average in mm/d

hydro = hist_dry.copy()

hydro  = hydro.rio.write_crs("EPSG:4326", inplace=True)

surface = hydro["surf_roff"]
subsurface = hydro["sub_surf_roff"]
melt = hydro["melt"]

# Set spatial dimensions for rioxarray
surface = surface.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
subsurface = subsurface.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
melt = melt.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import numpy as np

# Step 1: Open your sub-basin raster
src_path = "/home/jec18/HydroSheds/basins_agu/basin_133.tif"
with rasterio.open(src_path) as src:
    # Reproject to EPSG:3857 or UTM zone appropriate for your region
    dst_crs = "EPSG:3857"  # Example: Web Mercator in meters

    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds
    )

    kwargs = src.meta.copy()
    kwargs.update({
        'crs': dst_crs,
        'transform': transform,
        'width': width,
        'height': height
    })

    # Step 2: Create destination array
    dst = np.empty((src.count, height, width), dtype=np.float32)

    # Step 3: Reproject
    for i in range(1, src.count + 1):
        reproject(
            source=rasterio.band(src, i),
            destination=dst[i - 1],
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest
        )

# Step 4: Calculate pixel area in m²
res_x, res_y = transform[0], -transform[4]  # pixel width and height in meters
pixel_area_m2 = (res_x * res_y)
print(f"Pixel area: {pixel_area_m2:.6f} m²")  # Should be ~0.0009 km² for 30m data

surface_projected = surface.rio.reproject("EPSG:32718")  # Use appropriate UTM zone for your location
res = surface_projected.rio.resolution()
area_m2 = abs(res[0] * res[1])  # now in m²
print(f"Surface grid area: {area_m2} m2")

# Writing files watershed_area, basins_resampled and total_area_m2 for one scenario as it won't change!!!

# Paths
basin_dir = "/home/jec18/HydroSheds/basins_agu"
basin_paths = sorted(glob.glob(basin_dir + "/*.tif"))[:100]
resampled_dir = "/home/clara/NbS_paper/resampled_basins/"
os.makedirs(resampled_dir, exist_ok=True)

# Output metadata storage
total_area_m2_list = []
watershed_area_km2_list = []

for i, basin_path in enumerate(basin_paths):
    with rasterio.open(basin_path) as src:
        basins = src.read(1)
        basins = np.where(basins == src.nodata, 0, basins)
        src_transform = src.transform
        src_crs = src.crs

    # Replace NaNs with 0 for masking
    basins[np.isnan(basins)] = 0
    basins[basins == 255] = 0

    # Resample to model grid
    basins_resampled = np.zeros(surface.shape, dtype=np.float32)
    reproject(
        source=basins,
        destination=basins_resampled,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=surface.rio.transform(),
        dst_crs=surface.rio.crs,
        src_nodata=src.nodata,
        resampling=Resampling.max
    )
    basins_resampled[basins_resampled == 255] = 0
    basins_resampled[np.isnan(basins_resampled)] = 0

    # Compute total area
    total_area_m2 = np.nansum(basins_resampled * area_m2)
    watershed_area_km2 = total_area_m2 / 1e6

    total_area_m2_list.append(total_area_m2)
    watershed_area_km2_list.append(watershed_area_km2)

    # Save mask
    np.save(os.path.join(resampled_dir, f"basin_{i:03d}.npy"), basins_resampled)

# Save area metadata
np.save(os.path.join(resampled_dir, "total_area_m2.npy"), np.array(total_area_m2_list))
np.save(os.path.join(resampled_dir, "watershed_area_km2.npy"), np.array(watershed_area_km2_list))
