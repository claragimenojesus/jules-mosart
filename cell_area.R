# Clara Gimeno Jésus
# 8th August 2023

library(raster)
library(terra)

# Check WGS for our input grid
r2 <- rast("/home/clara/rahu_data/netcdf/jules_latlon_ESA_rahu.nc")
#class       : SpatRaster 
#dimensions  : 102, 2, 1  (nrow, ncol, nlyr)
#resolution  : 1, 0.03663889  (x, y)
#extent      : 0.5, 2.5, -73.88928, -70.15211  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 
#source      : jules_latlon_ESA_rahu.nc:lon_bnds 
#varname     : lon_bnds 
#name        : lon_bnds 
############## WGS84 - don't perform area calculation on this file as struggles to capture variable names


r <- rast("/home/clara/rahu_data/netcdf/jules_land_frac_ESA_rahu.nc")
a <- cellSize(r, unit="m")
crs(a) <- "epsg:4326"
a_df <- as.data.frame(a)
# mean(a[,1]) = 15.6146 km2

# Write nc file with grid areas
writeCDF(a, "/home/clara/rahu_data/netcdf/gridcell_area_rahu.nc", varname="area", unit="m2", prec="double", gridmap="", overwrite=TRUE)

# Check to see other ancil files give the same
#r3 <- rast("/home/clara/rahu_data/netcdf/jules_frac_9pft_2015_ESA_rahu.nc")
#a3 <- cellSize(r3, unit="km")
#a3_df <- as.data.frame(a3)
