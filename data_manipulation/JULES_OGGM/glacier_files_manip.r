# 06/12/2023 CGJ

# Code regridding and assembling JULES_OGGM glacier data to JULES-hydro grid. i.e. multiple glaciers contributing to one gridcell in JULES.

# NOTE: this code will need to be re-run once we obtain the rest of the glacier data for historical_00_18. As of
#pt 9th August 2023, there are 16 glaciers missing from the dataset. Check with Jon this data is mising.

# If data is missing, we just need to change the latlon_516.csv file and use the updated areas.csv and melt.csv files.

library(raster)
library(dplyr)
library(magrittr)
library(MASS)

# Read JULES_OGGM glacier center coordinates
pts <- read.csv("/home/clara/rahu_data/JULES_OGGM/historical_00_18/latlon_516.csv") 
pts2 <- pts[,2:3] 
pts2 <- relocate(pts2, 'Longitude')
pts2 <- as.matrix(pts2)

# Read JULES_OGGM glacier areas (516x19). Each glacier has an area for each year.
area_oggm <- read.csv('/home/clara/rahu_data/JULES_OGGM/historical_00_18/areas.csv')

# Read JULES-hydro gridcell area size
r <- raster("/home/clara/rahu_data/netcdf/gridcell_area_rahu.nc")

# Preallocate 3D dataset for glacier areas
DF <- array(data=NA, dim=c(75,102,19))
for (year in 1:19) {
    # Read each year of data
    area <- area_oggm[,year+1]

    # Transfer values associated with 'object' type spatial data (points, lines, polygons) to raster cells. If x represents points, each point is assigned to a grid cell. Points that fall on a border between cells are placed in the cell to the right and/or in the cell below. The value of a grid cell is determined by the values associated with the points and function fun.
    # Sum values of glacier areas to contributing JULES-hydro gridcell.
    rnew <- rasterize(pts2, r, area, fun=sum, na.rm=TRUE)
    df <- as.matrix(rnew, xy=TRUE)

    # Store in 3D array
    DF[,,year] <- df
}

# Now create % glacier area per gridcell area 3d array
perc <- array(data=NA, dim=c(75,102,19))
for (year in 1:19) {
    for (lat in 1:75) {
        for (lon in 1:102) {
            perc[lat,lon,year] = DF[lat,lon,year]/r[lat,lon,1] 
        }
    }
}

save(perc, file="/home/clara/rahu_data/JULES_OGGM/percentage_glacier_per_gridcell_area.rda")

# Let's do the same for glacier melt
pts <- read.csv("/home/clara/rahu_data/JULES_OGGM/historical_00_18/latlon_516.csv")
pts2 <- pts[,2:3]
pts2 <- relocate(pts2, 'Longitude')
pts2 <- as.matrix(pts2)

melt_oggm <- read.csv('/home/clara/rahu_data/JULES_OGGM/historical_00_18/yearly_melt.csv')
r <- raster("/home/clara/rahu_data/netcdf/gridcell_area_rahu.nc")

DF <- array(data=NA, dim=c(75,102,19))
for (year in 1:19) {
    melt <- melt_oggm[,year+1]
    rnew <- rasterize(pts2, r, melt, na.rm=TRUE, fun=sum)
    dfmelt <- as.matrix(rnew, xy=TRUE)
    DF[,,year] <- dfmelt
}

save(DF, file="/home/clara/rahu_data/JULES_OGGM/yearly_melt_gridded.rda")


############## Checks ##################################################
# count number of points in each gridcell after rasterize()

pts <- read.csv("/home/clara/rahu_data/JULES_OGGM/historical_00_18/latlon.csv")
pts2 <- pts[,2:3]
pts2 <- relocate(pts2, 'Longitude')
pts2 <- as.matrix(pts2)

area_oggm <- read.csv('/home/clara/rahu_data/JULES_OGGM/historical_00_18/areas.csv')
r <- raster("/home/clara/rahu_data/netcdf/gridcell_area_rahu.nc")

rcc <- rasterize(pts2, r, fun='count')
rcc_df <- as.matrix(rcc, xy=TRUE)
write.csv(rcc_df, "/home/clara/rahu_data/JULES_OGGM/number_glaciers_per_gridcell.csv", row.names=FALSE)
glacier_loc <- which(!is.na(rcc_df), TRUE) # x=lat (1st column) y=lon (2nd column) 

# Gridcells with only one glacier contributing:
#which(rcc_df==1)
# [1] 2427 2729 2804 3099 3101 3248 3774 5218 5292 5515 5593 5669 5831 5892 5985
#[16] 6116 6127 6274
#> which(rcc_df==1, TRUE)
#      row col
# [1,]  27  33
# [2,]  29  37
#[3,]  29  38
# [4,]  24  42
# [5,]  26  42
# [6,]  23  44
# [7,]  24  51
# [8,]  43  70
# [9,]  42  71
#[10,]  40  74
#[11,]  43  75
#[12,]  44  76
#[13,]  56  78
#[14,]  42  79
#[15,]  60  80
#[16,]  41  82
#[17,]  52  82
#[18,]  49  84

# Now we can check the DF at those locations is equal to the area of contributing glacier point

# Check performed for obvious locations in jupyter notebook /home/clara/rahu_data/JULES_OGGM/