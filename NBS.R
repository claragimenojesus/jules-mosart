#############################################
## NBS model for JULES
##
## project: RAHU
## Wouter Buytaert 17/08/2022
#############################################

## libraries. Install if not present

library(stars)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggthemes)

## Variables

# Each cocha and amuna consist of the following elements:
# - input time series
# - storage time series
# - output time series
# - coordinates
# - indices of the location in the netcdf matrix

cochas <- list()
amunas <- list()

# read JULES output data and prepare cochas and amunas

## read JULES output and plot for visual checking


f  <- read_ncdf("projects/JULES_vn6.1.S2.daily_hydrology.2018.2D.nc", var = "runoff")
f1 <- st_set_dimensions(f, 3, values = as.character(st_get_dimension_values(f, 3)))
f2 <- slice(f, index = 2, along = "time")

ggplot() +  
  geom_stars(data = f2[1], alpha = 0.8) + 
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(3, "cm"))

## 1. Qochas

# note: for st_extract to work, we need at least 2 points for some reason

#c <- st_point(c(-71,-13)) |> st_sfc(crs = st_crs(f))    # random point
#d <- st_point(c(-72,-14)) |> st_sfc(crs = st_crs(f))
#c <- c(c,d)

# or simply as matrix: 
#c <- rbind(c(-71,-13), c(-72, -14))

# but for now, just use random points (note, these may fall outside the catchment)
qochalocs = st_sample(st_as_sfc(st_bbox(runoff)), 20)

qocha_ts <- as.data.frame(st_extract(runoff, qochalocs))

for(i in unique(qocha_ts$geometry)) {
  ts <- qocha_ts[qocha_ts$geometry == i,]
  # apply cocha model
}

# return time series to JULES data



## 2. Amunas







## 3. restoration







## write JULES output data

