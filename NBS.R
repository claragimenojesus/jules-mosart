#############################################
## NBS model for JULES
##
## project: RAHU
## Contributors: Wouter Buytaert
##
## TODO: explore whether it makes sense to use R markdown for a project like this:
## https://rmarkdown.rstudio.com/index.html
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

## read JULES output data and prepare cochas and amunas

# read JULES output and plot for visual checking
# naming conventions: 
# Qs: surface runoff
# Qss: subsurface runoff

Qs  <- read_ncdf("projects/JULES_vn6.1.S2.daily_hydrology.2018.2D.nc", var = "runoff")
Qs <- st_set_dimensions(Qs, 3, values = as.character(st_get_dimension_values(Qs, 3)))

# do the same for Qss, and other relevant variables

# generate a single map as an example for visualization:

Qs2 <- slice(Qs, index = 2, along = "time")

ggplot() +  
  geom_stars(data = Qs2[1], alpha = 0.8) + 
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(3, "cm"))


## 1. Qochas

# note: for st_extract to work, we need at least 2 points for some reason
# TODO: figure out why st_extract does not work with 1 point (or what I am doing wrong)

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
  # apply cocha model. The basic steps are as follows:
  # 1. Decide on cocha volume Vc_max, and contributing area Ac
  # 2. Fill cocha: during the wet season the cocha fills at a rate determined by Ac, until it is full (Vc = Cv_max)
  # 3. Once Vc > 0, calculate recharge
  # 4. add recharge to Qss time series
  # 5. decide how to deal with Cocha overflow
  
}

# insert altered time series into the original JULES data (Qs, Qss, ...)



## 2. Amunas







## 3. restoration







## write altered JULES output data

