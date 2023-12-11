# Data manipulation for input for CC runs

library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)
library(stars)

# Objective:
# Convert future precipitation file to hourly data by extending the array and inputting the same values for hourly rainfall

# Format of old precip file to replicate
#stars object with 3 dimensions and 1 attribute
#attribute(s), summary of first 1e+05 cells:
#            Min. 1st Qu. Median         Mean      3rd Qu.        Max.
#precip [mm]    0       0      0 4.775065e-05 3.034025e-05 0.003136992
#dimension(s):
#     from     to         offset     delta  refsys x/y
#lon     1    102       -73.8893 0.0366519      NA [x]
#lat     1     75       -15.0211 0.0356078      NA [y]
#time    1 341880 1980-01-01 UTC   1 hours POSIXct  

# Read daily future precipitation file
future_precip <- read_ncdf('/mnt/homes/clara/rahu_data/new_climate_data/rcp85/pr_MIROC-ESM-CHEM_future_with_rcp85_r1i1p1BC_QDM_rolling_parallel_d03_JULES.nc', var="pr", proxy=FALSE)
# Rename variable name to match historical precipitation file format
future_precip <- setNames(future_precip, "precip [mm]")

# Convert to dataframe
precip_df <- raster::as.data.frame(future_precip,xy=TRUE)

# Convert to hourly data
precip_df[,4] <- precip_df[,4]/24

# Stack together daily data to generate hourly instances
stacked_df <- as.data.frame(lapply(precip_df, rep, each=24))

# Generate hourly time series
time<-seq(from=as.POSIXct("2019-01-01 00:00:00",tz='UTC'),to=as.POSIXct("2100-12-31 23:00:00",tz='UTC'),by="hour")

# Useful variables
latlon = 102*75
days = length(full_time)
hour = 24

# Preallocate space
full_times <- as.POSIXct(rep(NA,days*latlon))

for (j in 0:(days-1)) {
    full_times[(j*hour*latlon+1):(j*hour*latlon+hour*latlon)] <- rep(time[(j*hour+1):(j*hour+hour)], times=latlon)
}

colnames(stacked_df) <- c("lon","lat","time","precip..mm.") # change time column name to "time"

# Convert time column to hourly time series
stacked_df$time <- full_times

# Reconvert dataframe to nc file
# ??


###### Test on 10 days ####### works!
test10 <- future_precip[,,,1:10]
# this works! for 10 days but very long to run - second approach after
# df <- data.frame(lon = NA,lat=NA,time=array(NA,10*7650*24), precip=NA)
# for (lon in seq(0,101)) {
#   for (lat in seq(0, 74)) {
#     for (i in seq(0,9, by=1)) {
#       for (l in seq(0,23, by=1)) {
#        df[1+i*24+j,4] = precip_df[lon * 102 + lat * 75 + i + 1]
#       }
#     }
#   }
# }

# 1 st step reshape as a 3d array the daily data
test10_df <- raster::as.data.frame(test10,xy=TRUE)
test10_df[,4] <- test10_df[,4] / 24 # divide by 24 daily data for hourly tranformation

stacked_df10 <- as.data.frame(lapply(test10_df, rep, each=24))
time10<-seq(from=as.POSIXct("2019-01-01 00:00:00",tz='UTC'),to=as.POSIXct("2019-01-10 23:00:00",tz='UTC'),by="hour")

full_times10<- as.POSIXct(rep(NA,1836000))
#t <- as.POSIXct(rep(NA,1836000))
#t[1:183600] <- rep(time[1:24],times=7650) # this works for the first day

for (j in 0:9) {
    full_times10[(j*183600+1):(j*183600+183600)] <- rep(time10[(j*24+1):(j*24+24)], times=7650)
}

full_times10 <- as.data.frame(full_times10)

colnames(stacked_df10) <- c("lon","lat","time","precip..mm.") # change time column name to "time"

stacked_df10$time <- full_times10[,1]


# RECONVERT TO STARS/NC
# ??
