#!/usr/bin/env Rscript
library(ncdf4)

dirs <- list.dirs("/mnt/scratch/clara/CMIP5_rahu_data/rcp45") [-1] # output folders
#dirs85 <- list.dirs("/home/clara/rahu_data/new_climate_data/input_data/rcp85") [-1] # output folders
pr_dirs <- list.files("/home/clara/rahu_data/new_climate_data/rcp45", pattern="*.nc", full.names=TRUE)  # input folders for pr_CMIP5_future_with
#pr_dirs85 <- list.files("/home/clara/rahu_data/new_climate_data/rcp85", pattern="*.nc", full.names=TRUE)  # input folders for pr_CMIP5_future_with
t_dirs <- list.files("/home/clara/rahu_data/new_climate_data/rcp45/thourly", full.names=TRUE) # input folders for thourly thourly_CMIP5_with
#t_dirs85 <- list.files("/home/clara/rahu_data/new_climate_data/rcp85/thourly", full.names=TRUE) # input folders for thourly thourly_CMIP5_with

original <- nc_open("/home/clara/rahu_data/WRF_climate_fixed/precip.nc")
lon <- ncvar_get(original, "lon")
lat <- ncvar_get(original, "lat")
time <- as.array(seq(0,718775,1))

# Precipitation post-processing
for (i in 1:30) {
    future_precip <- nc_open(pr_dirs[i])
    future_precip2 <- ncvar_get(future_precip, var="pr")
    M <- array(NA, dim=c(102,75,29949*24))
    for (j in 0:29948) {
        M[,,(j*24+1):(j*24+24)]=future_precip2[,,(j+1)]/24
    }
    ncpath <- paste0(dirs[i],"/")
    ncname <- "precip"
    ncfname <- paste(ncpath, ncname, ".nc", sep="")
    dname <- "precip"
    londim <- ncdim_def("lon","degrees_east", as.double(lon))
    latdim <- ncdim_def("lat","degrees_north", as.double(lat))
    timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))
    prec_def <- ncvar_def("precip", "mm", list(londim, latdim, timedim), longname="hourly total precipitation",prec="double")
    ncout <- nc_create(ncfname, prec_def, force_v4=TRUE)
    ncvar_put(ncout , prec_def , M)
    ncout
    nc_close(ncout)
}

# Precipitation post-processing
#for (i in 5:30) {
#    future_precip <- nc_open(pr_dirs85[i])
#    future_precip2 <- ncvar_get(future_precip, var="pr")
#    M <- array(NA, dim=c(102,75,29949*24))
#    for (j in 0:29948) {
#        M[,,(j*24+1):(j*24+24)]=future_precip2[,,(j+1)]/24
#    }
#    ncpath <- paste0(dirs85[i],"/")
#    ncname <- "precip"
#    ncfname <- paste(ncpath, ncname, ".nc", sep="")
#    dname <- "precip"
#    londim <- ncdim_def("lon","degrees_east", as.double(lon))
#    latdim <- ncdim_def("lat","degrees_north", as.double(lat))
#    timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))
#    prec_def <- ncvar_def("precip", "mm", list(londim, latdim, timedim), longname="hourly total precipitation",prec="double")
#    ncout <- nc_create(ncfname, prec_def, force_v4=TRUE)
#    ncvar_put(ncout , prec_def , M)
#    ncout
#    nc_close(ncout)
#}

# Temperature post-processing
for (i in 1:30) {
future_temp <- nc_open(t_dirs[i])
temp <- ncvar_get(future_temp,"T2")
# convert units from degree C to Kelvin
temp <- temp + 273.15
ncpath <- paste0(dirs[i],"/")
ncname <- "t"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))
temp_def <- ncvar_def("t", "K", list(londim, latdim, timedim), longname="Air temperature at 2m",prec="double")
ncout <- nc_create(ncfname, temp_def, force_v4=TRUE)
ncvar_put(ncout , temp_def , temp)
ncout
nc_close(ncout)
}

#for (i in 3:30) {
#future_temp <- nc_open(t_dirs85[i])
#temp <- ncvar_get(future_temp,"T2")
# convert units from degree C to Kelvin
#temp <- temp + 273.15
#ncpath <- paste0(dirs85[i],"/")
#ncname <- "t"
#ncfname <- paste(ncpath, ncname, ".nc", sep="")
#londim <- ncdim_def("lon","degrees_east", as.double(lon))
#latdim <- ncdim_def("lat","degrees_north", as.double(lat))
#timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))
#temp_def <- ncvar_def("t", "K", list(londim, latdim, timedim), longname="Air temperature at 2m",prec="double")
#ncout <- nc_create(ncfname, temp_def, force_v4=TRUE)
#ncvar_put(ncout , temp_def , temp)
#ncout
#nc_close(ncout)
#}