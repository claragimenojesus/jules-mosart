#jeca@ichydro
library(stars)
library(RNetCDF)
library(tidyverse)
#creating source and destination directory. creating list for loop
sour_directory<-"/home/clara/rahu_data/WRF_climate_fixed/"
dest_dir<-"/home/clara/rahu_data/new_climate_data/meteo_input_data/"
list_weather<-c("wind.nc","rh.nc","pstar.nc","lw_down.nc","sw_down.nc","q.nc")


# a<-read_stars("/home/clara/rahu_data/WRF_climate_fixed/wind.nc") #testing

# dummy hourly timeseries to cover full climate scenarios required
dat1<-seq(from=as.POSIXct("2019-01-01",tz='UTC'),length.out=341880,by="hour")
dat2<-seq(from=as.POSIXct("2058-01-01",tz='UTC'),length.out=341880,by="hour")
dat3<-seq(from=as.POSIXct("2097-01-01",tz='UTC'),length.out=35016,by="hour")

# dat4<-c(dat1,dat2,dat3)

#Create loop
for (i in list_weather){
    a<-read_stars(paste0(sour_directory,i)) #load dataset
    b<-c(st_set_dimensions(a,'time',values=dat1),   #change time dimension
    st_set_dimensions(a,'time',values=dat2),    #change time dimension
    st_set_dimensions(a[,,,1:35016],'time',values=dat3),    #change time dim for shorter ncdf
    along="time")   #combining through time dimension
    # b[[1]]<-names(a)    #attribute name becomes messy, replacing name by original one
    write_stars(adrop(b[1]),paste0(dest_dir,i)) #generating file
}

# d<-seq(from=as.POSIXct("1980-01-01",tz='UTC'),length.out=341880,by="hour")
# b<-read_stars('/mnt/homes/clara/rahu_data/new_climate_data/rcp85/thourly/thourly_MIROC-ESM-CHEM_future_with_rcp85_r1i1p1BC_QDM_rolling_parallel_d03_JULES.nc'
# write_stars(obj=adrop(b[1]),dsn="") #generating file

# write_mdim(x=b, "wind.nc")
b2 <- slice(b, index = 2, along = time)

# Set dimensions to the same format as temperature and precipitation files
# wind
nc_path <- nc_open("/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data/wind.nc")
wind <- ncvar_get(nc_path,"wind")
lon <- ncvar_get(nc_path,"lon")
lat <- ncvar_get(nc_path, "lat")
time <- as.array(seq(0,718775,1))

ncpath <- "/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data_test/" 
ncname <- "wind"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))

wind_def <- ncvar_def("wind", "m/s", list(londim, latdim, timedim), longname="Wind speed at 2m",prec="double")

ncout <- nc_create(ncfname, wind_def, force_v4=TRUE)

ncvar_put(ncout , wind_def , wind)

ncout

nc_close(ncout)

# swdown
nc_path <- nc_open("/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data/sw_down.nc")
swdown <- ncvar_get(nc_path,"sw_down")

ncpath <- "/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data_test/" 
ncname <- "sw_down"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))

swdown_def <- ncvar_def("sw_down", "W/m2", list(londim, latdim,timedim), longname="Downward short wave flux at ground surface",prec="double")

ncout <- nc_create(ncfname, swdown_def, force_v4=TRUE)

ncvar_put(ncout , swdown_def , swdown)

ncout

nc_close(ncout)


# lwdown
nc_path <- nc_open("/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data/lw_down.nc")
lwdown <- ncvar_get(nc_path,"lw_down")

ncpath <- "/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data_test/" 
ncname <- "lw_down"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))

lwdown_def <- ncvar_def("lw_down", "W/m2", list(londim, latdim, timedim), longname="Downward long wave flux at ground surface",prec="double")

ncout <- nc_create(ncfname, lwdown_def, force_v4=TRUE)

ncvar_put(ncout , lwdown_def , lwdown)

ncout

nc_close(ncout)

# pstar
nc_path <- nc_open("/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data/pstar.nc")
pstar <- ncvar_get(nc_path,"pstar")

ncpath <- "/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data_test/" 
ncname <- "pstar"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))

pstar_def <- ncvar_def("pstar", "Pa", list(londim, latdim, timedim),longname="Surface pressure",prec="double")

ncout <- nc_create(ncfname, pstar_def, force_v4=TRUE)

ncvar_put(ncout , pstar_def , pstar)

ncout

nc_close(ncout)

# q
nc_path <- nc_open("/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data/q.nc")
q <- ncvar_get(nc_path,"q")

ncpath <- "/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data_test/" 
ncname <- "q"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))

q_def <- ncvar_def("q", "kg/kg", list(londim, latdim, timedim), longname="Specific humidity at 2m",prec="double")

ncout <- nc_create(ncfname, q_def, force_v4=TRUE)

ncvar_put(ncout , q_def , q)

ncout

nc_close(ncout)


# rh
nc_path <- nc_open("/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data/rh.nc")
rh <- ncvar_get(nc_path,"rh")

ncpath <- "/mnt/homes/clara/rahu_data/new_climate_data/meteo_input_data_test/" 
ncname <- "rh"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

londim <- ncdim_def("lon","degrees_east", as.double(lon))
latdim <- ncdim_def("lat","degrees_north", as.double(lat))
timedim <- ncdim_def("time", "hours since 2019-01-01 00:00:00", calendar="proleptic_gregorian", as.double(time))

rh_def <- ncvar_def("rh", "percent", list(londim, latdim, timedim), longname="Relative humidity at 2m",prec="double")

ncout <- nc_create(ncfname, rh_def, force_v4=TRUE)

ncvar_put(ncout , rh_def , rh)

ncout

nc_close(ncout)


