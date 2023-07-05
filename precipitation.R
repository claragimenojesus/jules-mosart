library(ncdf4)
library(ncdf4.helpers)
library(stars)
library(readr)
library(processNC)

precip_output <- nc_open("precip.nc")
precip <- ncvar_get(precip_output, "precip")
lon <- ncvar_get(precip_output, varid="lon")
lat <- ncvar_get(precip_output, varid="lat")
precip_output$dim$time$units
precip_time <- nc.get.time.series(precip_output, v="precip", time.dim.name="time")
head(precip_time)
dim(precip)

df <- aggregateNC("precip.nc","precip_year.nc", var="precip", startdate = "1980", enddate="2018", group_col=year)

daily_sum_precip <- precip_output %>% mutate(day=as.Date(time, format="%Y-%m-%d")) %>% group_by(day) %>% summarise(total_precip=sum(precip))



nlon <- dim(lon)
nlat <- dim(lat)
time <- ncvar_get(precip_output, "time")
nt <- dim(time)


precip_array <- ncvar_get(precip_output, "precip")
precip_vec_long <- as.vector(precip_array)
length(precip_vec_long)
precip_mat <- matrix(precip_vec_long, nrow=nlon*nlat, ncol=nt)
dim(precip_mat)


lonlat <- as.matrix(expand.grid(lon,lat))
precip_df02 <- data.frame(cbind(lonlat, precip_mat))



years <- seq(as.POSIXct("1980-01-01"), as.POSIXct("2018-12-31"), by="years")

precip <- read_ncdf("precip.nc", var="precip", proxy=FALSE)

precip_yearly <- aggregate(precip, by=years, FUN= sum, na.rm=TRUE )


years_2 <- seq(as.POSIXct("1980-01-01"), as.POSIXct("2018-12-31"), by="months")
precip_monthly <- aggregate(precip, by=years_2, FUN=sum, na.rm=TRUE)

by_t <- "years"
precip_yearly_2 <- aggregate(precip, by=by_t, FUN=sum, na.rm=TRUE)
