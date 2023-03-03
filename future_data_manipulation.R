# Data manipulation for input for CC runs

library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)

# Objective:
# Concatenate variables: wind, pstar, lw_down, q, rh, sw_down to convert from 1980 to 2018 timeframe to 2000 to 2100
# Convert future precipitation file to hourly data by extending the array and inputting the same values for hourly rainfall

# Original files - 1980-01-01 00:00:00 to 2018-12-31 12:59:59
hist_precip <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/precip.nc')
hist_lwd <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/lw_down.nc')
hist_pstar <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/pstar.nc')
hist_q <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/q.nc')
hist_swd <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/sw_down.nc')
hist_t <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/t.nc')
hist_wind <- nc_open('/mnt/homes/clara/rahu_data/WRF_climate_fixed/wind.nc')

old_time <- ncvar_get(hist_precip, "time")
old_precip <- ncvar_get(hist_precip, "precip")

# Future files - 2019-01-01 00:00:00 to 2100-12-31 12:59:59
future_tempfile <- nc_open('/mnt/homes/clara/rahu_data/new_climate_data/rcp85/thourly/thourly_MIROC-ESM-CHEM_future_with_rcp85_r1i1p1BC_QDM_rolling_parallel_d03_JULES.nc')
future_time <- ncvar_get(future_tempfile,"time2")
future_temp <- ncvar_get(future_tempfile, "T2")

# Daily precipitation - needs to be converted to hourly
pr <- nc_open('/mnt/homes/clara/rahu_data/new_climate_data/rcp85/pr_MIROC-ESM-CHEM_future_with_rcp85_r1i1p1BC_QDM_rolling_parallel_d03_JULES.nc')
future_precip_daily <- ncvar_get(pr,"pr")
counter = 1
future_precip_hourly <- array(numeric(),dim(c(length(lon),length(lat),length(future_time))))
for (i in seq(0,future_time[length(future_time)], by=24)) {
  future_precip_hourly[,,i:i+24] = future_precip_daily[,,counter]/24
  counter = counter + 1
}
