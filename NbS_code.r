#############################################
## NBS model for JULES
##
## project: RAHU
## Contributors: Jose Cuadros, Clara Gimeno Jésus, Wouter Buytaert
##
## TODO: explore whether it makes sense to use R markdown for a project like this:
## https://rmarkdown.rstudio.com/index.html
#############################################

## libraries. Install if not present

library(stars)
library(tidyverse)
library(viridis)
library(ggthemes)
library(ncdf4)
library(sf)
library(terra)

## Variables

# Path Clara
# list.files("/home/clara/project/rahu/jules-output",pattern="regrid")

# Each cocha and amuna consist of the following elements:
# - input time series
# - storage time series
# - output time series
# - coordinates
# - indices of the location in the netcdf matrix

#cochas <- list()
#amunas <- list()

## read JULES output data and prepare cochas and amunas
#gc() to check used space
# read JULES output and plot for visual checking
# naming conventions: 
# Qs: surface runoff
# Qss: subsurface runoff

#Load files

Jules  <- read_stars("/home/clara/project/rahu/jules-output/u-cj531/JULES_vn6.1.S2.daily_hydrology.2010.2D.nc", sub =c("sub_surf_roff","surf_roff","fao_et0"))
# fao_et0: FAO Penman-Monteith evapotranspiration for reference crop (kg m-2 s-1)

Jules<-st_set_crs(Jules,4326) #set 4326 epsg
#calculate pixel area
#grid_area <-cellSize(as(Jules[,,,1],"SpatRaster"))
grid_area <- rast("/home/clara/rahu_data/netcdf/gridcell_area_rahu.nc") # CG: 

# st_area(Jules) # Output: mean= 15672942 m^2
satcon<-read_stars("/home/clara/rahu_data/netcdf/jules_soil_props_2015_rosetta3_ESA_rahu.nc",sub="satcon") 

# Load GIS layers
qochas <- vect("/home/jec18/gis_layer/qochas/qochas.gpkg")  #This is a vector object
accum <- rast("/home/jec18/gis_layer/qochas/accum.tif") # accumulation properties of qochas
    #  CG: why do we load this if you have accum properties in qochas spatial vector?
dem <- rast("/home/jec18/gis_layer/qochas/dem.tif") # grid de 14400 14400
vub<-vect("/home/jec18/gis_layer/vub_catchment/vub_catchment/vub_catchment.shp")

gridcell <- rast("/home/clara/rahu_data/netcdf/gridcell_area_rahu.nc")


#### Calculate qochas properties
# Filtering qochas with most potential to accumulate water
qo1<-subset(qochas, qochas$accum_average>=quantile(qochas$accum_average,0.75,na.rm=T))
# Calculate qochas (vector) area
qo1$area<-expanse(qo1, unit="m")
# Number of qochas per pixel
qo2<-rasterize(centroids(qo1),as(Jules[,,,1],"SpatRaster"),fun=length) #matching jules resolution # CG: check sum does give 2057 qochas so well distirbuted
qo21<-rasterize(centroids(qo1),gridcell,fun='count') #matching jules resolution # CG: 
# Aggregated volume per pixel
Vmax<-qo2*5000 # can modify to match random distribution assignment. Assuming 5 m3 per qocha
# CG: maybe we should make the volume depend on the area of the geometry associated with each qocha?
# Qocha area per pixel
Aq <- qo2*2500 # Assuming 2500 m2 per qocha # CG: mean of qo1$area is 2831
# Contribution area per pixel
Ac<-rasterize(centroids(qo1),as(Jules[,,,1],"SpatRaster"),field="accum_average",fun=sum) #contribution area # CG: contribution area per pixel

# directory<- sort(list.files("projects",pattern="*.nc"))
# paste("projects/",sort(list.files("projects",pattern="*.nc")),sep="")
# do the same for Qss, and other relevant variables
# generate a single map as an example for visualization:
#Qs2 <- slice(Jules, index = 2, along = "time")

##### 1. Qochas modelling routine

### Before running this part, we should double check that all raster and netcdf have same spatial resolution
St<-matrix(0, dim(Jules)[1],dim(Jules)[2]) # Qochas starts empty # CG: is this the volume of water in the qocha at time t? yes
Qav= 0 #Dummy variable for maximum runoff available, CG: maximum Qin to the qocha
R=0 #Dummy for recharge
Ks=0 #Dummy for hydraulic conductivity

# CG: This for loop is for one year only 

for (t in 1:dim(Jules)[3]){ # Time loop
  for (i in 1:dim(Jules)[1]){ # Looping in x
    for (j in 1:dim(Jules)[2]){ # Looping in y 
      if (!is.na(qo2[j,i])){ # CG change in ij
        # Variables
        Ks = as.numeric(satcon[[1]][i,j,1]) # Superficial layer # why just this layer?
        Et = as.numeric(Jules[[3]][i,j,t]) # select third attribute
        Qs = as.numeric(Jules[[2]][i,j,t])
        Qss = as.numeric(Jules[[1]][i,j,t])
        vmax=Vmax[j,i] # CG change in ij
        ### Water balance at qocha. Units are in m3
        # Recharge = min (Ks, St_1)
        # R=min(Ks*Aq[j,i],St[i,j]) # CG : I don't agree with this, R is a rate not a volume
        # CG modified R formula

        R = min(Ks*86.4*Aq[j,i], St[i,j]) # Convert kg/m2/s to m3 infiltrating daily

        Qav = St[i,j] + Ac[j,i] * Qs * 86.4 - R - Et * Aq[j,i] * 86.4 # CG modified 
        # We should consider here using a better estimate of Et for open water, it should be higher than the one outputted in JULES for reference crop

        if (Qav>vmax){
          Qin=vmax - St[i,j] # modified by CG #### unsure
          St[i,j]=as.numeric(vmax)
          # Qof = Qav-Vmax # Overflow 
          # Qin=Ac[i,j]*Qs-Qof
        } else{
          St[i,j]=as.numeric(Qav)
          # Qof=0 # No overflow, new variable not needed
          Qin = Ac[j,i] * Qs * 86.4
        }
        ## Correct Qs and Qss
        Jules[[2]][i,j,t]<- as.numeric(Qs - Qin*Aq[j,i]/grid_area[j,i])# Correcting flows to be in kg/m2/s.
        Jules[[1]][i,j,t]<- as.numeric(Qss - R*Aq[j,i]/grid_area[j,i]) # Correcting flows to be in kg/m2/s.


        # other idea
        St[i,j] = St[i,j] + Ac[j,i] * Qs * 86.4 - R - Et * Aq[j,i] * 86.4

        if (St[i,j] > vmax) {
          Qin = vmax - St[i,j]
          St[i,j] = vmax
          Qout = St[i,j] - vmax
          Qs = Qs + Qout - Qin
          Qss = Qss + R
        } else {
          Qs = Qs - Ac*Qs - R - Et*Aq
          Qss = Qss + R
        }
      }
    }
  }
}


