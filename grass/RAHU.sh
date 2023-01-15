
###################################################
# GRASS GIS setup for the RAHU project.
#
# The grass database is installed in /var/data/GRASS/RAHU on cv-hydro2
# There is a (read-only) PERMANENT mapset accessible to anyone
# Everyone has their own mapset named after their username
# in which you can run this script without the risk of overwriting
# anything else
#
# Wouter Buytaert, 12/2022 - ...
###################################################

## START GRASS

# Notes:
# - Everyone works in their own mapset. If the mapset does not exist, it will be created.
# - The GRASS region (which defines the geographical extent and resolution in which any analysis will be performed) 
#   is defined as a rectangle covering the basin with some extra margin

MAPSET="$(whoami)"
grass /var/data/GRASS/latlong/$MAPSET -c
g.region n=13:00:00.5S s=15:30:00.5S e=70:00:00.5W w=72:30:00.5W res=00:00:01

## IMPORT DIGITAL ELEVATION MODEL

# Notes:
# - The SRTM DEM is available in the PERMANENT mapset. It was imported there using the following code:

#for LON in {69..82}; do
#    for LAT in {01..19}; do
#        r.in.gdal input=/var/data/srtm/USGS//s${LAT}_w0${LON}_1arc_v3.tif output=s${LAT}_w0${LON}_1arc_v3
#    done
#done
#MAPS=`g.list type=raster separator=comma pat="s*"`
#g.region raster=$MAPS -p
#r.patch input=$MAPS output=srtm_1arc_v3_peru
#g.remove type=rast pattern=*1arc_v3 -f

## CREATE DEPRESSIONLESS DEM

r.fill.dir input=srtm_1arc_v3_peru output=srtm_depless direction=srtm_flowdir

r.watershed input=srtm_depless drainage=srtm_drainage accumulation=srtm_accum

#3. r.map.calc to calculate difference between filled and not-filled DEM. This will identify depressions.
#4. assign null values to negative or no differences (i.e. no depression). This will allow to discard non-depression data and let adjacent pixels to vectorize in one unit.
#5. vectorize the raster difference map
#6. generate grid flow accumulation
#7. filter vectors by ecosystems, distance to river and elevation. Max distance to population centres?.
#8. extract maximum (or mean?) flow accumulation for each vector




