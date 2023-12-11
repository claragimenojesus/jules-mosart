# R script to convert latlon coordinates to km distance
# Clara Gimeno Jesus, Nov 2022

latlon_to_km <- function(lat1,lat2,lon1,lon2) {
  # Implement the Haversine formula
  
  # First, find the values of the lon/lat coordinates in radians
  lat1_r = lat1 / (180/pi)
  lat2_r = lat2 / (180/pi)
  lon1_r = lon1 / (180/pi)
  lon2_r = lon2 / (180/pi)
  
  R = 6371
  
  a = sin((lat2_r-lat1_r)/2) **2 + cos(lat1_r) * cos(lat2_r) * sin((lon2_r-lon1_r)/2) **2
  c = 2 * atan2(sqrt(a), sqrt(1-a))
  
  d = R*c
  
  return(d)
  
}