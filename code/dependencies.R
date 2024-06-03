source("R/set_paths.R")

library(pacman)
library(here)
library(tidyverse)
library(sf)
library(spatstat)
library(survey)
library(raster)
library(leaflet)
library(INLA)
library(malariaAtlas)
library(exactextractr)
library(ggsn)

## transforms sf file to Nepal Transverse Mercator projection (MUTM - Modified UTM)
transform_to_ntm <- function(my_obj){
  
  if ("sf" %in% class(my_obj)){
    rc <- st_transform(my_obj,
                       crs="+proj=tmerc +lat_0=0 +lon_0=84 +k=0.9999 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=293.17,726.18,245.36,0,0,0,0 +units=m +no_defs")
  } else if(class(my_obj) == "RasterLayer"){
    rc <- raster::projectRaster(my_obj,
                                crs="+proj=tmerc +lat_0=0 +lon_0=84 +k=0.9999 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=293.17,726.18,245.36,0,0,0,0 +units=m +no_defs",
                                method = "bilinear")
  } else {
    stop("can only transform rasterlayers and sf objects at the moment")
  }
  
}

