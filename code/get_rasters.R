# this script extracts and plots the spatial data used in the regression
# all datasets are publicaly available

# 1. Population density (world pop)
# 2. Elevation (world pop)
# 3. Travel time to nearest city (malaria atlas)


source("code/dependencies.R")
# read in shape file to boundary raster
npl_adm0 <- readRDS("data/shapefiles/npl_adm0.rds")%>%
  st_transform(., crs = "+proj=longlat +datum=WGS84 +no_defs") 

#########################
# 1. population density #
#########################

NPL_pop <- raster("data/spatial_covariates/npl_ppp_2020.tif")

# convert to easy format for ggplot
NPL_pop_spdf <- as(NPL_pop, "SpatialPixelsDataFrame")
NPL_pop_df <- as.data.frame(NPL_pop_spdf)
colnames(NPL_pop_df) <- c("value", "x", "y")

# lets have a quick look
p4 <- ggplot()+
  coord_sf(datum = NA, expand = FALSE)+
  geom_raster(data = NPL_pop_df, aes(x = x, y = y, fill = value)) +
  scale_fill_distiller(palette = "Purples", direction = 1,
                       name = "people per pixel")+
  labs(x = " ") + labs(y = " ") + theme_void()+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))


# lets have a quick look at it logged
p5 <- ggplot()+
  coord_sf(datum = NA, expand = FALSE)+
  geom_raster(data = NPL_pop_df, aes(x = x, y = y, fill = log(value))) +
  scale_fill_distiller(palette = "Purples", direction = 1,
                       name = "log people per pixel")+
  labs(x = " ") + labs(y = " ") + theme_void()+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))
# seem to be missing pixels indicating population was zero --> -Inf?
# let's check - answer is yes - THESE ARE NATIONAL PARKS

NPL_pop_df <- NPL_pop_df %>% 
  mutate(zero = case_when(value == 0 ~ 1,
                          value >0 ~ 0,
                          TRUE ~ as.numeric(NA)))
ggplot()+
  coord_sf(datum = NA, expand = FALSE)+
  geom_raster(data = NPL_pop_df, aes(x = x, y = y, fill = zero)) +
  labs(x = " ") + labs(y = " ") + theme_void()+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))


################ 
# 2. elevation #
################

# Source: WorldPop
# 100m resulation. 
# The value of each grid cell represents its elevation above the sea level (in meters).
# Data source: de Ferranti, J., 2017. 'Digital Elevation Data'. 
# Viewfinder Panoramas (www.viewfinderPanoramas.org/dem3.html); 
# based on NASA's Shuttle Radar Topography Mission (SRTM) data 
# (http://www2.jpl.nasa.gov/srtm/). Year = 2000


NPL_elev <- raster("data/spatial_covariates/npl_srtm_topo_100m.tif")

# convert to easy format for ggplot
NPL_elev_spdf <- as(NPL_elev, "SpatialPixelsDataFrame")
NPL_elev_df <- as.data.frame(NPL_elev_spdf)
colnames(NPL_elev_df) <- c("value", "x", "y")

# lets have a quick look
p6 <- ggplot()+
  coord_sf(datum = NA, expand = FALSE)+
  geom_raster(data = NPL_elev_df, aes(x = x, y = y, fill = value)) +
  scale_fill_distiller(palette = "Purples", direction = 1,
                       name = "elevation (metres above sea level)")+
  labs(x = " ") + labs(y = " ") + theme_void()+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))


#####################################
# 3. time to travel to nearest city #
#####################################

# get travel time to nearest city from malaria atlas
# get raster from malaria atlas
NPL_ttime_15 <- npl_adm0 %>%
  as_Spatial() %>%
  malariaAtlas::getRaster(
    surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015", 
    shp = .)

# convert to easy format for ggplot
NPL_ttime_15_spdf <- as(NPL_ttime_15, "SpatialPixelsDataFrame")
NPL_ttime_15_df <- as.data.frame(NPL_ttime_15_spdf)
colnames(NPL_ttime_15_df) <- c("value", "x", "y")

# coordinates for kathmandu for plotting
kathmandu <- data.frame(27.7172, 85.3240) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X27.7172, long=X85.324) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")

# lets have a quick look
p1 <- ggplot()+
  coord_sf(datum = NA, expand = FALSE)+
  geom_raster(data = NPL_ttime_15_df, aes(x = x, y = y, fill = value)) +
  scale_fill_distiller(palette = "Purples", n.breaks = 10, direction = 1,
                       name = "time to nearest city")+
  labs(x = " ") + labs(y = " ") + theme_void()+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))+
  geom_sf(data = kathmandu, size = 3, lty="solid", shape = 15, 
          fill = "black", color="black") +
  geom_sf_text(data = kathmandu, label = "Kathmandu", size = 3, 
               nudge_x = 0.4, nudge_y = 0)