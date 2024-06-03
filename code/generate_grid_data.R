# this script generates the dataset of covariates for spatial analyses

source("dependencies.R")

## set to true if need to rerun from scratch (otherwise loads existing grid)
rerun_all <- TRUE

# read in IVI survey data
df <- readRDS(paste0(my_path$data, "/merged.rds")) # this will not be reproducible as we cannot release the GPS data

# read in Nepal boundary shape file
npl_adm0 <- readRDS("data/shapefiles/npl_adm0.rds")%>%
  transform_to_ntm()
  
# read in the population (WorldPop)
npl_pop <- here("data", "spatial_covariates",
                "npl_ppp_2020.tif") %>%
  raster() %>% transform_to_ntm() # this transformation takes ages

# read in elevation (WorldPop)
npl_elev <- here("data", "spatial_covariates", "npl_srtm_topo_100m.tif") %>%
  raster() %>%  transform_to_ntm()

# get travel time to nearest city from malariaAtlas
# have to do this to match malariaAtlas, then transform to ntm
npl_adm0_sp <- readRDS("data/shapefiles/npl_adm0.rds")%>%
  st_transform(., crs = "+proj=longlat +datum=WGS84 +no_defs") 

npl_ttime_15 <- npl_adm0_sp %>%
  as_Spatial() %>%
  malariaAtlas::getRaster(
    surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015", 
    shp = .)

npl_ttime_15_ntm <- npl_ttime_15 %>%
  transform_to_ntm()

npl_ttime_15_ntm[npl_ttime_15_ntm==0] <- 0.0001  #for those in cities, setting this to a small number rather than having zeros


####################################################
## set up grid and bring in raster covariates etc ##
####################################################

## make a 5x5km grid over Nepal and then assigning ids
if(rerun_all){
  my_grid <- st_make_grid(npl_adm0, cellsize=c(5000,5000)) %>%
    st_sf(grid_id = 1:length(.))
  saveRDS(my_grid, here("data", "generated_data","five_by_five_grid.rds"))
}

my_grid <- readRDS(here("data","generated_data","five_by_five_grid.rds"))

## get grid centroids
grid_cents <- my_grid %>%
  group_by(grid_id) %>%
  st_centroid() %>%
  st_coordinates %>%
  data.frame

## get population in each grid cell
pop_grid <- exact_extract(npl_pop, my_grid,'sum')

## setting minimum to 10 people per 5km by 5km to help computation
pop_grid <- ifelse(pop_grid < 10, 10, pop_grid)

## adding population to grid cells
my_grid_pop <- my_grid %>%
  mutate(pop = pop_grid,
         logpop = log(pop_grid))

## get median elevation for each grid cell
elev_grid <- exact_extract(npl_elev, my_grid, 'median')

## get median travel time to major city for each grid cell
travel_time_grid <- exact_extract(npl_ttime_15_ntm, my_grid, 'median')

## add all covariates to grid cells
my_grid_covs <- my_grid_pop %>%
  mutate(travel_time = travel_time_grid,
         elevation = elev_grid)

## save grid with covariates
saveRDS(my_grid_covs,here("data","generated_data","grid_with_covs.rds"))


############################################
## match grid cells to sampled households ##
############################################

# please note this is not reproducible from the public repository 
# as the data has been deidentified

hhl_locs <- df %>%
  dplyr::select(HLongitude,HLatitude,H_ID) %>%
  st_as_sf(coords = c("HLongitude","HLatitude"),
           crs="+proj=longlat +datum=WGS84 +no_defs") %>%
  transform_to_ntm()

# add ntm northing and easting to survey data
hhl_locs_coords <- hhl_locs %>% 
  st_coordinates

df$H_easting <- hhl_locs_coords[,"X"]
df$H_northing <- hhl_locs_coords[,"Y"]

# grid id for each line in df = 69 grid cells sampled
grid_id <- st_within(hhl_locs,my_grid_covs) %>% unlist

# add grid_id and covariates to survey data
df$grid_id <- grid_id
gps_covs_data <- df %>%
  left_join(my_grid_covs %>% st_set_geometry(NULL)) %>%
  ## adding grid cell centroids
  left_join(grid_cents %>%
              mutate(grid_id = 1:nrow(grid_cents)) %>%
              rename(grid_easting=X,grid_northing=Y))

