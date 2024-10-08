# load shape files and save in rds format
source("code/dependencies.R")

adm0 <- sf::read_sf("npl_admbnda_nd_20201117_shp/npl_admbnda_adm0_nd_20201117.shp")
adm1 <- sf::read_sf("npl_admbnda_nd_20201117_shp/npl_admbnda_adm1_nd_20201117.shp")
adm2 <- sf::read_sf("npl_admbnda_nd_20201117_shp/npl_admbnda_adm2_nd_20201117.shp")
adm2_districts <- sf::read_sf("npl_admbnda_nd_20201117_shp/npl_admbnda_districts_nd_20201117.shp")

# add topogrpahical zones to districts
tmp <- data.frame(DIST_EN = adm2_districts$DIST_EN,
                  TOPO = NA)%>%
  mutate(TOPO = case_when(DIST_EN %in% c("Jhapa",
                                         "Siraha",
                                         "Saptari",
                                         "Morang",
                                         "Sunsari",
                                         "Dhanusha",
                                         "Mahottari",
                                         "Sarlahi",
                                         "Bara",
                                         "Parsa",
                                         "Rautahat",
                                         "Chitawan",
                                         "Kapilbastu",
                                         "Nawalparasi East",
                                         "Nawalparasi West",
                                         "Rupandehi",
                                         "Dang",
                                         "Banke",
                                         "Bardiya",
                                         "Kailali",
                                         "Kanchanpur") ~ "Terai",
                          DIST_EN %in% c("Taplejung",
                                         "Sankhuwasabha",
                                         "Solukhumbu",
                                         "Dolakha",
                                         "Sindhupalchok",
                                         "Rasuwa",
                                         "Manang",
                                         "Mustang",
                                         "Dolpa",
                                         "Mugu",
                                         "Jumla",
                                         "Kalikot",
                                         "Bajura",
                                         "Humla",
                                         "Bajhang",
                                         "Darchula") ~ "Mountain",
                          TRUE ~ "Hill"))

adm2_d <- adm2_districts %>% left_join(tmp)


saveRDS(adm0, file = "data/shapefiles/npl_adm0.rds")
saveRDS(adm1, file = "data/shapefiles/npl_adm1.rds")
saveRDS(adm2, file = "data/shapefiles/npl_adm2.rds")
saveRDS(adm2_d, file = "data/shapefiles/npl_adm2_districts.rds")

# check that the topography categories make sense - YEP
ggplot(adm2_d)+
  geom_sf(aes(fill = TOPO))+
  geom_sf_text(aes(label = DIST_EN))

# unite the districts into the 3 topographical regions
topo <- adm2_d %>%
  group_by(TOPO) %>% 
  summarise()
ggplot(topo)+
  geom_sf(aes(fill = TOPO))+
  geom_sf_text(aes(label = TOPO))

saveRDS(topo, file = "data/shapefiles/npl_topo.rds")
