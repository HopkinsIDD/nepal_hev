# grid-cell level spatial analysis of the HEV seroprevalence data from Nepal
# please note that in the absence of the household GPS data this code cannot be used to 
# reproduce the figures and estimates in the manuscript, only to show our methodology.

source("dependencies.R")
my_path <- set_paths()

###################################################################
## 1. read in survey data with spatial covariates and 5x5km grid ##
###################################################################

## take data from 1 row per individual to 1 row per grid cell
dat <- readRDS(here("data", "generated_data","survey_data_with_raster_covs.rds")) %>% 
  group_by(grid_id) %>%
  summarize(n_pos = sum(Sample_Status == "POSITIVE"),
            n_tested = n(),
            grid_id = grid_id[1],
            pop=pop[1],
            logpop=logpop[1],
            travel_time=travel_time[1],
            elevation=elevation[1],
            grid_easting=grid_easting[1],
            grid_northing=grid_northing[1],
            H_IDs=H_ID[1])

## bring in grid
my_grid <- readRDS(here("data", "generated_data","grid_with_covs.rds"))

# load shape file (national)
npl_adm0 <- readRDS("data/shapefiles/npl_adm0.rds")%>%
  transform_to_ntm()


###############################################################
## 2. set up SPDE (Stochastic Partial Differential Equation) ##
###############################################################

grid_locs_ll <- dat %>%
  dplyr::select(grid_easting, grid_northing, grid_id) %>% 
  st_as_sf(coords = c("grid_easting","grid_northing"),
           crs="+proj=tmerc +lat_0=0 +lon_0=84 +k=0.9999 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=293.17,726.18,245.36,0,0,0,0 +units=m +no_defs")

grid_locs_ntm <- grid_locs_ll %>% st_coordinates
dat$grid_easting <- grid_locs_ntm[,"X"]
dat$grid_northing <- grid_locs_ntm[,"Y"]

## make the convex hull for estimation around household coordinates
bnd <- inla.nonconvex.hull(grid_locs_ntm,convex=-0.1) ## inla mesh segment
mesh <- inla.mesh.2d(boundary = bnd ,
                     loc=grid_locs_ntm,
                     offset=c(-0.05, -0.05),
                     cutoff = 2000, ## if households are less than 2km apart, builds only a single vertex
                     max.edge=c(30000,50000)
)
## make a grid projector
grid_cents<- my_grid %>%
  group_by(grid_id) %>%
  st_centroid() %>%
  st_coordinates %>%
  data.frame %>%
  mutate(n = 100)

grid_projector <- inla.mesh.projector(mesh,loc=grid_cents %>% as.matrix)
A.pred <- grid_projector$proj$A ## get A matrix
plot(mesh,asp=1,main="hhl")
points(dat[,"grid_easting"] %>% unlist,dat[,"grid_northing"] %>% unlist,pch=19,cex=0.5,col="orange")

spde <- inla.spde2.matern(mesh=mesh, alpha=2) # alpha is Fractional operator order
s_index <- inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)

## observation/prediction matrix
A_est_multi <- inla.spde.make.A(mesh=mesh,
                                loc=dat %>%
                                  dplyr::select(grid_easting,grid_northing) %>%
                                  as.matrix)


################################
## 3. Fit spatial-only Models ##
################################

# First we will fit a model with only random spatial effects.

## model matrix using only spatial correlation
model_matrix <- model.matrix(~1,dat)[, -1]

## creating a stack called "est"
stack_est <- inla.stack(
  data=list(y=dat[,'n_pos'] %>% unlist,
            n=dat[,'n_tested'] %>% unlist), ## response
  A=list(A_est_multi,1), ## projector matrix
  effects=list(s_index,
               data.frame(Intercept=rep(1,nrow(dat)))),
  tag="est")

stack_pred <- inla.stack(
  data=list(y=NA,n=1),
  A=list(A.pred,1),
  effects=list(s_index,
               data.frame(Intercept=rep(1,nrow(my_grid)))),
  tag="pred")

est_and_pred_stack <- inla.stack(stack_pred,stack_est)

spatial_only_f <- y ~ -1 + Intercept + f(spatial.field, model = spde)

spatial_preds <- inla(spatial_only_f,
                      data=inla.stack.data(est_and_pred_stack),
                      family="binomial",
                      control.predictor=list(A=inla.stack.A(est_and_pred_stack),
                                             link=1,
                                             compute=TRUE),
                      control.compute=list(dic=FALSE,config=TRUE),
                      Ntrials=n)

## now we need to manipulate outputs to get estimates for posterior mean and sd into the grid

id.prd <- inla.stack.index(est_and_pred_stack, "pred")$data

## get posterior draws
sim <- inla.posterior.sample(1000, spatial_preds)
post_draws <- sapply(X=sim, FUN= function(x) { 1/(1+exp(-x$latent[id.prd])) }) # this sim takes a while
colnames(post_draws) <- paste0("pred_",1:ncol(post_draws))
post_draws <- post_draws %>% as.data.frame %>% mutate(grid_id = id.prd) %>% tibble

## number of people who have been infected in total
tmp = tidyr::pivot_longer(post_draws,cols=pred_1:pred_1000,names_to="draw",values_to="seroprev") %>% left_join(.,my_grid %>% dplyr::select(grid_id,pop) %>% st_set_geometry(NULL)) %>% mutate(infections=seroprev*pop)
tmp %>% group_by(draw) %>% summarize(infections=sum(infections)) %>% dplyr::select(infections) %>% unlist %>% quantile(.,c(.025,.5,.975))/sum(my_grid$pop)
tmp %>% group_by(draw) %>% summarize(infections=sum(infections)) %>% dplyr::select(infections) %>% unlist %>% quantile(.,c(.025,.5,.975))

spatial_means <- spatial_preds$summary.fitted.values$mean[id.prd]
spatial_sds <- spatial_preds$summary.fitted.values$sd[id.prd]

my_grid <- left_join(my_grid,
                     data.frame(grid_id = id.prd,
                                mean=spatial_means,
                                sd = spatial_sds,
                                rr = spatial_means/weighted.mean(spatial_means, w=exp(my_grid$logpop))                                        
                     ))

saveRDS(my_grid,here("data","generated_data","pred_spatial_only.rds"))

my_grid <- readRDS(here("data","generated_data","pred_spatial_only.rds"))

map_grid <-  st_intersection(my_grid,npl_adm0) 

#####################
# label coordinates #
#####################

# coordinates of Kathmandu
kathmandu <- data.frame(27.7172, 85.3240) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X27.7172, long=X85.324) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Kathmandu legend
kathmandu_leg <- data.frame(26.45, 83.4) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X26.45, long=X83.4) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Kathmandu text
kathmandu_txt <- data.frame(26.45, 84.15) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X26.45, long=X84.15) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Bagmati label so as not to obscur Kathmandu
Bagmati_txt <- data.frame(27.5172, 85.3240) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X27.5172, long=X85.324) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Province 2 label so as not to obscur sample sites
Province2_txt <- data.frame(26.97152, 85.6919) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X26.97152, long=X85.6919) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Province 1 label - allow change name to Koshi
Province1_txt <- data.frame(27.20986,87.27266) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X27.20986, long=X87.27266) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Sudurpaschim label so as not to obscur sample sites
Sudurpaschim_txt <- data.frame(29.42376, 80.93938) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X29.42376, long=X80.93938) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Terai label so as not to obscur sample sites
Terai_txt <- data.frame(27.68977, 83.00629) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X27.68977, long=X83.00629) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
# coordinates of Hill label so as not to obscur sample sites
Hill_txt <- data.frame(28.48977, 82.80629) %>% 
  as_tibble() %>% 
  dplyr::rename(lat=X28.48977, long=X82.80629) %>% 
  st_as_sf(.,coords = c("long","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")

# add province boundaries
npl_adm1 <- readRDS("data/shapefiles/npl_adm1.rds")%>%
  st_transform(., crs = "+proj=longlat +datum=WGS84 +no_defs")

# add ecological regions
npl_topo <- readRDS("data/shapefiles/npl_topo.rds")%>%
  st_transform(., crs = "+proj=longlat +datum=WGS84 +no_defs")

spatial_only_pred_plot <-
  ggplot() +
  geom_sf(data = map_grid , aes(fill = 100*mean), color=NA,lwd = 0) +
  geom_sf(data =  grid_locs_ll %>% transform_to_ntm(), color = 'lightpink', pch=3,size = 1.2) +
  geom_sf(data = npl_adm0, fill = NA, lwd = 0.1,color=alpha("white",0.5)) +
  geom_sf(data = npl_adm1, fill = NA, lwd = 0.1,color=alpha("white",0.5)) +
  geom_sf_text(data = npl_adm1 %>% filter(ADM1_EN!="Bagmati" & ADM1_EN!="Sudurpaschim" & ADM1_EN!="Province 2" & ADM1_EN!= "Province 1"), 
               stat = "sf_coordinates", aes(label = ADM1_EN), col = "white")+
  geom_sf_text(data = Bagmati_txt, label = "Bagmati", col = "white")+
  geom_sf_text(data = Province2_txt, label = "Madhesh", col = "white")+
  geom_sf_text(data = Province1_txt, label = "Koshi", col = "white")+
  geom_sf_text(data = Sudurpaschim_txt, label = "Sudurpaschim", col = "white")+
  geom_sf(data = kathmandu, size = 4, lty="solid", shape = 17, 
          fill = "red", color="red") +
  geom_sf(data = kathmandu_leg, size = 4, lty="solid", shape = 17, 
          fill = "red", color="red") +
  geom_sf_text(data = kathmandu_txt, label = "Kathmandu", size = 6)+
  scale_fill_viridis_c(direction=1,alpha=.9) +
  labs(fill="Seroprevalence (%)")    +
  coord_sf(datum = NA)  +
  ggsn::scalebar(npl_adm0, dist = 100, dist_unit = "km",transform=FALSE,location = "bottomleft",border.size = .3,
                 box.color = alpha("black",.6),box.fill = c(alpha("black",.6), "white")) + 
  labs(x = "") +
  labs(y = "") +
  theme_void()+
  theme(legend.position = c(0.8, 0.8))

ggsave(here("figs","sp_spatial_only.pdf"), width = 9, height = 6.2, units = "in")
ggsave(here("figs","sp_spatial_only.png"), width = 9, height = 6.2, units = "in")

## add topographical regions instead of provinces
spatial_only_pred_plot_topo <-
  ggplot() +
  geom_sf(data = map_grid , aes(fill = 100*mean), color=NA,lwd = 0) +
  geom_sf(data =  grid_locs_ll %>% transform_to_ntm(), color = 'lightpink', pch=3,size = 1.2) +
  geom_sf(data = npl_adm0, fill = NA, lwd = 0.1,color=alpha("white",0.5)) +
  geom_sf(data = npl_topo, fill = NA, lwd = 0.3,color=alpha("white",0.5)) +
  geom_sf_text(data = npl_topo %>% filter(GEO_REGION != "Terai" & GEO_REGION != "Hill"), 
               stat = "sf_coordinates", aes(label = GEO_REGION), col = "white", size = 5)+
  geom_sf_text(data = Terai_txt, label = "Terai", col = "white", size = 5)+
  geom_sf_text(data = Hill_txt, label = "Hill", col = "white", size = 5)+
  geom_sf(data = kathmandu, size = 4, lty="solid", shape = 17, 
          fill = "red", color="red") +
  geom_sf(data = kathmandu_leg, size = 4, lty="solid", shape = 17, 
          fill = "red", color="red") +
  geom_sf_text(data = kathmandu_txt, label = "Kathmandu", size = 6)+
  scale_fill_viridis_c(direction=1,alpha=.9) +
  labs(fill="Seroprevalence (%)")    +
  coord_sf(datum = NA)  +
  ggsn::scalebar(npl_adm0, dist = 100, dist_unit = "km",transform=FALSE,location = "bottomleft",border.size = .3,
                 box.color = alpha("black",.6),box.fill = c(alpha("black",.6), "white")) + 
  labs(x = "") +
  labs(y = "") +
  theme_void()+
  theme(legend.position = c(0.8, 0.8))

ggsave(here("figs","sp_spatial_only_topo.pdf"), width = 9, height = 6.2, units = "in")
ggsave(here("figs","sp_spatial_only_topo.png"), width = 9, height = 6.2, units = "in")
ggsave(here("figs","sp_spatial_only_topo.tif"), width = 9, height = 6.2, units = "in")

## Let's check how the population weighted mean from the map compares to the empirical estimate
map_grid$sp_pop <- map_grid$mean * map_grid$pop

sp_predict <- sum(map_grid$sp_pop) / sum(map_grid$pop)

###################################
## 4. Fit Models with covariates ##
###################################

## Covariate Models

#### Now let's add all the covariates (log-population, travel time and elevation):

model_matrix <- model.matrix(~logpop+
                              travel_time+                            
                              elevation,
                             data = dat)[, -1]

## creating a stack called "est"
stack_est <- inla.stack(
    data=list(y=dat[,'n_pos'] %>% unlist,
              n=dat[,'n_tested'] %>% unlist), ## response
    A=list(A_est_multi,1), ## projector matrix
    effects=list(s_index,
                 data.frame(Intercept=rep(1,nrow(dat)),model_matrix)),
    tag="est")

## mean imputation for missing covaraiates (due to clipping issues and imperfect alignment of rasters)
## the idea is to predict everywhere but then to exclude cells with missing covariate data from analyses
stack_pred <- inla.stack(
    data=list(y=NA,n=1),
    A=list(A.pred,1,1,1,1),
    effects=list(s_index,
                 data.frame(Intercept=rep(1,nrow(my_grid))),
                 list(logpop = ifelse(is.na(my_grid$logpop),mean(my_grid$logpop),my_grid$logpop)),
                 list(travel_time = ifelse(is.na(my_grid$travel_time),mean(my_grid$travel_time),my_grid$travel_time)),
                list(elevation = ifelse(is.na(my_grid$elevation),mean(my_grid$elevation),my_grid$elevation))),
    tag="pred")

est_and_pred_stack <- inla.stack(stack_pred,stack_est)

all_covs_f <- y ~ -1 + Intercept + logpop + travel_time + elevation + f(spatial.field, model = spde)

spatial_preds_all_covs <- inla(all_covs_f,
                               data=inla.stack.data(est_and_pred_stack),
                               family="binomial",
                               verbose=FALSE,
                               control.predictor=list(
                                   A=inla.stack.A(est_and_pred_stack),link=1,
                                   compute=TRUE),
                               control.compute=list(dic=FALSE,config=TRUE),
                               Ntrials=n)
saveRDS(file = "data/generated_data/all_covs.rds", summary(spatial_preds_all_covs))

## now we need to manipulate outputs to get estimates for posterior mean and sd into the grid
id.prd <- inla.stack.index(est_and_pred_stack, "pred")$data

## get posterior draws
sim <- inla.posterior.sample(1000, spatial_preds_all_covs)
post_draws <- sapply(X=sim, FUN= function(x) { 1/(1+exp(-x$latent[id.prd])) }) # this sim takes a while
colnames(post_draws) <- paste0("pred_",1:ncol(post_draws))
post_draws <- post_draws %>% as.data.frame %>% mutate(grid_id = id.prd) %>% tibble

## number of people who have been infected in total
tmp = tidyr::pivot_longer(post_draws,cols=pred_1:pred_1000,names_to="draw",values_to="seroprev") %>% left_join(.,my_grid %>% dplyr::select(grid_id,pop) %>% st_set_geometry(NULL)) %>% mutate(infections=seroprev*pop)
tmp %>% group_by(draw) %>% summarize(infections=sum(infections)) %>% dplyr::select(infections) %>% unlist %>% quantile(.,c(.025,.5,.975))/sum(my_grid$pop)

## get means
spatial_means2 <- spatial_preds_all_covs$summary.fitted.values$mean[id.prd]
spatial_sds2 <- spatial_preds_all_covs$summary.fitted.values$sd[id.prd]
my_grid <- left_join(my_grid,data.frame(
    grid_id = id.prd,
    mean2=spatial_means2,
    sd2 = spatial_sds2,
    rr2 = spatial_means2/
        weighted.mean(spatial_means2,w=exp(my_grid$logpop))            
))

saveRDS(my_grid,here("data","generated_data","preds_all.rds"))


## plot map of predicted seroprevalence
my_grid <-readRDS(here("generated_data","preds_all.rds"))

map_grid <-  st_intersection(my_grid,npl_adm0) 

covs_prev_plot <-
    ggplot() +
    geom_sf(data = map_grid , aes(fill = 100*mean2), color=NA,lwd = 0) +
    geom_sf(data = npl_adm0, fill = NA, lwd = 0.1,color=alpha("black",0.5)) +
    geom_sf(data =  grid_locs_ll %>% transform_to_ntm(), color = 'grey', pch=3,size = 1.2) +
    scale_fill_viridis_c(direction=1,alpha=.9) +
    labs(fill="Seroprevalence (%)")    +
    coord_sf(datum = NA)  +
    ggsn::scalebar(npl_adm0, dist = 100, dist_unit = "km",transform=FALSE,location = "bottomleft",border.size = .3,
                   box.color = alpha("black",.6),box.fill = c(alpha("black",.6), "white")) + 
    labs(x = "") +
    labs(y = "") +
    theme_void() +
  theme(legend.position = c(0.8, 0.8))

covs_prev_plot

ggsave(here("figs","sp_with_covs.pdf"), width = 9, height = 6.2, units = "in")
ggsave(here("figs","sp_with_covs.png"), width = 9, height = 6.2, units = "in")


#############################################################
## 5. Leave One Out Cross-Validation of the model (LOO-CV) ##
#############################################################

## leave one grid-id out of the fitting and predict on it
## save just the mean quantiles for each

rc = NULL
i=1

for(gc in dat$grid_id){
  print(i)
  
  tmp_dat <- dat %>%
    mutate(n_pos = ifelse(grid_id == gc,NA,n_pos))
  
  tmp_mod_matrix <- model.matrix(~logpop+travel_time+elevation,
                                 tmp_dat)[,-1]
  
  tmp_null_matrix <- model.matrix(~logpop,
                                  tmp_dat)[,-1]
  
  tmp_est_stack = inla.stack(
    data=list(y=tmp_dat[,'n_pos'] %>% unlist,
              n=dat[,'n_tested'] %>% unlist), ## response
    A=list(A_est_multi,1), ## projector matrix
    effects=list(s_index,
                 data.frame(Intercept=rep(1,nrow(tmp_dat)),tmp_mod_matrix)),
    tag="tmp")
  
  tmp_null_stack = inla.stack(
    data=list(y=tmp_dat[,'n_pos'] %>% unlist,
              n=dat[,'n_tested'] %>% unlist), ## response
    A=list(A_est_multi,1), ## projector matrix
    effects=list(s_index,
                 data.frame(Intercept=rep(1,nrow(tmp_dat)),tmp_null_matrix)),
    tag="tmp")
  
  all_covs_f <- y ~ -1 + Intercept + 
    logpop + 
    travel_time+                            
    elevation+
    f(spatial.field, model = spde)
  
  null_covs_f <- y ~ -1 + Intercept +f(spatial.field, model = spde)
  
  test <- inla(all_covs_f,
               data=inla.stack.data(tmp_est_stack),
               family="binomial",
               verbose=FALSE,
               control.predictor=list(
                 A=inla.stack.A(tmp_est_stack),link=1,
                 compute=TRUE),
               control.compute=list(dic=FALSE,config=TRUE),
               Ntrials=tmp_dat$n_tested)
  
  test_null <- inla(null_covs_f,
                    data=inla.stack.data(tmp_null_stack),
                    family="binomial",
                    verbose=FALSE,
                    control.predictor=list(
                      A=inla.stack.A(tmp_null_stack),link=1,
                      compute=TRUE),
                    control.compute=list(dic=FALSE,config=TRUE),
                    Ntrials=tmp_dat$n_tested)
  
  summary_hold_out = test$summary.fitted.values[which(is.na(tmp_dat$n_pos)), c("mean", "sd","0.025quant","0.5quant","0.975quant")]
  summary_others = test$summary.fitted.values[which(!is.na(tmp_dat$n_pos)), c("mean", "sd","0.025quant","0.5quant","0.975quant")]
  
  null_hold_out = test_null$summary.fitted.values[which(is.na(tmp_dat$n_pos)), c("mean", "sd","0.025quant","0.5quant","0.975quant")]
  null_others = test_null$summary.fitted.values[which(!is.na(tmp_dat$n_pos)), c("mean", "sd","0.025quant","0.5quant","0.975quant")]
  
  rc = rc %>% 
    bind_rows(
      data.frame(left_out_cell = gc,
                 full_mean = summary_hold_out$mean,
                 full_sd = summary_hold_out$sd,
                 full_median = summary_hold_out$`0.5quant`,
                 full_lb = summary_hold_out$`0.025quant`,
                 full_ub = summary_hold_out$`0.975quant`,
                 null_mean = null_hold_out$mean,
                 null_sd = null_hold_out$sd,
                 null_median = null_hold_out$`0.5quant`,
                 null_lb = null_hold_out$`0.025quant`,
                 null_ub = null_hold_out$`0.975quant`,
                 naive_mean = dat %>%
                   filter(grid_id!=gc) %>%
                   summarize(mean(n_pos/n_tested)) %>% unlist,
                 naive_median = dat %>%
                   filter(grid_id!=gc) %>%
                   summarize(median(n_pos/n_tested)) %>% unlist,
                 truth=dat %>%
                   filter(grid_id==gc) %>%
                   transmute(n_pos/n_tested) %>% unlist
      ))
  i <- i+1
}
## mean absolute error 
model_mae <- (rc$full_median - rc$truth) %>% abs %>% mean 
null_mae <- (rc$null_median - rc$truth) %>% abs %>% mean 
naive_mae <- (rc$naive_median - rc$truth) %>% abs %>% mean

# MAEs
model_mae
null_mae
naive_mae
## relative MAE
model_mae/naive_mae
model_mae/null_mae
null_mae/naive_mae
model_ae <- abs(rc$full_mean-rc$truth)
null_ae <- abs(rc$null_mean-rc$truth)
naive_ae <-abs(rc$naive_mean-rc$truth)
mean(model_ae < naive_ae)
mean(model_ae < null_ae)
mean(null_ae < naive_ae)

## bias 
model_bias <- (rc$mean - rc$truth) %>% mean
ggplot(rc,aes(x=truth)) + 
  geom_abline(intercept = 0, slope = 1,col='grey70') + xlim(0,.9) + ylim(0,.9) +
  geom_point(aes(y=full_mean),col='steelblue') +
  geom_smooth(aes(y=full_mean),method="lm",col='steelblue',lty=2,se=F)+ 
  geom_point(aes(y=null_mean),col='red3') +
  geom_smooth(aes(y=null_mean),method="lm",col='red3',lty=2,se=F)+ 
  xlab("proportion seropositive") +
  ylab("predicted proportion seropositive") +
  theme_bw()
ggsave(here("figs","cv_corplot.pdf"))
ggsave(here("figs","cv_corplot.png"))
saveRDS(rc, here("data","generated_data","loocv-results.rds"))


# for the model with just spatial random effects:

### Pearsons correlation of 
cor(rc$null_mean,rc$truth)  
###a mean absolute error of 
(rc$null_mean - rc$truth) %>% abs %>% mean
### and a bias of 
(rc$null_mean - rc$truth) %>% mean

# for the model with covariates:

## Does not improve on prediction of model with just spatial random effects
### Pearsons correlation of 
cor(rc$full_mean,rc$truth)  
###a mean absolute error of 
(rc$full_mean - rc$truth) %>% abs %>% mean
### and a bias of 
(rc$full_mean - rc$truth) %>% mean


