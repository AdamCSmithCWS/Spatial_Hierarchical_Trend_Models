### simulating fake population trajectories
# Download the Pacific Decadal Oscilation data and save to rds
# library(rsoi)
# pdo <- download_pdo() # Run on February 28, 2022
# saveRDS(pdo,file = "data/pdo.rds")
### pdo is monthly PDO values Jan 1900 - Jan 2022

library(tidyverse)
library(bbsBayes2) # https://github.com/bbsBayes/bbsBayes2


# select real BBS data for PAWR -------------------------------------------


species = "Pacific Wren"  
species_f <- gsub(species,pattern = " ",replacement = "_") #name without spaces

stratification = "bbs_usgs"

s <- stratify(by = stratification,
              species = species)

prep_data <- prepare_data(s)


prep_spat_data <- prepare_spatial(prep_data,
                      strata_map = load_map(stratification))
print(prep_spat_data$spatial_data$map) #visualise the spatial strata
  
 

# dataframe of all observations -------------------------------------------

real_data <- prep_spat_data$raw_data

min_year <- min(real_data$year)



# generate simulated trajectories -----------------------------------------



# set up spatial ----------------------------------------------------------
coord_f <- function(geometry, w_xy = "X"){ # function to extract single coord
  tt = sf::st_centroid(geometry) %>%
    sf::st_coordinates()%>%
    as.data.frame()
  return(tt[,w_xy])
}
base_map <- load_map(stratify_by = stratification) %>% 
  filter(strata_name %in% real_data$strata_name) %>% 
  select(strata_name,geom) %>% 
  mutate(x_centroid = coord_f(geom,"X"),
         y_centroid = coord_f(geom,"Y"),
         y_scale = as.numeric(scale(y_centroid)),
         north = y_scale - min(y_scale),
         south = max(y_scale) - y_scale)

vis_map <- ggplot(base_map)+
  geom_sf(aes(fill = south))+
  geom_sf_text(aes(x = x_centroid, y = y_centroid,
                   label = strata_name),
               size = 1)+
  scale_colour_viridis_c()
print(vis_map)






# a cyclical pattern similar to something climatic (~approximately 5 year cycle)
# shift the cycle intensity from east to west
# 
### using the a loess smooth of the realised time-series of the pacific decadal oscillation
### span of the smooth set to 10 years
loess_sm <- function(x,year,sp){
  sm <- loess(x~year,
             span = sp) 
  return(predict(sm))
}
pdo <- readRDS("data/pdo.rds")
pdo_ann <- pdo %>% 
  group_by(Year) %>% 
  summarise(annual_mean_pdo = mean(PDO)) %>% 
  filter(Year >= min_year,
         Year < 2022) %>% 
  mutate(loess_mean_pdo = loess_sm(annual_mean_pdo,Year,0.2)) %>% 
  rename(year = Year)

vis_pdo <- ggplot(data = pdo_ann,
                  aes(x = year,y = loess_mean_pdo))+
  geom_line()+
  geom_point(aes(x = year,y = annual_mean_pdo))
print(vis_pdo)


## Generate stratum smooths and intercepts ----------------------------------------
# latitudinal variation in the slopes
SLOPE_1 = 0.015# ~ 1.5%/year increase in north and stable in south

### use simple, smooth, x-y coordinate variation in betas and stratas

set.seed(2019)
random_year_effects <- rnorm(length(min_year:2021)*length(base_map$strata_name),0,0.05) 
strata_base_trajs <- base_map %>%
  sf::st_set_geometry(.,NULL) %>%
  expand_grid(.,year = min_year:2021) %>% 
  left_join(.,pdo_ann,by = "year") %>% 
  mutate(year_centered = year-1995,
         slope = SLOPE_1*sqrt(north),
         pdo_str = annual_mean_pdo*0.05 + random_year_effects,
         pdo_smooth = loess_mean_pdo*0.1*(north),
         intercept = (south*0.5),
         expected_strata = intercept + slope*year_centered + pdo_smooth  + pdo_str, 
         expected_count = exp(expected_strata))


vis_trajs <- ggplot(strata_base_trajs,
                    aes(x = year,y = expected_strata))+
  geom_line(aes(group = strata_name,colour = north))+
  scale_colour_viridis_c()

print(vis_trajs)

saveRDS(strata_base_trajs,
        "data/simulated_data_true_trajectories.rds")


#MAs <- round(log(c(1,5,10,20,50)),2)
MAs <- c(0.1,0.5,1,5,10)

# Loop Mean Abundance -----------------------------------------------------


for(ma in MAs){
  
  log_ma <- round(log(ma),2)
  ma_f <- gsub(as.character(ma),pattern = ".",replacement = "-",
               fixed = TRUE)
  
  sd_site <- 0.7
  set.seed(2019)
  site_effects <- rnorm(prep_spat_data$model_data$n_sites,0,sd_site)
  sd_obs <- 0.4
  set.seed(2019)
  observer_effects <- rnorm(prep_spat_data$model_data$n_observers,0,sd_obs)
  sd_noise <- 0.2
  set.seed(2019)
  noise_effects <- rnorm(prep_spat_data$model_data$n_counts,0,sd_noise)
  
  retrans <- 0.5*(sd_obs^2 + sd_noise^2)

  #adjust and select the strata mean log-scale expected values
    strata_base_trajs_2 <- strata_base_trajs %>% 
    mutate(expected_strata = log_ma + intercept + slope*year_centered + pdo_smooth  + pdo_str) %>% 
      select(strata_name,
             year,
             expected_strata)
  
    # # center the site and observer effects within each stratum - just to ensure that the variation in sampling
    ## through time has minimal effect on the stratum mean counts
     site_obs_df <- real_data %>% 
       select(strata_name,site,observer) %>% 
       distinct() %>% 
       group_by(strata_name) %>% 
       mutate(site_eff = as.numeric(scale(site_effects[site],scale = FALSE)),
               observer_eff = as.numeric(scale(observer_effects[observer],scale = FALSE))) #%>% 
       # group_by(strata_name) %>% 
       # summarise(test = mean(site_eff),
       #           test_obs = mean(observer_eff))
       
  # link to raw data
  sim_real_data <- real_data %>% 
    left_join(.,strata_base_trajs_2,
              by = c("strata_name","year")) %>% 
    left_join(.,site_obs_df,
              by = c("strata_name","site","observer")) %>% 
    mutate(expected_log_count = expected_strata + site_eff + observer_eff + noise_effects)
 
  mean_nuisance <-sim_real_data %>% 
    group_by(strata_name) %>% 
    summarise(mean_site = mean(site_effects[site]),
              mean_obs = mean(observer_effects[observer]),
              .groups = "keep")
  
  strata_base_trajs_2 <- strata_base_trajs_2 %>% 
    mutate(expected_mean_count = exp(expected_strata + retrans))
  
  saveRDS(strata_base_trajs_2,
          paste0("data/simulated_data_true_trajectories_",ma_f,".rds"))
  
  set.seed(2019)
  sim_counts <- rpois(prep_spat_data$model_data$n_counts,exp(sim_real_data$expected_log_count))
  
  prep_spat_data_new <- prep_spat_data
  prep_data_new <- prep_data
  

# replace counts with simulated counts ------------------------------------
  prep_spat_data_new$raw_data$count <- sim_counts
  prep_spat_data_new$model_data$count <- sim_counts
  
  saveRDS(prep_spat_data_new,paste0("data/simulated_spatial_data_mean",ma_f,".rds"))
  
  prep_data_new$raw_data$count <- sim_counts
  prep_data_new$model_data$count <- sim_counts
  saveRDS(prep_data_new,paste0("data/simulated_hier_data_mean",ma_f,".rds"))
  
  #save prepared data as rds for future modeling


}#end ma loop
