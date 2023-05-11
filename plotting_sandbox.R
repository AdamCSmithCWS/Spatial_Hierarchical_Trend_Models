# plotting sandbox




#setwd("C:/GitHub/Spatial_hierarchical_trend_models")
setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")


# load setup --------------------------------------------------------------


library(bbsBayes2)
library(tidyverse)
library(patchwork)
#library(HDInterval)
HDL <- function(x,int,upper = TRUE){
  b <- HDInterval::hdi(x,int)
  return(ifelse(upper,b[["upper"]],b[["lower"]]))
}
source("functions/posterior_summary_functions.R")
source("functions/map_trends.R")
source("functions/trend_functions.R")

species <- "Eastern Whip-poor-will"

stratification <- "bbs_usgs"


models = c("gamye","first_diff")

model_variants <- c("nonhier","hier","spatial")
inds_save <- NULL

for(model in models){
  for(model_variant in model_variants){
    
    if(model == "gamye" & model_variant == "nonhier"){next}
    
    fit <- readRDS(paste0("output/",paste(species,model,model_variant,sep = "_"),".rds"))
    
    
    # indices and trends ------------------------------------------------------
    
    inds <- generate_indices(fit)
    saveRDS(inds,paste0("output/",paste("indices",species,model,model_variant,sep = "_")))
    inds_out <- inds$indices %>% 
      mutate(species = species,
             model = model,
             model_variant = model_variant)
    
    inds_save <- bind_rows(inds_save,inds_out)
    

    
    if(model == "gamye"){
    inds_smooth <- generate_indices(fit,
                                    alternate_n = "n_smooth")
    saveRDS(inds_smooth,paste0("output/",paste("indices_smooth",species,model,model_variant,sep = "_")))
    
    
    inds_out <- inds_smooth$indices %>% 
      mutate(species = species,
             model = model,
             model_variant = model_variant)
    inds_save <- bind_rows(inds_save,inds_out)
    
 
    
    }
    

  }
}


stratification = "bbs_usgs"
# plot the trajectories ---------------------------------------------------
realized_strata_map <- load_map(stratify_by = stratification) %>% 
  filter(strata_name %in% inds_save$region) 

strat_grid <- geofacet::grid_auto(realized_strata_map,
                                  codes = "strata_name",
                                  names = "strata_name",
                                  seed = 2019)


base_map <- realized_strata_map <- load_map(stratify_by = stratification)

for(model in models){
  
  
  inds_hier <- readRDS(paste0("output/",paste("indices",species,model,"hier",sep = "_")))
  inds_spatial <- readRDS(paste0("output/",paste("indices",species,model,"spatial",sep = "_")))
  if(model == models[2]){
    inds_nonhier <- readRDS(paste0("output/",paste("indices",species,model,"nonhier",sep = "_")))
  }
  first_years <- c(rep(min(inds_hier$indices$year),2),seq(1971,2011,by = 5))

  pdf(paste("figures/trend_maps_",model,".pdf"),
      width = 11,
      height = 8.5)
  
  tmp1 <- inds_hier$indices %>% 
    filter(region_type == "stratum") %>% 
    mutate(model_variant = "hier")
  
  tmp2 <- inds_spatial$indices %>% 
    filter(region_type == "stratum") %>% 
    mutate(model_variant = "spatial")
  
   
  tmp <- bind_rows(tmp1,tmp2)
  
  if(model == models[2]){
    tmp3 <- inds_nonhier$indices %>% 
      filter(region_type == "stratum") %>% 
      mutate(model_variant = "nonhier")
    
    
    tmp <- bind_rows(tmp3,tmp)
    
  }
  geo_hier <- ggplot(data = tmp,
                     aes(x = year,y = index))+
    geom_line(aes(colour = model_variant))+
    geom_ribbon(aes(ymin = index_q_0.05,ymax = index_q_0.95,
                    fill = model_variant),
                alpha = 0.3)+
    scale_y_continuous(limits = c(0,NA))+
    scale_colour_viridis_d(end = 0.75,begin = 0.2,
                           aesthetics = c("colour","fill"))+
    #geom_point(aes(y = obs_mean,x = year),alpha = 0.3)+
    geofacet::facet_geo(~region,grid = strat_grid,scales = "free_y")+
    theme_minimal()
  
  print(geo_hier)
  
 tmpsimple <- tmp %>% 
   sf::st_drop_geometry()

   spaghetti1 <- ggplot(data = tmpsimple,
                        aes(x = year,y = index,
                            group = region))+
     geom_line(aes(colour = (n_routes_total)))+
     scale_colour_viridis_c()+
     scale_y_continuous(trans = "log10")+
     facet_wrap(vars(model_variant))
     
  print(spaghetti1)
  
   
   for(j in 1:length(first_years)){
     
    ys <- first_years[j]
    if(j != 1){
    ye <- ys+10
    }else{
      ye <- 2021
    }
    tt_hier <- generate_trends(inds_hier,
                          min_year = ys,
                          max_year = ye,
                          slope = TRUE)
    
    tt_spatial <- generate_trends(inds_spatial,
                               min_year = ys,
                               max_year = ye,
                               slope = TRUE)


    
  map_hier <- map_trends(tt_hier,
                         title = paste("Hierarchical",model))
    
  map_spatial <- map_trends(tt_spatial,
                            title = paste("Spatial",model))

        traj_hier <- plot_indices(inds_hier,
                                  title = FALSE,
                                  axis_title_size = 9,
                                  axis_text_size = 8)
        traj_spatial <- plot_indices(inds_spatial,
                                     title = FALSE,
                                     axis_title_size = 9,
                                     axis_text_size = 8)
        
        if(model == models[2]){
          
          tt_nonhier <- generate_trends(inds_nonhier,
                                        min_year = ys,
                                        max_year = ye,
                                        slope = TRUE)
          
          
          map_nonhier <- map_trends(tt_nonhier,
                                 title = paste("Non-Hierarchical",model))
          
          traj_nonhier <- plot_indices(inds_nonhier,
                                    title = FALSE,
                                    axis_title_size = 9,
                                    axis_text_size = 8)
          
          
          print(traj_nonhier[["continent"]] + traj_hier[["continent"]] + traj_spatial[["continent"]] +
                  map_nonhier + map_hier + map_spatial + 
                  plot_layout(ncol = 3,nrow = 2,
                              byrow = TRUE,
                              heights = c(1,4),
                              guides = "collect"))
        }else{
        
        
        
        print(traj_hier[["continent"]] + traj_spatial[["continent"]] +
                map_hier + map_spatial + 
                plot_layout(ncol = 2,nrow = 2,
                            byrow = TRUE,
                            heights = c(1,4),
                            guides = "collect"))
        
        }
        
  }
  dev.off()
  


  
} 
  
 
  




  


# CBC and Shorebirds ------------------------------------------------------

  data_sets <- c("CBC","Shorebird")
  models <- c("gamye","first_diff")
  model_variants <- c("hier","spatial")
  
## regional and overall trajectories
  
  for(data_set in data_sets){

    if(data_set == "Shorebird"){
      load(paste0("data/Red_Knot",data_set,"_data.RData"))
      data_prep <- data_1
      map <- realized_strata_map %>% 
        rename(strata_name = hex_name,
               stratum = stratn)
      
      rm(list = c("stan_data",
                  "neighbours",
                  "realized_strata_map",
                  "data_1"))
      
      strat_df <- data_prep %>% 
        select(hex_name,stratn) %>% 
        rename(stratum = stratn,
               strata_name = hex_name) %>% 
        distinct() 
        
      yrs <- data_prep %>% 
        select(yr,year) %>% 
        # rename(yr = year_vec,
        #        year = YearCollected) %>% 
        distinct() %>% 
        arrange(yr)
      
      
    }

   if(data_set == "CBC"){
     
     data_prep <- readRDS(paste0("output/dataframe_",data_set,"_spatial_gamye.rds"))
     
     strat_df <- data_prep %>% 
       select(strata_name,strata_vec,non_zero,area_sq_km) %>% 
       rename(stratum = strata_vec) %>% 
       distinct() %>% 
       mutate(area_weight = area_sq_km/sum(area_sq_km))
     
     yrs <- data_prep %>% 
       select(year_vec,count_year) %>% 
       rename(yr = year_vec,
              year = count_year) %>% 
       distinct()
     
     map <- strata_map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
       select(-area_sq_km) %>% 
      inner_join(.,strat_df,
                 by = "strata_name") %>% 
      arrange(stratum)
   }
    

# indices ------------------------------------------------------
    
    for(model in models){
      if(model == "first_diff" & data_set == "Shorebird"){next}
      for(model_variant in model_variants){
        
        fit <- readRDS(paste0("output/fit_",data_set,"_",model_variant,"_",model,".rds"))
        
        ind_samples <- posterior_samples(fit = fit,
                                         parm = "n",
                                         dims = c("stratum","yr")) %>% 
          inner_join(.,strat_df,
                     by = "stratum") %>% 
          inner_join(.,yrs,
                     by = "yr")
        
        saveRDS(ind_samples,paste0("output/",data_set,"_indices_",model,"_",model_variant,".rds"))
        
        if(data_set == "Shorebird"){
          N_samples <- posterior_samples(fit = fit,
                                         parm = "N",
                                         dims = c("yr")) %>% 
            inner_join(.,yrs,
                       by = "yr")
          
          saveRDS(N_samples,paste0("output/",data_set,"_continent_indices_",model,"_",model_variant,".rds"))
        }
        
        if(model == "gamye"){
          
          ind_samples_smooth <- posterior_samples(fit = fit,
                                           parm = "nsmooth",
                                           dims = c("stratum","yr")) %>% 
            inner_join(.,strat_df,
                       by = "stratum") %>% 
            inner_join(.,yrs,
                       by = "yr")
          
          saveRDS(ind_samples_smooth,paste0("output/",data_set,"_smooth_indices_",model,"_",model_variant,".rds"))
          
          if(data_set == "Shorebird"){
            N_samples <- posterior_samples(fit = fit,
                                             parm = "NSmooth",
                                             dims = c("yr")) %>% 
              inner_join(.,yrs,
                         by = "yr")
            
            saveRDS(N_samples,paste0("output/",data_set,"_continent_smooth_indices_",model,"_",model_variant,".rds"))
          }
          
          
        }
        
        
        
        
      }
      
    }
    
    
  }
  
  
  
  

# trajectories and trends -------------------------------------------------

  
  inds_out <- NULL
  trends_out <- NULL
  
  for(data_set in data_sets){
   

    if(data_set == "CBC"){
      data_prep <- readRDS(paste0("output/dataframe_",data_set,"_spatial_gamye.rds"))
      
      strat_df <- data_prep %>% 
        select(strata_name,strata_vec,non_zero,area_sq_km) %>% 
        rename(stratum = strata_vec) %>% 
        distinct() %>% 
        mutate(area_weight = area_sq_km/sum(area_sq_km))
      
      yrs <- data_prep %>% 
        select(year_vec,count_year) %>% 
        rename(yr = year_vec,
               year = count_year) %>% 
        distinct()
      

      #first_years <- c(min(yrs$year),floor(median(yrs$year)),seq(1970,2005,by = 5),max(yrs$year)-10)
      first_years <- c(rep(min(yrs$year),2),floor(median(yrs$year)),seq(1970,2005,by = 1),max(yrs$year)-10)
      
      
      map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
        select(-area_sq_km) %>% 
        inner_join(.,strat_df,
                   by = "strata_name") %>% 
        arrange(stratum)
    }
    
    if(data_set == "Shorebird"){
      load(paste0("data/Red_Knot",data_set,"_data.RData"))
      data_prep <- data_1
      map <- realized_strata_map %>% 
        rename(strata_name = hex_name,
               stratum = stratn)
      
      rm(list = c("stan_data",
                  "neighbours",
                  "realized_strata_map",
                  "data_1"))
      
      strat_df <- data_prep %>% 
        select(hex_name,stratn) %>% 
        rename(stratum = stratn,
               strata_name = hex_name) %>% 
        distinct() 
      
      yrs <- data_prep %>% 
        select(yr,year) %>% 
        # rename(yr = year_vec,
        #        year = YearCollected) %>% 
        distinct() %>% 
        arrange(yr)
      
      first_years <- c(rep(min(yrs$year),2),floor(median(yrs$year)),seq(1980,2005,by = 5),max(yrs$year)-10)
      
    }
    
    

    for(model in models){
      if(model == "first_diff" & data_set == "Shorebird"){next}
      
      for(model_variant in model_variants){
        
        ind_samples <- readRDS(paste0("output/",data_set,"_indices_",model,"_",model_variant,".rds"))
       
        
        
  inds_strat <- ind_samples %>%
    group_by(strata_name,year) %>% 
    summarise(index = median(.value),
              index_q_0.05 = HDL(.value,0.9,upper = FALSE),
              index_q_0.95 = HDL(.value,0.9,upper = TRUE),
              .groups = "keep") %>% 
    mutate(index_type = "full")
  
  
  

  
  if(model == "gamye"){
    ind_samples_smooth <- readRDS(paste0("output/",data_set,"_smooth_indices_",model,"_",model_variant,".rds"))
  
    inds_strat_smooth <- ind_samples_smooth %>%
      group_by(strata_name,year) %>% 
      summarise(index = median(.value),
                index_q_0.05 = HDL(.value,0.9,upper = FALSE),
                index_q_0.95 = HDL(.value,0.9,upper = TRUE),
                .groups = "keep") %>% 
      mutate(index_type = "smooth")
    
    inds_strat <- bind_rows(inds_strat,inds_strat_smooth)
    
  }
  
  if(model == "gamye"){
    ind_trend <- ind_samples_smooth
  }else{
    ind_trend <- ind_samples
  }

  
  trends_strata <- NULL
  
  for(j in 1:length(first_years)){
    
    ys <- first_years[j]
    
    if(j %in% c(1,3)){
      ye <- 2019
    }
    if(j == 2){
      ye <- first_years[3]
    }
    if(j > 3){
      ye <- ys+10
    }
    nyrs <- ye-ys
    
    trend_tmp <- ind_trend %>% 
    filter(year %in% c(ys,ye)) %>% 
    ungroup() %>% 
    select(-matches(match = "yr",ignore.case = FALSE)) %>% 
    pivot_wider(names_from = year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",ys),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",ye),.x,
                         fixed = TRUE))%>% 
    group_by(.draw,strata_name) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end/start,ny = nyrs),
              ch = chng(end/start),
              .groups = "keep") %>% 
    group_by(strata_name) %>% 
    summarise(trend = mean(t),
              lci = quantile(t,0.025,names = FALSE),
              uci = quantile(t,0.975,names = FALSE),
              width_CI = uci-lci) %>% 
    mutate(model = model,
           model_variant = model_variant,
           start_year = ys,
           end_year = ye,
           trend_length = nyrs,
           data_set = data_set)
    
  trends_strata <- bind_rows(trends_strata, trend_tmp)
  }
  
trends_out <- bind_rows(trends_out,trends_strata)
  
  if(data_set == "CBC"){
  inds_comp <- ind_samples %>%
    mutate(.value = .value*area_weight) %>% 
    group_by(.draw,year) %>% #draw-wise summary of area-weighted indices
    summarise(.value = sum(.value),
              .groups = "drop") %>% 
    group_by(year) %>% 
    summarise(index = median(.value),
              index_q_0.05 = HDL(.value,0.9,upper = FALSE),
              index_q_0.95 = HDL(.value,0.9,upper = TRUE),
              .groups = "keep") %>% 
    mutate(strata_name = "continent") %>% 
    mutate(index_type = "full")

  
  ind_comp_trend <- ind_samples %>%
    mutate(.value = .value*area_weight) %>% 
    group_by(.draw,year) %>% #draw-wise summary of area-weighted indices
    summarise(.value = sum(.value),
              .groups = "drop")
  
  
  if(model == "gamye"){
    
    
    inds_comp_smooth <- ind_samples_smooth %>%
      mutate(.value = .value*area_weight) %>% 
      group_by(.draw,year) %>% #draw-wise summary of area-weighted indices
      summarise(.value = sum(.value),
                .groups = "drop") %>% 
      group_by(year) %>% 
      summarise(index = median(.value),
                index_q_0.05 = HDL(.value,0.9,upper = FALSE),
                index_q_0.95 = HDL(.value,0.9,upper = TRUE),
                .groups = "keep") %>% 
      mutate(strata_name = "continent") %>% 
      mutate(index_type = "smooth")
    
    inds_comp <- bind_rows(inds_comp,inds_comp_smooth)
    
    ind_comp_trend <- ind_samples_smooth %>%
      mutate(.value = .value*area_weight) %>% 
      group_by(.draw,year) %>% #draw-wise summary of area-weighted indices
      summarise(.value = sum(.value),
                .groups = "drop")
  }
  
  }else{
    
    
    comp_samples <- readRDS(paste0("output/",data_set,"_continent_indices_",model,"_",model_variant,".rds"))
    
    
    
    inds_comp <- comp_samples %>%  
      group_by(year) %>% 
      summarise(index = median(.value),
                index_q_0.05 = HDL(.value,0.9,upper = FALSE),
                index_q_0.95 = HDL(.value,0.9,upper = TRUE),
                .groups = "keep") %>% 
      mutate(strata_name = "continent") %>% 
      mutate(index_type = "full")
    
 
      comp_samples <- readRDS(paste0("output/",data_set,"_continent_smooth_indices_",model,"_",model_variant,".rds"))
      
      
      
      inds_comp_smooth <- comp_samples %>%
        group_by(year) %>% 
        summarise(index = median(.value),
                  index_q_0.05 = HDL(.value,0.9,upper = FALSE),
                  index_q_0.95 = HDL(.value,0.9,upper = TRUE),
                  .groups = "keep") %>% 
        mutate(strata_name = "continent") %>% 
        mutate(index_type = "smooth")
    
      inds_comp <- bind_rows(inds_comp,inds_comp_smooth)
    
      ind_comp_trend <- comp_samples
      
  }
  inds_all <- bind_rows(inds_comp,
                        inds_strat) %>% 
    mutate(model = model,
           model_variant = model_variant,
           data_set = data_set)
    
  
  inds_out <- bind_rows(inds_out,inds_all)

  

# composite trends --------------------------------------------------------

  trends_comp <- NULL
  
  for(j in 1:length(first_years)){
    
    ys <- first_years[j]
    
    if(j %in% c(1,3)){
      ye <- 2019
    }
    if(j == 2){
      ye <- first_years[3]
    }
    if(j > 3){
      ye <- ys+10
    }
    
    nyrs <- ye-ys
    trend_tmp <- ind_comp_trend %>% 
      filter(year %in% c(ys,ye)) %>% 
      ungroup() %>% 
      select(-matches(match = "yr",ignore.case = FALSE)) %>% 
      pivot_wider(names_from = year,
                  values_from = .value,
                  names_prefix = "Y") %>% 
      rename_with(., ~gsub(replacement = "start",
                           pattern = paste0("Y",ys),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "end",
                           pattern = paste0("Y",ye),.x,
                           fixed = TRUE))%>% 
      group_by(.draw) %>% 
      summarise(end = sum(end),
                start = sum(start),
                t = texp(end/start,ny = nyrs),
                ch = chng(end/start),
                .groups = "drop") %>% 
      summarise(trend = mean(t),
                lci = quantile(t,0.025,names = FALSE),
                uci = quantile(t,0.975,names = FALSE),
                width_CI = uci-lci) %>% 
      mutate(model = model,
             model_variant = model_variant,
             start_year = ys,
             end_year = ye,
             trend_length = nyrs,
             data_set = data_set,
             strata_name = "continent")
    
    trends_comp <- bind_rows(trends_comp, trend_tmp)
  }
  
  trends_out <- bind_rows(trends_out,trends_comp)

      }
      
    }#model
    
   
    
  }#data_set
  
  saveRDS(inds_out,"output/All_CBC_and_shorebird_indices.rds")
  saveRDS(trends_out,"output/All_CBC_and_shorebird_trends.rds")
  
  
  

# index plots -------------------------------------------------------------

  
  inds_out <- readRDS("output/All_CBC_and_shorebird_indices.rds")  
tmp <- inds_out %>% 
  filter(strata_name == "continent")

tp = ggplot(data = tmp,
            aes(x = year, y = index))+
  geom_ribbon(aes(ymin = index_q_0.05,
                  ymax = index_q_0.95,
                  fill = paste(model,model_variant, index_type)),
              alpha = 0.1)+
  geom_line(aes(colour = paste(model,model_variant, index_type)))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))+
  facet_wrap(vars(data_set),
             scales = "free_y")

print(tp)
 

data_prep <- readRDS(paste0("output/dataframe_","CBC","_spatial_gamye.rds"))
obs_mean <- data_prep %>% 
  select(strata_name,count_year,how_many) %>% 
  group_by(strata_name,count_year) %>% 
  summarise(mean_obs = mean(how_many),
            .groups = "keep")

tmp <- inds_out %>% 
  filter(strata_name != "continent",
         data_set == "CBC",
         index_type == "full")

tp = ggplot(data = tmp,
            aes(x = year, y = index))+
  geom_ribbon(aes(ymin = index_q_0.05,
                  ymax = index_q_0.95,
                  fill = paste(model,model_variant)),
              alpha = 0.1)+
  geom_line(aes(colour = paste(model,model_variant)))+
  geom_point(data = obs_mean,
             aes(x = count_year,y = mean_obs),
             inherit.aes = FALSE,
             alpha = 0.2)+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))+
  coord_cartesian(ylim = c(0,40))+
  facet_wrap(vars(strata_name),
             scales = "free_y")

print(tp)





# Map Trends --------------------------------------------------------------
trends_out <- readRDS("output/All_CBC_and_shorebird_trends.rds")

pdf("Figures/all_CBC_Shorebird_trend_maps.pdf",
    width = 11,
    height = 8.5)

for(data_setname in data_sets){

if(data_setname == "Shorebird"){
  load(paste0("data/Red_Knot",data_setname,"_data.RData"))
  data_prep <- data_1
  map <- realized_strata_map %>% 
    rename(strata_name = hex_name,
           stratum = stratn)
  
  rm(list = c("stan_data",
              "neighbours",
              "realized_strata_map",
              "data_1"))
  
  strat_df <- data_prep %>% 
    select(hex_name,stratn) %>% 
    rename(stratum = stratn,
           strata_name = hex_name) %>% 
    distinct() 
  
  yrs <- data_prep %>% 
    select(yr,year) %>% 
    # rename(yr = year_vec,
    #        year = YearCollected) %>% 
    distinct() %>% 
    arrange(yr)
  
  
}

if(data_setname == "CBC"){
  
  data_prep <- readRDS(paste0("output/dataframe_",data_setname,"_spatial_gamye.rds"))
  
  strat_df <- data_prep %>% 
    select(strata_name,strata_vec,non_zero,area_sq_km) %>% 
    rename(stratum = strata_vec) %>% 
    distinct() %>% 
    mutate(area_weight = area_sq_km/sum(area_sq_km))
  
  yrs <- data_prep %>% 
    select(year_vec,count_year) %>% 
    rename(yr = year_vec,
           year = count_year) %>% 
    distinct()
  
  map <- strata_map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
    select(-area_sq_km) %>% 
    inner_join(.,strat_df,
               by = "strata_name") %>% 
    arrange(stratum)
}



for(modelsel in models){
  
  yrpairs <- trends_out %>% 
    filter(data_set == as.character(data_setname)) %>% 
    select(start_year,end_year) %>% 
    distinct() %>% 
    arrange(start_year)
  
  for(j in 1:nrow(yrpairs)){
    sy <- as.integer(yrpairs[j,"start_year"])
    ey <- as.integer(yrpairs[j,"end_year"])
    
 model_variantsel <- "hier"   
    
    trend_tmp1 <- trends_out %>% 
      filter(model == modelsel,
             model_variant == model_variantsel,
             data_set == data_setname,
             strata_name != "continent",
             start_year == sy,
             end_year == ey)
    m1 <- map_trends(trends = trend_tmp1,
               base_map_blank = map,
               title = paste(data_setname,model_variantsel,modelsel),
               region_name = "strata_name")
    
    model_variantsel <- "spatial"
    trend_tmp2 <- trends_out %>% 
      filter(model == modelsel,
             model_variant == model_variantsel,
             data_set == data_setname,
             strata_name != "continent",
             start_year == sy,
             end_year == ey)
    m2 <- map_trends(trend_tmp2,
                     base_map_blank = map,
                     title = paste(data_setname,model_variantsel,modelsel),
                     region_name = "strata_name")
    
    print(m1 + m2 + plot_layout(guides = "collect"))
  }
  }
  
}

dev.off()




























