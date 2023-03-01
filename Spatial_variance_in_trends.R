## exploring variance of spatial patterns


library(bbsBayes2)
library(tidyverse)
library(patchwork)

species <- "Eastern Whip-poor-will"

stratification <- "bbs_usgs"

base_map <- load_map(stratification)


# palette -----------------------------------------------------------------


breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
labls <- c(paste0("< ", breaks[1]),
           paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),
           paste0("> ",breaks[length(breaks)]))
labls <- paste0(labls, " %")



models = c("gamye","first_diff")

model_variants <- c("nonhier","hier","spatial")

# load saved model object from bbsBayes2 ----------------------------------

model <- models[1]
model_variant <- model_variants[3]

for(model in models){
  for(model_variant in model_variants){
    
    if(model == "gamye" & model_variant == "nonhier"){next}
    
    fit <- readRDS(paste0("output/",paste(species,model,model_variant,sep = "_"),".rds"))
    

    sdbeta <- fit$model_fit$draws(variables = "sdbeta",
                                           format = "draws_matrix")  
    
    
# identify the draws with high and low SDbeta
    sd_quartiles <- quantile(sdbeta,c(0.10,0.9))
    
    draws_in_left_tails <- which(sdbeta < sd_quartiles[[1]])
    draws_in_right_tails <- which(sdbeta > sd_quartiles[[2]])
    

# load indices from bbsBayes2 ---------------------------------------------
    
    
    inds <- generate_indices(fit,
                             regions = "stratum")
    
    inds_smooth <- generate_indices(fit,
                                    alternate_n = "n_smooth",
                                    regions = "stratum")
    
    
    left_draws <- sample(draws_in_left_tails,16)
    right_draws <- sample(draws_in_right_tails,16)
    

    
    
    plot_function <- function(indices = inds_smooth,
                              selected_draws = right_draws,
                              end_year = 2021,
                              start_year = 1966,
                              base_map_blank = base_map){
      
    max_year_num <- 1+(end_year-indices$meta_data$start_year)
    min_year_num <- 1+(start_year-indices$meta_data$start_year)
    nyears = max_year_num - min_year_num
    
    trends_plot <- purrr::map_dfr(indices$sample, 
                                       .f = ~ 100*(((.x[selected_draws,max_year_num] / .x[selected_draws,min_year_num])^(1/nyears))-1),
                                       .id = "strata_name") %>% 
      pivot_longer(.,cols = -strata_name,
                   names_to = ".draw") %>% 
      mutate(strata_name = gsub("stratum_","",
                               strata_name),
             t_plot = cut(value, breaks = c(-Inf, breaks, Inf),
                          labels = labls))
  
    
    pal <- setNames(
      c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
        "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"),
      levels(trends_plot$t_plot))
    
    
    trend_map <- base_map_blank %>% 
      right_join(.,trends_plot,
                 by = "strata_name",
                 multiple = "all")
    
  plot_left <- ggplot()+
    geom_sf(data = trend_map,
            aes(fill = t_plot))+
    scale_fill_manual(values = pal)+
    theme_minimal()+
    facet_wrap(vars(.draw))
  
  return(plot_left)
  
    }
    
    p_right = plot_function(start_year = 2011)
    p_left = plot_function(selected_draws = left_draws,
                           start_year = 2011)
    # select stratum indices from the identified SD draws

# calculate the trends from each of the draws

# map the trends from each of the draws

# panel or annimate the trends to visualise variation in spatial pattern








