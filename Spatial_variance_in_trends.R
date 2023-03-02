## exploring variance of spatial patterns


library(bbsBayes2)
library(tidyverse)
library(patchwork)
source("functions/variance_spatial_pattern.r")
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

pdf("figures/spatial_trend_variance.pdf",
    width = 11,
    height = 8.5)
for(model in models){
  for(model_variant in model_variants){
    
    if(model == "gamye" & model_variant == "nonhier"){next}
    
    fit <- readRDS(paste0("output/",paste(species,model,model_variant,sep = "_"),".rds"))
    

    sdbeta <- fit$model_fit$draws(variables = "sdbeta",
                                           format = "draws_matrix") %>% 
      rowMeans()
    
    
# identify the draws with high and low SDbeta
    sd_quartiles <- quantile(sdbeta,c(0.05,0.95))
    
    draws_in_left_tails <- which(sdbeta < sd_quartiles[[1]])
    draws_in_right_tails <- which(sdbeta > sd_quartiles[[2]])
    

# load indices from bbsBayes2 ---------------------------------------------
    
    if(model != "gamye"){
    inds <- generate_indices(fit,
                             regions = "stratum")
    }else{
    inds <- generate_indices(fit,
                                    alternate_n = "n_smooth",
                                    regions = "stratum")
    }
    
    left_draws <- sample(draws_in_left_tails,16)
    right_draws <- sample(draws_in_right_tails,16)
    

    p_right = plot_variance_spatial_pattern(indices = inds,
                                            start_year = 1966,
                                            end_year = 1976,
                            selected_draws = draws_in_right_tails,
                            title = paste(model,model_variant,"right-tail"))
    p_left = plot_variance_spatial_pattern(indices = inds,
                                           start_year = 1966,
                                           end_year = 1976,
                                           selected_draws = draws_in_left_tails,
                           title = paste(model,model_variant,"left-tail"))
    
    print(p_left + p_right + plot_layout(guides = "collect"))
 
    
  }
  
}

dev.off()

    # select stratum indices from the identified SD draws

# calculate the trends from each of the draws

# map the trends from each of the draws

# panel or annimate the trends to visualise variation in spatial pattern








