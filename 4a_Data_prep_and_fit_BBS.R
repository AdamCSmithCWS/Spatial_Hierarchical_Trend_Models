
setwd("C:/GitHub/Spatial_hierarchical_trend_models")
setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")

library(bbsBayes2)
library(tidyverse)

species <- "Eastern Whip-poor-will"
#species <- "Pacific Wren"

stratification <- "bbs_usgs"
models = c("gamye","first_diff")

model_variants <- c("nonhier","hier","spatial")


for(model in models){
  for(model_variant in model_variants){
    if(model == "gamye" & model_variant == "nonhier"){next}
s <- stratify(by = stratification,
              species = species)


p <- prepare_data(s)

if(model_variant == "spatial"){
ps <- prepare_spatial(p,
                      strata_map = load_map(stratification))
print(ps$spatial_data$map)

}else{
  ps <- p
}


pm <- prepare_model(ps,
                    model = model,
                    model_variant = model_variant)

fit <- run_model(pm,
                 refresh = 400,
                 adapt_delta = 0.8,
                 max_treedepth = 11,
                 output_dir = "output",
                 output_basename = paste(species,model,model_variant,sep = "_"))



summ <- get_summary(fit)

summ <- summ %>% 
  mutate(variable_type = stringr::str_extract(variable, "^\\w+"))

saveRDS(summ,paste0("output/",model,"_",model_variant,"_param_summary.rds"))


  }
}












# Model summaries ---------------------------------------------------------


rhat_ess_summ <- summ %>% 
  group_by(variable_type) %>% 
  summarise(min_ess = min(ess_bulk),
            max_rhat = max(rhat),
            med_ess = median(ess_bulk),
            med_rhat = median(rhat))

sds <- summ %>% filter(grepl("sd",variable))

beta_raw <- summ %>% filter(variable_type == "beta_raw")

ess_fail <- summ %>% filter(ess_bulk < 1000)

ess_fail_summ <- ess_fail %>% 
  group_by(variable_type) %>% 
  summarise(min_ess = min(ess_bulk),
            max_rhat = max(rhat))

rhat_fail <- summ %>% filter(rhat > 1.01)

rhat_fail_summ <- rhat_fail %>% 
  group_by(variable_type) %>% 
  summarise(min_ess = min(ess_bulk),
            max_rhat = max(rhat))



# Wood Thrush fine scale demonstration and eBird comparison ------------------------------------

setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")

library(bbsBayes2)
library(tidyverse)

species <- "Wood Thrush"
aou <- search_species(species)$aou #numerical species indicator from BBS database

stratification <- "latlong"
models = c("gam","first_diff","gamye")

model_variant <- "spatial"

for(model in models){
  
  
  
  s <- stratify(by = stratification,
                species = species)
  
  
  p <- prepare_data(s,
                    min_year = 2000)
  
    ps <- prepare_spatial(p,
                          strata_map = load_map(stratification))
    print(ps$spatial_data$map)
    
  
  
  
  pm <- prepare_model(ps,
                      model = model,
                      model_variant = model_variant)
  
  fit <- run_model(pm,
                   refresh = 400,
                   adapt_delta = 0.8,
                   iter_warmup = 1000,
                   iter_sampling = 4000,
                   thin = 2,
                   output_dir = "output",
                   output_basename = paste(species,model,model_variant,sep = "_"))
  
  
  
  summ <- get_summary(fit)
  
  summ <- summ %>% 
    mutate(variable_type = stringr::str_extract(variable, "^\\w+"))
  
  saveRDS(summ,paste0("output/",aou,"_",model,"_",model_variant,"_param_summary.rds"))
  
  
  
}


if(model == "gamye"){
  inds <- generate_indices(fit,
                           alternate_n = "n_smooth")
}else{
inds <- generate_indices(fit)
}

trends <- generate_trends(inds,
                          min_year = 2007,
                          max_year = 2021)

map <- plot_map(trends)

print(map)




# WOTH non-spatial latlong ------------------------------------------------



s <- stratify(by = stratification,
              species = species)


p <- prepare_data(s,
                  min_year = 2000)


model <- "gamye"
model <- "first_diff"

pm <- prepare_model(p,
                    model = model, 
                    model_variant = "hier")

fit <- run_model(pm,
                 refresh = 400,
                 adapt_delta = 0.8,
                 iter_warmup = 1000,
                 iter_sampling = 4000,
                 thin = 2,
                 output_dir = "output",
                 output_basename = paste(species,model,"hier",sep = "_"))



summ <- get_summary(fit)

summ <- summ %>% 
  mutate(variable_type = stringr::str_extract(variable, "^\\w+"))

saveRDS(summ,paste0("output/",aou,"_",model,"_","hier","_param_summary.rds"))



