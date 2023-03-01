
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



# indices and trends ------------------------------------------------------

inds <- generate_indices(fit)
inds_smooth <- generate_indices(fit,
                                alternate_n = "n_smooth")


trends <- generate_trends(inds_smooth)

pdf(paste0("figures/",species,"rolling_shorttrends_spatial.pdf"))
map <- plot_map(trends)
print(map)

for(ys in seq(1970,2016,by = 5)){
  tt <- generate_trends(inds_smooth,
                        min_year = ys,
                        max_year = ys + 10)
  print(plot_map(tt))
}
dev.off()






