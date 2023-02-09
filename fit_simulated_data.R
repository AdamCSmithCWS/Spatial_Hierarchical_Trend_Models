#fit simulated data


library(bbsBayes2)
library(tidyverse)

species = "simulated"

MAs <- c(0.1,0.5,1,5,10)

# Loop Mean Abundance -----------------------------------------------------
models <- c("gamye","first_diff")
model_variants <- c("hier","spatial")


for(ma in MAs[-5]){
  
  
  log_ma <- round(log(ma),2)
  ma_f <- gsub(as.character(ma),pattern = ".",replacement = "-",
               fixed = TRUE)
  
  preped_data <- readRDS(paste0("data/simulated_",model_variant,"_data_mean",ma_f,".rds"))

  preped_model <- prepare_model(preped_data,
                                  model = model,
                                model_variant = model_variant)

    fit <- run_model(preped_model,
                     refresh = 200,
                     adapt_delta = 0.9,
                     iter_warmup = 1500,
                     iter_sampling = 2000,
                     output_basename = paste(species,model,model_variant,ma_f,sep = "_"),
                     output_dir = "output")


}


inds <- generate_indices(fit,
                         regions = c("continent",
                                     "stratum",
                                     "prov_state"))
trajs <- plot_indices(inds)
trends <- generate_trends(inds,
                          min_year = 1995,
                          max_year = 2021)
map <- bbsBayes2::plot_map(trends)
print(map)

print(trajs[[1]])

geof <- plot_geofacet(inds,
                      add_observed_means = FALSE,
                      multiple = TRUE)
print(geof)
