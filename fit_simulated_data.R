#fit simulated data

setwd("C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")
library(bbsBayes2)
library(tidyverse)

species = "simulated"

MAs <- c(0.1,0.5,1,5,10)

# Loop Mean Abundance -----------------------------------------------------
models <- c("gamye","first_diff")
model_variants <- c("nonhier","hier","spatial")
model <- models[2]
model_variant <- model_variants[1]

for(ma in MAs[c(5,1,3,2,4)]){
  
  if(model == "gamye" & model_variant == "nonhier"){next}
  
  log_ma <- round(log(ma),2)
  ma_f <- gsub(as.character(ma),pattern = ".",replacement = "-",
               fixed = TRUE)
  
  if(model_variant == "nonhier"){
    preped_data <- readRDS(paste0("data/simulated_","hier","_data_mean",ma_f,".rds"))
    
  }else{
  preped_data <- readRDS(paste0("data/simulated_",model_variant,"_data_mean",ma_f,".rds"))
}
  preped_model <- prepare_model(preped_data,
                                  model = model,
                                model_variant = model_variant)

    fit <- run_model(preped_model,
                     refresh = 400,
                     adapt_delta = 0.9,
                     iter_warmup = 1500,
                     iter_sampling = 2000,
                     max_treedepth = 12,
                     output_basename = paste(species,model,model_variant,ma_f,sep = "_"),
                     output_dir = "output")


}


