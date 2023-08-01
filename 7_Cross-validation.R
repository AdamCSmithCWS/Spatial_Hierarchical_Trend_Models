## cross-validation for bbs models
#setwd("C:/GitHub/Spatial_hierarchical_trend_models")
setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")

library(bbsBayes2)
library(tidyverse)

species <- "Eastern Whip-poor-will"
sp_n <- search_species(species)$aou

stratification <- "bbs_usgs"
model = c("gamye")

model_variants <- c("hier","spatial")

s <- stratify(stratification,species)
p <- prepare_data(s, min_n_routes = 1)
map<-load_map(stratify_by = stratification)
sp<-prepare_spatial(p,map)


#variant <- model_variants[1]



  m_spatial <- prepare_model(sp,model,
                         model_variant = "spatial",
                         calculate_cv = TRUE) 
  saveRDS(m_spatial,paste0("data/cv_prepared_",model,"_spatial_",sp_n,".rds"))
  
  m_hier <- prepare_model(sp,model,
                             model_variant = "hier",
                             calculate_cv = TRUE) 
  saveRDS(m_spatial,paste0("data/cv_prepared_",model,"_hier_",sp_n,".rds"))
         
  
  
   
  

# fitting by variant in k separate r sessions -----------------------------

  setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")
  
  library(bbsBayes2)
  library(tidyverse)
  
  species <- "Eastern Whip-poor-will"
  sp_n <- search_species(species)$aou
  
  stratification <- "bbs_usgs"
  model = c("gamye")
  
  model_variants <- c("hier","spatial")
  
k <- # setting k to one of 1:10 in each separate R session
  
  for(variant in model_variants){

m_sel <- readRDS(paste0("data/cv_prepared_",model,"_",variant,"_",sp_n,".rds"))
  m_tmp <- run_model(m_sel,
                     refresh = 500,
                     k = k,
                     save_model = FALSE)

  sum_cv <- get_summary(m_tmp,variables = "log_lik_cv")
  
  sum_cv <- sum_cv %>%
    mutate(original_count_index = m_tmp$model_data$test)
  
  saveRDS(sum_cv,paste0("output/cv_sum_",sp_n,"_",model,"_",variant,"_",k,".rds"))
  
  print(paste(model,variant,k))
  
}






