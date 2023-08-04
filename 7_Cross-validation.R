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
  
  m_hier <- prepare_model(p,model,
                             model_variant = "hier",
                             calculate_cv = TRUE) 
  saveRDS(m_hier,paste0("data/cv_prepared_",model,"_hier_",sp_n,".rds"))
         
  
  
   
  

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
  
  for(variant in model_variants[1]){

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



# compile cv results ------------------------------------------------------

sum_cv <- NULL
for(variant in model_variants){
  m_sel <- readRDS(paste0("data/cv_prepared_",model,"_",variant,"_",sp_n,".rds"))
  raw_data <- m_sel$raw_data %>% 
    mutate(folds = m_sel$folds,
           original_count_index = row_number())
  sum_cv1 <- NULL
  for(k in 1:10){
    
    sum_cv_tmp <- readRDS(paste0("output/cv_sum_",sp_n,"_",model,"_",variant,"_",k,".rds"))
    
    sum_cv1 <- bind_rows(sum_cv1,sum_cv_tmp)
    
  }
  
  
    
  sum_cv1 <- sum_cv1 %>% 
    inner_join(.,raw_data,
               by = "original_count_index") %>% 
    mutate(model_variant = variant,
           model = model)

  cv_wide_t <- sum_cv1 %>% 
    select(mean,original_count_index,
           strata_name,year,route,count) %>% 
    rename_with(.,~gsub(x = .x,"mean",paste(variant,"lppd",sep = "_")))
  
  sum_cv <- bind_rows(sum_cv,sum_cv1)

  if(variant == model_variants[1]){
    cv_wide <- cv_wide_t
  }else{
  cv_wide <- left_join(cv_wide,cv_wide_t)
  }
  
}


totals <- sum_cv %>% 
  group_by(model_variant) %>% 
  summarise(sum_lppd = sum(mean),
            mean_lppd = mean(mean),
            sd_lppd = sd(mean))

diffs <- cv_wide %>%
  mutate(diff_lppd = spatial_lppd - hier_lppd)

diffs %>% 
  summarise(mean = mean(diff_lppd),
            sd = sd(diff_lppd),
            z = mean/sqrt(sd))



