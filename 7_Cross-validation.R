## cross-validation for bbs models
#setwd("C:/GitHub/Spatial_hierarchical_trend_models")
#setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")

library(bbsBayes2)
library(tidyverse)


comparisons <- data.frame(species = c("Wood Thrush",
                                      "Eastern Whip-poor-will",
                                      "Eastern Whip-poor-will",
                                      "Carolina Wren",
                                      "Baird's Sparrow",
                                      "Baird's Sparrow",
                                      "Scissor-tailed Flycatcher",
                                      "Scissor-tailed Flycatcher",
                                      "Rufous Hummingbird",
                                      "Rufous Hummingbird",
                                      "Bewick's Wren",
                                      "Mountain Bluebird"),
                          sp_n = c(7550,
                                   4171,
                                   4171,
                                   7180,
                                   5450,
                                   5450,
                                   4430,
                                   4430,
                                   4330,
                                   4330,
                                   7190,
                                   7680),
                          model = c("gamye",
                                    "gamye",
                                    "first_diff",
                                    "first_diff",
                                    "first_diff",
                                    "gamye",
                                    "first_diff",
                                    "gamye",
                                    "first_diff",
                                    "gamye",
                                    "first_diff",
                                    "gamye"),
                          stratification = c("bbs_usgs",
                                             "bbs_usgs",
                                             "bbs_usgs",
                                             "bbs_usgs",
                                             "latlong",
                                             "latlong",
                                             "latlong",
                                             "latlong",
                                             "latlong",
                                             "latlong",
                                             "bbs_usgs",
                                             "bbs_usgs"))



# Data preparation for Cross-validation -----------------------------
# wd setting may be necessary to run outside of RStudio for more stable behaviour
#setwd( "C:/Users/SmithAC/Documents/GitHub/Spatial_Hierarchical_Trend_Models")

model_variants <- c("hier","spatial")

for(j in 1:nrow(comparisons)){
  
  
species <- comparisons[j,"species"]
sp_n <- comparisons[j,"sp_n"]
stratification <- comparisons[j,"stratification"]
model <- comparisons[j,"model"]

s <- stratify(stratification,species,
              release = 2022)
p <- prepare_data(s, min_n_routes = 1)
map<-load_map(stratify_by = stratification)
sp<-prepare_spatial(p,map)


#variant <- model_variants[1]


# 
  m_spatial <- prepare_model(sp,model,
                         model_variant = "spatial",
                         calculate_cv = TRUE)
  saveRDS(m_spatial,paste0("data/cv_prepared_",model,"_spatial_",sp_n,".rds"))

  m_hier <- prepare_model(p,model,
                             model_variant = "hier",
                             calculate_cv = TRUE)
  saveRDS(m_hier,paste0("data/cv_prepared_",model,"_hier_",sp_n,".rds"))


  m_spatialfull <- prepare_model(sp,model,
                             model_variant = "spatial",
                             calculate_cv = FALSE) 
  saveRDS(m_spatialfull,paste0("data/cv_prepared_",model,"_spatial_full_",sp_n,".rds"))
  
  m_hierfull <- prepare_model(p,model,
                          model_variant = "hier",
                          calculate_cv = FALSE) 
  saveRDS(m_hierfull,paste0("data/cv_prepared_",model,"_hier_full_",sp_n,".rds"))
  
  
  
}
  
   
  

# Fitting loop by species model variant -----------------------------------



for(j in 1:nrow(comparisons)){
  
  
  species <- comparisons[j,"species"]
  sp_n <- comparisons[j,"sp_n"]
  stratification <- comparisons[j,"stratification"]
  model <- comparisons[j,"model"]
  

#k <- # setting k to one of 1:10 in each separate R session
 for(k in 1:10){ 
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

}

  
}


# compile cv results ------------------------------------------------------


# Loop over species and model comparisons ---------------------------------




model_variants <- c("hier","spatial")
sum_cv <- NULL
cv_wide <- NULL
for(i in 1:nrow(comparisons)){

  species <- comparisons[i,"species"]
  sp_n <- comparisons[i,"sp_n"]
  model <- comparisons[i,"model"]
  

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
           model = model,
           species = species,
           sp_n = sp_n)

  cv_wide_t <- sum_cv1 %>% 
    select(mean,original_count_index,
           strata_name,year,route,count) %>% 
    rename_with(.,~gsub(x = .x,"mean",paste(variant,"lppd",sep = "_"))) %>% 
    mutate(species = species,
           model = model,
           sp_n = sp_n)
  
  sum_cv <- bind_rows(sum_cv,sum_cv1)

  if(variant == model_variants[1]){
    cv_wide_s <- cv_wide_t
  }else{
  cv_wide_s <- left_join(cv_wide_s,cv_wide_t)
  }
  
}

cv_wide <- bind_rows(cv_wide,cv_wide_s)

 

}#end i loop



# summed lppd by species model variant
totals <- sum_cv %>% 
  group_by(model_variant,model,species) %>% 
  summarise(sum_lppd = sum(mean),
            mean_lppd = mean(mean),
            se_lppd = sd(mean)/sqrt(n())) %>% 
  arrange(model,species)

totals_plot <- ggplot(data = totals,
                      aes(x = species,y = mean_lppd, group = model, colour = model_variant))+
  geom_point(position = position_dodge(width = 0.4))+
  facet_wrap(vars(model))+
  coord_flip()

totals_plot

totals_wide <- totals %>% 
  select(species,model,model_variant,mean_lppd) %>% 
  mutate(mean_lppd = round(mean_lppd,3)) %>% 
  pivot_wider(names_from = model_variant,
              values_from = mean_lppd) %>% 
  arrange(species,model)

# point-wise differences spatial - hierarchical (so positive = higher lppd with spatial)
diffs <- cv_wide %>%
  mutate(diff_lppd = spatial_lppd - hier_lppd)
# simple summary of species and model differences including approximate z-score
# positive values = support for spatial model.
sp_diff_summary <- diffs %>% 
  group_by(species,model) %>% 
  summarise(mean = round(mean(diff_lppd),3),
            sd = round(sd(diff_lppd),3),
            n = n(),
            z = round(mean/(sd/sqrt(n)),3),
            .groups = "keep") %>% 
  select(species,model,z) %>% 
  left_join(.,comparisons,by = c("species","model")) %>% 
  select(-c(sp_n)) %>% 
  left_join(totals_wide,by = c("species","model")) %>% 
  arrange(stratification,species)


# Export data for Table S1 ------------------------------------------------
write.csv(sp_diff_summary,
          "Figures/cv_diff_summary.csv")



z_plot <- ggplot(data = sp_diff_summary,
                      aes(x = species,y = z))+
  geom_point(position = position_dodge(width = 0.4))+
  coord_flip(ylim = c(0,NA))+
  facet_wrap(vars(model), scales = "free")

z_plot



# Figure 7 -----------------------------------------

cv_diff<-sp_diff_summary#read.csv("second submission\\cv_diff_summary.csv")

#set order for species
cv_diff$species<-factor(
  cv_diff$species,levels=c("Scissor-tailed Flycatcher",
                           "Baird's Sparrow","Rufous Hummingbird","Carolina Wren",
                           "Bewick's Wren","Mountain Bluebird","Wood Thrush","Eastern Whip-poor-will"))
#bar graph
p<-ggplot()+coord_flip()+
  geom_bar(data=cv_diff,aes(x=species,y=z,fill=model),
           stat="identity",position=position_dodge2(width=0.9,preserve="single"))+
  theme_classic()+
  scale_fill_viridis_d(labels=c("First-difference","GAMYE"),
                       begin = 0.2,end = 0.9)+#,values=c("#440154FF","#FDE725FF"))+
  geom_vline(aes(xintercept=3.6),colour="grey40",linetype="dashed")+
  geom_hline(aes(yintercept = 0), colour = "grey40")+
  xlab("")+ylab("Z-score point-wise difference in lppd \n (spatial – non-spatial)")+
  guides(fill=guide_legend(title="Model", reverse = TRUE))+
  theme(legend.position = "right",
        legend.margin = margin(0,0,0,0),
        plot.margin = unit(c(1,4,1,1),"mm"),
        axis.title = element_text(size = 9))

p<-p+annotate(geom="text",x=c(8.4,3.4),y=c(29,29),
              label=c("Broad-grained Stratification","Fine-grained Stratification"),
              hjust = 1)


pdf("Figures/Figure_7.pdf",
    width = 7,
    height = 7)
print(p)
dev.off()


# Fit full datasets for visualisation -------------------------------------


for(j in 1:nrow(comparisons)){
  
  
  species <- comparisons[j,"species"]
  sp_n <- comparisons[j,"sp_n"]
  stratification <- comparisons[j,"stratification"]
  model <- comparisons[j,"model"]
  
  
  for(variant in model_variants){
    
    m_sel <- readRDS(paste0("data/cv_prepared_",model,"_",variant,"_full_",sp_n,".rds"))
    m_tmp <- run_model(m_sel,
                       refresh = 500,
                       output_dir = "output",
                       output_basename = paste0(sp_n,"_",model,"_",variant,"_full"))

    
    
    print(paste(model,variant))
    
  }
  
  
}  

