

# Simulation summary - predicted vs true

library(bbsBayes2)
library(tidyverse)
library(patchwork)

# load true parameter estimates -------------------------------------------

strata_base_trajs <- readRDS("data/simulated_data_true_trajectories.rds")


true_trajectories <- strata_base_trajs %>% 
  select(strata_name,
         year,
         expected_count)


#trend function from ratio of indices to %/year
texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}


#trend function from slope of log-linear regression to %/year
slope_trend <- function(x,y){
  n = length(x)
  x = log(x)
  sx = sum(x)
  sy = sum(y)
  ssy = sum(y^2)
  sxx = sum(x*y)
  b = (n*sxx - sy*sx)/(n*ssy - sy^2)
  return(((exp(b)-1)*100))
}



true_trends <- NULL
tlengths <- c(10,50)
for(ny in tlengths){
  y2s <- c(max(true_trajectories$year):(min(true_trajectories$year)+ny))
for(y2 in y2s){
  y1 <- y2 - ny

  
  tmpt <- true_trajectories %>% 
    filter(year %in% c(y1,y2)) %>% 
    pivot_wider(names_from = year,
                values_from = expected_count,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",y1),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",y2),.x,
                         fixed = TRUE)) %>% 
    mutate(true_trend = texp(end/start,ny),
           start_year = y1,
           end_year = y2,
           nyears = ny) %>% 
    select(strata_name,start_year,end_year,true_trend,nyears)
  
  
  tmpt2 <- true_trajectories %>% 
    filter(year %in% c(y1:y2)) %>% 
    group_by(strata_name) %>% 
    summarise(true_trend_slope = slope_trend(expected_count,year)) %>% 
    mutate(start_year = y1,
           end_year = y2,
           nyears = ny) %>% 
    select(strata_name,start_year,end_year,true_trend_slope,nyears)
  
  tmpt <- left_join(tmpt,tmpt2,
                    by = c("strata_name","start_year","end_year","nyears"))
  true_trends <- bind_rows(true_trends,tmpt)
  
}
}

tmpp <- ggplot(data = true_trends,
               aes(x = true_trend,y = true_trend_slope,colour = start_year))+
  geom_point(alpha = 0.3)+
  scale_colour_viridis_c(end = 0.8)+
  facet_wrap(vars(nyears))
  
print(tmpp)

species = "simulated"

MAs <- c(0.1,1,10)

models <- c("gamye","first_diff")
model_variants <- c("nonhier","hier","spatial")



# model <- models[1]
# model_variant <- model_variants[1]

# Explore predicted vs true trajectories and trends for simulations -----------------------
estimated_trends <- NULL
for(ma in MAs){
  for(model in models){
  
for(model_variant in model_variants){
  if(model == "gamye" & model_variant == "nonhier"){next}
  
 
  log_ma <- round(log(ma),2)
  ma_f <- gsub(as.character(ma),pattern = ".",replacement = "-",
               fixed = TRUE)
  
  

# load fitted estimates ---------------------------------------------------
  fit <- readRDS(paste0("output/",paste(species,model,model_variant,ma_f,sep = "_"),".rds"))
  

  

# estimate trends and trajectories ----------------------------------------

inds <- generate_indices(fit,
                         regions = "stratum")
 
  saveRDS(inds,file = paste0("output/","sim_indices_",ma_f,"_",model,"_",model_variant,".rds")) 

  if(model == "gamye"){
  inds <- generate_indices(fit,
                           regions = "stratum",
                           alternate_n = "n_smooth")
  
  saveRDS(inds,file = paste0("output/","sim_indices_smooth_",ma_f,"_",model,"_",model_variant,".rds")) 
  }
  
  for(ny in tlengths){
    y2s <- c(max(inds$indices$year):(min(inds$indices$year)+ny))
    for(y2 in y2s){
      y1 <- y2 - ny
      
      trend <- generate_trends(inds,
                               min_year = y1,
                               max_year = y2,
                               quantiles = c(0.05,0.95)) 
      
      ttmp <- trend$trends %>% 
        mutate(model = model,
               model_variant = model_variant,
               mean_true_abundance = ma)
      
  estimated_trends <- bind_rows(estimated_trends,
                                ttmp)
    
  
    }
    
  }
    
  
  
} # end model_variant

} # end ma
  
}# end model
saveRDS(estimated_trends,
        file = "output/estimated_trend_simulated_data.rds")


# grab the fitted model for the last model and model variant --------
# create table of the non_zero_weights to rescale the true trajectories

non_zero_w <- readRDS(paste0("output/",paste(species,model,model_variant,ma_f,sep = "_"),".rds"))
non_zero_w <- non_zero_w$raw_data %>% 
  select(strata_name,non_zero_weight) %>% 
  distinct()



# identify the peripheral strata -------------------------------------------
perif_strata <- strata_base_trajs %>% 
  select(y_scale,strata_name) %>% 
  distinct() %>% 
  filter(abs(y_scale) > 1) %>% 
  select(strata_name) %>% 
  unlist()

# merge estimated and true ------------------------------------------------

trend_comp <- estimated_trends %>% 
    left_join(.,true_trends,
              by = c("start_year","end_year",
                     "region" = "strata_name")) %>% 
  mutate(n_years = end_year-start_year,
         abs_dif_trend = abs(trend - true_trend),
         abs_dif_trend_slope = abs(trend - true_trend_slope)) %>% 
  pivot_wider(names_from = c(model, model_variant),
              values_from = c(abs_dif_trend,abs_dif_trend_slope),
              names_sep = "_",
              id_cols = c(start_year,end_year,
                         region,n_routes,mean_n_routes,
                         n_years,true_trend,mean_true_abundance)) %>% 
  filter(mean_true_abundance %in% c(0.1,1,10)) %>% 
  mutate(gamye_dif = abs_dif_trend_gamye_hier - abs_dif_trend_gamye_spatial,
         gamye_dif_slope = abs_dif_trend_slope_gamye_hier - abs_dif_trend_slope_gamye_spatial,
         first_diff_dif_nonhier = abs_dif_trend_first_diff_nonhier - abs_dif_trend_first_diff_hier,
         first_diff_dif_nonhier_slope = abs_dif_trend_slope_first_diff_nonhier - abs_dif_trend_slope_first_diff_hier,
         first_diff_dif = abs_dif_trend_first_diff_hier - abs_dif_trend_first_diff_spatial,
         first_diff_dif_slope = abs_dif_trend_slope_first_diff_hier - abs_dif_trend_slope_first_diff_spatial,
         strat_data = ifelse(n_routes > 11,"High","Low"),
         peripheral = ifelse(region %in% perif_strata,"Peripheral Strata","Core Strata"),
         term = ifelse(n_years == 10,"Short","Long"),
         true_abundance = factor(mean_true_abundance, ordered = TRUE))  
### add a north-south selection (edge and core distinction - similar and different from mean)



# 
# 
# trend_comp_plot = ggplot(data = trend_comp,
#                          aes(x = n_years,y = gamye_dif_slope))+
#   geom_point(aes(colour = mean_n_routes),alpha = 0.3,
#              position = position_jitter(width = 10))+
#   scale_colour_viridis_c(direction = -1)+
#   guides(colour = guide_legend(order = -1))+
#   facet_wrap(vars(mean_true_abundance))+
#   geom_hline(yintercept = 0,alpha = 0.5)
#   
# 
# trend_comp_plot2 = ggplot(data = trend_comp,
#                          aes(x = n_years,y = first_diff_dif_slope))+
#   geom_point(aes(colour = mean_n_routes),alpha = 0.3,
#              position = position_jitter(width = 10))+
#   scale_colour_viridis_c(direction = -1)+
#   guides(colour = guide_legend(order = -1))+
#   facet_wrap(vars(mean_true_abundance))+
#   geom_hline(yintercept = 0,alpha = 0.5)
# 
# 
# print(trend_comp_plot + trend_comp_plot2)
# 


trend_comp_plot = ggplot(data = trend_comp,
                         aes(x = term,y = gamye_dif,
                             colour = true_abundance))+
  geom_boxplot()+
  coord_cartesian(ylim = c(-2,2))+
  facet_wrap(vars(peripheral))+
  labs(title = "GAMYE")+
  guides( colour = guide_legend(title = "Simulated mean count"))+
  ylab("Difference in trend error \n (Hierarchical-Spatial)")+
  geom_hline(yintercept = 0,alpha = 0.5)



trend_comp_plot2 = ggplot(data = trend_comp,
                         aes(x = term,y = first_diff_dif,
                             colour = true_abundance))+
  geom_boxplot()+
  coord_cartesian(ylim = c(-2,2))+
  facet_wrap(vars(peripheral))+
  guides( colour = guide_legend(title = "Simulated mean count"))+
  labs(title = "First-Difference")+
  ylab("(Hierarchical - Spatial)")+
  geom_hline(yintercept = 0,alpha = 0.5)


trend_comp_plot3 = ggplot(data = trend_comp,
                         aes(x = term,y = first_diff_dif_nonhier,
                             colour = true_abundance))+
  geom_boxplot()+
  coord_cartesian(ylim = c(-2,2))+
  facet_wrap(vars(peripheral))+
  labs(title = "First-Difference")+
  guides( colour = guide_legend(title = "Simulated mean count"))+
  ylab("(Non-Hierarchical - Hierarchical)")+
  geom_hline(yintercept = 0,alpha = 0.5)



pdf("Figures/Figure_3.pdf",
    width = 7.5,
    height = 3.5)
print(trend_comp_plot + trend_comp_plot2 + trend_comp_plot3 +
        plot_layout(guides = "collect") &
        xlab("Trend Term") &
        guides( colour = guide_legend(title = "Simulated mean count")) &
        theme(text = element_text(size = 9),
              legend.position = "bottom"))
dev.off()


# plot the trajectories ---------------------------------------------------



inds_out <- NULL
stratification = "bbs_usgs"
realized_strata_map <- load_map(stratify_by = stratification) %>% 
  filter(strata_name %in% true_trends$strata_name) 

strat_grid <- geofacet::grid_auto(realized_strata_map,
                                  codes = "strata_name",
                                  names = "strata_name",
                                  seed = 2019)
   
for(ma in MAs){
  ma_f <- gsub(as.character(ma),pattern = ".",replacement = "-",
               fixed = TRUE)
  log_ma <- round(log(ma),6)
  
  true_trajectories <- readRDS(paste0("data/simulated_data_true_trajectories_",ma_f,".rds")) %>% 
    left_join(.,non_zero_w,
              by = "strata_name") %>% 
    mutate(expected_mean_count = expected_mean_count * non_zero_weight)
  
  for(model in models){
    ind_vars <- NULL
    
  for(model_variant in model_variants){
    
    if(model == "gamye" & model_variant == "nonhier"){next}
    
    
    inds <- readRDS(paste0("output/","sim_indices_",ma_f,"_",model,"_",model_variant,".rds"))
    
    inds <- inds$indices %>% 
      mutate(model = model,
             model_variant = model_variant,
             mean_true_abundance = ma,
             obs_mean = ifelse(n_routes == 0,NA,obs_mean))
    
    ind_vars <- bind_rows(ind_vars,
                          inds)
  }
 ind_vars <- ind_vars %>% 
   left_join(., true_trajectories,
             by = c("region" = "strata_name",
                    "year" = "year")) %>% 
   mutate(strata_name = region)
  
inds_out <- bind_rows(inds_out,ind_vars)
  
  g_inds <- suppressMessages(ggplot(data = ind_vars,aes(x = year,y = index))+
                               geom_line(aes(x = year,
                                             y = expected_mean_count),
                                         colour = "black",
                                         alpha = 0.9,
                                         linewidth = 1,
                                         inherit.aes = FALSE)+
                               geom_point(aes(x = year,y = obs_mean),alpha = 0.3,inherit.aes = FALSE)+
                               geom_ribbon(aes(ymin = index_q_0.05,
                                               ymax = index_q_0.95,
                                               fill = model_variant),
                                           alpha = 0.3)+
                               geom_line(aes(colour = model_variant))+
                               scale_colour_viridis_d(end = 0.5,begin = 0.2,
                                                      aesthetics = c("fill","colour"))+
                               #labs(title = paste("Simulated data mean Abundance",mean_ab))+
                               scale_y_continuous(limits = c(0,NA))+
                               # geofacet::facet_geo(~strata_name,grid = strat_grid,
                               #                     scales = "free_y")+
                               facet_wrap(vars(strata_name),
                                          scales = "free_y")+
                               theme_bw()+
                               xlab("Mean annual abundance")+
                               ylab("")+
                               theme(strip.text = element_text(size = 6),
                                     strip.background = element_blank(),
                                     panel.spacing = unit(0.1,"mm"),
                                     axis.text.x = element_text(size = 5))) 
  
  print(g_inds)
  
  pdf(paste0("Figures/Geofacet",model,"_",ma_f,"spatial_vs_non.pdf"),
      width = 8.5,
      height = 11)
  print(g_inds)
  dev.off()
  
}
}


saveRDS(inds_out,"output/all_saved_true_estimated_indices.rds")




# Simulation_trajectory ---------------------------------------------------


#### add the smooth componennt to the gam side of the graph (just the lines in the same 
####  model colours)
ind_smooth <- NULL
for(model_variant in model_variants[2:3]){
  inds <- readRDS(paste0("output/","sim_indices_smooth_",1,"_","gamye","_",model_variant,".rds"))
  
  inds <- inds$indices %>% 
    mutate(model = "gamye",
           model_variant = model_variant)
  
  ind_smooth <- bind_rows(ind_smooth,
                        inds)
  
}


model_variant_names = data.frame("model_variant" = model_variants,
                         variant_plot = factor(c("Non-hierarchical","Hierarchical",
                                                 "Spatial"),
                                               levels = c("Non-hierarchical","Hierarchical",
                                                          "Spatial"),
                                               ordered = TRUE))

# ind_smooth <- ind_smooth %>% 
#   filter(region %in% c("US-AK-5","CA-AB-10","US-CA-15")) %>% 
#   mutate(strata_name = factor(region,levels = c("US-AK-5","CA-AB-10","US-CA-15"),
#                               ordered = TRUE),
#          model_plot = ifelse(model == "first_diff","First Difference","GAMYE")) %>% 
#   left_join(.,model_variant_names,
#             by = "model_variant")

st_sel <- c("CA-ON-8","CA-MB-12","US-SC-27","US-VA-28","US-LA-26","US-LA-27")
inds_plot <- inds_out %>% 
  filter(strata_name %in% st_sel,
         mean_true_abundance == 1) %>% 
  mutate(strata_name = factor(strata_name,levels = st_sel,
                              ordered = TRUE),
         model_plot = ifelse(model == "first_diff","First Difference","GAMYE")) %>% 
  left_join(.,model_variant_names,
            by = "model_variant")
  


  
  inds_demo <- ggplot(data = inds_plot,
                      aes(x = year,y = index))+
  geom_line(aes(x = year,
                y = expected_mean_count),
            colour = "black",
            alpha = 0.9,
            linewidth = 1)+
    geom_line(aes(colour = variant_plot),
              linewidth = 1)+
    #geom_line(data = ind_smooth,
    #          aes(x = year, y = index,colour = variant_plot),
    #          linewidth = 0.8,
    #          linetype = 2)+
    geom_point(aes(x = year,y = obs_mean),alpha = 0.1,inherit.aes = FALSE)+
    geom_ribbon(aes(ymin = index_q_0.05,
                    ymax = index_q_0.95,
                    fill = variant_plot),
                alpha = 0.2)+
    scale_colour_viridis_d(end = 0.8,begin = 0.2,
                           aesthetics = c("colour","fill"),
                           direction = -1)+
    guides( colour = guide_legend(title = "Model Variant"),
            fill = guide_legend(title = "Model Variant"))+
    scale_y_continuous(limits = c(0,NA))+
    facet_grid(scales = "free_y",
               cols = vars(model_plot),
               rows = vars(strata_name))+
    theme_bw()+
    ylab("Mean annual abundance")+
    xlab("")+
    theme(strip.text = element_text(),
          strip.background = element_blank(),
          panel.spacing = unit(1.5,"mm"),
          axis.text.x = element_text(size = 5),
          legend.position = "bottom",
          legend.text = element_text(size = 7))

  print(inds_demo)    
  
  pdf("figures/Figure_2.pdf",
      width = 4,
      height = 7)
  print(inds_demo)
  dev.off()
  
  
