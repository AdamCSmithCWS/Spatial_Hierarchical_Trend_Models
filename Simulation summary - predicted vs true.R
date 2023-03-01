

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

MAs <- c(0.1,0.5,1,5,10)

models <- c("gamye","first_diff")
model_variants <- c("nonhier","hier","spatial")



# model <- models[1]
# model_variant <- model_variants[1]

# Explore predicted vs true trajectories and trends for simulations -----------------------
estimated_trends <- NULL
for(ma in MAs[c(5,1)]){
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
  mutate(gamye_dif = abs_dif_trend_gamye_hier - abs_dif_trend_gamye_spatial,
         gamye_dif_slope = abs_dif_trend_slope_gamye_hier - abs_dif_trend_slope_gamye_spatial,
         first_diff_dif_nonhier = abs_dif_trend_first_diff_nonhier - abs_dif_trend_first_diff_hier,
         first_diff_dif_nonhier_slope = abs_dif_trend_slope_first_diff_nonhier - abs_dif_trend_slope_first_diff_hier,
         first_diff_dif = abs_dif_trend_first_diff_hier - abs_dif_trend_first_diff_spatial,
         first_diff_dif_slope = abs_dif_trend_slope_first_diff_hier - abs_dif_trend_slope_first_diff_spatial,
         strat_data = ifelse(n_routes > 11,"High","Low"))  #%>% 
  #filter(start_year < 1975)
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
                         aes(x = n_years,y = gamye_dif_slope))+
  geom_boxplot(aes(group = n_years))+
  coord_cartesian(ylim = c(-1,1))+
  facet_wrap(vars(mean_true_abundance,strat_data),
             scales = "free_y")+
  geom_hline(yintercept = 0,alpha = 0.5)

trend_comp_plot2 = ggplot(data = trend_comp,
                         aes(x = n_years,y = first_diff_dif_slope))+
  geom_boxplot(aes(group = n_years))+
  coord_cartesian(ylim = c(-1,1))+
  facet_wrap(vars(mean_true_abundance,strat_data),
             scales = "free_y")+
  geom_hline(yintercept = 0,alpha = 0.5)

trend_comp_plot3 = ggplot(data = trend_comp,
                          aes(x = n_years,y = first_diff_dif_nonhier_slope))+
  geom_boxplot(aes(group = n_years))+
  facet_wrap(vars(mean_true_abundance,strat_data))+
  geom_hline(yintercept = 0,alpha = 0.5)

print(trend_comp_plot + trend_comp_plot2 + trend_comp_plot3)



stratification = "bbs_usgs"
# plot the trajectories ---------------------------------------------------
realized_strata_map <- load_map(stratify_by = stratification) %>% 
  filter(strata_name %in% true_trends$strata_name) 

strat_grid <- geofacet::grid_auto(realized_strata_map,
                                  codes = "strata_name",
                                  names = "strata_name",
                                  seed = 2019)
   
for(ma in MAs[c(5,1)]){
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

# Compare trend estimates -------------------------------------------------
# 
# trend_comparison <- NULL
# 
# for(ma in MAs){
#   mean_ab <- signif(exp(ma),2)
#   load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
#   realized_strata_map = strata_map
#   
#   
#   t_sel <- stratum_trends %>% 
#     filter(mean_abundance == mean_ab,
#            data == "BBS") %>% 
#     mutate(nyears = last_year - first_year)
#   t_sel$model <- c(rep("Spatial",nrow(t_sel)/2),
#                    rep("NonSpatial",nrow(t_sel)/2))
#   
#   t_sel <- t_sel %>% 
#     select(Stratum_Factored,first_year,last_year,
#            nyears, model,
#            trend,mean_abundance) %>%
#     distinct() %>% 
#     pivot_wider(.,names_from = model,
#                 values_from = trend)
#   
#   true_smooth <- log_true_traj %>% 
#     select(Stratum,Stratum_Factored,
#            Year, True_scaled_smooth)
#   tr_yrs <- t_sel %>% 
#     select(first_year,last_year) %>% 
#     distinct()
#   
#   true_trends <- NULL
#   for(j in 1:nrow(tr_yrs)){
#     y1 <- as.numeric(tr_yrs[j,1])
#     y2 <- as.numeric(tr_yrs[j,2])
#     ny = y2-y1
#     
#     tmpt <- true_smooth %>% 
#       filter(Year %in% c(y1,y2)) %>% 
#       pivot_wider(names_from = Year,
#                   values_from = True_scaled_smooth,
#                   names_prefix = "Y") %>% 
#       rename_with(., ~gsub(replacement = "start",
#                            pattern = paste0("Y",y1),.x,
#                            fixed = TRUE))%>% 
#       rename_with(., ~gsub(replacement = "end",
#                            pattern = paste0("Y",y2),.x,
#                            fixed = TRUE)) %>% 
#       mutate(true_trend = texp(end/start,ny),
#              first_year = y1,
#              last_year = y2) %>% 
#       select(Stratum,Stratum_Factored,first_year,last_year,true_trend)
#     
#     true_trends <- bind_rows(true_trends,tmpt)
#     
#   }
#   
#   t_sel <- left_join(t_sel,true_trends,
#                      by = c("Stratum_Factored",
#                             "first_year",
#                             "last_year"))
#   
#   t_sel <- t_sel %>% 
#     mutate(Spatial_error = abs(Spatial - true_trend),
#            NonSpatial_error = abs(`NonSpatial` - true_trend),
#            dif_error = Spatial_error - NonSpatial_error)
#   
#   trend_comparison <- bind_rows(trend_comparison,t_sel)
# }
# 
# 
# mean_dif <- trend_comparison %>% 
#   group_by(nyears,mean_abundance) %>% 
#   summarise(mean_dif = mean(dif_error),
#             sd_dif = sd(dif_error),
#             lq = quantile(dif_error,0.05),
#             uq = quantile(dif_error,0.95),
#             lci = mean_dif-(1.96*sd_dif),
#             uci = mean_dif+(1.96*sd_dif))
# 
# mns <- ggplot(data = mean_dif,aes(x = nyears,y = mean_dif))+
#   geom_errorbar(aes(ymin = lq,ymax = uq),
#                 alpha = 0.3,
#                 width = 0)+
#   geom_point()+
#   facet_wrap(vars(mean_abundance),
#              nrow = 2,
#              ncol = 3)
# 
# print(mns)
# trend_comparison$nyearsF <- factor(trend_comparison$nyears)
# box <- ggplot(data = trend_comparison,
#               aes(y = dif_error,x = nyearsF))+
#   geom_boxplot(alpha = 0.5)+
#   facet_wrap(vars(mean_abundance),
#              nrow = 2,
#              ncol = 3,
#              scales = "free_y")+
#   ylab("Difference absolute error in trends Spatial - Non-spatial")+
#   xlab("Length of trend (number of years)")+
#   geom_abline(slope = 0,intercept = 0,colour = "blue")+
#   theme_bw()
# pdf("Figures/Figure_5.pdf",
#     width = 7,
#     height = 8)
# print(box)
# dev.off()
# 
