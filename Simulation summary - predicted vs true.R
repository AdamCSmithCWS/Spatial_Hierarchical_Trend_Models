

# Simulation summary - predicted vs true

library(bbsBayes2)
library(tidyverse)







# Explore predicted vs true trajectories and trends for simulations -----------------------

for(ma in MAs){
  mean_ab <- signif(exp(ma),2)
  load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
  realized_strata_map = strata_map
  
  strat_grid <- geofacet::grid_auto(realized_strata_map,
                                    codes = "Stratum_Factored",
                                    names = "Stratum",
                                    seed = 2019)
  
  ind_sel <- indices_all_out %>% 
    filter(!is.na(True_scaled_smooth),
           mean_abundance == mean_ab,
           version == "smooth")
  ind_sel$model <- c(rep("Spatial",nrow(ind_sel)/2),
                     rep("Non-Spatial",nrow(ind_sel)/2))
  
  mean_obs <- indices_all_out %>% 
    filter(!is.na(True_scaled_smooth),
           mean_abundance == mean_ab,
           version == "full") %>% 
    select(true_year,mean_obs,zero,
           n_surveys,Stratum,Stratum_Factored,
           Year) %>% 
    distinct()
  
  
  g_inds <- suppressMessages(ggplot(data = ind_sel,aes(x = true_year,y = median))+
                               geom_line(aes(x = true_year,
                                             y = True_scaled_smooth),
                                         colour = "black",
                                         alpha = 0.9,
                                         size = 1,
                                         inherit.aes = FALSE)+
                               geom_point(data = mean_obs,
                                          aes(x = true_year,y = mean_obs*zero,
                                              alpha = n_surveys),
                                          size = 0.3,
                                          inherit.aes = FALSE)+
                               geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),
                                           alpha = 0.3)+
                               geom_line(aes(colour = model))+
                               scale_colour_viridis_d(end = 0.85,begin = 0.2,
                                                      aesthetics = c("fill","colour"))+
                               #labs(title = paste("Simulated data mean Abundance",mean_ab))+
                               scale_y_continuous(limits = c(0,NA))+
                               geofacet::facet_geo(~Stratum,grid = strat_grid,
                                                   scales = "free_y")+
                               theme_bw()+
                               xlab("Mean annual abundance")+
                               ylab("")+
                               theme(strip.text = element_text(size = 6),
                                     strip.background = element_blank(),
                                     panel.spacing = unit(0.1,"mm"),
                                     axis.text.x = element_text(size = 5))) 
  
  #print(g_inds)
  
  pdf(paste0("Figures/Geofacet",mean_ab,"spatial_vs_non.pdf"),
      width = 8.5,
      height = 11)
  print(g_inds)
  dev.off()
  
}


# Compare trend estimates -------------------------------------------------

trend_comparison <- NULL

for(ma in MAs){
  mean_ab <- signif(exp(ma),2)
  load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
  realized_strata_map = strata_map
  
  
  t_sel <- stratum_trends %>% 
    filter(mean_abundance == mean_ab,
           data == "BBS") %>% 
    mutate(nyears = last_year - first_year)
  t_sel$model <- c(rep("Spatial",nrow(t_sel)/2),
                   rep("NonSpatial",nrow(t_sel)/2))
  
  t_sel <- t_sel %>% 
    select(Stratum_Factored,first_year,last_year,
           nyears, model,
           trend,mean_abundance) %>%
    distinct() %>% 
    pivot_wider(.,names_from = model,
                values_from = trend)
  
  true_smooth <- log_true_traj %>% 
    select(Stratum,Stratum_Factored,
           Year, True_scaled_smooth)
  tr_yrs <- t_sel %>% 
    select(first_year,last_year) %>% 
    distinct()
  
  true_trends <- NULL
  for(j in 1:nrow(tr_yrs)){
    y1 <- as.numeric(tr_yrs[j,1])
    y2 <- as.numeric(tr_yrs[j,2])
    ny = y2-y1
    
    tmpt <- true_smooth %>% 
      filter(Year %in% c(y1,y2)) %>% 
      pivot_wider(names_from = Year,
                  values_from = True_scaled_smooth,
                  names_prefix = "Y") %>% 
      rename_with(., ~gsub(replacement = "start",
                           pattern = paste0("Y",y1),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "end",
                           pattern = paste0("Y",y2),.x,
                           fixed = TRUE)) %>% 
      mutate(true_trend = texp(end/start,ny),
             first_year = y1,
             last_year = y2) %>% 
      select(Stratum,Stratum_Factored,first_year,last_year,true_trend)
    
    true_trends <- bind_rows(true_trends,tmpt)
    
  }
  
  t_sel <- left_join(t_sel,true_trends,
                     by = c("Stratum_Factored",
                            "first_year",
                            "last_year"))
  
  t_sel <- t_sel %>% 
    mutate(Spatial_error = abs(Spatial - true_trend),
           NonSpatial_error = abs(`NonSpatial` - true_trend),
           dif_error = Spatial_error - NonSpatial_error)
  
  trend_comparison <- bind_rows(trend_comparison,t_sel)
}


mean_dif <- trend_comparison %>% 
  group_by(nyears,mean_abundance) %>% 
  summarise(mean_dif = mean(dif_error),
            sd_dif = sd(dif_error),
            lq = quantile(dif_error,0.05),
            uq = quantile(dif_error,0.95),
            lci = mean_dif-(1.96*sd_dif),
            uci = mean_dif+(1.96*sd_dif))

mns <- ggplot(data = mean_dif,aes(x = nyears,y = mean_dif))+
  geom_errorbar(aes(ymin = lq,ymax = uq),
                alpha = 0.3,
                width = 0)+
  geom_point()+
  facet_wrap(vars(mean_abundance),
             nrow = 2,
             ncol = 3)

print(mns)
trend_comparison$nyearsF <- factor(trend_comparison$nyears)
box <- ggplot(data = trend_comparison,
              aes(y = dif_error,x = nyearsF))+
  geom_boxplot(alpha = 0.5)+
  facet_wrap(vars(mean_abundance),
             nrow = 2,
             ncol = 3,
             scales = "free_y")+
  ylab("Difference absolute error in trends Spatial - Non-spatial")+
  xlab("Length of trend (number of years)")+
  geom_abline(slope = 0,intercept = 0,colour = "blue")+
  theme_bw()
pdf("Figures/Figure_5.pdf",
    width = 7,
    height = 8)
print(box)
dev.off()

