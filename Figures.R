#### Figures for publication
library(tidyverse)
library(bbsBayes2)
library(cmdstanr)
library(posterior)
library(sf)
library(patchwork)
library(geofacet)
library(ggrepel)
source("functions/neighbours_define.R")
 source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")



# 1 map of example strata connections ---------------------------------


species = "Pacific Wren"  
species_f <- gsub(species,pattern = " ",replacement = "_") #name without spaces

stratification = "bbs_usgs"

s <- stratify(by = stratification,
              species = species)

prep_data <- prepare_data(s)

all_strata_map <- load_map(stratification)

real_strata <- all_strata_map %>% 
  filter(strata_name %in% prep_data$meta_strata$strata_name)

real_strata_neighbours <- neighbours_define(real_strata,
                                     plot_dir = "maps/",
                                     plot_file = "simulated",
                                     save_plot_data = TRUE,
                                     strat_indicator = "strata_name",
                                     nn_fill = FALSE)

load("maps/simulated_data.RData")

xb = range(st_coordinates(real_strata_map)[,"X"])
yb = range(st_coordinates(real_strata_map)[,"Y"])

set.seed(1)
real_strata_map <- real_strata_map %>% 
  mutate(rand_strat = sample(1:nrow(real_strata_map)))

ggp <- ggplot()+ 
  geom_sf(data = real_strata_map,
          aes(fill = rand_strat),
          alpha = 1,
          colour = grey(0.8))+
  geom_segment(data=DA,
               aes(x = long, y = lat,
                   xend=long_to,yend=lat_to),
               inherit.aes = FALSE,
               linewidth=0.5,alpha=0.2) +
  geom_sf(data = centres, alpha = 0.9,colour = "white",
          size = 0.5) + 
  xlab("")+
  ylab("")+
  scale_fill_viridis_c()+
  theme_minimal() +
    coord_sf(xlim = xb,ylim = yb)+
    theme(legend.position = "none")

#print(ggp)


pdf(file = paste0("Figures/Figure_1.pdf"),
    width = 3.5,
    height = 4)
print(ggp)
dev.off()






# WOTH fine grain analysis ------------------------------------------------

species <- "Wood Thrush"
model <- "gam"
model_variant <- "Spatial"

fit <- readRDS(paste0("output/",paste(species,model,model_variant,sep = "_"),".rds"))

inds <- generate_indices(fit)

trends <- generate_trends(inds,
                          min_year = 2007)

map <- plot_map(trends,
                title = FALSE)
map_ext <- load_map("latlong") %>% 
  filter(strata_name %in% unique(inds$raw_data$strata_name)) %>% 
  sf::st_buffer(.,100000) %>% 
  sf::st_bbox()

map <- map + 
  coord_sf(xlim = map_ext[c("xmin","xmax")]+c(0,200000),
           ylim = map_ext[c("ymin","ymax")]+c(-200000,0))

pdf("figures/Figure_6.pdf",
    width = 4.5,
    height = 4)
map
dev.off()











# 2 simulated data trajectory geofacet ---------------------------------------------------


output_dir <- "output/"
#tp = "non_linear"
load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))
load(paste0("Data/",species_f,"BBS","_data.RData"))


  strat_grid <- geofacet::grid_auto(realized_strata_map,
                                    codes = "Stratum_Factored",
                                    seed = 2)

  
  fig4_geo = ggplot(data = nsmooth_comp2,aes(y = True_nsmooth,
                                         x = Year))+
    geom_ribbon(aes(ymin = lci,ymax = uci,
                    fill = version),alpha = 0.3)+
    geom_point(data = nsmooth_comp2,aes(x = Year,y = mean_count),
               alpha = 0.1,
               size = 0.2,
               inherit.aes = FALSE)+
    geom_line(aes(colour = version))+
    scale_colour_viridis_d(aesthetics = c("colour","fill"),
                           begin = 0,
                           end = 0.5,
                           direction = -1)+
    scale_y_continuous(limits = c(0,NA))+
    geofacet::facet_geo(~Stratum_Factored,grid = strat_grid,
                        scales = "free")+
    xlab("")+
    ylab("Mean annual smooth trajectory")+
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(size = 5))
  
  pdf(file = paste0("Figures/Figure_3.pdf"),
      width = 7,
      height = 10)
  print(fig4_geo)
  dev.off()


# 5 Trend comparisons ---------------------------------------------------
    MAs <- round(log(c(0.1,0.5,1,50)),2)# true mean abundances for different simulations
    
output_dir <- "output/"
tp = "non_linear"
load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))

strat_df <- as.data.frame(strata_mask) %>% 
  select(Stratum_Factored,masked)
sw_trends <- NULL
strat_trends <- NULL

for(sns in c("","nonSpatial_alt_")){#,"nonSpatial_"))
  for(ma in MAs){
    
    out_base_sim <- paste0("sim_",sns,tp,"_",ma,"_BBS")
    
    lbl <- "Spatial"
    if(sns == "nonSpatial_alt_"){
      lbl <- "NonSpatial"
    }
    lblm <- paste0(signif(exp(ma),1))
      
   
    
    load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
    sw_t <- all_trends %>% 
      filter(Region_type == "Survey_Wide_Mean",
             #last_year-first_year == 2019,
             last_year-first_year == 5) %>% 
      mutate(version = lbl,
             true_mean = lblm)
    sw_trends <- bind_rows(sw_trends,sw_t)
    
    strat_t <- all_trends %>% 
      filter(Region_type == "Stratum_Factored",
             #last_year-first_year == 2019,
             last_year-first_year == 5) %>% 
      mutate(version = lbl,
             true_mean = lblm)
    strat_trends <- bind_rows(strat_trends,strat_t)
  }
}



strat_trends_nm <- strat_trends %>% 
  left_join(.,strat_df,by = "Stratum_Factored") %>%
  #filter(first_year %in% c(1970:2009)) %>% 
  mutate(t_dif = true_trend - trend,
         t_abs_dif = abs(t_dif),
         t_se = ((uci-lci)/(1.96*2)),
         t_prec = 1/t_se^2)%>% 
  mutate(trend_time = factor(paste0(first_year,"-",last_year)))



mean_difs <- strat_trends_nm %>% 
  group_by(true_mean,
           version,
           trend_time) %>% 
  summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
            lci = quantile(t_abs_dif,0.05,na.rm = T),
            uci = quantile(t_abs_dif,0.95,na.rm = T))

sw_trends_ab <- sw_trends %>% 
  mutate(t_dif = true_trend - trend,
         t_abs_dif = abs(t_dif),
         t_se = ((uci-lci)/(1.96*2)),
         t_prec = 1/t_se^2)%>% 
  mutate(trend_time = factor(paste0(first_year,"-",last_year)))

mean_difs_sw <- sw_trends_ab %>% 
  group_by(true_mean,
           version,
           trend_time) %>% 
  summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
            lci = quantile(t_abs_dif,0.05,na.rm = T),
            uci = quantile(t_abs_dif,0.95,na.rm = T))

trends4means_plot <- ggplot(data = mean_difs,
                              aes(x = trend_time,
                                  y = mean_abs_dif,
                                  colour = version),
                              position = position_dodge(width = 0.5))+
  geom_errorbar(aes(colour = version,ymin = lci,
                    ymax = uci),
                alpha = 0.3,
                width = 0,
                position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5))+
  scale_y_continuous(limits = c(0,NA))+
  ylab("Absolute difference in trend (Estimated - True)")+
  xlab("Timespan of Trend")+
  theme_bw()+
  facet_wrap(vars(true_mean),
             nrow = 2,
             scales = "fixed")
print(trends4means_plot)

m1 = brm(t_abs_dif  ~ version*true_mean+trend_time + (1|Stratum_Factored),
        data = strat_trends_nm)
summary(m1)

m2 = brm(t_abs_dif  ~ version*true_mean+trend_time,
         data = sw_trends_ab)
summary(m2)


# trends4a_plot <- ggplot(data = strat_trends_nm,
#                        aes(x = first_year,y = t_abs_dif,fill = version))+
#   geom_boxplot(aes(fill = version),
#              position = position_dodge(width = 1))+
#   scale_y_continuous(limits = c(0,NA))+
#   ylab("Absolute difference in trend (Estimated - True)")+
#   xlab("Timespan of Trend")+
#   theme_bw()+
#   facet_wrap(vars(true_mean),
#              nrow = 2,
#              scales = "free")
# 
#   print(trends4a_plot)
# 
# 
#   
  # strat_trends_m <- strat_trends %>%
  #   filter(first_year %in% c(1970:2009),
  #          masked == TRUE,
  #          version %in% c("NonSpatial Masked","Spatial Masked")) %>%
  #   mutate(first_year = factor(paste0(first_year,"-2019")))

#   m2 = lm(t_abs_dif~version*true_mean+first_year,
#           data = strat_trends_m)
#   summary(m2)
#   
#   
#   trends4b_plot <- ggplot(data = strat_trends_m,
#                           aes(x = first_year,y = t_abs_dif,fill = version))+
#     geom_boxplot(aes(fill = version),
#                  position = position_dodge(width = 1),
#                  coef = 3)+
#     scale_y_continuous(limits = c(0,NA))+
#     ylab("Absolute difference in trend (Estimated - True)")+
#     xlab("Timespan of Trend")+
#     theme_bw()+
#     facet_wrap(vars(true_mean),
#                nrow = 2,
#                scales = "free")
#   
#   print(trends4b_plot)
 
  mean_difs_m <- strat_trends_m %>% 
    group_by(true_mean,
             version,
             first_year) %>% 
    summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
              lci = quantile(t_abs_dif,0.05,na.rm = T),
              uci = quantile(t_abs_dif,0.95,na.rm = T))
  
  trends4means_plot_m <- ggplot(data = mean_difs_m,
                          aes(x = first_year,
                              y = mean_abs_dif,
                              colour = version),
                          position = position_dodge(width = 0.5))+
    geom_errorbar(aes(colour = version,ymin = lci,
                      ymax = uci),
                  alpha = 0.3,
                  width = 0,
                  position = position_dodge(width = 0.5))+
    geom_point(position = position_dodge(width = 0.5))+
    scale_y_continuous(limits = c(0,NA))+
    ylab("Absolute difference in trend (Estimated - True)")+
    xlab("Timespan of Trend")+
    theme_bw()+
    facet_wrap(vars(true_mean),
               nrow = 2,
               scales = "free")
  
  pdf(file = paste0("Figures/Figure_5.pdf"),
      width = 7,
      height = 4)
  print(trends4means_plot_m)
  dev.off()
  
  
 
  # 6 real data overall trajectories for 3 species --------------------------

  
  load("output/real_data_summaries.RData")
  
  
  fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                  "Cinclus_mexicanus",
                                  "Red_Knot"),
                    species = c("Yellow-headed Blackbird",
                                "Cinclus_mexicanus",
                                "Red Knot"),
                    species_l = c("BBS - Yellow-headed Blackbird",
                                  "CBC - American Dipper",
                                  "Shorebird - Red Knot"))
  
  Indices_all_out <- Indices_all_out %>% 
    left_join(.,fls,by = "species")
Iplot <- ggplot(data = Indices_all_out,
                aes(x = true_year,y = median))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),
              alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"),
                         begin = 0,
                         end = 0.6,
                         direction = 1)+
  ylab("Survey-wide mean annual predictions")+
  xlab("")+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(vars(species_l),
             nrow = 3,
             ncol = 1,
             scales = "free_y")

pdf(file = "Figures/Figure_6.pdf",
    width = 3.5,
    height = 5)
print(Iplot)
dev.off()



# 7 long-term and short-term (3-gen) trend maps ---------------------------
# six panel, paired maps

load("output/real_data_summaries.RData")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                "Cinclus_mexicanus",
                                "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "American Dipper",
                              "Red Knot"),
                  data = c("BBS",
                           "CBC",
                           "Shorebird"),
                  out_base = c(paste0("Yellow-headed_Blackbird","_real_","BBS"),
                               paste0("Cinclus_mexicanus","_CBC_B"),
                               paste0("Red Knot","_Shorebird")),
                  y1 = c(1966,
                         1966,
                         1980),
                  strat_map_name = c("Stratum_Factored",
                                     "strata_vec",
                                     "stratn"))



t1 <- tt_map_list[[1]][["TY1966-2019"]]  +
  labs(subtitle = "BBS",
       title = "Long-term")
t2a <- tt_map_list[[1]][["TY1970-1980"]] +
  labs(title = "First ten years")
t2b <- tt_map_list[[1]][["TY2009-2019"]] +
  labs(title = "Last ten years")
t3 <- tt_map_list[[2]][["TY1966-2019"]] +
  labs(subtitle = "CBC")
t4a <- tt_map_list[[2]][["TY1970-1980"]]
t4b <- tt_map_list[[2]][["TY2009-2019"]]
t5 <- tt_map_list[[3]][["TY1980-2019"]] +
  labs(subtitle = "Shorebird")
t6a <- tt_map_list[[3]][["TY1980-1990"]]
t6b <- tt_map_list[[3]][["TY2009-2019"]]

tcomb = t1 + t2a +t2b + t3 + t4a +t4b + t5 + t6a +t6b +
  plot_layout(ncol = 3,byrow = TRUE,
              guides = "collect")

pdf(file = paste0("Figures/Figure_7.pdf"),
    width = 7,
    height = 8)


print(tcomb)
dev.off()




# 8 Spaghetti plot Red Knot -------------------------------------------------
load("output/real_data_summaries.RData")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                "Cinclus_mexicanus",
                                "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "Cinclus_mexicanus",
                              "Red Knot"),
                  species_l = c("BBS - Yellow-headed Blackbird",
                                "CBC - American Dipper",
                                "Shorebird - Red Knot"))

Indices_all_out <- Indices_all_out %>% 
  left_join(.,fls,by = "species")

I_RK <- Indices_all_out %>% 
  filter(species == "Red Knot")
Iplot <- ggplot(data = I_RK,
                aes(x = true_year,y = median))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),
              alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"),
                         begin = 0,
                         end = 0.6,
                         direction = 1)+
  ylab("Survey-wide mean annual predictions")+
  xlab("")+
  theme_classic()+
  theme(legend.position = "none")


load("output/Red_Knot_SW_indices.RData")
set.seed(4)
sel <- sample(1:max(sw_smooth$samples$.draw),200)
rk_inds <- sw_smooth$samples %>% 
  filter(.draw %in% sel)

start_val = rk_inds %>% 
  filter(true_year == 1980)%>% 
  mutate(start_value = log10(.value)) %>% 
  select(.draw,start_value) 
rk_inds <- rk_inds %>% 
  left_join(.,start_val,by = ".draw")

spagl <- ggplot(data = rk_inds,
                aes(x = true_year,
                    y = .value,
                    group = .draw,
                    colour = start_value))+
  geom_line(alpha = 0.3)+
  xlab("")+
  ylab("log(Survey-wide smooth)")+
  theme_classic()+
  scale_colour_viridis_c(end = 0.9)+
  scale_y_continuous(trans = "log10")+
  theme(legend.position = "none")

pdf(file = paste0("Figures/Figure_8.pdf"),
    width = 7,
    height = 4)

print(Iplot + spagl) 

dev.off()




