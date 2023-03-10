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
source("functions/map_trends.R")
# source("Functions/palettes.R")
species = "Pine Warbler"  
species_f <- gsub(species,pattern = " ",replacement = "_")



# 1 map of example strata connections ---------------------------------



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






# 6 WOTH fine grain analysis ------------------------------------------------

species <- "Wood Thrush"
model <- "first_diff"
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











  # 4 real data overall trajectories and maps for BBS --------------------------

model_variants <- c("non-hier","hier","spatial")

model_variant_names = data.frame("model_variant" = model_variants,
                                 variant_plot = factor(c("Non-hierarchical","Hierarchical",
                                                         "Spatial"),
                                                       levels = c("Non-hierarchical","Hierarchical",
                                                                  "Spatial"),
                                                       ordered = TRUE))

  species <- "Eastern Whip-poor-will"

stratification = "bbs_usgs"
base_map <- load_map(stratification)
model <- "gamye"
  

  inds_hier <- readRDS(paste0("output/",paste("indices",species,model,"hier",sep = "_")))
  inds_spatial <- readRDS(paste0("output/",paste("indices",species,model,"spatial",sep = "_")))
  inds_hier_smooth <- readRDS(paste0("output/",paste("indices_smooth",species,model,"hier",sep = "_")))
  inds_spatial_smooth <- readRDS(paste0("output/",paste("indices_smooth",species,model,"spatial",sep = "_")))
 
  
  inds_hier_plot <- inds_hier$indices %>% 
    filter(region == "continent") %>% 
    mutate(model_variant = "hier")
  inds_spatial_plot <- inds_spatial$indices %>% 
    filter(region == "continent") %>% 
    mutate(model_variant = "spatial")
  
  inds_plot <- bind_rows(inds_hier_plot,
                         inds_spatial_plot) %>% 
    left_join(.,model_variant_names,
              by = "model_variant")
  
  trajs <- ggplot(data = inds_plot,
                  aes(x = year, y = index))+
    geom_ribbon(aes(ymin = index_q_0.05, ymax = index_q_0.95,
                    fill = variant_plot),
                alpha = 0.3)+
    geom_line(aes(colour = variant_plot))+
    scale_colour_viridis_d(end = 0.8,begin = 0.2,
                           aesthetics = c("colour","fill"),
                           direction = -1)+
    guides( colour = guide_legend(title = "Model Variant"),
            fill = guide_legend(title = "Model Variant"))+
    scale_y_continuous(limits = c(0,NA))+
    ylab("Estimated annual relative abundance")+
    xlab("")+
    theme_classic()+
    theme(legend.position = "bottom",
          axis.title = element_text(size = 8))
   
  
   


  first_years <- c(1970,1995,2011)
  map_hier <- vector("list",length(first_years))
  map_spatial <- vector("list",length(first_years))
  tt_out <- NULL
  for(j in 1:length(first_years)){
    
    ys <- first_years[j]
    ye <- ys+10
    
    tt_hier <- generate_trends(inds_hier_smooth,
                               min_year = ys,
                               max_year = ye) 
    
    tt_hier <- tt_hier[["trends"]] %>% 
      filter(region_type == "stratum") %>% 
      mutate(model_variant = "hier")
    
    tt_spatial <- generate_trends(inds_spatial_smooth,
                                  min_year = ys,
                                  max_year = ye)
    tt_spatial <- tt_spatial[["trends"]] %>% 
      filter(region_type == "stratum") %>% 
      mutate(model_variant = "spatial")
    
   tt_out <- bind_rows(tt_out,
                       tt_hier)
   tt_out <- bind_rows(tt_out,
                       tt_spatial)
   

  }
  
  tt_out <- tt_out %>% 
    mutate(span = paste(start_year,end_year,sep = "-")) %>% 
    left_join(.,model_variant_names,
              by = "model_variant")
  tmap <- map_trends(tt_out,
                     facetgrid = TRUE,
                     base_map_blank = base_map,
                     facet_1 = "span",
                     facet_2 = "variant_plot",
                     legend_title = "Trend",
                     add_base = FALSE)+
    theme(legend.position = "bottom",
          legend.key.size = unit(3,"mm"),
          legend.text = element_text(size = 7))
    
  
  
  
  
  full <- trajs / tmap +
    plot_layout(design = c("
                           11
                           22
                           22
                           22
                           "))
    
  
  pdf("Figures/Figure_4.pdf",
      width = 5,
      height = 10)
  print(full)
  dev.off()
  

  
  
  
  

# 5 trajectories and trend maps for CBC and Shorebird ---------------------
  
  model_variants <- c("non-hier","hier","spatial")
  
  model_variant_names = data.frame("model_variant" = model_variants,
                                   variant_plot = factor(c("Non-hierarchical","Hierarchical",
                                                           "Spatial"),
                                                         levels = c("Non-hierarchical","Hierarchical",
                                                                    "Spatial"),
                                                         ordered = TRUE))
  
  inds_out <- readRDS("output/All_CBC_and_shorebird_indices.rds")  
  inds_cont1 <- inds_out %>% 
    filter(strata_name == "continent",
           index_type != "smooth",
           data_set == "CBC",
           model == "first_diff") %>% 
    left_join(.,model_variant_names,
              by = "model_variant")
  
  inds_cont2 <- inds_out %>% 
    filter(strata_name == "continent",
           index_type != "smooth",
           data_set == "Shorebird",
           model == "gamye") %>% 
    left_join(.,model_variant_names,
              by = "model_variant")
  
  
  
  traj1 <- ggplot(data = inds_cont1,
                  aes(x = year, y = index))+
    geom_ribbon(aes(ymin = index_q_0.05, ymax = index_q_0.95,
                    fill = variant_plot),
                alpha = 0.3)+
    geom_line(aes(colour = variant_plot))+
    scale_colour_viridis_d(end = 0.8,begin = 0.2,
                           aesthetics = c("colour","fill"),
                           direction = -1)+
    guides( colour = guide_legend(title = "Model Variant"),
            fill = guide_legend(title = "Model Variant"))+
    scale_y_continuous(limits = c(0,NA))+
    ylab("Estimated annual \n relative abundance")+
    xlab("")+
    theme_classic()+
    theme(legend.position = "bottom",
          axis.title = element_text(size = 8))
  
  traj2 <- ggplot(data = inds_cont2,
                  aes(x = year, y = index))+
    geom_ribbon(aes(ymin = index_q_0.05, ymax = index_q_0.95,
                    fill = variant_plot),
                alpha = 0.3)+
    geom_line(aes(colour = variant_plot))+
    scale_colour_viridis_d(end = 0.8,begin = 0.2,
                           aesthetics = c("colour","fill"),
                           direction = -1)+
    guides( colour = guide_legend(title = "Model Variant"),
            fill = guide_legend(title = "Model Variant"))+
    scale_y_continuous(limits = c(0,NA))+
    ylab("Estimated annual \n relative abundance")+
    xlab("")+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_text(size = 8))
  
  
  
  trends_out <- readRDS("output/All_CBC_and_shorebird_trends.rds")

  
  ### CBC map
  data_prep <- readRDS(paste0("output/dataframe_","CBC","_spatial_gamye.rds"))
  
  strat_df <- data_prep %>% 
    select(strata_name,strata_vec,non_zero,area_sq_km) %>% 
    rename(stratum = strata_vec) %>% 
    distinct() %>% 
    mutate(area_weight = area_sq_km/sum(area_sq_km))
  
 
  map_cbc<- strata_map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
    select(-area_sq_km) %>% 
    inner_join(.,strat_df,
               by = "strata_name") %>% 
    arrange(stratum)  
  
  rm(list = c("data_prep","strat_df"))

  ### Shorebird map
  load(paste0("data/Red_Knot","Shorebird","_data.RData"))
  map_shorebird <- realized_strata_map %>% 
    rename(strata_name = hex_name,
           stratum = stratn)
  
  rm(list = c("stan_data",
              "neighbours",
              "realized_strata_map",
              "data_1"))

  ### plot trends for selected periods
  
 
  trend_tmp1 <- trends_out %>% 
    filter(model == "first_diff",
           model_variant == "hier",
           data_set == "CBC",
           strata_name != "continent",
           start_year == 1975,
           end_year == 1985)
  m1 <- map_trends(trends = trend_tmp1,
                   base_map_blank = map_cbc,
                   title = "",
                   region_name = "strata_name")
  
  trend_tmp2 <- trends_out %>% 
    filter(model == "first_diff",
           model_variant == "spatial",
           data_set == "CBC",
           strata_name != "continent",
           start_year == 1975,
           end_year == 1985)
  m2 <- map_trends(trend_tmp2,
                   base_map_blank = map_cbc,
                   title = "",
                   region_name = "strata_name")
  

  
  
  trend_tmp3 <- trends_out %>% 
    filter(model == "gamye",
           model_variant == "hier",
           data_set == "Shorebird",
           strata_name != "continent",
           start_year == 1980,
           end_year == 1990)
  m3 <- map_trends(trends = trend_tmp3,
                   base_map_blank = map_shorebird,
                   title = "Hierarchical",
                   region_name = "strata_name")
  
  trend_tmp4 <- trends_out %>% 
    filter(model == "gamye",
           model_variant == "spatial",
           data_set == "Shorebird",
           strata_name != "continent",
           start_year == 1980,
           end_year == 1990)
  m4 <- map_trends(trend_tmp4,
                   base_map_blank = map_shorebird,
                   title = "Spatial",
                   region_name = "strata_name")
  
  print(m3 + m4 + plot_layout(guides = "collect") & 
          theme(legend.position = "none"))
  
  plot_shorebird <- traj2 + m3 + m4 + plot_layout(design = "
                                                  #11#
                                                  2233
                                                  ") & 
    theme(legend.position = "none")
  
  
  plot_cbc <- traj1 + m1 + m2 + plot_layout(guides = "collect",
                                            design = "
                                            #11#
                                            2233
                                            ") & 
    theme(legend.position = "bottom",
          legend.key.size = unit(3,"mm"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
  
  
  plot_both <- plot_shorebird / plot_cbc + plot_layout(nrow = 2,
                                                       ncol = 1,
                                                       heights = c(1,2)) &
    
    theme(plot.margin = unit(rep(1,4),"mm"))
  
  pdf("Figures/Figure_5.pdf",
      width = 7,
      height = 4)
print(plot_both)

dev.off()
