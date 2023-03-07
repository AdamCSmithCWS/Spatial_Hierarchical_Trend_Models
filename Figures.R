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











model_variants <- c("non-hier","hier","spatial")
  # 4 real data overall trajectories and maps for EAWP --------------------------
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





