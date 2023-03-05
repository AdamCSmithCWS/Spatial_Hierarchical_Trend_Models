map_trends <- function(trends = trends,
                       variable = "trend",
                       plot_trend = TRUE,
                          base_map_blank = base_map,
                          title = "",
                       legend_title = NULL,
                       zoom_out = 0.05){
  
  zoom_out <- zoom_out+1
  
  if(plot_trend){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls <- c(paste0("< ", breaks[1]),
             paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),
             paste0("> ",breaks[length(breaks)]))
  labls <- paste0(labls, " %")
  
  trends_plot <- trends$trends %>% 
    filter(region_type == "stratum") %>% 
    rename_with(.,~gsub(variable,"value",.x)) %>% 
    mutate(t_plot = cut(value, breaks = c(-Inf, breaks, Inf),
                        labels = labls,
                        ordered_result = TRUE))
  
  pal <- setNames(
    c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
      "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"),
    levels(trends_plot$t_plot))
  
  start_year <- min(trends_plot$start_year)
  end_year <- min(trends_plot$end_year)
  
  if(!is.null(legend_title)){
    l_title <- paste(legend_title,"\n",start_year,"-",end_year)
    
  }else{
  l_title <- paste("Trend \n",start_year,"-",end_year)
  }
  }else{
    

    trends_plot <- trends$trends %>% 
      filter(region_type == "stratum") %>% 
      rename_with(.,~gsub(variable,"t_plot",.x)) 
    
    pal <- scale_colour_viridis_c()
    start_year <- min(trends_plot$start_year)
    end_year <- min(trends_plot$end_year)
    
    if(!is.null(legend_title)){
      l_title <- paste(legend_title,"\n",start_year,"-",end_year)
      
    }else{
    l_title <- paste(variable,"\n",start_year,"-",end_year)
    }
    
  }
  
  
  
  

  
  trend_map <- base_map_blank %>% 
    right_join(.,trends_plot,
               by = c("strata_name" = "region"),
               multiple = "all")
  
  map_ext <- sf::st_bbox(trend_map)
  
  
  plot_out <- ggplot()+
    geom_sf(data = base_map_blank,
            alpha = 0)+
    geom_sf(data = trend_map,
            aes(fill = t_plot))+
    coord_sf(xlim = map_ext[c("xmin","xmax")]*zoom_out,
             ylim = map_ext[c("ymin","ymax")]*zoom_out)+
    labs(title = title)+
    theme_minimal()
  
  if(plot_trend){
    plot_out <- plot_out +
    scale_fill_manual(values = pal,
                      guide = guide_legend(title = l_title))
  }else{
    plot_out <- plot_out +
      scale_fill_viridis_c(guide = guide_legend(title = l_title))
  }
  
  return(plot_out)
  
}
