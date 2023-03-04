map_trends <- function(trends = trends,
                       variable = "trend",
                       plot_trend = TRUE,
                          base_map_blank = base_map,
                          title = ""){
  
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
  
  
  }else{
    

    trends_plot <- trends$trends %>% 
      filter(region_type == "stratum") %>% 
      rename_with(.,~gsub(variable,"t_plot",.x)) 
    
    pal <- scale_colour_viridis_c()
    
    
  }
  
  

  

  
  trend_map <- base_map_blank %>% 
    right_join(.,trends_plot,
               by = "strata_name",
               multiple = "all")
  
  plot_out <- ggplot()+
    geom_sf(data = trend_map,
            aes(fill = t_plot))+
    labs(title = title)+
    theme_minimal()
  
  if(plot_trend){
    plot_out <- plot_out +
    scale_fill_manual(values = pal)
  }else{
    plot_out <- plot_out +
      scale_fill_viridis_c()
  }
  
  return(plot_out)
  
}
