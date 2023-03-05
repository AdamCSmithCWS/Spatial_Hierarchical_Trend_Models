
# Trends function ---------------------------------------------------------
#tr_tmp <- trends_function(ind_list = ttinds)
trends_function <- function(samples = indices,
                            start_year = NULL,
                            end_year = NULL,
                            quant = 0.95,
                            parameter){
  
  samples = indices
  parameter = ind_list$parameter
  strat = ind_list$strat
  year = ind_list$year
  dims = ind_list$dims
  weights_df = ind_list$weights_df
  area = ind_list$area
  summary_regions = ind_list$summary_regions
  to_summarise = ind_list$to_summarise
  
  if(is.null(end_year)){
    end_year <- max(samples$true_year)
  }
  if(is.null(start_year)){
    start_year <- min(samples$true_year)
  }
  
  nyrs <- end_year-start_year
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)
  
  if(!is.null(weights_df) & to_summarise){
    
    
    indt <- samples %>% 
      filter(true_year %in% c(start_year,end_year)) %>% 
      ungroup() %>% 
      select(-matches(match = year,ignore.case = FALSE)) %>% 
      pivot_wider(names_from = true_year,
                  values_from = .value,
                  names_prefix = "Y") %>% 
      rename_with(., ~gsub(replacement = "start",
                           pattern = paste0("Y",start_year),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "end",
                           pattern = paste0("Y",end_year),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "stratum_trend",
                           pattern = summary_regions,.x,
                           fixed = TRUE))
    
    
    
    tt <- indt %>% 
      group_by(.draw,stratum_trend) %>% 
      summarise(end = sum(end),
                start = sum(start),
                t = texp(end/start,ny = nyrs),
                ch = chng(end/start),
                .groups = "keep") %>% 
      group_by(stratum_trend) %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lu,names = FALSE),
                uci = quantile(t,uu,names = FALSE),
                percent_change = median(ch),
                p_ch_lci = quantile(ch,lu,names = FALSE),
                p_ch_uci = quantile(ch,uu,names = FALSE),
                prob_decline = prob_dec(ch,0),
                prob_decline_GT30 = prob_dec(ch,-30),
                prob_decline_GT50 = prob_dec(ch,-50),
                prob_decline_GT70 = prob_dec(ch,-70))%>% 
      rename_with(., ~gsub(replacement = summary_regions,
                           pattern = "stratum_trend",.x,
                           fixed = TRUE))
    
    
    
    
  }else{ #else is.null weights_df
    
    if(!is.null(strat)){
      indt <- samples %>% 
        filter(true_year %in% c(start_year,end_year)) %>% 
        #ungroup() %>% 
        select(-matches(match = year,ignore.case = FALSE)) %>% 
        pivot_wider(names_from = true_year,
                    values_from = .value,
                    names_prefix = "Y") %>% 
        rename_with(., ~gsub(replacement = "start",
                             pattern = paste0("Y",start_year),.x,
                             fixed = TRUE))%>% 
        rename_with(., ~gsub(replacement = "end",
                             pattern = paste0("Y",end_year),.x,
                             fixed = TRUE))%>% 
        rename_with(., ~gsub(replacement = "stratum_trend",
                             pattern = strat,.x,
                             fixed = TRUE))
      
      
      
      tt <- indt %>% 
        group_by(.draw,stratum_trend) %>% 
        summarise(t = texp(end/start,ny = nyrs),
                  ch = chng(end/start),
                  .groups = "keep") %>% 
        group_by(stratum_trend) %>% 
        summarise(trend = mean(t),
                  lci = quantile(t,lu,names = FALSE),
                  uci = quantile(t,uu,names = FALSE),
                  percent_change = median(ch),
                  p_ch_lci = quantile(ch,lu,names = FALSE),
                  p_ch_uci = quantile(ch,uu,names = FALSE),
                  prob_decline = prob_dec(ch,0),
                  prob_decline_GT30 = prob_dec(ch,-30),
                  prob_decline_GT50 = prob_dec(ch,-50),
                  prob_decline_GT70 = prob_dec(ch,-70))%>% 
        rename_with(., ~gsub(replacement = strat,
                             pattern = "stratum_trend",.x,
                             fixed = TRUE))
      
    }else{
      
      indt <- samples %>% 
        filter(true_year %in% c(start_year,end_year)) %>% 
        #ungroup() %>% 
        select(-matches(match = year,ignore.case = FALSE)) %>% 
        pivot_wider(names_from = true_year,
                    values_from = .value,
                    names_prefix = "Y") %>% 
        rename_with(., ~gsub(replacement = "start",
                             pattern = paste0("Y",start_year),.x,
                             fixed = TRUE))%>% 
        rename_with(., ~gsub(replacement = "end",
                             pattern = paste0("Y",end_year),.x,
                             fixed = TRUE))
      
      
      
      tt <- indt %>% 
        group_by(.draw) %>% 
        summarise(t = texp(end/start,ny = nyrs),
                  ch = chng(end/start),
                  .groups = "keep") %>% 
        summarise(trend = mean(t),
                  lci = quantile(t,lu,names = FALSE),
                  uci = quantile(t,uu,names = FALSE),
                  percent_change = median(ch),
                  p_ch_lci = quantile(ch,lu,names = FALSE),
                  p_ch_uci = quantile(ch,uu,names = FALSE),
                  prob_decline = prob_dec(ch,0),
                  prob_decline_GT30 = prob_dec(ch,-30),
                  prob_decline_GT50 = prob_dec(ch,-50),
                  prob_decline_GT70 = prob_dec(ch,-70))
    }
    
    
  }
  
  tt <- tt %>% 
    mutate(Region_type = summary_regions)
  return(tt)
}


# helper functions --------------------------------------------------------

p_lt <- function(x,th){
  length(which(x < th))/length(x)
}


p_neg <- function(x){
  length(which(x < 0))/length(x)
}

texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}




chng <- function(x){
  (x-1)*100
}

prob_dec <- function(ch,thresh){
  
  length(which(ch < thresh))/length(ch)
}
