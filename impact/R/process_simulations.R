
# function to estimate proportion averted
prop_averted <- function(t.nv, t.v){
  
  # prop averted over stoch realisations
  diff.v <- (t.nv - t.v)/t.nv
  
  # set Inf values (divided by 0) to NA
  diff.v[which(is.infinite(diff.v))] <- NA
  diff.v[which(is.nan(diff.v))] <- NA
  
  res.diff <- apply(diff.v, c(1,3), sum)/dim(diff.v)[2]
  
  # yearly stats over param variability 
  res.diff.mean <- rowMeans(res.diff, na.rm = TRUE)
  res.diff.l95 <- rowQuantiles(res.diff, probs=0.025, na.rm = TRUE)
  res.diff.q50 <- rowQuantiles(res.diff, probs=0.5, na.rm = TRUE)
  res.diff.u95 <- rowQuantiles(res.diff, probs=0.975, na.rm = TRUE)
  
  return(list(res.diff.mean, res.diff.l95, res.diff.q50, res.diff.u95))
  
}




#  function to format cohort data to plot individual risk  / benefit 

format_cohort = function(data_cohort, col_names, sp9, metric){
  out = data_cohort %>% 
    select(any_of(col_names)) %>%
    mutate(measure = c(rep("mean", 9) , rep("low", 9), rep("up", 9)),
           sp9 = rep(sp9, 3),
           metric = metric) %>% 
    pivot_longer(cols = - c(measure,sp9,metric)) %>% 
    pivot_wider(names_from = measure, values_from = value)
  
  return(out)
}




