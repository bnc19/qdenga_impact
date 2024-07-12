# Functions relating to trial data formatting and plotting ---------------------


# 1: function to factorise VCD data when importing it at the start of any script
factor_VCD = function(VCD){
  VCD = VCD  %>%  
    mutate(age = factor(age, levels = c("4-5yrs", "6-11yrs", "12-16yrs", "all")), 
           serostatus = factor(serostatus,
                               levels = c("SN", "SP", "both"),
                               labels = c("seronegative", "seropositive", "both")),
           serotype = factor(serotype,levels = c("D1", "D2", "D3", "D4", "all"),
                             labels = c("DENV1", "DENV2", "DENV3", "DENV4", "all")),
           trial = factor(trial, levels = c("P", "V","both"), 
                          labels = c("placebo", "vaccine", "both")),
           outcome = factor(outcome, labels = c("symptomatic", "hospitalised"),
                            levels = c("symp", "hosp")))
  
  return(VCD)
}



# 2: function to factorise serology data when importing it at the start of any script
factor_serology= function(data){
  data = data  %>%  
    mutate(serostatus = factor(serostatus,
                               levels = c("SN", "SP"),
                               labels = c("seronegative", "seropositive")),
           serotype = factor(serotype),
           trial = factor(trial, levels = c("Placebo", "TAK"), 
                          labels = c("placebo", "vaccine")))
  
  return(data)
}

