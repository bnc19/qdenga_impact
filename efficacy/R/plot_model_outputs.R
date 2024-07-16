# function to add extra populations to VCD data and calculate AR --------------- 

calc_attack_rates = function(VCD, titre = F) {
  
# VCD by serostatus and month 
  VCD_BD = VCD %>%  
    filter(serostatus != "both", age == "all", serotype== "all") %>% 
    group_by(age, serotype, year,month,serostatus) %>%  
    summarise(Y= sum(Y),
              N = sum(N)) %>% 
    mutate(trial = "both")
    

# VCD by trial and month 
 VCD_VD = VCD %>%  
    filter(serostatus == "both", age == "all", serotype== "all") %>% 
    group_by(age, serotype,year, month,trial) %>%  
    summarise(Y= sum(Y),
              N = sum(N)) %>% 
   mutate(serostatus = "both")
 
# VCD by serostatus and serotype and trial  
 VCD_BVK = VCD %>%  
   filter(serostatus != "both", age == "all", serotype != "all") %>% 
   group_by(age, serotype,trial,serostatus) %>%  
   summarise(Y= sum(Y),
             N = mean(N)) %>% 
   mutate(month = 1000)
 
 # VCD by serostatus and age and trial  
 VCD_BVJ = VCD %>%  
   filter(serostatus != "both", age != "all", serotype == "all", month < 48) %>% 
   group_by(age, serotype,trial,serostatus) %>%  
   summarise(Y= sum(Y),
             N = mean(N)) %>% 
   mutate(month = 1000)

  # combine all data 
  
  out = VCD %>%  
    filter(!is.na(Y))  %>%  # remove year 4 / 4.5 age data as NA
    bind_rows(VCD_BD, VCD_VD,VCD_BVK,VCD_BVJ) %>% 
    mutate(outcome = "VCD") %>% 
    filter(!is.na(N)) %>%  # remove serotype age 
    mutate(age = factor(age)) %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
    return(out)
}

# function to calculate hospital attack rate -----------------------------------

calc_hosp_rates = function(hosp, VCD) {
  
  BVD_pop = VCD %>%  
    filter(age == "all", serotype == "all", serostatus != "both")  %>% 
    select(-Y, - X)
  
  # Hosp by serostatus, trial and month after months 
  BV_4 = hosp %>% 
    filter(age == "all", serotype != "all", month > 36) %>% 
    group_by (serostatus, trial, month) %>% 
    summarise(Y = sum(Y))  
  
  BVD = hosp %>%
    filter(age != "all", serotype == "all") %>% 
    group_by (serostatus, trial, month) %>% 
    summarise(Y=sum(Y)) %>% # BVD hosp upto 36 months  
    bind_rows(BV_4) %>% 
    left_join(BVD_pop)
  
  BD = BVD %>% 
    group_by(serostatus, month) %>% 
    summarise(Y = sum(Y),
              N= sum(N)) %>% 
    mutate(trial = "both", serotype = "all", age = "all") 
  
  VD = BVD %>% 
    group_by(trial, month) %>% 
    summarise(Y = sum(Y),
              N= sum(N)) %>% 
    mutate(serostatus = "both", serotype = "all", age = "all")
  
  # Hosp by age, serostatus, trial 
  
  BVJ_pop = VCD %>%  
    filter(age != "all", serotype == "all", serostatus != "both", month <= 36) %>% 
    group_by(age, trial, serostatus) %>% 
    summarise(N = mean(N))
  
  
  BVJ = hosp %>%  
    filter(age != "all", serotype == "all") %>% 
    group_by(serostatus, trial, age) %>% 
    summarise(Y = sum(Y)) %>% 
    left_join(BVJ_pop) %>% 
    mutate(month = 1000,
           serotype = "all")
  
  
  # Hosp by serotype, serostatus, trial 
  
  BVK_pop = VCD %>%  
    filter(age == "all", serotype == "all", serostatus != "both") %>% 
    group_by(trial, serostatus) %>% 
    summarise(N = mean(N))
  
    BVK = hosp %>%  
    filter(age == "all", serotype != "all") %>% 
    group_by(serostatus, trial, serotype) %>% 
    summarise(Y = sum(Y)) %>% 
    left_join(BVK_pop)%>% 
      mutate(month = 1000, 
             age = "all")
  
  # combine all data 
  
  out = bind_rows(VD, BD, BVD, BVJ, BVK) %>%
    mutate(outcome = "hosp",
           year = month / 12) %>% 
    mutate(age = factor(age)) %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  return(out)
}




# function to extract model results --------------------------------------------
extract_model_results = function(fit_ext){
  

    out =  fit_ext %>%  
      as.data.frame() %>%  
      mutate(ni = row_number()) %>%  
      pivot_longer(cols = - ni) %>%  
      group_by(name) %>% 
      mutate(value = value * 100) %>%  
      summarise(
        lower = quantile(value, 0.025, na.rm = T),
        mean = mean(value, na.rm = T),
        upper = quantile(value, 0.975, na.rm = T)
      )

 
      out2 = out %>%  
        separate(name, into = c("group", "name"), sep = "\\[") %>% 
        separate(name, into = c("name", NA), sep = "\\]") %>% 
        mutate(type = "model") 
    
        return(out2)
}



