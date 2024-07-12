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


# function to plot attack rates ------------------------------------------------

plot_attack_rate = function(VCD,
                            hosp,
                            file_path, 
                            hospital = F,
                            AR) {
  
  # add aggregated populations to data and calculate attack rates

  if (hospital == F) {
    fil = "1"
    AR_data = calc_attack_rates(VCD)
    y = "Symptomatic attack rate (%)"
  } else {
    fil = "2"
    AR_data = calc_hosp_rates(VCD = VCD, hosp = hosp)
    y = "Hospitalisation attack rate (%)"
    
  }

    AR_model = extract_model_results(AR)  



# plot attack rate by year ---------------------------------------------------
  
  AR_BD = AR_data %>%
    filter(serostatus != "both",
           age == "all",
           trial == "both",
           serotype == "all",
           month != 1000)
  
  
  AR_plot_BD =  AR_model %>%
    filter(group == "AR_BRD") %>%
    separate(name, into = c("serostatus","outcome", "time")) %>%
    filter(outcome == eval(rlang::parse_expr(fil))) %>% 
    mutate(serostatus= ifelse(serostatus == 1, "seronegative", "seropositive")) %>% 
    mutate(month = ifelse(time == 1, 12, ifelse(time == 2, 18, ifelse(time == 3, 
                  24, ifelse(time == 4, 36, ifelse( time == 5, 48, 54)))))) %>%  
    bind_rows(AR_BD) %>%
    ggplot(aes(x = as.numeric(month) , y = mean)) +
    geom_point(aes(color = serostatus,shape = type),
               position = position_dodge(width = 2),
               size = 3) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        color = serostatus,
        linetype = type
      ),
      position = position_dodge(width = 2),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = "Month", y = y) +
    scale_x_continuous(breaks = seq(0, 54, 6)) +
    scale_color_brewer(palette = "Paired")
    
  
# plot attack rate by serostauts and trial arm ---------------------------------
  
  pal_fill = scales::brewer_pal(palette = "Paired")(4)
  
  AR_VD = AR_data %>%
    filter(serostatus == "both",
           age == "all",
           trial != "both",
           serotype == "all",
           month != 1000)
  
  AR_plot_VD =  AR_model %>%
    filter(group == "AR_VRD") %>%
    separate(name, into = c("trial","outcome",  "time")) %>%
    filter(outcome == eval(rlang::parse_expr(fil))) %>%    
    mutate(trial= ifelse(trial == 1, "placebo", "TAK-003")) %>% 
    mutate(month = ifelse(time == 1, 12, ifelse(time == 2, 18, ifelse(time == 3, 
               24, ifelse(time == 4, 36, ifelse( time == 5, 48, 54)))))) %>%  
    bind_rows(AR_VD) %>%
    ggplot(aes(x = as.numeric(month) , y = mean)) +
    geom_point(aes(color = trial,shape = type),
               position = position_dodge(width = 2),
               size = 3) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        color = trial,
        linetype = type
      ),
      position = position_dodge(width = 2),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = "Month", y = "") +
    scale_x_continuous(breaks = seq(0, 54, 6)) +
    scale_color_manual(values = pal_fill[3:4]) +
    guides(shape = "none",
           linetype = "none")
  
  # plot attack rate by serotype, serostatus and trial arm -----------------------
  
  
  AR_BVK = AR_data %>%
    filter(serostatus != "both",
           age == "all",
           trial != "both",
           serotype != "all",
           month == 1000)
  
  AR_plot_BVK = AR_model %>%
    filter(group == "AR_BVKR") %>%
    separate(name, into = c("serostatus", "trial", "serotype", "outcome")) %>% 
    filter(outcome == eval(rlang::parse_expr(fil))) %>% 
    mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
           trial = ifelse(trial == 1, "placebo", "TAK-003"),
           serotype = ifelse(serotype == 1, "DENV-1",
                              ifelse(serotype == 2, "DENV-2",
                                     ifelse(serotype ==3, "DENV-3", 
                                            "DENV-4")))) %>% 
    bind_rows(AR_BVK) %>%
    unite(c(trial, serostatus), col = "x") %>%
    ggplot(aes(x = x, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = y) +
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(
      labels  = c(
        "Placebo \n Seronegative",
        "Placebo \n Seropositive",
        "TAK-003 \n Seronegative",
        "TAK-003 \n Seropositive"
      )
    ) +
    guides(shape = "none",
           linetype = "none")
  
  
# plot attack rate by age, serostatus and trial arm -----------------------
  
  
  AR_BVJ = AR_data %>%
    filter(serostatus != "both",
           age != "all",
           trial != "both",
           serotype == "all",
           month == 1000)
  
  AR_plot_BVJ = AR_model %>%
    filter(group == "AR_BVJR") %>%
    separate(name, into = c("serostatus", "trial", "age", "outcome")) %>% 
    filter(outcome == eval(rlang::parse_expr(fil))) %>% 
    mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
           trial = ifelse(trial == 1, "placebo", "TAK-003"),
           age = factor(age,
                        labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                        levels = c(1,2,3))) %>% 
    bind_rows(AR_BVJ) %>%
    unite(c(trial, serostatus), col = "x") %>%
    ggplot(aes(x = x, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = age,
        group = interaction(type, age)
      ),
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, age),
        linetype = type,
        color = age
      ),
      position = position_dodge(width =  0.5),
      width =  0.4,
      linewidth = 1
    ) +
    labs(x = " ", y = "") +
    scale_color_brewer(palette = "Accent") +
    scale_x_discrete(
      labels  = c(
        "Placebo \n Seronegative",
        "Placebo \n Seropositive",
        "TAK-003 \n Seronegative",
        "TAK-003 \n Seropositive"
      )
    )+
    guides(shape = "none",
           linetype = "none")
  

  # save grid plot of AR figures -----------------------------------------------
  
  AR_grid_plot = plot_grid(
    AR_plot_BD,
    AR_plot_VD,
    AR_plot_BVK,
    AR_plot_BVJ,
    ncol = 2,
    labels = c("a", "b", "c", "d"),
    align = "vh"
  )
  
  ggsave(
    plot = AR_grid_plot,
    filename = paste0(file_path, "/AR_grid_plot_outcome_", fil, ".png"),
    height = 42,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.9
  )
  
# plot at finest resolution
  
if(hospital == F){

AR_BVJD = AR_data %>%
  filter(serostatus != "both",
           age != "all",
           trial != "both",
           serotype == "all",
           month != 1000)

AR_BVKD = AR_data %>%
  filter(serostatus != "both",
         age == "all",
         trial != "both",
         serotype != "all",
         month != 1000)  

  
AR_plot_BVJD =  AR_model %>%
    filter(group == "AR_BVJRD") %>%
    separate(name, into = c("serostatus", "trial", "age", "outcome", "time")) %>% 
    filter(outcome == eval(rlang::parse_expr(fil))) %>% 
    mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
           trial = ifelse(trial == 1, "placebo", "TAK-003"),
           age = factor(age,
                        labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                        levels = c(1,2,3)),) %>% 
    mutate(month = ifelse(time == 1, 12, ifelse(time == 2, 18, 
                   ifelse(time == 3,24,  ifelse(time == 4, 36, 
                   ifelse( time == 5, 48, 54)))))) %>%
  filter(month <= 36) %>% 
  bind_rows(AR_BVJD) %>%
    ggplot(aes(x = month, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = age,
        group = interaction(type, age)
      ),
      position = position_dodge(width = 2),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, age),
        linetype = type,
        color = age
      ),
      position = position_dodge(width =  2),
      width =  0.4,
      linewidth = 1
    ) +
    facet_grid(serostatus ~ trial) +
    labs(x = "Month", y = y) +
    scale_color_brewer(palette = "Accent") 

AR_plot_BVKD = AR_model %>%
  filter(group == "AR_BVKRD") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome", "time")) %>% 
  filter(outcome == eval(rlang::parse_expr(fil))) %>% 
  mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
         trial = ifelse(trial == 1, "placebo", "TAK-003"),
         serotype = ifelse(serotype == 1, "DENV-1",
                           ifelse(serotype == 2, "DENV-2",
                                  ifelse(serotype ==3, "DENV-3", 
                                         "DENV-4")))) %>% 
  mutate(month = ifelse(time == 1, 12, ifelse(time == 2, 18, 
                 ifelse(time == 3,24,  ifelse(time == 4, 36, 
                 ifelse( time == 5, 48, 54)))))) %>%
  bind_rows(AR_BVKD) %>%
  ggplot(aes(x = month, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width = 4),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width = 4),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(serostatus ~ trial) +
  labs(x = "Month", y = y) +
  scale_color_brewer(palette = "Set2") 


  
ggsave(
  plot = AR_plot_BVJD,
  filename = paste0(file_path, "/AR_BVJD.png"),
  height = 42,
  width = 42,
  units = "cm",
  dpi = 300,
  scale = 0.9
)

ggsave(
  plot = AR_plot_BVKD,
  filename = paste0(file_path, "/AR_BVKD.png"),
  height = 42,
  width = 60,
  units = "cm",
  dpi = 300,
  scale = 0.9
)  
  
}
  
  
}

# function to plot vaccine efficacy --------------------------------------------

plot_VE  = function(file_path,
                    hospital = F,
                    serostatus,
                    VE = NULL
) {
  
    VE_model = extract_model_results(VE)  

 
  VE_model$lower = ifelse( VE_model$lower < -200, -200,  VE_model$lower )
  
  if(hospital ==F){ 
    fil = "1"
  y = "Efficacy against symptomatic disease (%)"
  } else {
    fil= "2" 
    y = "Efficacy against hospitalisation (%)"
    }
  

  
# plot VE by serostatus for each serotype --------------------------------------
  
  VE_BKT =  VE_model %>%
    filter(group == "VE_BKRT") %>%
    separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>% 
    filter(outcome == eval(rlang::parse_expr(fil))) %>% 
    mutate(serotype = factor(serotype, labels = c("DENV-1","DENV-2", "DENV-3", "DENV-4"))) %>% 
    mutate(month = as.numeric(month)) 
  
  VE_BKT$serostatus = factor(VE_BKT$serostatus, labels = serostatus)

  
  VE_BKT_plot = VE_BKT %>% 
    ggplot(aes(x = month , y = mean)) +
    geom_line() +
    geom_ribbon(
      aes(
        ymin = lower ,
        ymax = upper ,
      ), alpha = 0.5
    ) +
    labs(x = "", y = y) +
    scale_x_continuous(breaks = seq(0, 54, 6)) +
    geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
    facet_grid(serostatus~serotype) 
  
  ggsave(
    plot = VE_BKT_plot,
    filename = paste0(file_path, "/VE_BKT_plot_outcome_", fil, ".png"),
    height = 21,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.8
  )
  
# plot VE  age, serostatus for each serotype  ----------------------------------
  
  VE_BKJT =  VE_model %>%
    filter(group == "VE") %>%
    separate(name, into = c("serostatus", "serotype", "age","outcome", "month")) %>%
    filter(outcome == eval(rlang::parse_expr(fil))) %>% 
    mutate(
      age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs")),
      serotype = factor(serotype, labels = c("DENV-1", "DENV-2", "DENV-3", "DENV-4")),
      month = as.numeric(month))
 
  
  VE_BKJT$serostatus = factor(VE_BKJT$serostatus, labels = serostatus)
  
  VE_BKJT_plot = VE_BKJT %>% 
    filter(serostatus != "multitypic") %>% 
    ggplot(aes(x = month , y = mean)) +
    geom_line(
      aes(
        color = serostatus, 
      )
    ) +
    geom_ribbon(
      aes(
        ymin = lower ,
        ymax = upper ,
        fill = serostatus),
      alpha = 0.5
    ) +
    labs(x = "", y = y) +
    scale_x_continuous(breaks = seq(0, 54, 6)) +
    geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired") +
    facet_grid(serotype ~ age) 
  
  ggsave(
    plot = VE_BKJT_plot,
    filename = paste0(file_path,"/VE_BKJT_plot_outcome_", fil, ".png"),
    height = 41,
    width = 42,
    units = "cm",
    dpi = 300,
    scale = 0.8
  )
}


# function to plot titres ----------------------------------------------------------------------------------------
plot_titres = function(file_path, n = NULL){
  

  # format model output and add data 
  model_titres = n %>%   
    as.data.frame() %>%  
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = log(value)) %>%  
    summarise(
      lower = quantile(value, 0.025),
      mean = mean(value),
      upper = quantile(value, 0.975)
    ) %>%  
    filter(grepl("n\\[", name)) %>% 
    separate(name, into = c(NA, "name"), sep = "\\[") %>% 
    separate(name, into = c("name", NA), sep = "\\]") %>% 
    separate(name, into = c("serostatus", "serotype", "time")) %>%  
    mutate(time= as.numeric(time),
           serostatus = ifelse(serostatus == 1, "seronegative", 
                               "seropositive"),
           serotype = paste0("DENV-", serotype))
    
  
 plot_titres = model_titres %>% 
    ggplot(aes(x = time, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper,
                    ), alpha = 0.2) +
    theme(legend.position ="top") +
    facet_grid(serostatus~ serotype) +
    labs(x = "Month", y = "log titre")
  
  
  ggsave(plot_titres, filename = paste0(file_path,"/plot_titres.jpg"),
         width = 40, height = 20, unit ="cm")
  
  
}


# function to plot everything --------------------------------------------------
plot_output = function(VCD=NULL,
                       hosp=NULL,
                       titre_model=T,
                       file_path,
                       AR = NULL,
                       VE = NULL,
                       n = NULL,
                       serostatus= c("seronegative", "monotypic", "multitypic")) {
  library(tidyverse)
  library(Hmisc)
  library(cowplot)
  
  theme_set(
    theme_bw() +
      theme(
        text = element_text(size = 16),
        legend.position ="top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)
      ))
  
  # data
if(is.null(VCD)) {
  VCD =  read.csv("data/processed/vcd_data.csv")
  VCD = factor_VCD(VCD)
} 

if(is.null(hosp)) {
  hosp = read.csv("data/processed/hosp_data.csv")
  hosp = factor_VCD(hosp)
}
  
if (is.null(AR)){  
 AR = readRDS(paste0(file_path, "/AR.RDS"))
}
if (is.null(VE)){  
 VE = readRDS(paste0(file_path, "/VE.RDS"))
}
if(is.null(n)){
 n = readRDS(paste0(file_path, "/n.RDS"))
}

# plot attack rates 
  plot_attack_rate(
    VCD =  VCD,
    hosp = hosp,
    file_path = file_path,
    AR = AR
  )
  
  # plot hosp attack rates 
  plot_attack_rate(
    VCD =  VCD,
    hosp = hosp,
    file_path = file_path,
    hospital = T,
    AR = AR
  )
  
  # plot VE
  plot_VE(file_path = file_path,
          serostatus = serostatus,
          VE = VE)
  
  # plot VE against hosp
  plot_VE(file_path = file_path,
          hospital = T,
          serostatus = serostatus,
          VE=VE)
  
  
# if titre plot plot imputed titres 
  if (titre_model == T) {
    plot_titres(file_path = file_path,
                n=n) }
}
