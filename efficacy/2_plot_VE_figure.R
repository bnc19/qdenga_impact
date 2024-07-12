# Script to plot the main figure including model fit and the vaccine efficacy 
# estimate (Fig. 2) 

# load functions 
library(tidyverse)
library(Hmisc)

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 14),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) 
  
)

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
age_fill2 = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")
trial_fill = c("#C51B8A", "#99CC99")


# source files 
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
path = "output/final/M30_FINAL/"

# thesis path 
# path = "output/final/M32_FINAL/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
VCD =  read.csv("data/processed/vcd_data.csv")
hosp = read.csv("data/processed/hosp_data.csv")

# Tidy data 
VCD = factor_VCD(VCD)
hosp = factor_VCD(hosp)

# summarise over iterations 
VE_model = extract_model_results(VE)  
AR_model = extract_model_results(AR)  

# calculate attack rates from data 
AR_data = calc_attack_rates(VCD)
HR_data = calc_hosp_rates(VCD=VCD, hosp = hosp)
data= bind_rows(AR_data, HR_data)

# plot VE 
VE_plot =  VE_model %>%
  filter(group == "VE") %>%
  separate(name, into = c("serostatus", "serotype", "age","outcome", "month")) %>%
   filter(serostatus != 3) %>% # low power to est. multitypic so don't plot here 
  mutate(
    age = ifelse(age ==1, 1, 2), # age groups 2 and 3 have same VE 
    age = factor(age, labels = c("4-5yrs","6-16yrs")),
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")),
    serotype = factor(serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, labels = c("seronegative", "monotypic"))) %>% 
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = age), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age), alpha = 0.6) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 54,12)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  facet_grid(outcome +  serostatus ~ serotype, scales = "free_y" ) +
  theme(legend.position = c(0.09,0.32),
        legend.title = element_blank()) +
  scale_fill_manual(values = age_fill2) +
  scale_colour_manual(values = age_fill2) 

#Plot AR  across time 
AR_VRD = data %>%
  filter(serostatus == "both",
         age == "all",
         trial != "both",
         serotype == "all") %>%
  select(-serotype,-age,-serostatus,-X, -year) %>%
  mutate(outcome = factor(
    outcome, levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")))

time_plot =  AR_model %>%
  filter(group == "AR_VRD") %>%
  separate(name, into = c("trial", "outcome", "month")) %>%
  mutate(
    trial = factor(trial, labels = c("placebo", "vaccine")),
    month = ifelse(month == 1, 12,
                   ifelse(
                     month == 2, 18,
                     ifelse(month == 3, 24,
                            ifelse(month == 4, 36,
                                   ifelse(month == 5, 48, 54)))
                   )),
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))
  ) %>%
  bind_rows(AR_VRD) %>%
  mutate(month = as.numeric(month)) %>%
  ggplot(aes(x = month, y = mean)) +
  geom_point(aes(color = trial, shape = type),
             position = position_dodge(width = 5),
             size = 2) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      color = trial,
      linetype = type
    ),
    position = position_dodge(width = 5),
    width =  0.4,
    linewidth = .8
  ) +
  labs(x = "Month", y = "Attack rate (%)") +
  scale_colour_manual(values = trial_fill) +
  scale_x_continuous(breaks = c(12,18,24,36,48,54)) +
  facet_wrap(~outcome, ncol = 1, strip.position="right") +
  theme(legend.position =  c(0.73,0.31),
        legend.title = element_blank()) 

# Plot AR by age and trial and serostatus 

AR_BVJR = data %>%
  filter(serotype == "all" |
           age != "all", serostatus != "both", month == 1000) %>%
  select(-month,-X, -year,-serotype) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  ))

age_plot =  AR_model %>%
    filter(group == "AR_BVJR") %>%
    separate(name, into = c("serostatus", "trial", "age", "outcome")) %>% 
    mutate(outcome = factor(outcome,labels = c("symptomatic", "hospitalised"))) %>% 
    mutate(serostatus = factor(serostatus, labels=c("seronegative", "seropositive")),
           trial = factor(trial, labels=c("placebo", "vaccine")),
           age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                        levels = c(1,2,3))) %>% 
    bind_rows(AR_BVJR) %>%
    unite(c(trial, serostatus), col = "x") %>%
    ggplot(aes(x = x, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = age,
        group = interaction(type, age)
      ),
      position = position_dodge(width = 0.7),
      size = 2
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, age),
        linetype = type,
        color = age
      ),
      position = position_dodge(width =  0.7),
      width =  0.4,
      linewidth = .8
    ) +
    labs(x = " ", y ="" ) +
  scale_colour_manual(values = age_fill) +
  scale_x_discrete(
      labels  = c(
        "P\nSN",
        "P\nSP",
        "V\nSN",
        "V\nSP"
      )
    ) +
    guides(shape = "none",
           linetype = "none") +
  facet_wrap(~outcome, ncol = 1, strip.position="right") +
  theme(legend.position =  c(0.7,0.35),
        legend.title = element_blank()) 
    

# Plot AR by serotype and trial and serostatus 
  
AR_BVKR = data %>%
  filter(serotype != "all" |
           age == "all", serostatus != "both", month == 1000) %>%
  select(-month,-X, -year,-age) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  ))

  
serotype_plot = AR_model %>%
  filter(group == "AR_BVKR") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome")) %>%
  mutate(outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))) %>%
  mutate(
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    trial = factor(trial, labels = c("placebo", "vaccine")),
    serotype = factor(serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4"))) %>%
  bind_rows(AR_BVKR) %>% 
    unite(c(trial, serostatus), col = "x") %>%
    ggplot(aes(x = x, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 0.7),
      size = 2
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  0.7),
      width =  0.4,
      linewidth = .8
    ) +
    labs(x = " ", y =" " ) +
  scale_colour_manual(values = serotype_fill) +
  scale_x_discrete(
      labels  = c(
        "P\nSN",
        "P\nSP",
        "V\nSN",
        "V\nSP"
      )
    ) +
    guides(shape = "none",
           linetype = "none") +
  facet_wrap(~ outcome, ncol = 1, strip.position="right") +
  theme(legend.position = c(0.72,0.31),
        legend.title = element_blank()) 


# combine all plots 
g1 = cowplot::plot_grid(time_plot, NULL, age_plot,  NULL,
                        serotype_plot, ncol=5, rel_widths = c(1,-0.002,1,-0.002,1),
                        axis = "tblr", align = "h",
                        labels = c("a","", "b","", "c"))

g2 = cowplot::plot_grid(g1, VE_plot, rel_heights = c(1,1.4),
                        ncol =1, labels = c("", "d"))  


ggsave(
  plot = g2,
  filename =  "output/figures/main_fit_ve_fig.png",
  height = 21,
  width = 18,
  units = "cm",
  dpi = 600
)


