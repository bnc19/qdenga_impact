################################################################################
# Script to plot example individual impact using demo code 
################################################################################


library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

pathin = "demo/processed" 
pathout = "demo/figures"
sp9 = c(20, 40, 60, 80) ## DEMO ONLY RUNS 4 TRANSMISSION SETTINGS

source("R/process_simulations.R")


# plotting  
col = c("#67A9CF", "#C51B8A", "#99CC99")
col_sero = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")
pd = position_dodge(0.9)

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      legend.background = element_rect(fill = NA, colour = NA),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
)

# Read in data 
cohort_file = "cohort_0.8vc_prop_averted/csv/cohort_prop_averted_"
BRA_df.cohort = read.csv(paste0(pathin, "/VS_D15/", cohort_file, "overall_uncert_BRA_6_95perc.csv"))
BRA_df.cohort.m = read.csv(paste0(pathin, "/VS_D15/", cohort_file, "mean_uncert_BRA_6_95perc.csv"))


# serostatus impact 
col_names_b = c(
  "out_cohort_dis_all_vacc",
  "out_cohort_dis_all_vacc_pos",
  "out_cohort_dis_all_vacc_neg",
  "out_cohort_sdis_all_vacc",
  "out_cohort_sdis_all_vacc_pos",
  "out_cohort_sdis_all_vacc_neg"
)

################################################################################
################# Individual overall impact and by serostatus ##################
################################################################################

# BRA 
VS_over_df_B_b = format_cohort(BRA_df.cohort, col_names_b, rev(1:9 * 10) , "o")
VS_mean_df_B_b = format_cohort(BRA_df.cohort.m, col_names_b, rev(1:9 * 10) , "m")

# VS  
VS_both_B_b = VS_mean_df_B_b %>% 
  bind_rows(VS_over_df_B_b) %>%
  filter(sp9 %in% c(20, 40, 60, 80)) %>%  ### DEMO ONLY RUNS SP9 20, 40, 60, 80
  mutate(type = rep(c(rep("symptomatic", 3), rep("hospitalised", 3)), 8),
         baseline = rep(c("all", "seropositive", "seronegative"), 16),
         vacc_model = factor("VS_D15")) %>% 
  select(-name) %>% 
  pivot_wider(names_from = metric, values_from = mean:up)

  
BRA_VS_D15_plot_b = VS_both_B_b %>% 
  ggplot() + 
  geom_line(aes(x=sp9, y = mean_o*100)) +
  geom_ribbon(aes(x=sp9, ymin=low_o*100, ymax=up_o*100), alpha=0.2) +
  geom_ribbon(aes(x=sp9, ymin=low_m*100, ymax=up_m*100), alpha=0.4) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_grid(rows = vars(baseline), cols = vars(type), 
             scales="free_y", switch = 'y') +
  ylab("Proportion averted (%)")+
  xlab("Average seroprevalence \nin 9-year-olds") +
  ylim(c(-50,100))+
  theme(legend.position = "none") +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  scale_x_continuous(breaks = c(20, 40, 60, 80))
  

ggsave(filename =  paste0(pathout, "/BRA_prop_averted_cohort_VS_D15.png"), 
       plot = BRA_VS_D15_plot_b, 
       device = "png",
       width = 8.8,
       height = 10,
       units = "cm",
       dpi = 600)

