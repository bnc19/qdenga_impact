################################################################################
# Script to plot example population impact using demo code 
################################################################################


library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

pathout = "demo/processed"
dir.create("demo/figures")

serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")
col = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
)



# scenarios 
countries = c("BRA")
vacc.age = 6
vacc.cov = 0.8 
cumy = c(10) # show impact over 10 years 
cv.y = 0.8


################################################################################
###################### Overall averted population  #############################
################################################################################

df.out.cum = read.csv(paste0(pathout, "/VS_D15/", "1.2956_0.6833_1.0553/","/csv/Summary_cumulative_impacts_all_serotypes.csv"))
df.out.cum$vacc_mode = "VS_D15"

BRA_pop_impact = df.out.cum %>% 
  filter(year == 10, pop == "pop", base  == "all") %>%
  separate(serop9, into = c(NA, "serop9")) %>% 
  mutate(serop9 = as.numeric(serop9),
         type = ifelse(type == "hospitalisation", 
                       "hospitalised", "symptomatic")) %>% 
  ggplot() +
  geom_line(aes(x=serop9, y=mean*100)) +
  geom_ribbon(aes(x=serop9, ymin=low*100, ymax=upper*100), alpha = 0.3) +
  facet_wrap(~ type, ncol = 2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ylab(paste0("Proportion averted \nover ", cumy, " years (%)")) +
  xlab("Average seroprevalence in 9-year-olds")

ggsave(filename = "demo/figures/BRA_pop_impact_VS_D15.png",
       plot = BRA_pop_impact,
       device = "png",
       width = 18,
       height = 10,
       units = "cm",
       dpi = 600)
