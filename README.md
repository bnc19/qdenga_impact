# Qdenga impact


# About this repository
This repository contains all the code needed to reproduce the results presented in the manuscript: Efficacy, public health impact and optimal use of the Takeda dengue vaccine. 

All data used for model calibration is available in the repository under *Efficacy/data*. 

All required packages and their versions are available in the DESCRIPTION file. Instructions to install and download cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan).

# Repository structure

This repository is divided into two sections: 

* *Efficacy* - code to fit the cohort Bayesian survival model to published clinical incidence data in Stan to estimate Qdenga's efficacy. **Note that the code is currently set up to run the final model locally for 1000 iterations with the first 500 discarded as burn in. This is to demo the models quickly (~15 minutes on a standard computer). To obtain the results provided in the manuscript, four chains were run for 10,000 iterations, discarding the first 5000 iterations as burnin.**
  
* *Impact* - code to sample from the posterior distributions obtained in *efficacy* and simulate impact using a stochastic compartmental dengue transmission model. 

## Efficacy

**This folder contains**

### Scripts
* *1_VE_model_fitting.R*: Main script to calibrate the stan model to the publicly available data. The script is set up to run the final model (M30) using 1000 iterations. All other model variants presented in the manuscript can be run below. All model variants are run using the same stan model (*final_model.csv*)  with different parameters turned on and off, using the flags. 
* *2_plot_VE_figure.R*: Main script to plot the attack rates and vaccine efficacy estimated in *1_VE_model_fitting.R*, to reconstruct Figure 1 in the manuscript. 

### Data

* *hosp_data.csv*: hospital VCD case data extracted from the published phase III clinical trial [1-5]. 
* *vcd_data.csv*: symptomatic VCD case extracted from the published phase III clinical trial [1-5]. 
* *n0_new.csv*: initial neutralising antibody titres induced by Qdenga vaccination, fitted from the published phase III clinical trial [1-5]. 
* *seropositive_by_age_baseline.csv*: number of baseline seropositive individuals by age group, extracted from the published phase III clinical trial [1-5].

### Models

* *final_model.csv*: Rstan model to fit a the Bayesian survival model to case data, assuming binomial and multinomial likelihoods. Flags in the data block turn on and off different parameters, used to test the 31 model variants presented in the manuscript. 


### R
* *factor_data.R*: functions to format the data.
* *format_stan_data.R*: function to format all data input to the Stan model used for model calibration. Outputs data in a list format.
* *run_model_fitting.R*: function to calibrate the stan model to the publically available data.
* *plot_model_outputs.R*: functions to calculate the observed attack rates and format model outputs for plotting. 

# References 

1.	Biswal, S. et al. Efficacy of a tetravalent dengue vaccine in healthy children and adolescents. N. Engl. J. Med. 381, 2009–2019 (2019).
2.	Biswal, S. et al. Efficacy of a tetravalent dengue vaccine in healthy children aged 4–16 years: a randomised, placebo-controlled, phase 3 trial. The Lancet 395, 1423–1433 (2020).
3.	López-Medina, E. et al. Efficacy of a dengue vaccine candidate (TAK-003) in healthy children and adolescents 2 years after vaccination. J. Infect. Dis. 225, 1521–1532 (2022).
4.	Rivera, L. et al. Three-year efficacy and safety of Takeda’s dengue vaccine candidate (TAK-003). Clin. Infect. Dis. 75, 107–117 (2022).
5.	Tricou, V. et al. Long-term efficacy and safety of a tetravalent dengue vaccine (TAK-003): 4·5-year results from a phase 3, randomised, double-blind, placebo-controlled trial. Lancet Glob. Health 12, e257–e270 (2024).

