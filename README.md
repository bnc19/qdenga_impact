# Qdenga impact


# About this repository
This repository contains the code needed to reproduce the main results presented in the manuscript: Efficacy, public health impact and optimal use of the Takeda dengue vaccine. 

All data used for model calibration is available in the repository under **Efficacy/data**. 

All required packages and their versions are available in the DESCRIPTION file. Instructions to download and install cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan). Instructions to download and install Rtools can be found [here](https://cran.r-project.org/bin/windows/Rtools/).

# Repository structure

This repository is divided into two sections: 

* *Efficacy* - code to fit the cohort Bayesian survival model to published clinical incidence data in Stan to estimate Qdenga's efficacy. **Note that the code is currently set up to run the final model locally for 1000 iterations, with the first 500 discarded as burn-in. This is to demo the models quickly (~15 minutes on a standard computer). To obtain the results provided in the manuscript, four chains were run for 10,000 iterations, discarding the first 5000 iterations as burn-in.**
  
* *Impact* - code to sample from the posterior distribution obtained in *efficacy* and simulate impact using a stochastic compartmental dengue transmission model. **Note that the code is currently set up to run four vaccine scenarios using 10 simulations for each of 50 posterior samples. To obtain the results provided in the manuscript, 50 simulations were run for each of 200 posterior samples across 4032 scenarios.**


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

* *final_model.csv*: Rstan model to fit the Bayesian survival model to case data, assuming binomial and multinomial likelihoods. Flags in the data block turn on and off different parameters used to test the 31 model variants presented in the manuscript. 


### R
* *factor_data.R*: functions to format the data.
* *format_stan_data.R*: function to format all data input to the Stan model used for model calibration. Outputs data in a list format.
* *run_model_fitting.R*: function to calibrate the stan model to the publically available data.
* *plot_model_outputs.R*: functions to calculate the observed attack rates and format model outputs for plotting. 


## Impact

For the main manuscript, we ran the following scenarios: 

  - Assumed demography: Brazil or Philippines
  - Vaccine mechanism of protection: only against disease or also against infection
  - Duration of vaccine decay: 15 years or 5 years
  - Age of vaccination: 6 to 12 years
  - Transmission intensity setting (average 9 years seroprevalence): 10% to 90%
  - Vaccine coverage: 20%, 40%, 60%, 80%
  - Pre-vaccination screening to identify seropositive individuals: yes/no

The code in this repo shows how to run four scenarios (Brazil, the vaccine protects disease only for 15 years; vaccination of age 6 with 80% coverage, without pre-vaccination screening; in transmission settings with 20%, 40%, 60%, and 80% 9-year-old seropositivity).

For these scenarios, the demo samples the posterior VE estimates 50 times and runs 10 stochastic simulations per sample.

**This folder contains**

* *ps_new24.csv*: 200 posterior samples of the survival model, obtained by running the code in the *efficacy* section of this repo. 

### Scripts
* *1.run_impact_simualtions.R*: Main script to run the transmission model with and without vaccination to estimate the impact of Qdenga vaccination. The script is set up to run the model *model_stock_tak_2.R*, parameterised using 50 posterior samples from *ps_new24.csv*. For each posterior sample, 10 model simulations are run to demonstrate the code. To recreate the results from the paper, 50 simulations were run for each of 200 posterior samples. 
* *2.process_pop_results.R*: Process output of model simulations to estimate the proportion of symptomatic and hospitalised cases averted in the entire population over 10 years of vaccination (population impact), using the demo scenarios. 
* *3.plot_pop_impact.R*: Plot population impact to reconstruct Figure 2a (VS_D15, sp9 20, 40, 60, 80) in the manuscript. 
* *4.process_ind_results.R*: Process output of model simulations to estimate the proportion of symptomatic and hospitalised cases averted in the first vaccinated cohort, over 10 years of vaccination (individual impact), using the demo scenarios. 
* *5.plot_pop_impact.R*: Plot individual impact by serostatus to reconstruct Figure 4 (VS_D15, sp9 20, 40, 60, 80) in the manuscript. 


### Models
 
* *model_stoch_tak_2.R*:  The transmission model used to estimate the population- and individual-level impact of Qdenga.


### R
* *run_who_cluster.R*: functions to equilibrate the model and run the model with and without vaccination.
* *process_demog.R*: function to import the demography files. 
* *config_params.R*: function to define all transmission model parameters, including the VE posterior draws.
* *process_simulations.R*: functions to process simulation and calculate impact. 

### demog
* age, births, initial population, life expectancy, and mortality demography data for Brazil and the Philippines, used to parameterise the births and deaths rate of the simulation model *model_stoch_tak_2.R*

# References 

1.	Biswal, S. et al. Efficacy of a tetravalent dengue vaccine in healthy children and adolescents. N. Engl. J. Med. 381, 2009–2019 (2019).
2.	Biswal, S. et al. Efficacy of a tetravalent dengue vaccine in healthy children aged 4–16 years: a randomised, placebo-controlled, phase 3 trial. The Lancet 395, 1423–1433 (2020).
3.	López-Medina, E. et al. Efficacy of a dengue vaccine candidate (TAK-003) in healthy children and adolescents 2 years after vaccination. J. Infect. Dis. 225, 1521–1532 (2022).
4.	Rivera, L. et al. Three-year efficacy and safety of Takeda’s dengue vaccine candidate (TAK-003). Clin. Infect. Dis. 75, 107–117 (2022).
5.	Tricou, V. et al. Long-term efficacy and safety of a tetravalent dengue vaccine (TAK-003): 4·5-year results from a phase 3, randomised, double-blind, placebo-controlled trial. Lancet Glob. Health 12, e257–e270 (2024).

