# Qdenga impact


# About this repository
This repository contains all the code needed to reproduce the results presented in the manuscript: Efficacy, public health impact and optimal use of the Takeda dengue vaccine. 

All data used for model calibration is available in the repository under *efficacy/data*. 

All required packages and their versions are available in the DESCRIPTION file. Instructions to install and download cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan).

# Repository structure

This repository is divided into two sections: 

* *efficacy* - code to fit the cohort Bayesian survival model to published clinical incidence data in Stan to estimate Qdenga's efficacy. **Note that the code is currently set up to run the final model locally for 1000 iterations with the first 500 discarded as burn in. This is to demo the models quickly (~4 hours on a standard computer). To obtain the results provided in the manuscript, four chains were run for 10,000 iterations, discarding the first 5000 iterations as burnin.**
  
* *impact* - code to sample from the posterior distributions obtained in *efficacy* and simulate impact using a stochastic compartmental dengue transmission model. 

## efficacy

**This folder contains**

### Data

* *hosp_data.csv*: hospital VCD case data extracted from published phase III clinical trial (1-5)
* *vcd_data.csv*: symptomatic VCD case extracted from published phase III clinical trial (1-5)
* *n0_new.csv*: initial neutralising antibody titres induced by Qdenga vaccination
* *seropositive_by_age_baseline.csv*: number of baseline seropositive individuals by age group 
