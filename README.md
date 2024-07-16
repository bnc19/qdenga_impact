# Qdenga impact


# About this repository
This repository contains all the code needed to reproduce the results presented in the manuscript: Efficacy, public health impact and optimal use of the Takeda dengue vaccine. 

All data used for model calibration is available in the repository under _efficacy/data_. 

The model fitting requires the package cmdStan. Instructions to install and download cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan).

# Repository structure

This repository is divided into two sections: 

* *efficacy* - code to fit the cohort Bayesian survival model to published clinical incidence data in Stan to estimate Qdenga's efficacy.
* *impact* - code to sample from the posterior distributions obtained in *efficacy* and simulate impact using a stochastic compartmental dengue transmission model. 
