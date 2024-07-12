# Qdenga impact


# About this repository
This repository contains all the code needed to reproduce the results presented in Efficacy, public health impact and optimal use of the Takeda dengue vaccine, in full. All data used is made available in the repository.

All required packages and their versions are available in the DESCRIPTION file. Instructions to install and download cmdStan can be found [here](https://mc-stan.org/users/interfaces/cmdstan).

# Repository structure

This repository is divided into two sections: 

* *efficacy* - code to fit the cohort Bayesian survival model to published clinical incidence data in Stan to estimate Qdenga's efficacy.
* *impact* - code to sample from the posterior distributions obtained in *efficacy* and simulate impact using a stochastic compartmental dengue transmission model. 
