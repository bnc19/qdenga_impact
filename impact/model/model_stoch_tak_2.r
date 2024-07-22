## Neil's dengue model transcribed from BM to odin
##

  PI_C <- 3.14159265359 ## PI
  DT  <-  1.0  ## timestep in days - keep at 1
  
  TIME <- step*DT
  
  YL <- user() ## year length in days
  
  YEAR <- TIME/YL ## current model time in years
  
  NYO <- 20 ## number of year output in post-vacc output vectors
  NUM_YEAR_ACCUM <- 4 ## time to accumulate incidence for a couple of output variables
  
  age_per <- YL ## how frequently aging occurs - needs to be year length for precise match to demography
  # N_age <- 28 ## number of age classes
  # N_age_p1 <- N_age + 1
  # agec[1:20] <- 1 ## first 20 are 1 year wide
  # agec[21:28] <- 10 ## then next 8 are 10 years wide
  N_age <- 40 ## number of age classes
  N_age_p1 <- N_age + 1
  agec[1:33] <- 1 ## first 33 are 1 year wide
  agec[34] <-7 ### then 7 years
  agec[35:40] <- 10 ## then next 8 are 10 years wide
  ## lower boundary of each age class. Extra entry gives upper boundary of last age class 
  ageb[1] <- 0
  ageb[2:(N_age_p1)] <- ageb[i-1]+agec[i-1]
  
  N_sim <- user() # simulate 5 million people (initial pop size)
  
  age_removal_d[,] <- user()
  age_rate_d[,] <- user()
  pop_size_d[] <- user()
  births_d[,] <- user()
  life_expec_d[,] <- user()
  coverage_d[,] <- user()
  
  CALIB_YEAR <- user() ## calendar year used for R0 calibration against demography (for age-dependent exposure)
  VACC_YEAR <- user() ## calendar year vaccination starts
  VACC_END_YEAR <- user() ## calendar year vaccination ends
  FIRST_YEAR <- user() ## first calendar year in demography data
  LAST_YEAR <- user() ## last calendar year in demography data
  max_rel_year <- LAST_YEAR-FIRST_YEAR ## number of years covered by demography data
  
  EQUILIB_YEARS <- user() ## years to run model before FIRST_YEAR reached
  
  YEAR_STEP <- 1 ## time between population size estimates given in demography files
  
  VaccOn <- user() ## 1 if vaccination, 0 if not
  VaccRoutineAge <- user()
  ZeroVE <- user() ## 1 for zero VE, 1 for standard VE model
  DoVEInf <- user() ## 1 for VE against inf, 0 for against disease
  
  vca <- VaccRoutineAge ## age at routine vacc (<20, assuming current age structure)
  
  vacc_child_age <- if(vca == 0) N_age_p1 else vca ## VaccRoutineAge==0 means no routine vacc
  ## simulation year that vaccination starts (i.e. when calendar year reaches VACC_YEAR)
  vacc_child_starttime <- EQUILIB_YEARS+VACC_YEAR-FIRST_YEAR
  ## simulation year that vaccination stops (i.e. when calendar year reaches VACC_END_YEAR)
  vacc_child_stoptime <- vacc_child_starttime+VACC_END_YEAR-VACC_YEAR
  ## how many years have passed since vaccination started
  YEARS_POST_VACC <- YEAR-vacc_child_starttime

  ## index to coverage table row
  cov_year <- floor(YEARS_POST_VACC)+1
  ## current coverage (routine only)
  gavi_cov <- if( VaccOn == 0 || YEARS_POST_VACC<0 || YEARS_POST_VACC>VACC_END_YEAR-VACC_YEAR) 0 else coverage_d[as.integer(cov_year),2]
  
  ## new vaccine efficacy model params
  nc0[,] <- user() # Ab at t=0 (by serotype and serostatus) - includes 1/(exp(-ts[j]*hs[j])+exp(-ts[j]*hl[j])) factor
  hs[] <- user() # Ab fast decay rate (by serostatus) - rate not half life
  hl[] <- user() # Ab slow decay rate (by serostatus)
  ts[] <- user() # Ab decay fast-slow transition time (by serostatus)
  nc50_inf[,] <- user() # Ab n50 for dis (by serotype and serostatus)  
  nc50_dis[,] <- user() # Ab n50 for dis (by serotype and serostatus)
  nc50_sdis[,] <- user() # Ab n50 for sdis (by serotype and serostatus)
  L_dis[,] <- user() # enhancement for symp dis (by serotype and serostatus)
  L_sdis[,] <- user() # enhancement for severe dis or hosp (by serotype and serostatus)
  ws[] <- user() # Ab dose response slope param (by serotype)
  nc50_age_under6 <- user()
  nc50_age[1:N_age] <- if(i<6) exp(nc50_age_under6) else 1
  max_ve_decay_time <- user() # how long (days) to let VE decline
  
  ## facility exists to have different coverage by serostatus (i.e. to model testing)
  ## variable controlling whether vaccination can vary by serostatus (0=no, 1=yes)
  vacc_by_serostatus <- user()
  sero_test_sens <-0.896 ## WHO value
  sero_test_spec <- 0.947 ## WHO vaue
  vacc_child_coverage <- if(vacc_by_serostatus==1) gavi_cov*sero_test_sens else gavi_cov
  vacc_child_coverage_S <- if(vacc_by_serostatus==1) gavi_cov*(1.0-sero_test_spec) else gavi_cov
  
  ## catch-up age range. If over N_age then no catch-up
  vacc_cu_minage <- N_age_p1
  vacc_cu_maxage <- N_age_p1
  ## catch-up coverage, again by serostatus
  vacc_cu_coverage <- 0.0
  vacc_cu_coverage_S <- 0.0
  ## year of catch-up campaign - normally on the same year that routine starts
  vacc_cu_time <- vacc_child_starttime
  ## rounded version of vacc_cu_time so that it happens during a time step when aging happens
  ## note that in current implementation, all catch-up doses are given in a single time step
  ## +1 in line below ensures catch-up happens day after aging, not the same day
  vacc_cu_rndtime <- floor(vacc_cu_time*YL/age_per)*age_per+1
  
  #### dynamic demog
  
  ## current model year relative to FIRST_YEAR
  rel_year <- YEAR - EQUILIB_YEARS  
  
  ## row in age_removal demog file corresponding to current year
  year_row <- 1+(if(rel_year<0) 0 else if(rel_year>max_rel_year) max_rel_year else floor(rel_year))
  
  ## year used to calibrate R0 (for age-dependent exposure)
  year_calib <- 1+CALIB_YEAR-FIRST_YEAR
  
  death[1:N_age] <- age_removal_d[as.integer(year_row),1+i]
  ## correction for first year to allow for continuous births
  death[1] <- 0.4691*death[1]*death[1]+1.9686*death[1]
  ## change time units to days
  deathrt[1:N_age] <- death[i]/YL*DT
  
  ## not really the mean age of people in each age class, but biassed below midpoint to allow for deaths
  mean_age[1:N_age] <- (0.75*ageb[i]+0.25*ageb[i+1])
  
  ## proportion of each age class which age to the next each year
  cur_age_rate[1:N_age] <- age_rate_d[as.integer(year_row),1+i]
  ## age rate is only non-zero once per age_per
  agerts <- if(floor(TIME/age_per) == TIME/age_per) age_per/YL else 0
  ## agerts <- DT/YL
  agert[1:N_age] <- cur_age_rate[i]*agerts  ## /agec[i]
  
  ## initial value used for initialising human population 
  N_init_age0[1:N_age] <- pop_size_d[i+1]
  N_init <- sum(N_init_age0)
  pop_scale <- if(N_sim<=0) 1.0 else N_sim/N_init
  N_init_age[1:N_age] <- pop_size_d[i+1]*pop_scale

  ## currently used for initialising mosquito model
  N_eq <- sum(N_init_age)
  
  ## Life expectancy for each age class - used to calculate YLLs
  life_expec[1:N_age] <- life_expec_d[as.integer(year_row),1+i]
  ## life expectancy used for R0 calibration
  init_lifespan <- life_expec_d[as.integer(year_calib),2]
  init_life_expec[1:N_age] <- life_expec_d[as.integer(year_calib),1+i]
  
  ## births happen continuously
  Nb <- births_d[as.integer(year_row),2]
  births <- pop_scale*Nb*DT/YL
  
  ## disease severity
  cfr <- user()  ## disease related death not modelled explicitly - this only used for YLL calculation
  daly_dis <- user()  ## DALYs associated with symptomatic disease
  daly_sdis <- user() ## DALYs associated with severe disease (on top of mild disease)
  
  ## economics
  disc_fact <- user() ## discounting factor (p.a.) used for DALYs
  disc_fact_vacc <- user() ## discounting for vaccine doses
  cum_disc <- if(YEARS_POST_VACC<0) 0 else 1.0/(disc_fact^YEARS_POST_VACC)
  cum_disc_vacc <- if(YEARS_POST_VACC<0) 0 else 1.0/(disc_fact_vacc^YEARS_POST_VACC)
  

  ## proportion of infections causing symptomatic disease in unvaccinated
  
  dis_pri[] <- user() ## prop of primary infectious which are symptomatic
  dis_sec[] <- user()
  dis_tert[] <- user()
  dis_quart[] <- user()
  
  ## proportion of infections causing severe disease
  sdis_pri[] <- user() ## prop of primary infectious which are severe/hosp
  sdis_sec[] <- user()
  sdis_tert[] <- user()
  sdis_quart[] <- user()
  
  
  ## mosquito params
  
  omega <- user()  ## density dep power
  delta <- user()  ## adult mos death rate
  epsilon <- user()  ## larval maturation
  Rm <- user() ##  mosq rep no.
  ##  Mwt = adult wild type female density per person
  ##  this param can also be changed to give required R0 rather than beta_hm
  Mwt <- user() 
  sigma <- user() ##  low density limit of larval death rate
  ## fecundity - computed to give required Rm
  gamma <- Rm*delta*(epsilon+sigma)/epsilon
  ## carrying capacity - computed to give required Mwt
  Kc_mean <- Mwt*delta*( (epsilon*(gamma-delta)/(delta*sigma)-1) ^(-1/omega) )/epsilon
  Kc_season <- user() ## seasomal variation in carrying capacity
  Kc <- Kc_mean*(1+Kc_season*cos(2*PI_C*TIME/YL))
  kappa <- user()  ## biting rate/day
  
  Beta_mh_mean <- user() ## mosquito to human transmission prob (assumed 1)
  Beta_mh_season <- user() ## seasonal variation in transmission prob
  Beta_mh <- Beta_mh_mean*(1+Beta_mh_season*cos(2*PI_C*TIME/YL))
  ## external force of infection on people
  extInf <- user()
  
  ## Facility to vary human exposure by age
  Acrit <- 3
  FOIagescale <- 1.0
  FOIas[1:N_age] <- if(i>Acrit) FOIagescale else 1
  suscinitpop[1:N_age] <- FOIas[i]*init_life_expec[i]
  R0agescale <- sum(init_life_expec)/sum(suscinitpop)
  
  R0_req <- user() ## required R0 (not precise)
  Rel_R01 <- user() ## R0 scaling for DENV1
  Rel_R02 <- user() ## R0 scaling for DENV2
  Rel_R03 <- user() ## R0 scaling for DENV3
  Rel_R04 <- user() ## R0 scaling for DENV4
  
  Beta_hm <- R0agescale*R0_req/(kappa*kappa*Mwt*inf_per*Beta_mh_mean/(1+delta*eip)/delta)
  Beta_hm_1 <- Beta_hm*Rel_R01
  Beta_hm_2 <- Beta_hm*Rel_R02
  Beta_hm_3 <- Beta_hm*Rel_R03
  Beta_hm_4 <- Beta_hm*Rel_R04
  
  ## Wolbachia params - move to R later
  Wb_cyto <- 1
  Wb_mat <- 0.95
  Wb_fM <- 0.95
  Wb_fF <- 0.85
  Wb_eff <- 0.6
  Wb_relsusc1 <- 1-Wb_eff
  Wb_relsusc2 <- 1-Wb_eff
  Wb_relsusc3 <- 1-Wb_eff
  Wb_relsusc4 <- 1-Wb_eff
  Wb_relinf1 <- 1
  Wb_relinf2 <- 1
  Wb_relinf3 <- 1
  Wb_relinf4 <- 1
  Wb_introtime <- 10000  ## years
  Wb_introlevel <- 0.0
  Wb_introduration <- 60
  Wb_introrate <- Wb_introlevel*Mwt*N_eq*DT/Wb_introduration
  delta_wb <- delta/Wb_fM
  
  incub <- user() ## incubation period in humans
  eip <- user() ## EIP
  inf_per <- user()  ## duration of human infectiousness
  dur_cross_prot <- user() ## duration of human cross-protection after infection
  nu <- DT/dur_cross_prot ##  waning rate of cross-protection 
  
  R0_1 <- kappa*kappa*Mwt*Beta_hm_1*inf_per*Beta_mh_mean/(1+delta*eip)/delta/R0agescale ##  eqn correct for exp distrib eip
  R0_2 <- kappa*kappa*Mwt*Beta_hm_2*inf_per*Beta_mh_mean/(1+delta*eip)/delta/R0agescale
  R0_3 <- kappa*kappa*Mwt*Beta_hm_3*inf_per*Beta_mh_mean/(1+delta*eip)/delta/R0agescale
  R0_4 <- kappa*kappa*Mwt*Beta_hm_4*inf_per*Beta_mh_mean/(1+delta*eip)/delta/R0agescale
  eq_FOI1 <- R0_1/init_lifespan ## in time units of years
  eq_FOI2 <- R0_2/init_lifespan
  eq_FOI3 <- R0_3/init_lifespan
  eq_FOI4 <- R0_4/init_lifespan
  
  ## catch-up can be weighted by age in theory
  vacc_cu_age_weight[1:N_age_p1] <- 1
  
  ## vacc_cu_age_weight[0:1] <- 0
  ## vacc_cu_age_weight[2:5] <- 1   ## 0.86
  ## vacc_cu_age_weight[6:11] <- 1  ## 0.59
  ## vacc_cu_age_weight[12:14] <- 1
  ## vacc_cu_age_weight[15:N_age] <- 0
  
  
  ## modelling vaccination as part of aging, with arbitrary waning
  
  ## proportion NOT vaccinated as they age out of an age class - seropositives first
  ## j=1 corresponds to unvaccinated states
  vacc_noncov[1:N_age_p1,1] <- (if((YEARS_POST_VACC>=0) && (YEAR<vacc_child_stoptime) && (i-1 == vacc_child_age)) (1-vacc_child_coverage) else 1)
  ## j = 2,3 are vaccinated groups, so vacc_noncov==1 (already vaccinated)
  vacc_noncov[1:N_age_p1,2:3] <- 1
  ## now campaigns - happen on a single day
  vcu_noncov[1:N_age,1] <- (if((TIME == vacc_cu_rndtime) && (i >= vacc_cu_minage) && (i <= vacc_cu_maxage)) (1-vacc_cu_coverage*vacc_cu_age_weight[i]) else 1)
  vcu_noncov[1:N_age,2:3] <- 1
  
  ## now seronegatives
  vacc_noncov_S[1:N_age_p1,1] <- (if((YEARS_POST_VACC>=0) && (YEAR<vacc_child_stoptime) && (i-1 == vacc_child_age)) (1-vacc_child_coverage_S) else 1)
  vacc_noncov_S[1:N_age_p1,2:3] <- 1
  vcu_noncov_S[1:N_age,1] <- (if((TIME == vacc_cu_rndtime) && (i >= vacc_cu_minage) && (i <= vacc_cu_maxage)) (1-vacc_cu_coverage_S)*vacc_cu_age_weight[i] else 1)
  vcu_noncov_S[1:N_age,2:3] <- 1
  
  
  ## arbitrary time waning of VE against infection
  ## model would need to be extended to let VE against disease vary similarly
  ## ve_cu_time and ve_child_time are 1 immediately after vaccination, but are multiplied by ve_pri_param etc

  ## complicated if clause represents 2 doses at 0 and 3 months (for Takeda)
  cu_time_from_last_dose_base <- if(TIME<vacc_cu_rndtime+YL/4) (TIME-vacc_cu_rndtime) else if(TIME>=vacc_cu_rndtime+YL/4) (TIME-vacc_cu_rndtime-YL/4) else 0
  cu_time_from_last_dose <- if(cu_time_from_last_dose_base>max_ve_decay_time) max_ve_decay_time else cu_time_from_last_dose_base
  
  youngest_cu_age <- if(vacc_cu_minage<=N_age) ageb[as.integer(vacc_cu_minage+1)]+floor(YEAR-vacc_cu_time) else 1000
  oldest_cu_age <- if(vacc_cu_minage<=N_age) ageb[as.integer(vacc_cu_maxage+1)]+floor(YEAR-vacc_cu_time) else 1000
  
  # vacc_ct[1:N_age] <- if((i<as.integer(vca)+1) || (YEARS_POST_VACC<ageb[i]-ageb[as.integer(vca)+1])) 1e8 else floor(YEAR)-ageb[i]+ageb[as.integer(vca)+1]
  vacc_ct[1:N_age] <- if(i<as.integer(vca)+1) 1e8 else if(i<=20) floor(YEAR)-ageb[i]+ageb[as.integer(vca)+1] else floor(YEAR)-30+ageb[as.integer(vca)+1]
  
  
  ## complicated if clause represents 2 doses at 0 and 3 months (for Takeda)
  time_from_last_dose_base[1:N_age] <- if((YEAR>=vacc_ct[i]) && (YEAR<vacc_ct[i]+0.25)) YL*(YEAR-vacc_ct[i]) else if(YEAR>=vacc_ct[i]+0.25) YL*(YEAR-vacc_ct[i]-0.25) else 0
  time_from_last_dose[1:N_age] <- if(time_from_last_dose_base[i]>max_ve_decay_time) max_ve_decay_time else time_from_last_dose_base[i]
    
  initial(out_time_from_last_dose[1:N_age]) <- 0
  update(out_time_from_last_dose[1:N_age]) <- time_from_last_dose[i]
  dim(out_time_from_last_dose) <- N_age

  ### time varying VE (RR)
  
  ####  following only works if cu campaign targets children older than routine campaign
  

  nc[1:4,1:3,1:N_age] <- nc0[i,j] * (exp(-time_from_last_dose[k]*hs[j]-ts[j]*hl[j])+exp(-time_from_last_dose[k]*hl[j]-ts[j]*hs[j]))
  nc_cu[1:4,1:3] <- nc0[i,j] *  (exp(-cu_time_from_last_dose*hs[j]-ts[j]*hl[j])+exp(-cu_time_from_last_dose*hl[j]-ts[j]*hs[j]))
  
  initial(out_nc[1:4,1:3,1:N_age]) <- 0
  update(out_nc[1:4,1:3,1:N_age]) <- nc[i,j,k]
  dim(out_nc) <- c(4,3,N_age)
  
  ## relative risks of infection, disease and severe disease in vaccinated people
  ## note that RR_dis = RR(dis|inf)
  ## If DoVEInf=1, assume that VE(inf) = sqrt of total protection, with all enhancement in VE(dis|inf)
  ## This is approximately consistent with findings of https://pubmed.ncbi.nlm.nih.gov/27418050/
  
  RR_inf_vc[1:4,1,1:N_age] <- if((ZeroVE==0) && (DoVEInf)) 1/(1+(nc[i,j,k]/(nc50_inf[i,j]*nc50_age[k]))^ws[i]) else 1
  RR_inf_cu[1:4,1,1:N_age] <- if((ZeroVE==0) && (DoVEInf) && (ageb[k]>=youngest_cu_age) && (ageb[k]<=oldest_cu_age)) 1/(1+(nc_cu[i,j]/(nc50_inf[i,j]*nc50_age[k]))^ws[i]) else 1
  
  RR_inf_vc[1:4,2:3,1:N_age] <- if((ZeroVE==0) && (DoVEInf)) 1/(1+(nc[i,j,k]/(nc50_inf[i,j]*nc50_age[k]))^ws[i]) else 1
  RR_inf_cu[1:4,2:3,1:N_age] <- if((ZeroVE==0) && (DoVEInf) && (ageb[k]>=youngest_cu_age) && (ageb[k]<=oldest_cu_age)) 1/(1+(nc_cu[i,j]/(nc50_inf[i,j]*nc50_age[k]))^ws[i]) else 1
  
  RR_dis_vc[1:4,1:3,1:N_age] <- if(ZeroVE==1) 1 else (L_dis[i,j]/(1+(nc[i,j,k]/(nc50_dis[i,j]*nc50_age[k]))^ws[i]))/RR_inf_vc[i,j,k] 
  RR_sdis_vc[1:4,1:3,1:N_age] <- if(ZeroVE==1) 1 else L_sdis[i,j]/(1+(nc[i,j,k]/(nc50_sdis[i,j]*nc50_age[k]))^ws[i])/RR_inf_vc[i,j,k]
  
  RR_dis_cu[1:4,1:3,1:N_age] <- if((ZeroVE==1) || (ageb[k]<youngest_cu_age) || (ageb[k]>oldest_cu_age)) 1 else L_dis[i,j]/(1+(nc_cu[i,j]/(nc50_dis[i,j]*nc50_age[k]))^ws[i])/RR_inf_cu[i,j,k]
  RR_sdis_cu[1:4,1:3,1:N_age] <- if((ZeroVE==1) || (ageb[k]<youngest_cu_age) || (ageb[k]>oldest_cu_age)) 1 else  L_sdis[i,j]/(1+(nc_cu[i,j]/(nc50_sdis[i,j]*nc50_age[k]))^ws[i])/RR_inf_cu[i,j,k]
  
  RR_inf[1:4,1:3,1:N_age] <- RR_inf_vc[i,j,k]*RR_inf_cu[i,j,k]
  RR_dis[1:4,1:3,1:N_age] <- RR_dis_vc[i,j,k]*RR_dis_cu[i,j,k] 
  RR_sdis[1:4,1:3,1:N_age] <- RR_sdis_vc[i,j,k]*RR_sdis_cu[i,j,k]
  
  initial(out_RR[1:4,1:3,1:N_age]) <- 0
  update(out_RR[1:4,1:3,1:N_age]) <- RR_dis[i,j,k]
  dim(out_RR) <- c(4,3,N_age)
  
  
  ## rho params are susceptibility to infection
  ## assume no immunity beyond initial R state following each infection
  
  rho_pri <- 1
  rho_sec <- 1
  rho_tert <- 1
  rho_quart <- 1
  
  rho1[1:N_age,1] <- rho_pri
  rho2[1:N_age,1] <- rho_pri
  rho3[1:N_age,1] <- rho_pri
  rho4[1:N_age,1] <- rho_pri
  rho21[1:N_age,1] <- rho_sec
  rho31[1:N_age,1] <- rho_sec
  rho41[1:N_age,1] <- rho_sec
  rho12[1:N_age,1] <- rho_sec
  rho32[1:N_age,1] <- rho_sec
  rho42[1:N_age,1] <- rho_sec
  rho13[1:N_age,1] <- rho_sec
  rho23[1:N_age,1] <- rho_sec
  rho43[1:N_age,1] <- rho_sec
  rho14[1:N_age,1] <- rho_sec
  rho24[1:N_age,1] <- rho_sec
  rho34[1:N_age,1] <- rho_sec
  rho231[1:N_age,1] <- rho_tert
  rho241[1:N_age,1] <- rho_tert
  rho341[1:N_age,1] <- rho_tert
  rho132[1:N_age,1] <- rho_tert
  rho142[1:N_age,1] <- rho_tert
  rho342[1:N_age,1] <- rho_tert
  rho123[1:N_age,1] <- rho_tert
  rho143[1:N_age,1] <- rho_tert
  rho243[1:N_age,1] <- rho_tert
  rho124[1:N_age,1] <- rho_tert
  rho134[1:N_age,1] <- rho_tert
  rho234[1:N_age,1] <- rho_tert
  rho2341[1:N_age,1] <- rho_quart
  rho1342[1:N_age,1] <- rho_quart
  rho1243[1:N_age,1] <- rho_quart
  rho1234[1:N_age,1] <- rho_quart

  rho1[1:N_age,2:3] <- rho_pri*RR_inf[1,1,i]
  rho2[1:N_age,2:3] <- rho_pri*RR_inf[2,1,i]
  rho3[1:N_age,2:3] <- rho_pri*RR_inf[3,1,i]
  rho4[1:N_age,2:3] <- rho_pri*RR_inf[4,1,i]
  rho21[1:N_age,2:3] <- rho_sec*RR_inf[1,2,i]
  rho31[1:N_age,2:3] <- rho_sec*RR_inf[1,2,i]
  rho41[1:N_age,2:3] <- rho_sec*RR_inf[1,2,i]
  rho12[1:N_age,2:3] <- rho_sec*RR_inf[2,2,i]
  rho32[1:N_age,2:3] <- rho_sec*RR_inf[2,2,i]
  rho42[1:N_age,2:3] <- rho_sec*RR_inf[2,2,i]
  rho13[1:N_age,2:3] <- rho_sec*RR_inf[3,2,i]
  rho23[1:N_age,2:3] <- rho_sec*RR_inf[3,2,i]
  rho43[1:N_age,2:3] <- rho_sec*RR_inf[3,2,i]
  rho14[1:N_age,2:3] <- rho_sec*RR_inf[4,2,i]
  rho24[1:N_age,2:3] <- rho_sec*RR_inf[4,2,i]
  rho34[1:N_age,2:3] <- rho_sec*RR_inf[4,2,i]
  rho231[1:N_age,2:3] <- rho_tert*RR_inf[1,3,i]
  rho241[1:N_age,2:3] <- rho_tert*RR_inf[1,3,i]
  rho341[1:N_age,2:3] <- rho_tert*RR_inf[1,3,i]
  rho132[1:N_age,2:3] <- rho_tert*RR_inf[2,3,i]
  rho142[1:N_age,2:3] <- rho_tert*RR_inf[2,3,i]
  rho342[1:N_age,2:3] <- rho_tert*RR_inf[2,3,i]
  rho123[1:N_age,2:3] <- rho_tert*RR_inf[3,3,i]
  rho143[1:N_age,2:3] <- rho_tert*RR_inf[3,3,i]
  rho243[1:N_age,2:3] <- rho_tert*RR_inf[3,3,i]
  rho124[1:N_age,2:3] <- rho_tert*RR_inf[4,3,i]
  rho134[1:N_age,2:3] <- rho_tert*RR_inf[4,3,i]
  rho234[1:N_age,2:3] <- rho_tert*RR_inf[4,3,i]
  rho2341[1:N_age,2:3] <- rho_quart*RR_inf[1,3,i]
  rho1342[1:N_age,2:3] <- rho_quart*RR_inf[2,3,i]
  rho1243[1:N_age,2:3] <- rho_quart*RR_inf[3,3,i]
  rho1234[1:N_age,2:3] <- rho_quart*RR_inf[4,3,i]
  

  ## phi params relate to infectiousness
  
  ##  infectiousness of symptomatic cases relative to asymp
  phi_dis_enhance <- user()
  phi_ed <- phi_dis_enhance-1
  phi_scale[1:4] <- 1/(1+dis_pri[i]*phi_ed)  ##  to keep R0 definition valid
  phi_pri[1:4] <- phi_scale[i]
  phi_sec[1:4] <- phi_scale[i]
  phi_tert[1:4] <- phi_scale[i]
  phi_quart[1:4] <- phi_scale[i]
  

  ## p9p <- 1-((1-p9)^2)
  ## phi_scale <- (1+0.45*phi_ed+p9p*(2+0.95*phi_ed))/(1+dis_pri*phi_ed+p9p*(2+(dis_sec+dis_tert)*phi_ed))/(1+0.45*phi_ed)  ##  to keep p9 approx accurate
  
 
  
  ################ Main model initialisation ##############
  
  ## Mosquito model first
  ## mosquito model is deterministic currently
  
  initial(Lwt) <- Mwt*N_eq*delta/epsilon
  initial(Mwt_S) <- Mwt*N_eq*0.9997
  initial(Mwt_E1) <- 0
  initial(Mwt_E2) <- 0
  initial(Mwt_E3) <- 0
  initial(Mwt_E4) <- 0
  initial(Mwt_I1) <- floor(Mwt*N_eq*0.00005)
  initial(Mwt_I2) <- floor(Mwt*N_eq*0.00006)
  initial(Mwt_I3) <- floor(Mwt*N_eq*0.00007)
  initial(Mwt_I4) <- floor(Mwt*N_eq*0.00008)
  initial(Lwb) <- 0
  initial(Mwb_S) <- 0
  initial(Mwb_E1) <- 0
  initial(Mwb_E2) <- 0
  initial(Mwb_E3) <- 0
  initial(Mwb_E4) <- 0
  initial(Mwb_I1) <- 0
  initial(Mwb_I2) <- 0
  initial(Mwb_I3) <- 0
  initial(Mwb_I4) <- 0
  
  Mwt_tot <- Mwt_S+Mwt_E1+Mwt_E2+Mwt_E3+Mwt_E4+Mwt_I1+Mwt_I2+Mwt_I3+Mwt_I4
  Mwb_tot <- Mwb_S+Mwb_E1+Mwb_E2+Mwb_E3+Mwb_E4+Mwb_I1+Mwb_I2+Mwb_I3+Mwb_I4
  R0t_1 <- kappa*kappa*(Mwt_tot+Wb_relsusc1*Wb_relinf1*Mwb_tot)*Beta_hm_1*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  R0t_2 <- kappa*kappa*(Mwt_tot+Wb_relsusc2*Wb_relinf2*Mwb_tot)*Beta_hm_2*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  R0t_3 <- kappa*kappa*(Mwt_tot+Wb_relsusc3*Wb_relinf3*Mwb_tot)*Beta_hm_3*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  R0t_4 <- kappa*kappa*(Mwt_tot+Wb_relsusc4*Wb_relinf4*Mwb_tot)*Beta_hm_4*inf_per*Beta_mh/(1+delta*eip)/delta/NT/R0agescale
  
  Lwt_birth <- (DT*gamma*(Mwt_tot*(Mwt_tot+(1-Wb_cyto)*Mwb_tot)/(Mwt_tot+Mwb_tot)+Wb_fF*(1-Wb_mat)*Mwb_tot))
  Lwb_birth <- (DT*gamma*Wb_fF*Wb_mat*Mwb_tot)
  L_deathrt <- DT*sigma*((1+((Lwt+Lwb)/(Kc*NT))^omega))
  O_Lwt <- (DT*epsilon+L_deathrt)*(Lwt)
  Lwt_mature <- (DT*epsilon/(DT*epsilon+L_deathrt))*(O_Lwt)
  update(Lwt) <- Lwt_birth+Lwt-O_Lwt
  O_Lwb <- (DT*epsilon+L_deathrt)*(Lwb)
  Lwb_mature <- (DT*epsilon/(DT*epsilon+L_deathrt))*(O_Lwb)
  update(Lwb) <- Lwb_birth+Lwb-O_Lwb
  
  Mwt_FOI1 <- DT*Beta_hm_1*kappa*infectious1
  Mwt_FOI2 <- DT*Beta_hm_2*kappa*infectious2
  Mwt_FOI3 <- DT*Beta_hm_3*kappa*infectious3
  Mwt_FOI4 <- DT*Beta_hm_4*kappa*infectious4
  O_Mwt_S <- (DT*delta+Mwt_FOI1+Mwt_FOI2+Mwt_FOI3+Mwt_FOI4)*(Mwt_S)
  Mwt_inf1 <- (Mwt_FOI1/(DT*delta+Mwt_FOI1+Mwt_FOI2+Mwt_FOI3+Mwt_FOI4))*(O_Mwt_S)
  Mwt_inf2 <- (Mwt_FOI2/(DT*delta+Mwt_FOI2+Mwt_FOI3+Mwt_FOI4))*(O_Mwt_S-Mwt_inf1)
  Mwt_inf3 <- (Mwt_FOI3/(DT*delta+Mwt_FOI3+Mwt_FOI4))*(O_Mwt_S-Mwt_inf1-Mwt_inf2)
  Mwt_inf4 <- (Mwt_FOI4/(DT*delta+Mwt_FOI4))*(O_Mwt_S-Mwt_inf1-Mwt_inf2-Mwt_inf3)
  update(Mwt_S) <- Lwt_mature+Mwt_S-O_Mwt_S
  O_Mwt_E1 <- (DT*(delta+1/eip))*(Mwt_E1)
  Mwt_E1_incub <- (1/(delta*eip+1))*(O_Mwt_E1)
  update(Mwt_E1) <- Mwt_inf1+Mwt_E1-O_Mwt_E1
  O_Mwt_E2 <- (DT*(delta+1/eip))*(Mwt_E2)
  Mwt_E2_incub <- (1/(delta*eip+1))*(O_Mwt_E2)
  update(Mwt_E2) <- Mwt_inf2+Mwt_E2-O_Mwt_E2
  O_Mwt_E3 <- (DT*(delta+1/eip))*(Mwt_E3)
  Mwt_E3_incub <- (1/(delta*eip+1))*(O_Mwt_E3)
  update(Mwt_E3) <- Mwt_inf3+Mwt_E3-O_Mwt_E3
  O_Mwt_E4 <- (DT*(delta+1/eip))*(Mwt_E4)
  Mwt_E4_incub <- (1/(delta*eip+1))*(O_Mwt_E4)
  update(Mwt_E4) <- Mwt_inf4+Mwt_E4-O_Mwt_E4
  update(Mwt_I1) <- Mwt_E1_incub+Mwt_I1-(DT*delta)*(Mwt_I1)
  update(Mwt_I2) <- Mwt_E2_incub+Mwt_I2-(DT*delta)*(Mwt_I2)
  update(Mwt_I3) <- Mwt_E3_incub+Mwt_I3-(DT*delta)*(Mwt_I3)
  update(Mwt_I4) <- Mwt_E4_incub+Mwt_I4-(DT*delta)*(Mwt_I4)
  
  Mwb_FOI1 <- Wb_relsusc1*Mwt_FOI1
  Mwb_FOI2 <- Wb_relsusc2*Mwt_FOI2
  Mwb_FOI3 <- Wb_relsusc3*Mwt_FOI3
  Mwb_FOI4 <- Wb_relsusc4*Mwt_FOI4
  O_Mwb_S <- (DT*delta_wb+Mwb_FOI1+Mwb_FOI2+Mwb_FOI3+Mwb_FOI4)*(Mwb_S)
  Mwb_inf1 <- (Mwb_FOI1/(DT*delta_wb+Mwb_FOI1+Mwb_FOI2+Mwb_FOI3+Mwb_FOI4))*(O_Mwb_S)
  Mwb_inf2 <- (Mwb_FOI2/(DT*delta_wb+Mwb_FOI2+Mwb_FOI3+Mwb_FOI4))*(O_Mwb_S-Mwb_inf1)
  Mwb_inf3 <- (Mwb_FOI3/(DT*delta_wb+Mwb_FOI3+Mwb_FOI4))*(O_Mwb_S-Mwb_inf1-Mwb_inf2)
  Mwb_inf4 <- (Mwb_FOI4/(DT*delta_wb+Mwb_FOI4))*(O_Mwb_S-Mwb_inf1-Mwb_inf2-Mwb_inf3)
  Mwb_intro <- if((TIME>=Wb_introtime*YL) && (TIME<Wb_introtime*YL+Wb_introduration)) (Wb_introrate) else 0
  update(Mwb_S) <- Lwb_mature+Mwb_intro+Mwb_S-O_Mwb_S
  O_Mwb_E1 <- (DT*(delta_wb+1/eip))*(Mwb_E1)
  Mwb_E1_incub <- (1/(delta*eip+1))*(O_Mwb_E1)
  update(Mwb_E1) <- Mwb_inf1+Mwb_E1-O_Mwb_E1
  O_Mwb_E2 <- (DT*(delta_wb+1/eip))*(Mwb_E2)
  Mwb_E2_incub <- (1/(delta_wb*eip+1))*(O_Mwb_E2)
  update(Mwb_E2) <- Mwb_inf2+Mwb_E2-O_Mwb_E2
  O_Mwb_E3 <- (DT*(delta_wb+1/eip))*(Mwb_E3)
  Mwb_E3_incub <- (1/(delta_wb*eip+1))*(O_Mwb_E3)
  update(Mwb_E3) <- Mwb_inf3+Mwb_E3-O_Mwb_E3
  O_Mwb_E4 <- (DT*(delta_wb+1/eip))*(Mwb_E4)
  Mwb_E4_incub <- (1/(delta_wb*eip+1))*(O_Mwb_E4)
  update(Mwb_E4) <- Mwb_inf4+Mwb_E4-O_Mwb_E4
  update(Mwb_I1) <- Mwb_E1_incub+Mwb_I1-(DT*delta_wb)*(Mwb_I1)
  update(Mwb_I2) <- Mwb_E2_incub+Mwb_I2-(DT*delta_wb)*(Mwb_I2)
  update(Mwb_I3) <- Mwb_E3_incub+Mwb_I3-(DT*delta_wb)*(Mwb_I3)
  update(Mwb_I4) <- Mwb_E4_incub+Mwb_I4-(DT*delta_wb)*(Mwb_I4)
  
  #### Human model ###
  
  ## state variables are matrices of dimensions N_age x 3
  ## var[,1] = unvaccinated
  ## var[,2] = vaccinated when seropositive 
  ## var[,3] = vaccinated when seronegative

  initial(S[1:N_age,1]) <- floor(0.5+N_init_age[i]*exp(-(eq_FOI1+eq_FOI2+eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(S[1:N_age,2:3]) <- 0
  initial(I1[1:N_age,1:3]) <- 0
  initial(I2[1:N_age,1:3]) <- 0
  initial(I3[1:N_age,1:3]) <- 0
  initial(I4[1:N_age,1:3]) <- 0
  
  initial(R1[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*exp(-(eq_FOI2+eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(R2[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*exp(-(eq_FOI1+eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(R3[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI3*mean_age[i]))*exp(-(eq_FOI1+eq_FOI2+eq_FOI4)*mean_age[i]))
  initial(R4[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI1+eq_FOI2+eq_FOI3)*mean_age[i]))
  initial(R1[1:N_age,2:3]) <- 0
  initial(R2[1:N_age,2:3]) <- 0
  initial(R3[1:N_age,2:3]) <- 0
  initial(R4[1:N_age,2:3]) <- 0
  
  initial(I21[1:N_age,1:3]) <- 0
  initial(I31[1:N_age,1:3]) <- 0
  initial(I41[1:N_age,1:3]) <- 0
  initial(I12[1:N_age,1:3]) <- 0
  initial(I32[1:N_age,1:3]) <- 0
  initial(I42[1:N_age,1:3]) <- 0
  initial(I13[1:N_age,1:3]) <- 0
  initial(I23[1:N_age,1:3]) <- 0
  initial(I43[1:N_age,1:3]) <- 0
  initial(I14[1:N_age,1:3]) <- 0
  initial(I24[1:N_age,1:3]) <- 0
  initial(I34[1:N_age,1:3]) <- 0
  
  initial(R12[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*exp(-(eq_FOI3+eq_FOI4)*mean_age[i]))
  initial(R13[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*exp(-(eq_FOI2+eq_FOI4)*mean_age[i]))
  initial(R14[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI2+eq_FOI3)*mean_age[i]))
  initial(R23[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*exp(-(eq_FOI1+eq_FOI4)*mean_age[i]))
  initial(R24[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI1+eq_FOI3)*mean_age[i]))
  initial(R34[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-(eq_FOI1+eq_FOI2)*mean_age[i]))
  initial(R12[1:N_age,2:3]) <- 0
  initial(R13[1:N_age,2:3]) <- 0
  initial(R14[1:N_age,2:3]) <- 0
  initial(R23[1:N_age,2:3]) <- 0
  initial(R24[1:N_age,2:3]) <- 0
  initial(R34[1:N_age,2:3]) <- 0
  
  initial(I231[1:N_age,1:3]) <- 0
  initial(I241[1:N_age,1:3]) <- 0
  initial(I341[1:N_age,1:3]) <- 0
  initial(I132[1:N_age,1:3]) <- 0
  initial(I142[1:N_age,1:3]) <- 0
  initial(I342[1:N_age,1:3]) <- 0
  initial(I123[1:N_age,1:3]) <- 0
  initial(I143[1:N_age,1:3]) <- 0
  initial(I243[1:N_age,1:3]) <- 0
  initial(I124[1:N_age,1:3]) <- 0
  initial(I134[1:N_age,1:3]) <- 0
  initial(I234[1:N_age,1:3]) <- 0
  
  initial(R123[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*exp(-eq_FOI4*mean_age[i]))
  initial(R124[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-eq_FOI3*mean_age[i]))
  initial(R134[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-eq_FOI2*mean_age[i]))
  initial(R234[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i]))*exp(-eq_FOI1*mean_age[i]))
  initial(R123[1:N_age,2:3]) <- 0
  initial(R124[1:N_age,2:3]) <- 0
  initial(R134[1:N_age,2:3]) <- 0
  initial(R234[1:N_age,2:3]) <- 0
  
  initial(I2341[1:N_age,1:3]) <- 0
  initial(I1342[1:N_age,1:3]) <- 0
  initial(I1243[1:N_age,1:3]) <- 0
  initial(I1234[1:N_age,1:3]) <- 0
  
  initial(R1234[1:N_age,1]) <- floor(0.5+N_init_age[i]*(1-exp(-eq_FOI1*mean_age[i]))*(1-exp(-eq_FOI2*mean_age[i]))*(1-exp(-eq_FOI3*mean_age[i]))*(1-exp(-eq_FOI4*mean_age[i])))
  initial(R1234[1:N_age,2:3]) <- 0
  
  Ntotal[1:N_age,1:3] <- S[i,j]+I1[i,j]+I2[i,j]+I3[i,j]+I4[i,j]+R1[i,j]+R2[i,j]+R3[i,j]+R4[i,j]+I21[i,j]+I31[i,j]+I41[i,j]+I12[i,j]+I32[i,j]+I42[i,j]+I13[i,j]+I23[i,j]+I43[i,j]+I14[i,j]+I24[i,j]+I34[i,j]+R12[i,j]+R13[i,j]+R14[i,j]+R23[i,j]+R24[i,j]+R34[i,j]+I231[i,j]+I241[i,j]+I341[i,j]+I132[i,j]+I142[i,j]+I342[i,j]+I123[i,j]+I143[i,j]+I243[i,j]+I124[i,j]+I134[i,j]+I234[i,j]+R123[i,j]+R124[i,j]+R134[i,j]+R234[i,j]+I2341[i,j]+I1342[i,j]+I1243[i,j]+I1234[i,j]+R1234[i,j]
  NT <- sum(Ntotal)
  Ntotal_nv[1:N_age] <- Ntotal[i,1]
  Ntotal_nvS[1:N_age] <- S[i,1]
  Ntotal_v[1:N_age] <- Ntotal[i,2]+Ntotal[i,3]
  Ntotal_vS[1:N_age] <- Ntotal[i,3]
  Ntotal_all[1:N_age] <- Ntotal[i,1]+Ntotal[i,2]+Ntotal[i,3]
  NTv <- sum(Ntotal_v)
  NTnv <- sum(Ntotal_nv)
  NTvS <- sum(Ntotal_vS)
  NTvE <- NTv-NTvS
  


  ## these are just derived variables
  
  initial(Ntotal_out[1:N_age]) <- N_init_age[i]
  update(Ntotal_out[1:N_age]) <- Ntotal_all[i]
  initial(NT_out) <- sum(N_init_age)
  update(NT_out) <- NT
  initial(NTv_out) <- 0
  update(NTv_out) <- NTv
  initial(NTnv_out) <- sum(N_init_age)
  update(NTnv_out) <- NTnv
  initial(NTvS_out) <- 0
  update(NTvS_out) <- NTvS
  initial(NTvE_out) <- 0
  update(NTvE_out) <- NTvE
  initial(prop_seroneg_at_9_out) <- 0
  update(prop_seroneg_at_9_out) <- prop_seroneg_at_9
  initial(prop_seroneg_at_vca_out) <- 0
  update(prop_seroneg_at_vca_out) <- prop_seroneg_at_vca
  
  
  
  ## infectiousness weighted aggregated infectious compartments
  
  
  Y1[1:N_age,1] <- phi_pri[1]*(1+dis_pri[1]*phi_ed)*inf_1[i,j]+phi_sec[1]*(1+dis_sec[1]*phi_ed)*(inf_21[i,j]+inf_31[i,j]+inf_41[i,j])+phi_tert[1]*(1+dis_tert[1]*phi_ed)*(inf_231[i,j]+inf_241[i,j]+inf_341[i,j])+phi_quart[1]*(1+dis_quart[1]*phi_ed)*inf_2341[i,j]
  Y1[1:N_age,2:3] <- phi_pri[1]*(1+dis_pri[1]*RR_dis[1,1,i]*phi_ed)*inf_1[i,j]+phi_sec[1]*(1+dis_sec[1]*RR_dis[1,2,i]*phi_ed)*(inf_21[i,j]+inf_31[i,j]+inf_41[i,j])+phi_tert[1]*(1+dis_tert[1]*RR_dis[1,3,i]*phi_ed)*(inf_231[i,j]+inf_241[i,j]+inf_341[i,j])+phi_quart[1]*(1+dis_quart[1]*RR_dis[1,3,i]*phi_ed)*inf_2341[i,j]
  
  Y2[1:N_age,1] <- phi_pri[2]*(1+dis_pri[2]*phi_ed)*inf_2[i,j]+phi_sec[2]*(1+dis_sec[2]*phi_ed)*(inf_12[i,j]+inf_32[i,j]+inf_42[i,j])+phi_tert[2]*(1+dis_tert[2]*phi_ed)*(inf_132[i,j]+inf_142[i,j]+inf_342[i,j])+phi_quart[2]*(1+dis_quart[2]*phi_ed)*inf_1342[i,j]
  Y2[1:N_age,2:3] <- phi_pri[2]*(1+dis_pri[2]*RR_dis[2,1,i]*phi_ed)*inf_2[i,j]+phi_sec[2]*(1+dis_sec[2]*RR_dis[2,2,i]*phi_ed)*(inf_12[i,j]+inf_32[i,j]+inf_42[i,j])+phi_tert[2]*(1+dis_tert[2]*RR_dis[2,3,i]*phi_ed)*(inf_132[i,j]+inf_142[i,j]+inf_342[i,j])+phi_quart[2]*(1+dis_quart[2]*RR_dis[2,3,i]*phi_ed)*inf_1342[i,j]
  
  Y3[1:N_age,1] <- phi_pri[3]*(1+dis_pri[3]*phi_ed)*inf_3[i,j]+phi_sec[3]*(1+dis_sec[3]*phi_ed)*(inf_13[i,j]+inf_23[i,j]+inf_43[i,j])+phi_tert[3]*(1+dis_tert[3]*phi_ed)*(inf_123[i,j]+inf_143[i,j]+inf_243[i,j])+phi_quart[3]*(1+dis_quart[3]*phi_ed)*inf_1243[i,j]
  Y3[1:N_age,2:3] <- phi_pri[3]*(1+dis_pri[3]*RR_dis[3,1,i]*phi_ed)*inf_3[i,j]+phi_sec[3]*(1+dis_sec[3]*RR_dis[3,2,i]*phi_ed)*(inf_13[i,j]+inf_23[i,j]+inf_43[i,j])+phi_tert[3]*(1+dis_tert[3]*RR_dis[3,3,i]*phi_ed)*(inf_123[i,j]+inf_143[i,j]+inf_243[i,j])+phi_quart[3]*(1+dis_quart[3]*RR_dis[3,3,i]*phi_ed)*inf_1243[i,j]
  
  Y4[1:N_age,1] <- phi_pri[4]*(1+dis_pri[4]*phi_ed)*inf_4[i,j]+phi_sec[4]*(1+dis_sec[4]*phi_ed)*(inf_14[i,j]+inf_24[i,j]+inf_34[i,j])+phi_tert[4]*(1+dis_tert[4]*phi_ed)*(inf_124[i,j]+inf_134[i,j]+inf_234[i,j])+phi_quart[4]*(1+dis_quart[4]*phi_ed)*inf_1234[i,j]
  Y4[1:N_age,2:3] <- phi_pri[4]*(1+dis_pri[4]*RR_dis[4,1,i]*phi_ed)*inf_4[i,j]+phi_sec[4]*(1+dis_sec[4]*RR_dis[4,2,i]*phi_ed)*(inf_14[i,j]+inf_24[i,j]+inf_34[i,j])+phi_tert[4]*(1+dis_tert[4]*RR_dis[4,3,i]*phi_ed)*(inf_124[i,j]+inf_134[i,j]+inf_234[i,j])+phi_quart[4]*(1+dis_quart[4]*RR_dis[4,3,i]*phi_ed)*inf_1234[i,j]
  
  
  Y1T <- sum(Y1)/NT
  Y2T <- sum(Y2)/NT
  Y3T <- sum(Y3)/NT
  Y4T <- sum(Y4)/NT
  
  initial(Y1T_out) <- 0
  update(Y1T_out) <- Y1T
  initial(Y2T_out) <- 0
  update(Y2T_out) <- Y2T
  initial(Y3T_out) <- 0
  update(Y3T_out) <- Y3T
  initial(Y4T_out) <- 0
  update(Y4T_out) <- Y4T
  
  ## background force of infection (all serotypes) to keep infection going
  
  extInfRand <- extInf ## random(0,2*extInf)
  
  ## model human incubation and infectious period with aggregate serotype specific overlapping compartments
  ## rather than for every single infected state
  
  initial(exposed1) <- 0
  initial(exposed2) <- 0
  initial(exposed3) <- 0
  initial(exposed4) <- 0
  update(exposed1) <- exposed1 + Y1T - exposed1*DT/incub
  update(exposed2) <- exposed2 + Y2T - exposed2*DT/incub
  update(exposed3) <- exposed3 + Y3T - exposed3*DT/incub
  update(exposed4) <- exposed4 + Y4T - exposed4*DT/incub
  
  initial(infectious1) <- 0
  initial(infectious2) <- 0
  initial(infectious3) <- 0
  initial(infectious4) <- 0
  update(infectious1) <- infectious1*(1-DT/inf_per) + exposed1*DT/incub
  update(infectious2) <- infectious2*(1-DT/inf_per) + exposed2*DT/incub
  update(infectious3) <- infectious3*(1-DT/inf_per) + exposed3*DT/incub
  update(infectious4) <- infectious4*(1-DT/inf_per) + exposed4*DT/incub
  
  ## FOI on people (multiplied by DT)

  FOI1 <- DT*Beta_mh*kappa*(Mwt_I1+Mwb_I1*Wb_relinf1)/NT+extInfRand
  FOI2 <- DT*Beta_mh*kappa*(Mwt_I2+Mwb_I2*Wb_relinf2)/NT+extInfRand
  FOI3 <- DT*Beta_mh*kappa*(Mwt_I3+Mwb_I3*Wb_relinf3)/NT+extInfRand
  FOI4 <- DT*Beta_mh*kappa*(Mwt_I4+Mwb_I4*Wb_relinf4)/NT+extInfRand
  
  ## allow for age-related exposure
  
  FOI1a[1:N_age] <- FOIas[i]*FOI1
  FOI2a[1:N_age] <- FOIas[i]*FOI2
  FOI3a[1:N_age] <- FOIas[i]*FOI3
  FOI4a[1:N_age] <- FOIas[i]*FOI4
  
  #### derived variables to output
  
  ## output cumulative FOI
  
  initial(inc_FOI1) <- 0
  initial(inc_FOI2) <- 0
  initial(inc_FOI3) <- 0
  initial(inc_FOI4) <- 0
  update(inc_FOI1) <- inc_FOI1+FOI1
  update(inc_FOI2) <- inc_FOI2+FOI2
  update(inc_FOI3) <- inc_FOI3+FOI3
  update(inc_FOI4) <- inc_FOI4+FOI4
  initial(tot_inc_FOI) <- 0
  update(tot_inc_FOI) <- (if(NUM_YEAR_ACCUM*floor(YEAR/NUM_YEAR_ACCUM)==YEAR) 0 else tot_inc_FOI)+FOI1+FOI2+FOI3+FOI4

  ## track and output cumulative annual human infection incidence (primary, secondary, tert+quart)
  
  infection_pri[1:N_age,1:3] <- inf_1[i,j]+inf_2[i,j]+inf_3[i,j]+inf_4[i,j]
  infection_sec[1:N_age,1:3] <- inf_21[i,j]+inf_31[i,j]+inf_41[i,j]+inf_12[i,j]+inf_32[i,j]+inf_42[i,j]+inf_13[i,j]+inf_23[i,j]+inf_43[i,j]+inf_14[i,j]+inf_24[i,j]+inf_34[i,j]
  infection_tq[1:N_age,1:3] <- inf_231[i,j]+inf_241[i,j]+inf_341[i,j]+inf_132[i,j]+inf_142[i,j]+inf_342[i,j]+inf_123[i,j]+inf_143[i,j]+inf_243[i,j]+inf_124[i,j]+inf_134[i,j]+inf_234[i,j]+inf_2341[i,j]+inf_1342[i,j]+inf_1243[i,j]+inf_1234[i,j]

  initial(cum_infection_pri[1:N_age]) <- 0
  update(cum_infection_pri[1:N_age]) <- cum_infection_pri[i] + infection_pri[i,1]+infection_pri[i,2]+infection_pri[i,3]
  initial(cum_infection_priV[1:N_age]) <- 0
  update(cum_infection_priV[1:N_age]) <- cum_infection_priV[i] + infection_pri[i,2]+infection_pri[i,3]
  
  initial(cum_infection_sec[1:N_age]) <- 0
  update(cum_infection_sec[1:N_age]) <- cum_infection_sec[i] + infection_sec[i,1]+infection_sec[i,2]+infection_sec[i,3]
  initial(cum_infection_secV[1:N_age]) <- 0
  update(cum_infection_secV[1:N_age]) <- cum_infection_secV[i] + infection_sec[i,2]+infection_sec[i,3]
  
  initial(cum_infection_tq[1:N_age]) <- 0
  update(cum_infection_tq[1:N_age]) <- cum_infection_tq[i] + infection_tq[i,1]+infection_tq[i,2]+infection_tq[i,3]
  initial(cum_infection_tqV[1:N_age]) <- 0
  update(cum_infection_tqV[1:N_age]) <- cum_infection_tqV[i] + infection_tq[i,2]+infection_tq[i,3]

  
  ## Takeda WHO outputs:
  ##   - symp cases, hosp cases, YLLs
  ##   - total and by serotype
  ##   - for each of the first 20 years of vaccination
  ##   - in (a) whole population, (b) vaccinated population, (c) seroneg vacc, (d) seropos vacc
  ## plus: number vaccinated per year
  ## plus most of above (other than serotype specific outputs) calculated with discounting (for health economic calcs)
  
  
  out_update_switch[1:NYO] <- if(YEARS_POST_VACC>=(i-1) && YEARS_POST_VACC<i) 1 else 0
  
  ## VE over time
  
  initial(VE_seroneg[1:NYO]) <-0
  initial(VE_seropos_mono[1:NYO]) <-0
  initial(VE_seropos_multi[1:NYO]) <-0
  
  update(VE_seroneg[1:NYO]) <- if(out_update_switch[i]==0||(vca+i-1)>N_age) VE_seroneg[i] else RR_dis[1,1,as.integer(vca)+i-1]
  update(VE_seropos_mono[1:NYO]) <- if(out_update_switch[i]==0||(vca+i-1)>N_age) VE_seropos_mono[i] else RR_dis[1,2,as.integer(vca)+i-1]
  update(VE_seropos_multi[1:NYO]) <- if(out_update_switch[i]==0||(vca+i-1)>N_age) VE_seropos_multi[i] else RR_dis[1,3,as.integer(vca)+i-1]
  
  ## seropositive at 9 and vacc age- outputs cumulative average over first n (1-20) years of vacc period
  
  prop_seroneg_at_9 <- (S[10,1]+S[10,2]+S[10,3])/(Ntotal[10,1]+Ntotal[10,2]+Ntotal[10,3])
  prop_seroneg_at_vca <- (S[as.integer(vca)+1,1]+S[as.integer(vca)+1,2]+S[as.integer(vca)+1,3])/(Ntotal[as.integer(vca)+1,1]+Ntotal[as.integer(vca)+1,2]+Ntotal[as.integer(vca)+1,3])
  initial(cum_prop_seroneg_at_9) <- 0
  update(cum_prop_seroneg_at_9) <- cum_prop_seroneg_at_9 + (if(YEARS_POST_VACC<0) 0 else DT*prop_seroneg_at_9)
  av_prop_seroneg_at_9 <- if(YEARS_POST_VACC<=0) 0 else cum_prop_seroneg_at_9/YEARS_POST_VACC/YL
  
  initial(out_prop_seroneg_at_9[1:NYO]) <-0
  update(out_prop_seroneg_at_9[1:NYO]) <- if(out_update_switch[i]==0) out_prop_seroneg_at_9[i] else av_prop_seroneg_at_9
  
  initial(cum_prop_seroneg_at_vca) <- 0
  update(cum_prop_seroneg_at_vca) <- cum_prop_seroneg_at_vca + (if(YEARS_POST_VACC<0) 0 else DT*prop_seroneg_at_vca)
  av_prop_seroneg_at_vca <- if(YEARS_POST_VACC<=0) 0 else cum_prop_seroneg_at_vca/YEARS_POST_VACC/YL
  
  initial(out_prop_seroneg_at_vca[1:NYO]) <-0
  update(out_prop_seroneg_at_vca[1:NYO]) <- if(out_update_switch[i]==0) out_prop_seroneg_at_vca[i] else av_prop_seroneg_at_vca
  
  ## track symptomatic disease incidence and output cumulative totals
  disease_sero[1:N_age,1,1] <- dis_pri[1]*inf_1[i,j]+dis_sec[1]*(inf_21[i,j]+inf_31[i,j]+inf_41[i,j])+dis_tert[1]*(inf_231[i,j]+inf_241[i,j]+inf_341[i,j])+dis_quart[1]*inf_2341[i,j]
  disease_sero[1:N_age,2:3,1] <- RR_dis[1,1,i]*dis_pri[1]*inf_1[i,j]+RR_dis[1,2,i]*dis_sec[1]*(inf_21[i,j]+inf_31[i,j]+inf_41[i,j])+RR_dis[1,3,i]*(dis_tert[1]*(inf_231[i,j]+inf_241[i,j]+inf_341[i,j])+dis_quart[1]*inf_2341[i,j])
  
  disease_sero[1:N_age,1,2] <- dis_pri[2]*inf_2[i,j]+dis_sec[2]*(inf_12[i,j]+inf_32[i,j]+inf_42[i,j])+dis_tert[2]*(inf_132[i,j]+inf_142[i,j]+inf_342[i,j])+dis_quart[2]*inf_1342[i,j]
  disease_sero[1:N_age,2:3,2] <- RR_dis[2,1,i]*dis_pri[2]*inf_2[i,j]+RR_dis[2,2,i]*dis_sec[2]*(inf_12[i,j]+inf_32[i,j]+inf_42[i,j])+RR_dis[2,3,i]*(dis_tert[2]*(inf_132[i,j]+inf_142[i,j]+inf_342[i,j])+dis_quart[2]*inf_1342[i,j])
  
  disease_sero[1:N_age,1,3] <- dis_pri[3]*inf_3[i,j]+dis_sec[3]*(inf_13[i,j]+inf_23[i,j]+inf_43[i,j])+dis_tert[3]*(inf_123[i,j]+inf_143[i,j]+inf_243[i,j])+dis_quart[3]*inf_1243[i,j]
  disease_sero[1:N_age,2:3,3] <- RR_dis[3,1,i]*dis_pri[3]*inf_3[i,j]+RR_dis[3,2,i]*dis_sec[3]*(inf_13[i,j]+inf_23[i,j]+inf_43[i,j])+RR_dis[3,3,i]*(dis_tert[3]*(inf_123[i,j]+inf_143[i,j]+inf_243[i,j])+dis_quart[3]*inf_1243[i,j])
  
  disease_sero[1:N_age,1,4] <- dis_pri[4]*inf_4[i,j]+dis_sec[4]*(inf_14[i,j]+inf_24[i,j]+inf_34[i,j])+dis_tert[4]*(inf_124[i,j]+inf_134[i,j]+inf_234[i,j])+dis_quart[4]*inf_1234[i,j]
  disease_sero[1:N_age,2:3,4] <- RR_dis[4,1,i]*dis_pri[4]*inf_4[i,j]+RR_dis[4,2,i]*dis_sec[4]*(inf_14[i,j]+inf_24[i,j]+inf_34[i,j])+RR_dis[4,3,i]*(dis_tert[4]*(inf_124[i,j]+inf_134[i,j]+inf_234[i,j])+dis_quart[4]*inf_1234[i,j])
  
  disease_sero_vacc_pri[1:N_age,1] <- dis_pri[1]*RR_dis[1,1,i]*(inf_1[i,2]+inf_1[i,3])
  disease_sero_vacc_pri[1:N_age,2] <- dis_pri[2]*RR_dis[2,1,i]*(inf_2[i,2]+inf_2[i,3])
  disease_sero_vacc_pri[1:N_age,3] <- dis_pri[3]*RR_dis[3,1,i]*(inf_3[i,2]+inf_3[i,3])
  disease_sero_vacc_pri[1:N_age,4] <- dis_pri[4]*RR_dis[4,1,i]*(inf_4[i,2]+inf_4[i,3])
  
  dis_sero_unvacc[1:4] <- sum(disease_sero[,1,i])
  dis_sero_vacc[1:4] <- sum(disease_sero[,2,i])+sum(disease_sero[,3,i])
  dis_sero_pop[1:4] <- dis_sero_unvacc[i]+dis_sero_vacc[i]
  dis_sero_vacc_pri[1:4] <- sum(disease_sero_vacc_pri[,i])
  dis_sero_vacc_secp[1:4] <- dis_sero_vacc[i]-dis_sero_vacc_pri[i] ## corrected after outputs3
  dis_sero_vacc_neg[1:4] <- sum(disease_sero[,3,i])
  dis_sero_vacc_pos[1:4] <- sum(disease_sero[,2,i])
  dis_all_pop <- sum(dis_sero_pop)
  dis_all_vacc <- sum(dis_sero_vacc)
  dis_all_unvacc <- sum(dis_sero_unvacc)
  dis_all_vacc_neg <- sum(dis_sero_vacc_neg)  
  dis_all_vacc_pos <- sum(dis_sero_vacc_pos)
  dis_all_vacc_pri <- sum(dis_sero_vacc_pri)  
  dis_all_vacc_secp <- sum(dis_sero_vacc_secp)
  disc_dis_all_pop <- dis_all_pop*cum_disc
  disc_dis_all_vacc <- dis_all_vacc*cum_disc

  initial(out_dis_all_pop[1:NYO]) <-0
  initial(out_dis_sero_pop[1:NYO,1:4]) <-0
  initial(out_dis_all_unvacc[1:NYO]) <-0
  initial(out_dis_sero_unvacc[1:NYO,1:4]) <-0
  initial(out_dis_all_vacc[1:NYO]) <-0
  initial(out_dis_all_vacc_neg[1:NYO]) <-0
  initial(out_dis_all_vacc_pos[1:NYO]) <-0
  initial(out_dis_all_vacc_pri[1:NYO]) <-0
  initial(out_dis_all_vacc_secp[1:NYO]) <-0
  initial(out_dis_sero_vacc[1:NYO,1:4]) <-0
  initial(out_dis_sero_vacc_neg[1:NYO,1:4]) <-0
  initial(out_dis_sero_vacc_pos[1:NYO,1:4]) <-0
  initial(out_dis_sero_vacc_pri[1:NYO,1:4]) <-0
  initial(out_dis_sero_vacc_secp[1:NYO,1:4]) <-0
  initial(out_disc_dis_all_pop[1:NYO]) <-0
  initial(out_disc_dis_all_vacc[1:NYO]) <-0
  
  update(out_dis_all_pop[1:NYO]) <- out_dis_all_pop[i] + (if(out_update_switch[i]==0) 0 else dis_all_pop)
  update(out_dis_sero_pop[1:NYO,1:4]) <- out_dis_sero_pop[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_pop[j])
  update(out_dis_all_unvacc[1:NYO]) <- out_dis_all_unvacc[i] + (if(out_update_switch[i]==0) 0 else dis_all_unvacc)
  update(out_dis_sero_unvacc[1:NYO,1:4]) <- out_dis_sero_unvacc[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_unvacc[j])
  update(out_dis_all_vacc[1:NYO]) <- out_dis_all_vacc[i] + (if(out_update_switch[i]==0) 0 else dis_all_vacc)
  update(out_dis_all_vacc_neg[1:NYO]) <- out_dis_all_vacc_neg[i] + (if(out_update_switch[i]==0) 0 else dis_all_vacc_neg)
  update(out_dis_all_vacc_pos[1:NYO]) <- out_dis_all_vacc_pos[i] + (if(out_update_switch[i]==0) 0 else dis_all_vacc_pos)
  update(out_dis_all_vacc_pri[1:NYO]) <- out_dis_all_vacc_pri[i] + (if(out_update_switch[i]==0) 0 else dis_all_vacc_pri)
  update(out_dis_all_vacc_secp[1:NYO]) <- out_dis_all_vacc_secp[i] + (if(out_update_switch[i]==0) 0 else dis_all_vacc_secp)
  update(out_dis_sero_vacc[1:NYO,1:4]) <- out_dis_sero_vacc[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_vacc[j])
  update(out_dis_sero_vacc_neg[1:NYO,1:4]) <- out_dis_sero_vacc_neg[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_vacc_neg[j])
  update(out_dis_sero_vacc_pos[1:NYO,1:4]) <- out_dis_sero_vacc_pos[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_vacc_pos[j])
  update(out_dis_sero_vacc_pri[1:NYO,1:4]) <- out_dis_sero_vacc_pri[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_vacc_pri[j])
  update(out_dis_sero_vacc_secp[1:NYO,1:4]) <- out_dis_sero_vacc_secp[i,j] + (if(out_update_switch[i]==0) 0 else dis_sero_vacc_secp[j])
  update(out_disc_dis_all_pop[1:NYO]) <- out_disc_dis_all_pop[i] + (if(out_update_switch[i]==0) 0 else disc_dis_all_pop)
  update(out_disc_dis_all_vacc[1:NYO]) <- out_disc_dis_all_vacc[i] + (if(out_update_switch[i]==0) 0 else disc_dis_all_vacc)
  
  initial(tot_inc_dis) <- 0
  update(tot_inc_dis) <- (if(NUM_YEAR_ACCUM*floor(YEAR/NUM_YEAR_ACCUM)==YEAR) 0 else tot_inc_dis)+1e5*dis_all_pop/NT
  
  intYPV <- floor(YEARS_POST_VACC)+1 # add 2 since want age group vca+1
  
  cohort_nvacc_neg <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else Ntotal[as.integer(intYPV+vca),3]
  cohort_nvacc_pos <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else Ntotal[as.integer(intYPV+vca),2]
  dim(out_cohort_nvacc_neg) <- NYO
  dim(out_cohort_nvacc_pos) <- NYO
  initial(out_cohort_nvacc_neg[1:NYO]) <-0
  initial(out_cohort_nvacc_pos[1:NYO]) <-0
  update(out_cohort_nvacc_neg[1:NYO]) <- if(out_update_switch[i]==0) out_cohort_nvacc_neg[i] else cohort_nvacc_neg
  update(out_cohort_nvacc_pos[1:NYO]) <- if(out_update_switch[i]==0) out_cohort_nvacc_pos[i] else cohort_nvacc_pos
  
  dim(cohort_dis_sero_unvacc) <- 4
  dim(cohort_dis_sero_vacc) <- 4
  dim(cohort_dis_sero_vacc_pri) <- 4
  dim(cohort_dis_sero_vacc_secp) <- 4
  dim(cohort_dis_sero_vacc_neg) <- 4
  dim(cohort_dis_sero_vacc_pos) <- 4
  
  cohort_dis_sero_unvacc[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else disease_sero[as.integer(intYPV+vca),1,i]
  cohort_dis_sero_vacc[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else sum(disease_sero[as.integer(intYPV+vca),2:3,i])
  cohort_dis_sero_vacc_pri[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else disease_sero_vacc_pri[as.integer(intYPV+vca),i]
  cohort_dis_sero_vacc_secp[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else cohort_dis_sero_vacc[i]-cohort_dis_sero_vacc_pri[i]
  cohort_dis_sero_vacc_neg[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else disease_sero[as.integer(intYPV+vca),3,i]
  cohort_dis_sero_vacc_pos[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else disease_sero[as.integer(intYPV+vca),2,i]

  dim(out_cohort_dis_sero_unvacc) <- c(NYO,4)
  dim(out_cohort_dis_sero_vacc) <- c(NYO,4)
  dim(out_cohort_dis_sero_vacc_pri) <- c(NYO,4)
  dim(out_cohort_dis_sero_vacc_secp) <- c(NYO,4)
  dim(out_cohort_dis_sero_vacc_neg) <- c(NYO,4)
  dim(out_cohort_dis_sero_vacc_pos) <- c(NYO,4)
  
  initial(out_cohort_dis_sero_unvacc[1:NYO,1:4]) <-0
  initial(out_cohort_dis_sero_vacc[1:NYO,1:4]) <-0
  initial(out_cohort_dis_sero_vacc_pri[1:NYO,1:4]) <-0
  initial(out_cohort_dis_sero_vacc_secp[1:NYO,1:4]) <-0
  initial(out_cohort_dis_sero_vacc_neg[1:NYO,1:4]) <-0
  initial(out_cohort_dis_sero_vacc_pos[1:NYO,1:4]) <-0
  
  update(out_cohort_dis_sero_unvacc[1:NYO,1:4]) <- out_cohort_dis_sero_unvacc[i,j] + (if(out_update_switch[i]==0) 0 else cohort_dis_sero_unvacc[j])
  update(out_cohort_dis_sero_vacc[1:NYO,1:4]) <- out_cohort_dis_sero_vacc[i,j] + (if(out_update_switch[i]==0) 0 else cohort_dis_sero_vacc[j])
  update(out_cohort_dis_sero_vacc_pri[1:NYO,1:4]) <- out_cohort_dis_sero_vacc_pri[i,j] + (if(out_update_switch[i]==0) 0 else cohort_dis_sero_vacc_pri[j])
  update(out_cohort_dis_sero_vacc_secp[1:NYO,1:4]) <- out_cohort_dis_sero_vacc_secp[i,j] + (if(out_update_switch[i]==0) 0 else cohort_dis_sero_vacc_secp[j])
  update(out_cohort_dis_sero_vacc_neg[1:NYO,1:4]) <- out_cohort_dis_sero_vacc_neg[i,j] + (if(out_update_switch[i]==0) 0 else cohort_dis_sero_vacc_neg[j])
  update(out_cohort_dis_sero_vacc_pos[1:NYO,1:4]) <- out_cohort_dis_sero_vacc_pos[i,j] + (if(out_update_switch[i]==0) 0 else cohort_dis_sero_vacc_pos[j])


  ## severe disease
  sdisease_sero[1:N_age,1,1] <- sdis_pri[1]*inf_1[i,j]+sdis_sec[1]*(inf_21[i,j]+inf_31[i,j]+inf_41[i,j])+sdis_tert[1]*(inf_231[i,j]+inf_241[i,j]+inf_341[i,j])+sdis_quart[1]*inf_2341[i,j]
  sdisease_sero[1:N_age,2:3,1] <- RR_sdis[1,1,i]*sdis_pri[1]*inf_1[i,j]+RR_sdis[1,2,i]*sdis_sec[1]*(inf_21[i,j]+inf_31[i,j]+inf_41[i,j])+RR_sdis[1,3,i]*(sdis_tert[1]*(inf_231[i,j]+inf_241[i,j]+inf_341[i,j])+sdis_quart[1]*inf_2341[i,j])
  
  sdisease_sero[1:N_age,1,2] <- sdis_pri[2]*inf_2[i,j]+sdis_sec[2]*(inf_12[i,j]+inf_32[i,j]+inf_42[i,j])+sdis_tert[2]*(inf_132[i,j]+inf_142[i,j]+inf_342[i,j])+sdis_quart[2]*inf_1342[i,j]
  sdisease_sero[1:N_age,2:3,2] <- RR_sdis[2,1,i]*sdis_pri[2]*inf_2[i,j]+RR_sdis[2,2,i]*sdis_sec[2]*(inf_12[i,j]+inf_32[i,j]+inf_42[i,j])+RR_sdis[2,3,i]*(sdis_tert[2]*(inf_132[i,j]+inf_142[i,j]+inf_342[i,j])+sdis_quart[2]*inf_1342[i,j])
  
  sdisease_sero[1:N_age,1,3] <- sdis_pri[3]*inf_3[i,j]+sdis_sec[3]*(inf_13[i,j]+inf_23[i,j]+inf_43[i,j])+sdis_tert[3]*(inf_123[i,j]+inf_143[i,j]+inf_243[i,j])+sdis_quart[3]*inf_1243[i,j]
  sdisease_sero[1:N_age,2:3,3] <- RR_sdis[3,1,i]*sdis_pri[3]*inf_3[i,j]+RR_sdis[3,2,i]*sdis_sec[3]*(inf_13[i,j]+inf_23[i,j]+inf_43[i,j])+RR_sdis[3,3,i]*(sdis_tert[3]*(inf_123[i,j]+inf_143[i,j]+inf_243[i,j])+sdis_quart[3]*inf_1243[i,j])
  
  sdisease_sero[1:N_age,1,4] <- sdis_pri[4]*inf_4[i,j]+sdis_sec[4]*(inf_14[i,j]+inf_24[i,j]+inf_34[i,j])+sdis_tert[4]*(inf_124[i,j]+inf_134[i,j]+inf_234[i,j])+sdis_quart[4]*inf_1234[i,j]
  sdisease_sero[1:N_age,2:3,4] <- RR_sdis[4,1,i]*sdis_pri[4]*inf_4[i,j]+RR_sdis[4,2,i]*sdis_sec[4]*(inf_14[i,j]+inf_24[i,j]+inf_34[i,j])+RR_sdis[4,3,i]*(sdis_tert[4]*(inf_124[i,j]+inf_134[i,j]+inf_234[i,j])+sdis_quart[4]*inf_1234[i,j])
  
  sdisease_sero_vacc_pri[1:N_age,1] <- sdis_pri[1]*RR_sdis[1,1,i]*(inf_1[i,2]+inf_1[i,3])
  sdisease_sero_vacc_pri[1:N_age,2] <- sdis_pri[2]*RR_sdis[2,1,i]*(inf_2[i,2]+inf_2[i,3])
  sdisease_sero_vacc_pri[1:N_age,3] <- sdis_pri[3]*RR_sdis[3,1,i]*(inf_3[i,2]+inf_3[i,3])
  sdisease_sero_vacc_pri[1:N_age,4] <- sdis_pri[4]*RR_sdis[4,1,i]*(inf_4[i,2]+inf_4[i,3])
  
  sdis_sero_unvacc[1:4] <- sum(sdisease_sero[,1,i])
  sdis_sero_vacc[1:4] <- sum(sdisease_sero[,2,i])+sum(sdisease_sero[,3,i])
  sdis_sero_pop[1:4] <- sdis_sero_unvacc[i]+sdis_sero_vacc[i]
  sdis_sero_vacc_pri[1:4] <- sum(sdisease_sero_vacc_pri[,i])
  sdis_sero_vacc_secp[1:4] <- sdis_sero_vacc[i]-sdis_sero_vacc_pri[i] # corrected after outputs3
  sdis_sero_vacc_neg[1:4] <- sum(sdisease_sero[,3,i])
  sdis_sero_vacc_pos[1:4] <- sum(sdisease_sero[,2,i])
  sdis_all_pop <- sum(sdis_sero_pop)
  sdis_all_vacc <- sum(sdis_sero_vacc)
  sdis_all_unvacc <- sum(sdis_sero_unvacc)
  sdis_all_vacc_neg <- sum(sdis_sero_vacc_neg)  
  sdis_all_vacc_pos <- sum(sdis_sero_vacc_pos)
  sdis_all_vacc_pri <- sum(sdis_sero_vacc_pri)  
  sdis_all_vacc_secp <- sum(sdis_sero_vacc_secp)
  disc_sdis_all_pop <- sdis_all_pop*cum_disc
  disc_sdis_all_vacc <- sdis_all_vacc*cum_disc
  
  initial(out_sdis_all_pop[1:NYO]) <-0
  initial(out_sdis_sero_pop[1:NYO,1:4]) <-0
  initial(out_sdis_all_unvacc[1:NYO]) <-0
  initial(out_sdis_sero_unvacc[1:NYO,1:4]) <-0
  initial(out_sdis_all_vacc[1:NYO]) <-0
  initial(out_sdis_all_vacc_neg[1:NYO]) <-0
  initial(out_sdis_all_vacc_pos[1:NYO]) <-0
  initial(out_sdis_all_vacc_pri[1:NYO]) <-0
  initial(out_sdis_all_vacc_secp[1:NYO]) <-0
  initial(out_sdis_sero_vacc[1:NYO,1:4]) <-0
  initial(out_sdis_sero_vacc_neg[1:NYO,1:4]) <-0
  initial(out_sdis_sero_vacc_pos[1:NYO,1:4]) <-0
  initial(out_sdis_sero_vacc_pri[1:NYO,1:4]) <-0
  initial(out_sdis_sero_vacc_secp[1:NYO,1:4]) <-0
  initial(out_disc_sdis_all_pop[1:NYO]) <-0
  initial(out_disc_sdis_all_vacc[1:NYO]) <-0
  
  update(out_sdis_all_pop[1:NYO]) <- out_sdis_all_pop[i] + (if(out_update_switch[i]==0) 0 else sdis_all_pop)
  update(out_sdis_sero_pop[1:NYO,1:4]) <- out_sdis_sero_pop[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_pop[j])
  update(out_sdis_all_unvacc[1:NYO]) <- out_sdis_all_unvacc[i] + (if(out_update_switch[i]==0) 0 else sdis_all_unvacc)
  update(out_sdis_sero_unvacc[1:NYO,1:4]) <- out_sdis_sero_unvacc[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_unvacc[j])
  update(out_sdis_all_vacc[1:NYO]) <- out_sdis_all_vacc[i] + (if(out_update_switch[i]==0) 0 else sdis_all_vacc)
  update(out_sdis_all_vacc_neg[1:NYO]) <- out_sdis_all_vacc_neg[i] + (if(out_update_switch[i]==0) 0 else sdis_all_vacc_neg)
  update(out_sdis_all_vacc_pos[1:NYO]) <- out_sdis_all_vacc_pos[i] + (if(out_update_switch[i]==0) 0 else sdis_all_vacc_pos)
  update(out_sdis_all_vacc_pri[1:NYO]) <- out_sdis_all_vacc_pri[i] + (if(out_update_switch[i]==0) 0 else sdis_all_vacc_pri)
  update(out_sdis_all_vacc_secp[1:NYO]) <- out_sdis_all_vacc_secp[i] + (if(out_update_switch[i]==0) 0 else sdis_all_vacc_secp)
  update(out_sdis_sero_vacc[1:NYO,1:4]) <- out_sdis_sero_vacc[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_vacc[j])
  update(out_sdis_sero_vacc_neg[1:NYO,1:4]) <- out_sdis_sero_vacc_neg[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_vacc_neg[j])
  update(out_sdis_sero_vacc_pos[1:NYO,1:4]) <- out_sdis_sero_vacc_pos[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_vacc_pos[j])
  update(out_sdis_sero_vacc_pri[1:NYO,1:4]) <- out_sdis_sero_vacc_pri[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_vacc_pri[j])
  update(out_sdis_sero_vacc_secp[1:NYO,1:4]) <- out_sdis_sero_vacc_secp[i,j] + (if(out_update_switch[i]==0) 0 else sdis_sero_vacc_secp[j])
  update(out_disc_sdis_all_pop[1:NYO]) <- out_disc_sdis_all_pop[i] + (if(out_update_switch[i]==0) 0 else disc_sdis_all_pop)
  update(out_disc_sdis_all_vacc[1:NYO]) <- out_disc_sdis_all_vacc[i] + (if(out_update_switch[i]==0) 0 else disc_sdis_all_vacc)
  
  dim(cohort_sdis_sero_unvacc) <- 4
  dim(cohort_sdis_sero_vacc) <- 4
  dim(cohort_sdis_sero_vacc_pri) <- 4
  dim(cohort_sdis_sero_vacc_secp) <- 4
  dim(cohort_sdis_sero_vacc_neg) <- 4
  dim(cohort_sdis_sero_vacc_pos) <- 4
  
  cohort_sdis_sero_unvacc[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else sdisease_sero[as.integer(intYPV+vca),1,i]
  cohort_sdis_sero_vacc[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else sum(sdisease_sero[as.integer(intYPV+vca),2:3,i])
  cohort_sdis_sero_vacc_pri[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else sdisease_sero_vacc_pri[as.integer(intYPV+vca),i]
  cohort_sdis_sero_vacc_secp[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else cohort_sdis_sero_vacc[i]-cohort_sdis_sero_vacc_pri[i]
  cohort_sdis_sero_vacc_neg[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else sdisease_sero[as.integer(intYPV+vca),3,i]
  cohort_sdis_sero_vacc_pos[1:4] <- if((YEARS_POST_VACC<0)||((YEARS_POST_VACC)>=20)) 0 else sdisease_sero[as.integer(intYPV+vca),2,i]
  
  dim(out_cohort_sdis_sero_unvacc) <- c(NYO,4)
  dim(out_cohort_sdis_sero_vacc) <- c(NYO,4)
  dim(out_cohort_sdis_sero_vacc_pri) <- c(NYO,4)
  dim(out_cohort_sdis_sero_vacc_secp) <- c(NYO,4)
  dim(out_cohort_sdis_sero_vacc_neg) <- c(NYO,4)
  dim(out_cohort_sdis_sero_vacc_pos) <- c(NYO,4)
  
  initial(out_cohort_sdis_sero_unvacc[1:NYO,1:4]) <-0
  initial(out_cohort_sdis_sero_vacc[1:NYO,1:4]) <-0
  initial(out_cohort_sdis_sero_vacc_pri[1:NYO,1:4]) <-0
  initial(out_cohort_sdis_sero_vacc_secp[1:NYO,1:4]) <-0
  initial(out_cohort_sdis_sero_vacc_neg[1:NYO,1:4]) <-0
  initial(out_cohort_sdis_sero_vacc_pos[1:NYO,1:4]) <-0
  
  update(out_cohort_sdis_sero_unvacc[1:NYO,1:4]) <- out_cohort_sdis_sero_unvacc[i,j] + (if(out_update_switch[i]==0) 0 else cohort_sdis_sero_unvacc[j])
  update(out_cohort_sdis_sero_vacc[1:NYO,1:4]) <- out_cohort_sdis_sero_vacc[i,j] + (if(out_update_switch[i]==0) 0 else cohort_sdis_sero_vacc[j])
  update(out_cohort_sdis_sero_vacc_pri[1:NYO,1:4]) <- out_cohort_sdis_sero_vacc_pri[i,j] + (if(out_update_switch[i]==0) 0 else cohort_sdis_sero_vacc_pri[j])
  update(out_cohort_sdis_sero_vacc_secp[1:NYO,1:4]) <- out_cohort_sdis_sero_vacc_secp[i,j] + (if(out_update_switch[i]==0) 0 else cohort_sdis_sero_vacc_secp[j])
  update(out_cohort_sdis_sero_vacc_neg[1:NYO,1:4]) <- out_cohort_sdis_sero_vacc_neg[i,j] + (if(out_update_switch[i]==0) 0 else cohort_sdis_sero_vacc_neg[j])
  update(out_cohort_sdis_sero_vacc_pos[1:NYO,1:4]) <- out_cohort_sdis_sero_vacc_pos[i,j] + (if(out_update_switch[i]==0) 0 else cohort_sdis_sero_vacc_pos[j])
  
  
  ## years of life lost
  
  yll_sero[1:N_age,1:3,1:4] <- sdisease_sero[i,j,k]*cfr*life_expec[i]
  yll_sero_vacc_pri_age[1:N_age,1:4] <- sdisease_sero_vacc_pri[i,j]*cfr*life_expec[i]
  disc_yll_all[1:N_age,1:3] <- sum(sdisease_sero[i,j,])*cfr*(1.0-1.0/(disc_fact^life_expec[i]))/log(disc_fact)
  
  yll_sero_unvacc[1:4] <- sum(yll_sero[,1,i])
  yll_sero_vacc[1:4] <- sum(yll_sero[,2,i])+sum(yll_sero[,3,i])
  yll_sero_pop[1:4] <- yll_sero_unvacc[i]+yll_sero_vacc[i]
  yll_sero_vacc_pri[1:4] <- sum(yll_sero_vacc_pri_age[,i])
  yll_sero_vacc_secp[1:4] <- yll_sero_vacc[i]-yll_sero_vacc_neg[i]
  yll_sero_vacc_neg[1:4] <- sum(yll_sero[,3,i])
  yll_sero_vacc_pos[1:4] <- sum(yll_sero[,2,i])
  yll_all_pop <- sum(yll_sero_pop)
  yll_all_vacc <- sum(yll_sero_vacc)
  yll_all_unvacc <- sum(yll_sero_unvacc)
  yll_all_vacc_neg <- sum(yll_sero_vacc_neg)  
  yll_all_vacc_pos <- sum(yll_sero_vacc_pos)
  yll_all_vacc_pri <- sum(yll_sero_vacc_pri)  
  yll_all_vacc_secp <- sum(yll_sero_vacc_secp)
  disc_yll_all_pop <- cum_disc*sum(disc_yll_all)
  disc_yll_all_vacc <- cum_disc*(sum(disc_yll_all[,2])+sum(disc_yll_all[,3]))

  
  initial(out_yll_all_pop[1:NYO]) <-0
  initial(out_yll_sero_pop[1:NYO,1:4]) <-0
  initial(out_yll_all_unvacc[1:NYO]) <-0
  initial(out_yll_sero_unvacc[1:NYO,1:4]) <-0
  initial(out_yll_all_vacc[1:NYO]) <-0
  initial(out_yll_all_vacc_neg[1:NYO]) <-0
  initial(out_yll_all_vacc_pos[1:NYO]) <-0
  initial(out_yll_all_vacc_pri[1:NYO]) <-0
  initial(out_yll_all_vacc_secp[1:NYO]) <-0
  initial(out_yll_sero_vacc[1:NYO,1:4]) <-0
  initial(out_yll_sero_vacc_neg[1:NYO,1:4]) <-0
  initial(out_yll_sero_vacc_pos[1:NYO,1:4]) <-0
  initial(out_yll_sero_vacc_pri[1:NYO,1:4]) <-0
  initial(out_yll_sero_vacc_secp[1:NYO,1:4]) <-0
  initial(out_disc_yll_all_pop[1:NYO]) <-0
  initial(out_disc_yll_all_vacc[1:NYO]) <-0
  
  update(out_yll_all_pop[1:NYO]) <- out_yll_all_pop[i] + (if(out_update_switch[i]==0) 0 else yll_all_pop)
  update(out_yll_sero_pop[1:NYO,1:4]) <- out_yll_sero_pop[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_pop[j])
  update(out_yll_all_unvacc[1:NYO]) <- out_yll_all_unvacc[i] + (if(out_update_switch[i]==0) 0 else yll_all_unvacc)
  update(out_yll_sero_unvacc[1:NYO,1:4]) <- out_yll_sero_unvacc[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_unvacc[j])
  update(out_yll_all_vacc[1:NYO]) <- out_yll_all_vacc[i] + (if(out_update_switch[i]==0) 0 else yll_all_vacc)
  update(out_yll_all_vacc_neg[1:NYO]) <- out_yll_all_vacc_neg[i] + (if(out_update_switch[i]==0) 0 else yll_all_vacc_neg)
  update(out_yll_all_vacc_pos[1:NYO]) <- out_yll_all_vacc_pos[i] + (if(out_update_switch[i]==0) 0 else yll_all_vacc_pos)
  update(out_yll_all_vacc_pri[1:NYO]) <- out_yll_all_vacc_pri[i] + (if(out_update_switch[i]==0) 0 else yll_all_vacc_pri)
  update(out_yll_all_vacc_secp[1:NYO]) <- out_yll_all_vacc_secp[i] + (if(out_update_switch[i]==0) 0 else yll_all_vacc_secp)
  update(out_yll_sero_vacc[1:NYO,1:4]) <- out_yll_sero_vacc[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_vacc[j])
  update(out_yll_sero_vacc_neg[1:NYO,1:4]) <- out_yll_sero_vacc_neg[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_vacc_neg[j])
  update(out_yll_sero_vacc_pos[1:NYO,1:4]) <- out_yll_sero_vacc_pos[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_vacc_pos[j])
  update(out_yll_sero_vacc_pri[1:NYO,1:4]) <- out_yll_sero_vacc_pri[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_vacc_pri[j])
  update(out_yll_sero_vacc_secp[1:NYO,1:4]) <- out_yll_sero_vacc_secp[i,j] + (if(out_update_switch[i]==0) 0 else yll_sero_vacc_secp[j])
  update(out_disc_yll_all_pop[1:NYO]) <- out_disc_yll_all_pop[i] + (if(out_update_switch[i]==0) 0 else disc_yll_all_pop)
  update(out_disc_yll_all_vacc[1:NYO]) <- out_disc_yll_all_vacc[i] + (if(out_update_switch[i]==0) 0 else disc_yll_all_vacc)
  
  
  ## count vaccinated
  num_child_vacc_age[1:N_age] <- (if((YEARS_POST_VACC>=0) && (YEAR<vacc_child_stoptime) && (i == vacc_child_age)) (vacc_child_coverage) else 0)*Ntotal_nv[i]*agert[i]
  num_child_vacc_age_neg[1:N_age] <- (if((YEARS_POST_VACC>=0) && (YEAR<vacc_child_stoptime) && (i == vacc_child_age)) (vacc_child_coverage) else 0)*Ntotal_nvS[i]*agert[i]
  num_cu_vacc_age[1:N_age] <- (if((TIME == vacc_cu_rndtime) && (ageb[i]>=vacc_cu_minage) && (ageb[i]<=vacc_cu_maxage)) (vacc_cu_coverage*vacc_cu_age_weight[i]) else 0)*Ntotal_nv[i]
  num_cu_vacc_age_neg[1:N_age] <- (if((TIME == vacc_cu_rndtime) && (ageb[i]>=vacc_cu_minage) && (ageb[i]<=vacc_cu_maxage)) (vacc_cu_coverage*vacc_cu_age_weight[i]) else 0)*Ntotal_nvS[i]
  
  num_child_vacc_sum <- (sum(num_child_vacc_age)+sum(num_cu_vacc_age))
  num_child_vacc_sum_neg <- (sum(num_child_vacc_age_neg)+sum(num_cu_vacc_age_neg))
  num_child_vacc_sum_pos <- num_child_vacc_sum-num_child_vacc_sum_neg
  num_child_vacc_per100k <- num_child_vacc_sum*1e5/NT
  
  initial(out_nvacc_all_pop[1:NYO]) <-0
  initial(out_nvacc_all_neg[1:NYO]) <-0
  initial(out_nvacc_all_pos[1:NYO]) <-0
  initial(out_disc_nvacc_all_pop[1:NYO]) <-0
  
  update(out_nvacc_all_pop[1:NYO]) <- out_nvacc_all_pop[i] + (if(out_update_switch[i]==0) 0 else num_child_vacc_sum)
  update(out_nvacc_all_neg[1:NYO]) <- out_nvacc_all_neg[i] + (if(out_update_switch[i]==0) 0 else num_child_vacc_sum_neg)
  update(out_nvacc_all_pos[1:NYO]) <- out_nvacc_all_pos[i] + (if(out_update_switch[i]==0) 0 else num_child_vacc_sum_pos)
  update(out_disc_nvacc_all_pop[1:NYO]) <- out_disc_nvacc_all_pop[i] + (if(out_update_switch[i]==0) 0 else num_child_vacc_sum*cum_disc_vacc)
  
  
  ###
  ### update equations for human state variables
  ###
  
  ## helpful to define these to simplify later updates
  nveai[1] <- 1
  nnveai[1] <- 0
  nveai[2:3] <- 1
  nnveai[2:3] <- 0
  

  O_S[1:N_age,1:3] <- rbinom((S[i,j]), (rho1[i,j]*FOI1a[i]+rho2[i,j]*FOI2a[i]+rho3[i,j]*FOI3a[i]+rho4[i,j]*FOI4a[i]+deathrt[i]))
  inf_1[1:N_age,1:3] <- rbinom((O_S[i,j]), (rho1[i,j]*FOI1a[i]/(rho1[i,j]*FOI1a[i]+rho2[i,j]*FOI2a[i]+rho3[i,j]*FOI3a[i]+rho4[i,j]*FOI4a[i]+deathrt[i])))
  inf_2[1:N_age,1:3] <- rbinom((O_S[i,j]-inf_1[i,j]), (rho2[i,j]*FOI2a[i]/(rho2[i,j]*FOI2a[i]+rho3[i,j]*FOI3a[i]+rho4[i,j]*FOI4a[i]+deathrt[i])))
  inf_3[1:N_age,1:3] <- rbinom((O_S[i,j]-inf_1[i,j]-inf_2[i,j]), (rho3[i,j]*FOI3a[i]/(rho3[i,j]*FOI3a[i]+rho4[i,j]*FOI4a[i]+deathrt[i])))
  inf_4[1:N_age,1:3] <- rbinom((O_S[i,j]-inf_1[i,j]-inf_2[i,j]-inf_3[i,j]), (rho4[i,j]*FOI4a[i]/(rho4[i,j]*FOI4a[i]+deathrt[i])))
  age_S[2:N_age_p1,1:3] <- rbinom(S[i-1,j]-O_S[i-1,j], agert[i-1])
  age_S[1,1] <- rbinom(10000000, births/10000000) # rpois(births)
  age_S[1,2:3] <- 0
  nvacc_S[2:N_age,1] <- rbinom(age_S[i,j],vacc_noncov_S[i,j]) 
  nvacc_S[1,1] <-  age_S[1,1]
  nvacc_S[1:N_age,2] <- age_S[i,j]
  nvacc_S[1:N_age,3] <- age_S[i,j]+age_S[i,1]-nvacc_S[i,1]
  n_S[1:N_age,1:3] <-  nvacc_S[i,j]+S[i,j]-O_S[i,j]-age_S[i+1,j]


  O_I1[1:N_age,1:3] <- rbinom((I1[i,j]), (nu+deathrt[i]))
  recov_1[1:N_age,1:3] <- rbinom((O_I1[i,j]), (nu/(nu+deathrt[i])))
  infV_1[1:N_age,1] <- inf_1[i,j]+nnveai[j]*inf_1[i,2]
  infV_1[1:N_age,2] <- nveai[j]*inf_1[i,j]
  infV_1[1:N_age,3] <- inf_1[i,j]
  age_I1[2:N_age_p1,1:3] <- rbinom((I1[i-1,j]-O_I1[i-1,j]+infV_1[i-1,j]), agert[i-1])
  age_I1[1,1:3] <- 0
  nvacc_I1[2:N_age,1] <- rbinom(age_I1[i,j],vacc_noncov[i,j]) 
  nvacc_I1[1,1] <-  0
  nvacc_I1[1:N_age,2] <- age_I1[i,j]+age_I1[i,1]-nvacc_I1[i,1]
  nvacc_I1[1:N_age,3] <- age_I1[i,j]
  n_I1[1:N_age,1:3] <- nvacc_I1[i,j]+I1[i,j]-O_I1[i,j]+infV_1[i,j]-age_I1[i+1,j]
  
  O_I2[1:N_age,1:3] <- rbinom((I2[i,j]), (nu+deathrt[i]))
  recov_2[1:N_age,1:3] <- rbinom((O_I2[i,j]), (nu/(nu+deathrt[i])))
  infV_2[1:N_age,1] <- inf_2[i,j]+nnveai[j]*inf_2[i,2]
  infV_2[1:N_age,2] <- nveai[j]*inf_2[i,j]
  infV_2[1:N_age,3] <- inf_2[i,j]
  age_I2[2:N_age_p1,1:3] <- rbinom((I2[i-1,j]-O_I2[i-1,j]+infV_2[i-1,j]), agert[i-1])
  age_I2[1,1:3] <- 0
  nvacc_I2[2:N_age,1] <- rbinom(age_I2[i,j],vacc_noncov[i,j]) 
  nvacc_I2[1,1] <-  0
  nvacc_I2[1:N_age,2] <- age_I2[i,j]+age_I2[i,1]-nvacc_I2[i,1]
  nvacc_I2[1:N_age,3] <- age_I2[i,j]
  n_I2[1:N_age,1:3] <- nvacc_I2[i,j]+I2[i,j]-O_I2[i,j]+infV_2[i,j]-age_I2[i+1,j]
  
  O_I3[1:N_age,1:3] <- rbinom((I3[i,j]), (nu+deathrt[i]))
  recov_3[1:N_age,1:3] <- rbinom((O_I3[i,j]), (nu/(nu+deathrt[i])))
  infV_3[1:N_age,1] <- inf_3[i,j]+nnveai[j]*inf_3[i,2]
  infV_3[1:N_age,2] <- nveai[j]*inf_3[i,j]
  infV_3[1:N_age,3] <- inf_3[i,j]
  age_I3[2:N_age_p1,1:3] <- rbinom((I3[i-1,j]-O_I3[i-1,j]+infV_3[i-1,j]), agert[i-1])
  age_I3[1,1:3] <- 0
  nvacc_I3[2:N_age,1] <- rbinom(age_I3[i,j],vacc_noncov[i,j]) 
  nvacc_I3[1,1] <-  0
  nvacc_I3[1:N_age,2] <- age_I3[i,j]+age_I3[i,1]-nvacc_I3[i,1]
  nvacc_I3[1:N_age,3] <- age_I3[i,j]
  n_I3[1:N_age,1:3] <- nvacc_I3[i,j]+I3[i,j]-O_I3[i,j]+infV_3[i,j]-age_I3[i+1,j]
  
  O_I4[1:N_age,1:3] <- rbinom((I4[i,j]), (nu+deathrt[i]))
  recov_4[1:N_age,1:3] <- rbinom((O_I4[i,j]), (nu/(nu+deathrt[i])))
  infV_4[1:N_age,1] <- inf_4[i,j]+nnveai[j]*inf_4[i,2]
  infV_4[1:N_age,2] <- nveai[j]*inf_4[i,j]
  infV_4[1:N_age,3] <- inf_4[i,j]
  age_I4[2:N_age_p1,1:3] <- rbinom((I4[i-1,j]-O_I4[i-1,j]+infV_4[i-1,j]), agert[i-1])
  age_I4[1,1:3] <- 0
  nvacc_I4[2:N_age,1] <- rbinom(age_I4[i,j],vacc_noncov[i,j]) 
  nvacc_I4[1,1] <-  0
  nvacc_I4[1:N_age,2] <- age_I4[i,j]+age_I4[i,1]-nvacc_I4[i,1]
  nvacc_I4[1:N_age,3] <- age_I4[i,j]
  n_I4[1:N_age,1:3] <- nvacc_I4[i,j]+I4[i,j]-O_I4[i,j]+infV_4[i,j]-age_I4[i+1,j]
  
  O_R1[1:N_age,1:3] <- rbinom((R1[i,j]), (rho12[i,j]*FOI2a[i]+rho13[i,j]*FOI3a[i]+rho14[i,j]*FOI4a[i]+deathrt[i]))
  inf_12[1:N_age,1:3] <- rbinom((O_R1[i,j]), (rho12[i,j]*FOI2a[i]/(rho12[i,j]*FOI2a[i]+rho13[i,j]*FOI3a[i]+rho14[i,j]*FOI4a[i]+deathrt[i])))
  inf_13[1:N_age,1:3] <- rbinom((O_R1[i,j]-inf_12[i,j]), (rho13[i,j]*FOI3a[i]/(rho13[i,j]*FOI3a[i]+rho14[i,j]*FOI4a[i]+deathrt[i])))
  inf_14[1:N_age,1:3] <- rbinom((O_R1[i,j]-inf_12[i,j]-inf_13[i,j]), (rho14[i,j]*FOI4a[i]/(rho14[i,j]*FOI4a[i]+deathrt[i])))
  age_R1[2:N_age_p1,1:3] <- rbinom((R1[i-1,j]-O_R1[i-1,j]+recov_1[i-1,j]), agert[i-1])
  age_R1[1,1:3] <- 0
  nvacc_R1[2:N_age,1] <- rbinom(age_R1[i,j],vacc_noncov[i,j]) 
  nvacc_R1[1,1] <-  0
  nvacc_R1[1:N_age,2] <- age_R1[i,j]+age_R1[i,1]-nvacc_R1[i,1]
  nvacc_R1[1:N_age,3] <- age_R1[i,j]  
  n_R1[1:N_age,1:3] <- nvacc_R1[i,j]+R1[i,j]-O_R1[i,j]+recov_1[i,j]-age_R1[i+1,j]
  
  O_R2[1:N_age,1:3] <- rbinom((R2[i,j]), (rho21[i,j]*FOI1a[i]+rho23[i,j]*FOI3a[i]+rho24[i,j]*FOI4a[i]+deathrt[i]))
  inf_21[1:N_age,1:3] <- rbinom((O_R2[i,j]), (rho21[i,j]*FOI1a[i]/(rho21[i,j]*FOI1a[i]+rho23[i,j]*FOI3a[i]+rho24[i,j]*FOI4a[i]+deathrt[i])))
  inf_23[1:N_age,1:3] <- rbinom((O_R2[i,j]-inf_21[i,j]), (rho23[i,j]*FOI3a[i]/(rho23[i,j]*FOI3a[i]+rho24[i,j]*FOI4a[i]+deathrt[i])))
  inf_24[1:N_age,1:3] <- rbinom((O_R2[i,j]-inf_21[i,j]-inf_23[i,j]), (rho24[i,j]*FOI4a[i]/(rho24[i,j]*FOI4a[i]+deathrt[i])))
  age_R2[2:N_age_p1,1:3] <- rbinom((R2[i-1,j]-O_R2[i-1,j]+recov_2[i-1,j]), agert[i-1])
  age_R2[1,1:3] <- 0
  nvacc_R2[2:N_age,1] <- rbinom(age_R2[i,j],vacc_noncov[i,j]) 
  nvacc_R2[1,1] <-  0
  nvacc_R2[1:N_age,2] <- age_R2[i,j]+age_R2[i,1]-nvacc_R2[i,1]
  nvacc_R2[1:N_age,3] <- age_R2[i,j]  
  n_R2[1:N_age,1:3] <- nvacc_R2[i,j]+R2[i,j]-O_R2[i,j]+recov_2[i,j]-age_R2[i+1,j]
  
  O_R3[1:N_age,1:3] <- rbinom((R3[i,j]), (rho31[i,j]*FOI1a[i]+rho32[i,j]*FOI2a[i]+rho34[i,j]*FOI4a[i]+deathrt[i]))
  inf_31[1:N_age,1:3] <- rbinom((O_R3[i,j]), (rho31[i,j]*FOI1a[i]/(rho31[i,j]*FOI1a[i]+rho32[i,j]*FOI2a[i]+rho34[i,j]*FOI4a[i]+deathrt[i])))
  inf_32[1:N_age,1:3] <- rbinom((O_R3[i,j]-inf_31[i,j]), (rho32[i,j]*FOI2a[i]/(rho32[i,j]*FOI2a[i]+rho34[i,j]*FOI4a[i]+deathrt[i])))
  inf_34[1:N_age,1:3] <- rbinom((O_R3[i,j]-inf_31[i,j]-inf_32[i,j]), (rho34[i,j]*FOI4a[i]/(rho34[i,j]*FOI4a[i]+deathrt[i])))
  age_R3[2:N_age_p1,1:3] <- rbinom((R3[i-1,j]-O_R3[i-1,j]+recov_3[i-1,j]), agert[i-1])
  age_R3[1,1:3] <- 0
  nvacc_R3[2:N_age,1] <- rbinom(age_R3[i,j],vacc_noncov[i,j]) 
  nvacc_R3[1,1] <-  0
  nvacc_R3[1:N_age,2] <- age_R3[i,j]+age_R3[i,1]-nvacc_R3[i,1]
  nvacc_R3[1:N_age,3] <- age_R3[i,j]  
  n_R3[1:N_age,1:3] <- nvacc_R3[i,j]+R3[i,j]-O_R3[i,j]+recov_3[i,j]-age_R3[i+1,j]
  
  O_R4[1:N_age,1:3] <- rbinom((R4[i,j]), (rho41[i,j]*FOI1a[i]+rho42[i,j]*FOI2a[i]+rho43[i,j]*FOI3a[i]+deathrt[i]))
  inf_41[1:N_age,1:3] <- rbinom((O_R4[i,j]), (rho41[i,j]*FOI1a[i]/(rho41[i,j]*FOI1a[i]+rho42[i,j]*FOI2a[i]+rho43[i,j]*FOI3a[i]+deathrt[i])))
  inf_42[1:N_age,1:3] <- rbinom((O_R4[i,j]-inf_41[i,j]), (rho42[i,j]*FOI2a[i]/(rho42[i,j]*FOI2a[i]+rho43[i,j]*FOI3a[i]+deathrt[i])))
  inf_43[1:N_age,1:3] <- rbinom((O_R4[i,j]-inf_41[i,j]-inf_42[i,j]), (rho43[i,j]*FOI3a[i]/(rho43[i,j]*FOI3a[i]+deathrt[i])))
  age_R4[2:N_age_p1,1:3] <- rbinom((R4[i-1,j]-O_R4[i-1,j]+recov_4[i-1,j]), agert[i-1])
  age_R4[1,1:3] <- 0
  nvacc_R4[2:N_age,1] <- rbinom(age_R4[i,j],vacc_noncov[i,j]) 
  nvacc_R4[1,1] <-  0
  nvacc_R4[1:N_age,2] <- age_R4[i,j]+age_R4[i,1]-nvacc_R4[i,1]
  nvacc_R4[1:N_age,3] <- age_R4[i,j]  
  n_R4[1:N_age,1:3] <- nvacc_R4[i,j]+R4[i,j]-O_R4[i,j]+recov_4[i,j]-age_R4[i+1,j]

  O_I12[1:N_age,1:3] <- rbinom((I12[i,j]), (nu+deathrt[i]))
  recov_12[1:N_age,1:3] <- rbinom((O_I12[i,j]), (nu/(nu+deathrt[i])))
  age_I12[2:N_age_p1,1:3] <- rbinom((I12[i-1,j]-O_I12[i-1,j]+inf_12[i-1,j]), agert[i-1])
  age_I12[1,1:3] <- 0
  nvacc_I12[2:N_age,1] <- rbinom(age_I12[i,j],vacc_noncov[i,j]) 
  nvacc_I12[1,1] <-  0
  nvacc_I12[1:N_age,2] <- age_I12[i,j]+age_I12[i,1]-nvacc_I12[i,1]
  nvacc_I12[1:N_age,3] <- age_I12[i,j]  
  n_I12[1:N_age,1:3] <- nvacc_I12[i,j]+I12[i,j]-O_I12[i,j]+inf_12[i,j]-age_I12[i+1,j]
  
  O_I13[1:N_age,1:3] <- rbinom((I13[i,j]), (nu+deathrt[i]))
  recov_13[1:N_age,1:3] <- rbinom((O_I13[i,j]), (nu/(nu+deathrt[i])))
  age_I13[2:N_age_p1,1:3] <- rbinom((I13[i-1,j]-O_I13[i-1,j]+inf_13[i-1,j]), agert[i-1])
  age_I13[1,1:3] <- 0
  nvacc_I13[2:N_age,1] <- rbinom(age_I13[i,j],vacc_noncov[i,j]) 
  nvacc_I13[1,1] <-  0
  nvacc_I13[1:N_age,2] <- age_I13[i,j]+age_I13[i,1]-nvacc_I13[i,1]
  nvacc_I13[1:N_age,3] <- age_I13[i,j]  
  n_I13[1:N_age,1:3] <- nvacc_I13[i,j]+I13[i,j]-O_I13[i,j]+inf_13[i,j]-age_I13[i+1,j]
  
  O_I14[1:N_age,1:3] <- rbinom((I14[i,j]), (nu+deathrt[i]))
  recov_14[1:N_age,1:3] <- rbinom((O_I14[i,j]), (nu/(nu+deathrt[i])))
  age_I14[2:N_age_p1,1:3] <- rbinom((I14[i-1,j]-O_I14[i-1,j]+inf_14[i-1,j]), agert[i-1])
  age_I14[1,1:3] <- 0
  nvacc_I14[2:N_age,1] <- rbinom(age_I14[i,j],vacc_noncov[i,j]) 
  nvacc_I14[1,1] <-  0
  nvacc_I14[1:N_age,2] <- age_I14[i,j]+age_I14[i,1]-nvacc_I14[i,1]
  nvacc_I14[1:N_age,3] <- age_I14[i,j]  
  n_I14[1:N_age,1:3] <- nvacc_I14[i,j]+I14[i,j]-O_I14[i,j]+inf_14[i,j]-age_I14[i+1,j]
  
  O_I21[1:N_age,1:3] <- rbinom((I21[i,j]), (nu+deathrt[i]))
  recov_21[1:N_age,1:3] <- rbinom((O_I21[i,j]), (nu/(nu+deathrt[i])))
  age_I21[2:N_age_p1,1:3] <- rbinom((I21[i-1,j]-O_I21[i-1,j]+inf_21[i-1,j]), agert[i-1])
  age_I21[1,1:3] <- 0
  nvacc_I21[2:N_age,1] <- rbinom(age_I21[i,j],vacc_noncov[i,j]) 
  nvacc_I21[1,1] <-  0
  nvacc_I21[1:N_age,2] <- age_I21[i,j]+age_I21[i,1]-nvacc_I21[i,1]
  nvacc_I21[1:N_age,3] <- age_I21[i,j]  
  n_I21[1:N_age,1:3] <- nvacc_I21[i,j]+I21[i,j]-O_I21[i,j]+inf_21[i,j]-age_I21[i+1,j]
  
  O_I23[1:N_age,1:3] <- rbinom((I23[i,j]), (nu+deathrt[i]))
  recov_23[1:N_age,1:3] <- rbinom((O_I23[i,j]), (nu/(nu+deathrt[i])))
  age_I23[2:N_age_p1,1:3] <- rbinom((I23[i-1,j]-O_I23[i-1,j]+inf_23[i-1,j]), agert[i-1])
  age_I23[1,1:3] <- 0
  nvacc_I23[2:N_age,1] <- rbinom(age_I23[i,j],vacc_noncov[i,j]) 
  nvacc_I23[1,1] <-  0
  nvacc_I23[1:N_age,2] <- age_I23[i,j]+age_I23[i,1]-nvacc_I23[i,1]
  nvacc_I23[1:N_age,3] <- age_I23[i,j]  
  n_I23[1:N_age,1:3] <- nvacc_I23[i,j]+I23[i,j]-O_I23[i,j]+inf_23[i,j]-age_I23[i+1,j]
  
  O_I24[1:N_age,1:3] <- rbinom((I24[i,j]), (nu+deathrt[i]))
  recov_24[1:N_age,1:3] <- rbinom((O_I24[i,j]), (nu/(nu+deathrt[i])))
  age_I24[2:N_age_p1,1:3] <- rbinom((I24[i-1,j]-O_I24[i-1,j]+inf_24[i-1,j]), agert[i-1])
  age_I24[1,1:3] <- 0
  nvacc_I24[2:N_age,1] <- rbinom(age_I24[i,j],vacc_noncov[i,j]) 
  nvacc_I24[1,1] <-  0
  nvacc_I24[1:N_age,2] <- age_I24[i,j]+age_I24[i,1]-nvacc_I24[i,1]
  nvacc_I24[1:N_age,3] <- age_I24[i,j]  
  n_I24[1:N_age,1:3] <- nvacc_I24[i,j]+I24[i,j]-O_I24[i,j]+inf_24[i,j]-age_I24[i+1,j]
  
  O_I31[1:N_age,1:3] <- rbinom((I31[i,j]), (nu+deathrt[i]))
  recov_31[1:N_age,1:3] <- rbinom((O_I31[i,j]), (nu/(nu+deathrt[i])))
  age_I31[2:N_age_p1,1:3] <- rbinom((I31[i-1,j]-O_I31[i-1,j]+inf_31[i-1,j]), agert[i-1])
  age_I31[1,1:3] <- 0
  nvacc_I31[2:N_age,1] <- rbinom(age_I31[i,j],vacc_noncov[i,j]) 
  nvacc_I31[1,1] <-  0
  nvacc_I31[1:N_age,2] <- age_I31[i,j]+age_I31[i,1]-nvacc_I31[i,1]
  nvacc_I31[1:N_age,3] <- age_I31[i,j]  
  n_I31[1:N_age,1:3] <- nvacc_I31[i,j]+I31[i,j]-O_I31[i,j]+inf_31[i,j]-age_I31[i+1,j]
  
  O_I32[1:N_age,1:3] <- rbinom((I32[i,j]), (nu+deathrt[i]))
  recov_32[1:N_age,1:3] <- rbinom((O_I32[i,j]), (nu/(nu+deathrt[i])))
  age_I32[2:N_age_p1,1:3] <- rbinom((I32[i-1,j]-O_I32[i-1,j]+inf_32[i-1,j]), agert[i-1])
  age_I32[1,1:3] <- 0
  nvacc_I32[2:N_age,1] <- rbinom(age_I32[i,j],vacc_noncov[i,j]) 
  nvacc_I32[1,1] <-  0
  nvacc_I32[1:N_age,2] <- age_I32[i,j]+age_I32[i,1]-nvacc_I32[i,1]
  nvacc_I32[1:N_age,3] <- age_I32[i,j]  
  n_I32[1:N_age,1:3] <- nvacc_I32[i,j]+I32[i,j]-O_I32[i,j]+inf_32[i,j]-age_I32[i+1,j]
  
  O_I34[1:N_age,1:3] <- rbinom((I34[i,j]), (nu+deathrt[i]))
  recov_34[1:N_age,1:3] <- rbinom((O_I34[i,j]), (nu/(nu+deathrt[i])))
  age_I34[2:N_age_p1,1:3] <- rbinom((I34[i-1,j]-O_I34[i-1,j]+inf_34[i-1,j]), agert[i-1])
  age_I34[1,1:3] <- 0
  nvacc_I34[2:N_age,1] <- rbinom(age_I34[i,j],vacc_noncov[i,j]) 
  nvacc_I34[1,1] <-  0
  nvacc_I34[1:N_age,2] <- age_I34[i,j]+age_I34[i,1]-nvacc_I34[i,1]
  nvacc_I34[1:N_age,3] <- age_I34[i,j]  
  n_I34[1:N_age,1:3] <- nvacc_I34[i,j]+I34[i,j]-O_I34[i,j]+inf_34[i,j]-age_I34[i+1,j]
  
  O_I41[1:N_age,1:3] <- rbinom((I41[i,j]), (nu+deathrt[i]))
  recov_41[1:N_age,1:3] <- rbinom((O_I41[i,j]), (nu/(nu+deathrt[i])))
  age_I41[2:N_age_p1,1:3] <- rbinom((I41[i-1,j]-O_I41[i-1,j]+inf_41[i-1,j]), agert[i-1])
  age_I41[1,1:3] <- 0
  nvacc_I41[2:N_age,1] <- rbinom(age_I41[i,j],vacc_noncov[i,j]) 
  nvacc_I41[1,1] <-  0
  nvacc_I41[1:N_age,2] <- age_I41[i,j]+age_I41[i,1]-nvacc_I41[i,1]
  nvacc_I41[1:N_age,3] <- age_I41[i,j]  
  n_I41[1:N_age,1:3] <- nvacc_I41[i,j]+I41[i,j]-O_I41[i,j]+inf_41[i,j]-age_I41[i+1,j]
  
  O_I42[1:N_age,1:3] <- rbinom((I42[i,j]), (nu+deathrt[i]))
  recov_42[1:N_age,1:3] <- rbinom((O_I42[i,j]), (nu/(nu+deathrt[i])))
  age_I42[2:N_age_p1,1:3] <- rbinom((I42[i-1,j]-O_I42[i-1,j]+inf_42[i-1,j]), agert[i-1])
  age_I42[1,1:3] <- 0
  nvacc_I42[2:N_age,1] <- rbinom(age_I42[i,j],vacc_noncov[i,j]) 
  nvacc_I42[1,1] <-  0
  nvacc_I42[1:N_age,2] <- age_I42[i,j]+age_I42[i,1]-nvacc_I42[i,1]
  nvacc_I42[1:N_age,3] <- age_I42[i,j]  
  n_I42[1:N_age,1:3] <- nvacc_I42[i,j]+I42[i,j]-O_I42[i,j]+inf_42[i,j]-age_I42[i+1,j]
  
  O_I43[1:N_age,1:3] <- rbinom((I43[i,j]), (nu+deathrt[i]))
  recov_43[1:N_age,1:3] <- rbinom((O_I43[i,j]), (nu/(nu+deathrt[i])))
  age_I43[2:N_age_p1,1:3] <- rbinom((I43[i-1,j]-O_I43[i-1,j]+inf_43[i-1,j]), agert[i-1])
  age_I43[1,1:3] <- 0
  nvacc_I43[2:N_age,1] <- rbinom(age_I43[i,j],vacc_noncov[i,j]) 
  nvacc_I43[1,1] <-  0
  nvacc_I43[1:N_age,2] <- age_I43[i,j]+age_I43[i,1]-nvacc_I43[i,1]
  nvacc_I43[1:N_age,3] <- age_I43[i,j]  
  n_I43[1:N_age,1:3] <- nvacc_I43[i,j]+I43[i,j]-O_I43[i,j]+inf_43[i,j]-age_I43[i+1,j]
  
  O_R12[1:N_age,1:3] <- rbinom((R12[i,j]), (rho123[i,j]*FOI3a[i]+rho124[i,j]*FOI4a[i]+deathrt[i]))
  inf_123[1:N_age,1:3] <- rbinom((O_R12[i,j]), (rho123[i,j]*FOI3a[i]/(rho123[i,j]*FOI3a[i]+rho124[i,j]*FOI4a[i]+deathrt[i])))
  inf_124[1:N_age,1:3] <- rbinom((O_R12[i,j]-inf_123[i,j]), (rho124[i,j]*FOI4a[i]/(rho124[i,j]*FOI4a[i]+deathrt[i])))
  age_R12[2:N_age_p1,1:3] <- rbinom((R12[i-1,j]-O_R12[i-1,j]+recov_12[i-1,j]+recov_21[i-1,j]), agert[i-1])
  age_R12[1,1:3] <- 0
  nvacc_R12[2:N_age,1] <- rbinom(age_R12[i,j],vacc_noncov[i,j]) 
  nvacc_R12[1,1] <-  0
  nvacc_R12[1:N_age,2] <- age_R12[i,j]+age_R12[i,1]-nvacc_R12[i,1]
  nvacc_R12[1:N_age,3] <- age_R12[i,j] 
  n_R12[1:N_age,1:3] <- nvacc_R12[i,j]+R12[i,j]-O_R12[i,j]+recov_12[i,j]+recov_21[i,j]-age_R12[i+1,j]

  O_R13[1:N_age,1:3] <- rbinom((R13[i,j]), (rho132[i,j]*FOI2a[i]+rho134[i,j]*FOI4a[i]+deathrt[i]))
  inf_132[1:N_age,1:3] <- rbinom((O_R13[i,j]), (rho132[i,j]*FOI2a[i]/(rho132[i,j]*FOI2a[i]+rho134[i,j]*FOI4a[i]+deathrt[i])))
  inf_134[1:N_age,1:3] <- rbinom((O_R13[i,j]-inf_132[i,j]), (rho134[i,j]*FOI4a[i]/(rho134[i,j]*FOI4a[i]+deathrt[i])))
  age_R13[2:N_age_p1,1:3] <- rbinom((R13[i-1,j]-O_R13[i-1,j]+recov_13[i-1,j]+recov_31[i-1,j]), agert[i-1])
  age_R13[1,1:3] <- 0
  nvacc_R13[2:N_age,1] <- rbinom(age_R13[i,j],vacc_noncov[i,j]) 
  nvacc_R13[1,1] <-  0
  nvacc_R13[1:N_age,2] <- age_R13[i,j]+age_R13[i,1]-nvacc_R13[i,1]
  nvacc_R13[1:N_age,3] <- age_R13[i,j] 
  n_R13[1:N_age,1:3] <- nvacc_R13[i,j]+R13[i,j]-O_R13[i,j]+recov_13[i,j]+recov_31[i,j]-age_R13[i+1,j]
  
  O_R14[1:N_age,1:3] <- rbinom((R14[i,j]), (rho142[i,j]*FOI2a[i]+rho143[i,j]*FOI3a[i]+deathrt[i]))
  inf_142[1:N_age,1:3] <- rbinom((O_R14[i,j]), (rho142[i,j]*FOI2a[i]/(rho142[i,j]*FOI2a[i]+rho143[i,j]*FOI3a[i]+deathrt[i])))
  inf_143[1:N_age,1:3] <- rbinom((O_R14[i,j]-inf_142[i,j]), (rho143[i,j]*FOI3a[i]/(rho143[i,j]*FOI3a[i]+deathrt[i])))
  age_R14[2:N_age_p1,1:3] <- rbinom((R14[i-1,j]-O_R14[i-1,j]+recov_14[i-1,j]+recov_41[i-1,j]), agert[i-1])
  age_R14[1,1:3] <- 0
  nvacc_R14[2:N_age,1] <- rbinom(age_R14[i,j],vacc_noncov[i,j]) 
  nvacc_R14[1,1] <-  0
  nvacc_R14[1:N_age,2] <- age_R14[i,j]+age_R14[i,1]-nvacc_R14[i,1]
  nvacc_R14[1:N_age,3] <- age_R14[i,j] 
  n_R14[1:N_age,1:3] <- nvacc_R14[i,j]+R14[i,j]-O_R14[i,j]+recov_14[i,j]+recov_41[i,j]-age_R14[i+1,j]
  
  O_R23[1:N_age,1:3] <- rbinom((R23[i,j]), (rho231[i,j]*FOI1a[i]+rho234[i,j]*FOI4a[i]+deathrt[i]))
  inf_231[1:N_age,1:3] <- rbinom((O_R23[i,j]), (rho231[i,j]*FOI1a[i]/(rho231[i,j]*FOI1a[i]+rho234[i,j]*FOI4a[i]+deathrt[i])))
  inf_234[1:N_age,1:3] <- rbinom((O_R23[i,j]-inf_231[i,j]), (rho234[i,j]*FOI4a[i]/(rho234[i,j]*FOI4a[i]+deathrt[i])))
  age_R23[2:N_age_p1,1:3] <- rbinom((R23[i-1,j]-O_R23[i-1,j]+recov_23[i-1,j]+recov_32[i-1,j]), agert[i-1])
  age_R23[1,1:3] <- 0
  nvacc_R23[2:N_age,1] <- rbinom(age_R23[i,j],vacc_noncov[i,j]) 
  nvacc_R23[1,1] <-  0
  nvacc_R23[1:N_age,2] <- age_R23[i,j]+age_R23[i,1]-nvacc_R23[i,1]
  nvacc_R23[1:N_age,3] <- age_R23[i,j] 
  n_R23[1:N_age,1:3] <- nvacc_R23[i,j]+R23[i,j]-O_R23[i,j]+recov_23[i,j]+recov_32[i,j]-age_R23[i+1,j]
  
  O_R24[1:N_age,1:3] <- rbinom((R24[i,j]), (rho241[i,j]*FOI1a[i]+rho243[i,j]*FOI3a[i]+deathrt[i]))
  inf_241[1:N_age,1:3] <- rbinom((O_R24[i,j]), (rho241[i,j]*FOI1a[i]/(rho241[i,j]*FOI1a[i]+rho243[i,j]*FOI3a[i]+deathrt[i])))
  inf_243[1:N_age,1:3] <- rbinom((O_R24[i,j]-inf_241[i,j]), (rho243[i,j]*FOI3a[i]/(rho243[i,j]*FOI3a[i]+deathrt[i])))
  age_R24[2:N_age_p1,1:3] <- rbinom((R24[i-1,j]-O_R24[i-1,j]+recov_24[i-1,j]+recov_42[i-1,j]), agert[i-1])
  age_R24[1,1:3] <- 0
  nvacc_R24[2:N_age,1] <- rbinom(age_R24[i,j],vacc_noncov[i,j]) 
  nvacc_R24[1,1] <-  0
  nvacc_R24[1:N_age,2] <- age_R24[i,j]+age_R24[i,1]-nvacc_R24[i,1]
  nvacc_R24[1:N_age,3] <- age_R24[i,j] 
  n_R24[1:N_age,1:3] <- nvacc_R24[i,j]+R24[i,j]-O_R24[i,j]+recov_24[i,j]+recov_42[i,j]-age_R24[i+1,j]
  
  O_R34[1:N_age,1:3] <- rbinom((R34[i,j]), (rho341[i,j]*FOI1a[i]+rho342[i,j]*FOI2a[i]+deathrt[i]))
  inf_341[1:N_age,1:3] <- rbinom((O_R34[i,j]), (rho341[i,j]*FOI1a[i]/(rho341[i,j]*FOI1a[i]+rho342[i,j]*FOI2a[i]+deathrt[i])))
  inf_342[1:N_age,1:3] <- rbinom((O_R34[i,j]-inf_341[i,j]), (rho342[i,j]*FOI2a[i]/(rho342[i,j]*FOI2a[i]+deathrt[i])))
  age_R34[2:N_age_p1,1:3] <- rbinom((R34[i-1,j]-O_R34[i-1,j]+recov_34[i-1,j]+recov_43[i-1,j]), agert[i-1])
  age_R34[1,1:3] <- 0
  nvacc_R34[2:N_age,1] <- rbinom(age_R34[i,j],vacc_noncov[i,j]) 
  nvacc_R34[1,1] <-  0
  nvacc_R34[1:N_age,2] <- age_R34[i,j]+age_R34[i,1]-nvacc_R34[i,1]
  nvacc_R34[1:N_age,3] <- age_R34[i,j] 
  n_R34[1:N_age,1:3] <- nvacc_R34[i,j]+R34[i,j]-O_R34[i,j]+recov_34[i,j]+recov_43[i,j]-age_R34[i+1,j]
 
  O_I123[1:N_age,1:3] <- rbinom((I123[i,j]), (nu+deathrt[i]))
  recov_123[1:N_age,1:3] <- rbinom((O_I123[i,j]), (nu/(nu+deathrt[i])))
  age_I123[2:N_age_p1,1:3] <- rbinom((I123[i-1,j]-O_I123[i-1,j]+inf_123[i-1,j]), agert[i-1])
  age_I123[1,1:3] <- 0
  nvacc_I123[2:N_age,1] <- rbinom(age_I123[i,j],vacc_noncov[i,j]) 
  nvacc_I123[1,1] <-  0
  nvacc_I123[1:N_age,2] <- age_I123[i,j]+age_I123[i,1]-nvacc_I123[i,1]
  nvacc_I123[1:N_age,3] <- age_I123[i,j]  
  n_I123[1:N_age,1:3] <- nvacc_I123[i,j]+I123[i,j]-O_I123[i,j]+inf_123[i,j]-age_I123[i+1,j]
  
  O_I124[1:N_age,1:3] <- rbinom((I124[i,j]), (nu+deathrt[i]))
  recov_124[1:N_age,1:3] <- rbinom((O_I124[i,j]), (nu/(nu+deathrt[i])))
  age_I124[2:N_age_p1,1:3] <- rbinom((I124[i-1,j]-O_I124[i-1,j]+inf_124[i-1,j]), agert[i-1])
  age_I124[1,1:3] <- 0
  nvacc_I124[2:N_age,1] <- rbinom(age_I124[i,j],vacc_noncov[i,j]) 
  nvacc_I124[1,1] <-  0
  nvacc_I124[1:N_age,2] <- age_I124[i,j]+age_I124[i,1]-nvacc_I124[i,1]
  nvacc_I124[1:N_age,3] <- age_I124[i,j]  
  n_I124[1:N_age,1:3] <- nvacc_I124[i,j]+I124[i,j]-O_I124[i,j]+inf_124[i,j]-age_I124[i+1,j]
  
  O_I132[1:N_age,1:3] <- rbinom((I132[i,j]), (nu+deathrt[i]))
  recov_132[1:N_age,1:3] <- rbinom((O_I132[i,j]), (nu/(nu+deathrt[i])))
  age_I132[2:N_age_p1,1:3] <- rbinom((I132[i-1,j]-O_I132[i-1,j]+inf_132[i-1,j]), agert[i-1])
  age_I132[1,1:3] <- 0
  nvacc_I132[2:N_age,1] <- rbinom(age_I132[i,j],vacc_noncov[i,j]) 
  nvacc_I132[1,1] <-  0
  nvacc_I132[1:N_age,2] <- age_I132[i,j]+age_I132[i,1]-nvacc_I132[i,1]
  nvacc_I132[1:N_age,3] <- age_I132[i,j]  
  n_I132[1:N_age,1:3] <- nvacc_I132[i,j]+I132[i,j]-O_I132[i,j]+inf_132[i,j]-age_I132[i+1,j]
  
  O_I134[1:N_age,1:3] <- rbinom((I134[i,j]), (nu+deathrt[i]))
  recov_134[1:N_age,1:3] <- rbinom((O_I134[i,j]), (nu/(nu+deathrt[i])))
  age_I134[2:N_age_p1,1:3] <- rbinom((I134[i-1,j]-O_I134[i-1,j]+inf_134[i-1,j]), agert[i-1])
  age_I134[1,1:3] <- 0
  nvacc_I134[2:N_age,1] <- rbinom(age_I134[i,j],vacc_noncov[i,j]) 
  nvacc_I134[1,1] <-  0
  nvacc_I134[1:N_age,2] <- age_I134[i,j]+age_I134[i,1]-nvacc_I134[i,1]
  nvacc_I134[1:N_age,3] <- age_I134[i,j]  
  n_I134[1:N_age,1:3] <- nvacc_I134[i,j]+I134[i,j]-O_I134[i,j]+inf_134[i,j]-age_I134[i+1,j]
  
  O_I142[1:N_age,1:3] <- rbinom((I142[i,j]), (nu+deathrt[i]))
  recov_142[1:N_age,1:3] <- rbinom((O_I142[i,j]), (nu/(nu+deathrt[i])))
  age_I142[2:N_age_p1,1:3] <- rbinom((I142[i-1,j]-O_I142[i-1,j]+inf_142[i-1,j]), agert[i-1])
  age_I142[1,1:3] <- 0
  nvacc_I142[2:N_age,1] <- rbinom(age_I142[i,j],vacc_noncov[i,j]) 
  nvacc_I142[1,1] <-  0
  nvacc_I142[1:N_age,2] <- age_I142[i,j]+age_I142[i,1]-nvacc_I142[i,1]
  nvacc_I142[1:N_age,3] <- age_I142[i,j]  
  n_I142[1:N_age,1:3] <- nvacc_I142[i,j]+I142[i,j]-O_I142[i,j]+inf_142[i,j]-age_I142[i+1,j]
  
  O_I143[1:N_age,1:3] <- rbinom((I143[i,j]), (nu+deathrt[i]))
  recov_143[1:N_age,1:3] <- rbinom((O_I143[i,j]), (nu/(nu+deathrt[i])))
  age_I143[2:N_age_p1,1:3] <- rbinom((I143[i-1,j]-O_I143[i-1,j]+inf_143[i-1,j]), agert[i-1])
  age_I143[1,1:3] <- 0
  nvacc_I143[2:N_age,1] <- rbinom(age_I143[i,j],vacc_noncov[i,j]) 
  nvacc_I143[1,1] <-  0
  nvacc_I143[1:N_age,2] <- age_I143[i,j]+age_I143[i,1]-nvacc_I143[i,1]
  nvacc_I143[1:N_age,3] <- age_I143[i,j]  
  n_I143[1:N_age,1:3] <- nvacc_I143[i,j]+I143[i,j]-O_I143[i,j]+inf_143[i,j]-age_I143[i+1,j]
  
  O_I231[1:N_age,1:3] <- rbinom((I231[i,j]), (nu+deathrt[i]))
  recov_231[1:N_age,1:3] <- rbinom((O_I231[i,j]), (nu/(nu+deathrt[i])))
  age_I231[2:N_age_p1,1:3] <- rbinom((I231[i-1,j]-O_I231[i-1,j]+inf_231[i-1,j]), agert[i-1])
  age_I231[1,1:3] <- 0
  nvacc_I231[2:N_age,1] <- rbinom(age_I231[i,j],vacc_noncov[i,j]) 
  nvacc_I231[1,1] <-  0
  nvacc_I231[1:N_age,2] <- age_I231[i,j]+age_I231[i,1]-nvacc_I231[i,1]
  nvacc_I231[1:N_age,3] <- age_I231[i,j]  
  n_I231[1:N_age,1:3] <- nvacc_I231[i,j]+I231[i,j]-O_I231[i,j]+inf_231[i,j]-age_I231[i+1,j]
  
  O_I234[1:N_age,1:3] <- rbinom((I234[i,j]), (nu+deathrt[i]))
  recov_234[1:N_age,1:3] <- rbinom((O_I234[i,j]), (nu/(nu+deathrt[i])))
  age_I234[2:N_age_p1,1:3] <- rbinom((I234[i-1,j]-O_I234[i-1,j]+inf_234[i-1,j]), agert[i-1])
  age_I234[1,1:3] <- 0
  nvacc_I234[2:N_age,1] <- rbinom(age_I234[i,j],vacc_noncov[i,j]) 
  nvacc_I234[1,1] <-  0
  nvacc_I234[1:N_age,2] <- age_I234[i,j]+age_I234[i,1]-nvacc_I234[i,1]
  nvacc_I234[1:N_age,3] <- age_I234[i,j]  
  n_I234[1:N_age,1:3] <- nvacc_I234[i,j]+I234[i,j]-O_I234[i,j]+inf_234[i,j]-age_I234[i+1,j]
  
  O_I241[1:N_age,1:3] <- rbinom((I241[i,j]), (nu+deathrt[i]))
  recov_241[1:N_age,1:3] <- rbinom((O_I241[i,j]), (nu/(nu+deathrt[i])))
  age_I241[2:N_age_p1,1:3] <- rbinom((I241[i-1,j]-O_I241[i-1,j]+inf_241[i-1,j]), agert[i-1])
  age_I241[1,1:3] <- 0
  nvacc_I241[2:N_age,1] <- rbinom(age_I241[i,j],vacc_noncov[i,j]) 
  nvacc_I241[1,1] <-  0
  nvacc_I241[1:N_age,2] <- age_I241[i,j]+age_I241[i,1]-nvacc_I241[i,1]
  nvacc_I241[1:N_age,3] <- age_I241[i,j]  
  n_I241[1:N_age,1:3] <- nvacc_I241[i,j]+I241[i,j]-O_I241[i,j]+inf_241[i,j]-age_I241[i+1,j]
  
  O_I243[1:N_age,1:3] <- rbinom((I243[i,j]), (nu+deathrt[i]))
  recov_243[1:N_age,1:3] <- rbinom((O_I243[i,j]), (nu/(nu+deathrt[i])))
  age_I243[2:N_age_p1,1:3] <- rbinom((I243[i-1,j]-O_I243[i-1,j]+inf_243[i-1,j]), agert[i-1])
  age_I243[1,1:3] <- 0
  nvacc_I243[2:N_age,1] <- rbinom(age_I243[i,j],vacc_noncov[i,j]) 
  nvacc_I243[1,1] <-  0
  nvacc_I243[1:N_age,2] <- age_I243[i,j]+age_I243[i,1]-nvacc_I243[i,1]
  nvacc_I243[1:N_age,3] <- age_I243[i,j]  
  n_I243[1:N_age,1:3] <- nvacc_I243[i,j]+I243[i,j]-O_I243[i,j]+inf_243[i,j]-age_I243[i+1,j]
  
  O_I341[1:N_age,1:3] <- rbinom((I341[i,j]), (nu+deathrt[i]))
  recov_341[1:N_age,1:3] <- rbinom((O_I341[i,j]), (nu/(nu+deathrt[i])))
  age_I341[2:N_age_p1,1:3] <- rbinom((I341[i-1,j]-O_I341[i-1,j]+inf_341[i-1,j]), agert[i-1])
  age_I341[1,1:3] <- 0
  nvacc_I341[2:N_age,1] <- rbinom(age_I341[i,j],vacc_noncov[i,j]) 
  nvacc_I341[1,1] <-  0
  nvacc_I341[1:N_age,2] <- age_I341[i,j]+age_I341[i,1]-nvacc_I341[i,1]
  nvacc_I341[1:N_age,3] <- age_I341[i,j]  
  n_I341[1:N_age,1:3] <- nvacc_I341[i,j]+I341[i,j]-O_I341[i,j]+inf_341[i,j]-age_I341[i+1,j]
  
  O_I342[1:N_age,1:3] <- rbinom((I342[i,j]), (nu+deathrt[i]))
  recov_342[1:N_age,1:3] <- rbinom((O_I342[i,j]), (nu/(nu+deathrt[i])))
  age_I342[2:N_age_p1,1:3] <- rbinom((I342[i-1,j]-O_I342[i-1,j]+inf_342[i-1,j]), agert[i-1])
  age_I342[1,1:3] <- 0
  nvacc_I342[2:N_age,1] <- rbinom(age_I342[i,j],vacc_noncov[i,j]) 
  nvacc_I342[1,1] <-  0
  nvacc_I342[1:N_age,2] <- age_I342[i,j]+age_I342[i,1]-nvacc_I342[i,1]
  nvacc_I342[1:N_age,3] <- age_I342[i,j]  
  n_I342[1:N_age,1:3] <- nvacc_I342[i,j]+I342[i,j]-O_I342[i,j]+inf_342[i,j]-age_I342[i+1,j]
  
  O_R123[1:N_age,1:3] <- rbinom((R123[i,j]), (rho1234[i,j]*FOI4a[i]+deathrt[i]))
  inf_1234[1:N_age,1:3] <- rbinom((O_R123[i,j]), (rho1234[i,j]*FOI4a[i]/(rho1234[i,j]*FOI4a[i]+deathrt[i])))
  age_R123[2:N_age_p1,1:3] <- rbinom((R123[i-1,j]-O_R123[i-1,j]+recov_123[i-1,j]+recov_132[i-1,j]+recov_231[i-1,j]), agert[i-1])
  age_R123[1,1:3] <- 0
  nvacc_R123[2:N_age,1] <- rbinom(age_R123[i,j],vacc_noncov[i,j]) 
  nvacc_R123[1,1] <-  0
  nvacc_R123[1:N_age,2] <- age_R123[i,j]+age_R123[i,1]-nvacc_R123[i,1]
  nvacc_R123[1:N_age,3] <- age_R123[i,j] 
  n_R123[1:N_age,1:3] <- nvacc_R123[i,j]+R123[i,j]-O_R123[i,j]+recov_123[i,j]+recov_132[i,j]+recov_231[i,j]-age_R123[i+1,j]
  
  O_R124[1:N_age,1:3] <- rbinom((R124[i,j]), (rho1243[i,j]*FOI3a[i]+deathrt[i]))
  inf_1243[1:N_age,1:3] <- rbinom((O_R124[i,j]), (rho1243[i,j]*FOI3a[i]/(rho1243[i,j]*FOI3a[i]+deathrt[i])))
  age_R124[2:N_age_p1,1:3] <- rbinom((R124[i-1,j]-O_R124[i-1,j]+recov_124[i-1,j]+recov_142[i-1,j]+recov_241[i-1,j]), agert[i-1])
  age_R124[1,1:3] <- 0
  nvacc_R124[2:N_age,1] <- rbinom(age_R124[i,j],vacc_noncov[i,j]) 
  nvacc_R124[1,1] <-  0
  nvacc_R124[1:N_age,2] <- age_R124[i,j]+age_R124[i,1]-nvacc_R124[i,1]
  nvacc_R124[1:N_age,3] <- age_R124[i,j] 
  n_R124[1:N_age,1:3] <- nvacc_R124[i,j]+R124[i,j]-O_R124[i,j]+recov_124[i,j]+recov_142[i,j]+recov_241[i,j]-age_R124[i+1,j]
  
  O_R134[1:N_age,1:3] <- rbinom((R134[i,j]), (rho1342[i,j]*FOI2a[i]+deathrt[i]))
  inf_1342[1:N_age,1:3] <- rbinom((O_R134[i,j]), (rho1342[i,j]*FOI2a[i]/(rho1342[i,j]*FOI2a[i]+deathrt[i])))
  age_R134[2:N_age_p1,1:3] <- rbinom((R134[i-1,j]-O_R134[i-1,j]+recov_134[i-1,j]+recov_143[i-1,j]+recov_341[i-1,j]), agert[i-1])
  age_R134[1,1:3] <- 0
  nvacc_R134[2:N_age,1] <- rbinom(age_R134[i,j],vacc_noncov[i,j]) 
  nvacc_R134[1,1] <-  0
  nvacc_R134[1:N_age,2] <- age_R134[i,j]+age_R134[i,1]-nvacc_R134[i,1]
  nvacc_R134[1:N_age,3] <- age_R134[i,j] 
  n_R134[1:N_age,1:3] <- nvacc_R134[i,j]+R134[i,j]-O_R134[i,j]+recov_134[i,j]+recov_143[i,j]+recov_341[i,j]-age_R134[i+1,j]
  
  O_R234[1:N_age,1:3] <- rbinom((R234[i,j]), (rho2341[i,j]*FOI1a[i]+deathrt[i]))
  inf_2341[1:N_age,1:3] <- rbinom((O_R234[i,j]), (rho2341[i,j]*FOI1a[i]/(rho2341[i,j]*FOI1a[i]+deathrt[i])))
  age_R234[2:N_age_p1,1:3] <- rbinom((R234[i-1,j]-O_R234[i-1,j]+recov_234[i-1,j]+recov_243[i-1,j]+recov_342[i-1,j]), agert[i-1])
  age_R234[1,1:3] <- 0
  nvacc_R234[2:N_age,1] <- rbinom(age_R234[i,j],vacc_noncov[i,j]) 
  nvacc_R234[1,1] <-  0
  nvacc_R234[1:N_age,2] <- age_R234[i,j]+age_R234[i,1]-nvacc_R234[i,1]
  nvacc_R234[1:N_age,3] <- age_R234[i,j] 
  n_R234[1:N_age,1:3] <- nvacc_R234[i,j]+R234[i,j]-O_R234[i,j]+recov_234[i,j]+recov_243[i,j]+recov_342[i,j]-age_R234[i+1,j]

  O_I1234[1:N_age,1:3] <- rbinom((I1234[i,j]), (nu+deathrt[i]))
  recov_1234[1:N_age,1:3] <- rbinom((O_I1234[i,j]), (nu/(nu+deathrt[i])))
  age_I1234[2:N_age_p1,1:3] <- rbinom((I1234[i-1,j]-O_I1234[i-1,j]+inf_1234[i-1,j]), agert[i-1])
  age_I1234[1,1:3] <- 0
  nvacc_I1234[2:N_age,1] <- rbinom(age_I1234[i,j],vacc_noncov[i,j]) 
  nvacc_I1234[1,1] <-  0
  nvacc_I1234[1:N_age,2] <- age_I1234[i,j]+age_I1234[i,1]-nvacc_I1234[i,1]
  nvacc_I1234[1:N_age,3] <- age_I1234[i,j]  
  n_I1234[1:N_age,1:3] <- nvacc_I1234[i,j]+I1234[i,j]-O_I1234[i,j]+inf_1234[i,j]-age_I1234[i+1,j]
  
  O_I1243[1:N_age,1:3] <- rbinom((I1243[i,j]), (nu+deathrt[i]))
  recov_1243[1:N_age,1:3] <- rbinom((O_I1243[i,j]), (nu/(nu+deathrt[i])))
  age_I1243[2:N_age_p1,1:3] <- rbinom((I1243[i-1,j]-O_I1243[i-1,j]+inf_1243[i-1,j]), agert[i-1])
  age_I1243[1,1:3] <- 0
  nvacc_I1243[2:N_age,1] <- rbinom(age_I1243[i,j],vacc_noncov[i,j]) 
  nvacc_I1243[1,1] <-  0
  nvacc_I1243[1:N_age,2] <- age_I1243[i,j]+age_I1243[i,1]-nvacc_I1243[i,1]
  nvacc_I1243[1:N_age,3] <- age_I1243[i,j]  
  n_I1243[1:N_age,1:3] <- nvacc_I1243[i,j]+I1243[i,j]-O_I1243[i,j]+inf_1243[i,j]-age_I1243[i+1,j]
  
  O_I1342[1:N_age,1:3] <- rbinom((I1342[i,j]), (nu+deathrt[i]))
  recov_1342[1:N_age,1:3] <- rbinom((O_I1342[i,j]), (nu/(nu+deathrt[i])))
  age_I1342[2:N_age_p1,1:3] <- rbinom((I1342[i-1,j]-O_I1342[i-1,j]+inf_1342[i-1,j]), agert[i-1])
  age_I1342[1,1:3] <- 0
  nvacc_I1342[2:N_age,1] <- rbinom(age_I1342[i,j],vacc_noncov[i,j]) 
  nvacc_I1342[1,1] <-  0
  nvacc_I1342[1:N_age,2] <- age_I1342[i,j]+age_I1342[i,1]-nvacc_I1342[i,1]
  nvacc_I1342[1:N_age,3] <- age_I1342[i,j]  
  n_I1342[1:N_age,1:3] <- nvacc_I1342[i,j]+I1342[i,j]-O_I1342[i,j]+inf_1342[i,j]-age_I1342[i+1,j]
  
  O_I2341[1:N_age,1:3] <- rbinom((I2341[i,j]), (nu+deathrt[i]))
  recov_2341[1:N_age,1:3] <- rbinom((O_I2341[i,j]), (nu/(nu+deathrt[i])))
  age_I2341[2:N_age_p1,1:3] <- rbinom((I2341[i-1,j]-O_I2341[i-1,j]+inf_2341[i-1,j]), agert[i-1])
  age_I2341[1,1:3] <- 0
  nvacc_I2341[2:N_age,1] <- rbinom(age_I2341[i,j],vacc_noncov[i,j]) 
  nvacc_I2341[1,1] <-  0
  nvacc_I2341[1:N_age,2] <- age_I2341[i,j]+age_I2341[i,1]-nvacc_I2341[i,1]
  nvacc_I2341[1:N_age,3] <- age_I2341[i,j]  
  n_I2341[1:N_age,1:3] <- nvacc_I2341[i,j]+I2341[i,j]-O_I2341[i,j]+inf_2341[i,j]-age_I2341[i+1,j]
  
  O_R1234[1:N_age,1:3] <- rbinom((R1234[i,j]), (deathrt[i]))
  age_R1234[2:N_age_p1,1:3] <- rbinom((R1234[i-1,j]-O_R1234[i-1,j]+recov_1234[i-1,j]+recov_1243[i-1,j]+recov_1342[i-1,j]+recov_2341[i-1,j]), agert[i-1])
  age_R1234[1,1:3] <- 0
  nvacc_R1234[2:N_age,1] <- rbinom(age_R1234[i,j],vacc_noncov[i,j]) 
  nvacc_R1234[1,1] <-  0
  nvacc_R1234[1:N_age,2] <- age_R1234[i,j]+age_R1234[i,1]-nvacc_R1234[i,1]
  nvacc_R1234[1:N_age,3] <- age_R1234[i,j] 
  n_R1234[1:N_age,1:3] <- nvacc_R1234[i,j]+R1234[i,j]-O_R1234[i,j]+recov_1234[i,j]+recov_1243[i,j]+recov_1342[i,j]+recov_2341[i,j]-age_R1234[i+1,j]
  

  ## to make the following stochastic, need to
  ##     a. floor() the deterministic catch-up vaccination flows [easy but tedious] - also check they occur after aging does...
  ##     b. make waning stochastic - requres another binomial sampling step
  ## a and b should never clash, as only unvaccinated people are vaccinated, and only vaccinated people wane
  
  ## need to update to remove stochastic decay of vacc at some point. rates set to 0 for now
  vacc_decay <-0
  vaccS_decay <-0
  
  ncu_S[1:N_age] <- rbinom(n_S[i,1],vcu_noncov_S[i,1])
  vw_S[1:N_age] <- rbinom(n_S[i,2],vaccS_decay)
  update(S[1:N_age,1]) <- ncu_S[i]
  update(S[1:N_age,2]) <- n_S[i,2]-vw_S[i]
  update(S[1:N_age,3]) <- n_S[i,1]-ncu_S[i]+n_S[i,3]+vw_S[i]
  
  ncu_I1[1:N_age] <- rbinom(n_I1[i,1],vcu_noncov[i,1])
  vw_I1[1:N_age] <- rbinom(n_I1[i,2],vacc_decay)
  update(I1[1:N_age,1]) <- ncu_I1[i]
  update(I1[1:N_age,2]) <- n_I1[i,j]+n_I1[i,1]-ncu_I1[i]-vw_I1[i]
  update(I1[1:N_age,3]) <- n_I1[i,j]+vw_I1[i]
  ncu_I2[1:N_age] <- rbinom(n_I2[i,1],vcu_noncov[i,1])
  vw_I2[1:N_age] <- rbinom(n_I2[i,2],vacc_decay)
  update(I2[1:N_age,1]) <- ncu_I2[i]
  update(I2[1:N_age,2]) <- n_I2[i,j]+n_I2[i,1]-ncu_I2[i]-vw_I2[i]
  update(I2[1:N_age,3]) <- n_I2[i,j]+vw_I2[i]
  ncu_I3[1:N_age] <- rbinom(n_I3[i,1],vcu_noncov[i,1])
  vw_I3[1:N_age] <- rbinom(n_I3[i,2],vacc_decay)
  update(I3[1:N_age,1]) <- ncu_I3[i]
  update(I3[1:N_age,2]) <- n_I3[i,j]+n_I3[i,1]-ncu_I3[i]-vw_I3[i]
  update(I3[1:N_age,3]) <- n_I3[i,j]+vw_I3[i]
  ncu_I4[1:N_age] <- rbinom(n_I4[i,1],vcu_noncov[i,1])
  vw_I4[1:N_age] <- rbinom(n_I4[i,2],vacc_decay)
  update(I4[1:N_age,1]) <- ncu_I4[i]
  update(I4[1:N_age,2]) <- n_I4[i,j]+n_I4[i,1]-ncu_I4[i]-vw_I4[i]
  update(I4[1:N_age,3]) <- n_I4[i,j]+vw_I4[i]
  
  ncu_R1[1:N_age] <- rbinom(n_R1[i,1],vcu_noncov[i,1])
  vw_R1[1:N_age] <- rbinom(n_R1[i,2],vacc_decay)
  update(R1[1:N_age,1]) <- ncu_R1[i]
  update(R1[1:N_age,2]) <- n_R1[i,j]+n_R1[i,1]-ncu_R1[i]-vw_R1[i]
  update(R1[1:N_age,3]) <- n_R1[i,j]+vw_R1[i]
  ncu_R2[1:N_age] <- rbinom(n_R2[i,1],vcu_noncov[i,1])
  vw_R2[1:N_age] <- rbinom(n_R2[i,2],vacc_decay)
  update(R2[1:N_age,1]) <- ncu_R2[i]
  update(R2[1:N_age,2]) <- n_R2[i,j]+n_R2[i,1]-ncu_R2[i]-vw_R2[i]
  update(R2[1:N_age,3]) <- n_R2[i,j]+vw_R2[i]
  ncu_R3[1:N_age] <- rbinom(n_R3[i,1],vcu_noncov[i,1])
  vw_R3[1:N_age] <- rbinom(n_R3[i,2],vacc_decay)
  update(R3[1:N_age,1]) <- ncu_R3[i]
  update(R3[1:N_age,2]) <- n_R3[i,j]+n_R3[i,1]-ncu_R3[i]-vw_R3[i]
  update(R3[1:N_age,3]) <- n_R3[i,j]+vw_R3[i]
  ncu_R4[1:N_age] <- rbinom(n_R4[i,1],vcu_noncov[i,1])
  vw_R4[1:N_age] <- rbinom(n_R4[i,2],vacc_decay)
  update(R4[1:N_age,1]) <- ncu_R4[i]
  update(R4[1:N_age,2]) <- n_R4[i,j]+n_R4[i,1]-ncu_R4[i]-vw_R4[i]
  update(R4[1:N_age,3]) <- n_R4[i,j]+vw_R4[i]
  
  ncu_I12[1:N_age] <- rbinom(n_I12[i,1],vcu_noncov[i,1])
  vw_I12[1:N_age] <- rbinom(n_I12[i,2],vacc_decay)
  update(I12[1:N_age,1]) <- ncu_I12[i]
  update(I12[1:N_age,2]) <- n_I12[i,j]+n_I12[i,1]-ncu_I12[i]-vw_I12[i]
  update(I12[1:N_age,3]) <- n_I12[i,j]+vw_I12[i]
  ncu_I13[1:N_age] <- rbinom(n_I13[i,1],vcu_noncov[i,1])
  vw_I13[1:N_age] <- rbinom(n_I13[i,2],vacc_decay)
  update(I13[1:N_age,1]) <- ncu_I13[i]
  update(I13[1:N_age,2]) <- n_I13[i,j]+n_I13[i,1]-ncu_I13[i]-vw_I13[i]
  update(I13[1:N_age,3]) <- n_I13[i,j]+vw_I13[i]
  ncu_I14[1:N_age] <- rbinom(n_I14[i,1],vcu_noncov[i,1])
  vw_I14[1:N_age] <- rbinom(n_I14[i,2],vacc_decay)
  update(I14[1:N_age,1]) <- ncu_I14[i]
  update(I14[1:N_age,2]) <- n_I14[i,j]+n_I14[i,1]-ncu_I14[i]-vw_I14[i]
  update(I14[1:N_age,3]) <- n_I14[i,j]+vw_I14[i]
  ncu_I21[1:N_age] <- rbinom(n_I21[i,1],vcu_noncov[i,1])
  vw_I21[1:N_age] <- rbinom(n_I21[i,2],vacc_decay)
  update(I21[1:N_age,1]) <- ncu_I21[i]
  update(I21[1:N_age,2]) <- n_I21[i,j]+n_I21[i,1]-ncu_I21[i]-vw_I21[i]
  update(I21[1:N_age,3]) <- n_I21[i,j]+vw_I21[i]
  ncu_I23[1:N_age] <- rbinom(n_I23[i,1],vcu_noncov[i,1])
  vw_I23[1:N_age] <- rbinom(n_I23[i,2],vacc_decay)
  update(I23[1:N_age,1]) <- ncu_I23[i]
  update(I23[1:N_age,2]) <- n_I23[i,j]+n_I23[i,1]-ncu_I23[i]-vw_I23[i]
  update(I23[1:N_age,3]) <- n_I23[i,j]+vw_I23[i]
  ncu_I24[1:N_age] <- rbinom(n_I24[i,1],vcu_noncov[i,1])
  vw_I24[1:N_age] <- rbinom(n_I24[i,2],vacc_decay)
  update(I24[1:N_age,1]) <- ncu_I24[i]
  update(I24[1:N_age,2]) <- n_I24[i,j]+n_I24[i,1]-ncu_I24[i]-vw_I24[i]
  update(I24[1:N_age,3]) <- n_I24[i,j]+vw_I24[i]
  ncu_I31[1:N_age] <- rbinom(n_I31[i,1],vcu_noncov[i,1])
  vw_I31[1:N_age] <- rbinom(n_I31[i,2],vacc_decay)
  update(I31[1:N_age,1]) <- ncu_I31[i]
  update(I31[1:N_age,2]) <- n_I31[i,j]+n_I31[i,1]-ncu_I31[i]-vw_I31[i]
  update(I31[1:N_age,3]) <- n_I31[i,j]+vw_I31[i]
  ncu_I32[1:N_age] <- rbinom(n_I32[i,1],vcu_noncov[i,1])
  vw_I32[1:N_age] <- rbinom(n_I32[i,2],vacc_decay)
  update(I32[1:N_age,1]) <- ncu_I32[i]
  update(I32[1:N_age,2]) <- n_I32[i,j]+n_I32[i,1]-ncu_I32[i]-vw_I32[i]
  update(I32[1:N_age,3]) <- n_I32[i,j]+vw_I32[i]
  ncu_I34[1:N_age] <- rbinom(n_I34[i,1],vcu_noncov[i,1])
  vw_I34[1:N_age] <- rbinom(n_I34[i,2],vacc_decay)
  update(I34[1:N_age,1]) <- ncu_I34[i]
  update(I34[1:N_age,2]) <- n_I34[i,j]+n_I34[i,1]-ncu_I34[i]-vw_I34[i]
  update(I34[1:N_age,3]) <- n_I34[i,j]+vw_I34[i]
  ncu_I41[1:N_age] <- rbinom(n_I41[i,1],vcu_noncov[i,1])
  vw_I41[1:N_age] <- rbinom(n_I41[i,2],vacc_decay)
  update(I41[1:N_age,1]) <- ncu_I41[i]
  update(I41[1:N_age,2]) <- n_I41[i,j]+n_I41[i,1]-ncu_I41[i]-vw_I41[i]
  update(I41[1:N_age,3]) <- n_I41[i,j]+vw_I41[i]
  ncu_I42[1:N_age] <- rbinom(n_I42[i,1],vcu_noncov[i,1])
  vw_I42[1:N_age] <- rbinom(n_I42[i,2],vacc_decay)
  update(I42[1:N_age,1]) <- ncu_I42[i]
  update(I42[1:N_age,2]) <- n_I42[i,j]+n_I42[i,1]-ncu_I42[i]-vw_I42[i]
  update(I42[1:N_age,3]) <- n_I42[i,j]+vw_I42[i]
  ncu_I43[1:N_age] <- rbinom(n_I43[i,1],vcu_noncov[i,1])
  vw_I43[1:N_age] <- rbinom(n_I43[i,2],vacc_decay)
  update(I43[1:N_age,1]) <- ncu_I43[i]
  update(I43[1:N_age,2]) <- n_I43[i,j]+n_I43[i,1]-ncu_I43[i]-vw_I43[i]
  update(I43[1:N_age,3]) <- n_I43[i,j]+vw_I43[i]
  
  ncu_R12[1:N_age] <- rbinom(n_R12[i,1],vcu_noncov[i,1])
  vw_R12[1:N_age] <- rbinom(n_R12[i,2],vacc_decay)
  update(R12[1:N_age,1]) <- ncu_R12[i]
  update(R12[1:N_age,2]) <- n_R12[i,j]+n_R12[i,1]-ncu_R12[i]-vw_R12[i]
  update(R12[1:N_age,3]) <- n_R12[i,j]+vw_R12[i]
  ncu_R13[1:N_age] <- rbinom(n_R13[i,1],vcu_noncov[i,1])
  vw_R13[1:N_age] <- rbinom(n_R13[i,2],vacc_decay)
  update(R13[1:N_age,1]) <- ncu_R13[i]
  update(R13[1:N_age,2]) <- n_R13[i,j]+n_R13[i,1]-ncu_R13[i]-vw_R13[i]
  update(R13[1:N_age,3]) <- n_R13[i,j]+vw_R13[i]
  ncu_R14[1:N_age] <- rbinom(n_R14[i,1],vcu_noncov[i,1])
  vw_R14[1:N_age] <- rbinom(n_R14[i,2],vacc_decay)
  update(R14[1:N_age,1]) <- ncu_R14[i]
  update(R14[1:N_age,2]) <- n_R14[i,j]+n_R14[i,1]-ncu_R14[i]-vw_R14[i]
  update(R14[1:N_age,3]) <- n_R14[i,j]+vw_R14[i] 
  ncu_R23[1:N_age] <- rbinom(n_R23[i,1],vcu_noncov[i,1])
  vw_R23[1:N_age] <- rbinom(n_R23[i,2],vacc_decay)
  update(R23[1:N_age,1]) <- ncu_R23[i]
  update(R23[1:N_age,2]) <- n_R23[i,j]+n_R23[i,1]-ncu_R23[i]-vw_R23[i]
  update(R23[1:N_age,3]) <- n_R23[i,j]+vw_R23[i]
  ncu_R24[1:N_age] <- rbinom(n_R24[i,1],vcu_noncov[i,1])
  vw_R24[1:N_age] <- rbinom(n_R24[i,2],vacc_decay)
  update(R24[1:N_age,1]) <- ncu_R24[i]
  update(R24[1:N_age,2]) <- n_R24[i,j]+n_R24[i,1]-ncu_R24[i]-vw_R24[i]
  update(R24[1:N_age,3]) <- n_R24[i,j]+vw_R24[i] 
  ncu_R34[1:N_age] <- rbinom(n_R34[i,1],vcu_noncov[i,1])
  vw_R34[1:N_age] <- rbinom(n_R34[i,2],vacc_decay)
  update(R34[1:N_age,1]) <- ncu_R34[i]
  update(R34[1:N_age,2]) <- n_R34[i,j]+n_R34[i,1]-ncu_R34[i]-vw_R34[i]
  update(R34[1:N_age,3]) <- n_R34[i,j]+vw_R34[i]
  
  ncu_I123[1:N_age] <- rbinom(n_I123[i,1],vcu_noncov[i,1])
  vw_I123[1:N_age] <- rbinom(n_I123[i,2],vacc_decay)
  update(I123[1:N_age,1]) <- ncu_I123[i]
  update(I123[1:N_age,2]) <- n_I123[i,j]+n_I123[i,1]-ncu_I123[i]-vw_I123[i]
  update(I123[1:N_age,3]) <- n_I123[i,j]+vw_I123[i]
  ncu_I124[1:N_age] <- rbinom(n_I124[i,1],vcu_noncov[i,1])
  vw_I124[1:N_age] <- rbinom(n_I124[i,2],vacc_decay)
  update(I124[1:N_age,1]) <- ncu_I124[i]
  update(I124[1:N_age,2]) <- n_I124[i,j]+n_I124[i,1]-ncu_I124[i]-vw_I124[i]
  update(I124[1:N_age,3]) <- n_I124[i,j]+vw_I124[i]
  ncu_I132[1:N_age] <- rbinom(n_I132[i,1],vcu_noncov[i,1])
  vw_I132[1:N_age] <- rbinom(n_I132[i,2],vacc_decay)
  update(I132[1:N_age,1]) <- ncu_I132[i]
  update(I132[1:N_age,2]) <- n_I132[i,j]+n_I132[i,1]-ncu_I132[i]-vw_I132[i]
  update(I132[1:N_age,3]) <- n_I132[i,j]+vw_I132[i]
  ncu_I134[1:N_age] <- rbinom(n_I134[i,1],vcu_noncov[i,1])
  vw_I134[1:N_age] <- rbinom(n_I134[i,2],vacc_decay)
  update(I134[1:N_age,1]) <- ncu_I134[i]
  update(I134[1:N_age,2]) <- n_I134[i,j]+n_I134[i,1]-ncu_I134[i]-vw_I134[i]
  update(I134[1:N_age,3]) <- n_I134[i,j]+vw_I134[i]
  ncu_I142[1:N_age] <- rbinom(n_I142[i,1],vcu_noncov[i,1])
  vw_I142[1:N_age] <- rbinom(n_I142[i,2],vacc_decay)
  update(I142[1:N_age,1]) <- ncu_I142[i]
  update(I142[1:N_age,2]) <- n_I142[i,j]+n_I142[i,1]-ncu_I142[i]-vw_I142[i]
  update(I142[1:N_age,3]) <- n_I142[i,j]+vw_I142[i]
  ncu_I143[1:N_age] <- rbinom(n_I143[i,1],vcu_noncov[i,1])
  vw_I143[1:N_age] <- rbinom(n_I143[i,2],vacc_decay)
  update(I143[1:N_age,1]) <- ncu_I143[i]
  update(I143[1:N_age,2]) <- n_I143[i,j]+n_I143[i,1]-ncu_I143[i]-vw_I143[i]
  update(I143[1:N_age,3]) <- n_I143[i,j]+vw_I143[i]
  ncu_I231[1:N_age] <- rbinom(n_I231[i,1],vcu_noncov[i,1])
  vw_I231[1:N_age] <- rbinom(n_I231[i,2],vacc_decay)
  update(I231[1:N_age,1]) <- ncu_I231[i]
  update(I231[1:N_age,2]) <- n_I231[i,j]+n_I231[i,1]-ncu_I231[i]-vw_I231[i]
  update(I231[1:N_age,3]) <- n_I231[i,j]+vw_I231[i]
  ncu_I234[1:N_age] <- rbinom(n_I234[i,1],vcu_noncov[i,1])
  vw_I234[1:N_age] <- rbinom(n_I234[i,2],vacc_decay)
  update(I234[1:N_age,1]) <- ncu_I234[i]
  update(I234[1:N_age,2]) <- n_I234[i,j]+n_I234[i,1]-ncu_I234[i]-vw_I234[i]
  update(I234[1:N_age,3]) <- n_I234[i,j]+vw_I234[i]
  ncu_I241[1:N_age] <- rbinom(n_I241[i,1],vcu_noncov[i,1])
  vw_I241[1:N_age] <- rbinom(n_I241[i,2],vacc_decay)
  update(I241[1:N_age,1]) <- ncu_I241[i]
  update(I241[1:N_age,2]) <- n_I241[i,j]+n_I241[i,1]-ncu_I241[i]-vw_I241[i]
  update(I241[1:N_age,3]) <- n_I241[i,j]+vw_I241[i]
  ncu_I243[1:N_age] <- rbinom(n_I243[i,1],vcu_noncov[i,1])
  vw_I243[1:N_age] <- rbinom(n_I243[i,2],vacc_decay)
  update(I243[1:N_age,1]) <- ncu_I243[i]
  update(I243[1:N_age,2]) <- n_I243[i,j]+n_I243[i,1]-ncu_I243[i]-vw_I243[i]
  update(I243[1:N_age,3]) <- n_I243[i,j]+vw_I243[i]
  ncu_I341[1:N_age] <- rbinom(n_I341[i,1],vcu_noncov[i,1])
  vw_I341[1:N_age] <- rbinom(n_I341[i,2],vacc_decay)
  update(I341[1:N_age,1]) <- ncu_I341[i]
  update(I341[1:N_age,2]) <- n_I341[i,j]+n_I341[i,1]-ncu_I341[i]-vw_I341[i]
  update(I341[1:N_age,3]) <- n_I341[i,j]+vw_I341[i]
  ncu_I342[1:N_age] <- rbinom(n_I342[i,1],vcu_noncov[i,1])
  vw_I342[1:N_age] <- rbinom(n_I342[i,2],vacc_decay)
  update(I342[1:N_age,1]) <- ncu_I342[i]
  update(I342[1:N_age,2]) <- n_I342[i,j]+n_I342[i,1]-ncu_I342[i]-vw_I342[i]
  update(I342[1:N_age,3]) <- n_I342[i,j]+vw_I342[i]
  
  ncu_R123[1:N_age] <- rbinom(n_R123[i,1],vcu_noncov[i,1])
  vw_R123[1:N_age] <- rbinom(n_R123[i,2],vacc_decay)
  update(R123[1:N_age,1]) <- ncu_R123[i]
  update(R123[1:N_age,2]) <- n_R123[i,j]+n_R123[i,1]-ncu_R123[i]-vw_R123[i]
  update(R123[1:N_age,3]) <- n_R123[i,j]+vw_R123[i]
  ncu_R124[1:N_age] <- rbinom(n_R124[i,1],vcu_noncov[i,1])
  vw_R124[1:N_age] <- rbinom(n_R124[i,2],vacc_decay)
  update(R124[1:N_age,1]) <- ncu_R124[i]
  update(R124[1:N_age,2]) <- n_R124[i,j]+n_R124[i,1]-ncu_R124[i]-vw_R124[i]
  update(R124[1:N_age,3]) <- n_R124[i,j]+vw_R124[i]
  ncu_R134[1:N_age] <- rbinom(n_R134[i,1],vcu_noncov[i,1])
  vw_R134[1:N_age] <- rbinom(n_R134[i,2],vacc_decay)
  update(R134[1:N_age,1]) <- ncu_R134[i]
  update(R134[1:N_age,2]) <- n_R134[i,j]+n_R134[i,1]-ncu_R134[i]-vw_R134[i]
  update(R134[1:N_age,3]) <- n_R134[i,j]+vw_R134[i] 
  ncu_R234[1:N_age] <- rbinom(n_R234[i,1],vcu_noncov[i,1])
  vw_R234[1:N_age] <- rbinom(n_R234[i,2],vacc_decay)
  update(R234[1:N_age,1]) <- ncu_R234[i]
  update(R234[1:N_age,2]) <- n_R234[i,j]+n_R234[i,1]-ncu_R234[i]-vw_R234[i]
  update(R234[1:N_age,3]) <- n_R234[i,j]+vw_R234[i]
  
  ncu_I1234[1:N_age] <- rbinom(n_I1234[i,1],vcu_noncov[i,1])
  vw_I1234[1:N_age] <- rbinom(n_I1234[i,2],vacc_decay)
  update(I1234[1:N_age,1]) <- ncu_I1234[i]
  update(I1234[1:N_age,2]) <- n_I1234[i,j]+n_I1234[i,1]-ncu_I1234[i]-vw_I1234[i]
  update(I1234[1:N_age,3]) <- n_I1234[i,j]+vw_I1234[i]
  ncu_I1243[1:N_age] <- rbinom(n_I1243[i,1],vcu_noncov[i,1])
  vw_I1243[1:N_age] <- rbinom(n_I1243[i,2],vacc_decay)
  update(I1243[1:N_age,1]) <- ncu_I1243[i]
  update(I1243[1:N_age,2]) <- n_I1243[i,j]+n_I1243[i,1]-ncu_I1243[i]-vw_I1243[i]
  update(I1243[1:N_age,3]) <- n_I1243[i,j]+vw_I1243[i]
  ncu_I1342[1:N_age] <- rbinom(n_I1342[i,1],vcu_noncov[i,1])
  vw_I1342[1:N_age] <- rbinom(n_I1342[i,2],vacc_decay)
  update(I1342[1:N_age,1]) <- ncu_I1342[i]
  update(I1342[1:N_age,2]) <- n_I1342[i,j]+n_I1342[i,1]-ncu_I1342[i]-vw_I1342[i]
  update(I1342[1:N_age,3]) <- n_I1342[i,j]+vw_I1342[i]
  ncu_I2341[1:N_age] <- rbinom(n_I2341[i,1],vcu_noncov[i,1])
  vw_I2341[1:N_age] <- rbinom(n_I2341[i,2],vacc_decay)
  update(I2341[1:N_age,1]) <- ncu_I2341[i]
  update(I2341[1:N_age,2]) <- n_I2341[i,j]+n_I2341[i,1]-ncu_I2341[i]-vw_I2341[i]
  update(I2341[1:N_age,3]) <- n_I2341[i,j]+vw_I2341[i]
  
  ncu_R1234[1:N_age] <- rbinom(n_R1234[i,1],vcu_noncov[i,1])
  vw_R1234[1:N_age] <- rbinom(n_R1234[i,2],vacc_decay)
  update(R1234[1:N_age,1]) <- ncu_R1234[i]
  update(R1234[1:N_age,2]) <- n_R1234[i,j]+n_R1234[i,1]-ncu_R1234[i]-vw_R1234[i]
  update(R1234[1:N_age,3]) <- n_R1234[i,j]+vw_R1234[i]
  
  
## array dimension
  

  dim(age_removal_d) <- c(151,N_age_p1)
  dim(age_rate_d) <-  c(151,N_age_p1)
  dim(pop_size_d) <-  N_age_p1
  dim(births_d) <- c(151,2)
  dim(life_expec_d) <-  c(151,N_age_p1)
  dim(coverage_d) <- c(77,2)
  
  dim(nc0) <- c(4,3)
  dim(hs) <- 3
  dim(hl) <- 3
  dim(ts) <- 3
  dim(nc50_inf) <- c(4,3)
  dim(nc50_dis) <- c(4,3)
  dim(nc50_sdis) <- c(4,3)
  dim(L_dis) <- c(4,3)
  dim(L_sdis) <- c(4,3)
  dim(ws) <- 4
  dim(nc50_age) <- N_age
  
  dim(dis_pri) <- 4
  dim(dis_sec) <- 4
  dim(dis_tert) <- 4
  dim(dis_quart) <- 4
  
  dim(sdis_pri) <- 4
  dim(sdis_sec) <- 4
  dim(sdis_tert) <- 4
  dim(sdis_quart) <- 4
  
  dim(agec) <- N_age
  
  dim(death) <- N_age
  dim(ageb) <- N_age_p1
  dim(mean_age) <- N_age
  dim(init_life_expec) <- N_age
  dim(deathrt) <- N_age
  dim(cur_age_rate) <- N_age
  dim(agert) <- N_age
  dim(N_init_age0) <- N_age
  dim(N_init_age) <- N_age
  dim(life_expec) <- N_age
  dim(FOIas) <- N_age
  dim(Ntotal_out) <- N_age
  dim(suscinitpop) <- N_age
  dim(vacc_cu_age_weight) <- N_age_p1
  dim(vacc_noncov) <- c(N_age_p1,3)
  dim(vacc_noncov_S) <- c(N_age_p1,3)
  dim(vcu_noncov) <- c(N_age,3)
  dim(vcu_noncov_S) <- c(N_age,3)
  
  dim(nc) <- c(4,3,N_age)
  dim(RR_inf_vc) <- c(4,3,N_age) 
  dim(RR_dis_vc) <- c(4,3,N_age)  
  dim(RR_sdis_vc) <- c(4,3,N_age)  
  
  dim(nc_cu) <- c(4,3)
  dim(RR_inf_cu) <- c(4,3,N_age)  
  dim(RR_dis_cu) <- c(4,3,N_age)  
  dim(RR_sdis_cu) <- c(4,3,N_age)  
  
  dim(RR_inf) <- c(4,3,N_age)  
  dim(RR_dis) <- c(4,3,N_age)  
  dim(RR_sdis) <- c(4,3,N_age)  
  
  dim(out_update_switch) <- NYO
  
  dim(VE_seroneg) <- NYO
  dim(VE_seropos_mono) <- NYO
  dim(VE_seropos_multi) <- NYO
  
  dim(out_prop_seroneg_at_9) <- NYO
  dim(out_prop_seroneg_at_vca) <- NYO
  
  dim(disease_sero) <- c(N_age, 3, 4)
  dim(disease_sero_vacc_pri) <- c(N_age, 4)
  dim(dis_sero_unvacc) <- 4
  dim(dis_sero_vacc) <- 4
  dim(dis_sero_pop) <- 4
  dim(dis_sero_vacc_neg) <- 4
  dim(dis_sero_vacc_pos) <- 4
  dim(dis_sero_vacc_pri) <- 4
  dim(dis_sero_vacc_secp) <- 4
  
  dim(out_dis_all_pop) <- NYO
  dim(out_dis_sero_pop) <- c(NYO,4)
  dim(out_dis_all_unvacc) <- NYO
  dim(out_dis_sero_unvacc) <- c(NYO,4)
  dim(out_dis_all_vacc) <- NYO
  dim(out_dis_all_vacc_neg) <- NYO
  dim(out_dis_all_vacc_pos) <- NYO
  dim(out_dis_all_vacc_pri) <- NYO
  dim(out_dis_all_vacc_secp) <- NYO
  dim(out_dis_sero_vacc) <- c(NYO,4)
  dim(out_dis_sero_vacc_neg) <- c(NYO,4)
  dim(out_dis_sero_vacc_pos) <- c(NYO,4)
  dim(out_dis_sero_vacc_pri) <- c(NYO,4)
  dim(out_dis_sero_vacc_secp) <- c(NYO,4)
  dim(out_disc_dis_all_pop) <- NYO
  dim(out_disc_dis_all_vacc) <- NYO
  
  dim(sdisease_sero) <- c(N_age, 3, 4)
  dim(sdisease_sero_vacc_pri) <- c(N_age, 4)
  dim(sdis_sero_unvacc) <- 4
  dim(sdis_sero_vacc) <- 4
  dim(sdis_sero_pop) <- 4
  dim(sdis_sero_vacc_neg) <- 4
  dim(sdis_sero_vacc_pos) <- 4
  dim(sdis_sero_vacc_pri) <- 4
  dim(sdis_sero_vacc_secp) <- 4
  
  dim(out_sdis_all_pop) <- NYO
  dim(out_sdis_sero_pop) <- c(NYO,4)
  dim(out_sdis_all_unvacc) <- NYO
  dim(out_sdis_sero_unvacc) <- c(NYO,4)
  dim(out_sdis_all_vacc) <- NYO
  dim(out_sdis_all_vacc_neg) <- NYO
  dim(out_sdis_all_vacc_pos) <- NYO
  dim(out_sdis_all_vacc_pri) <- NYO
  dim(out_sdis_all_vacc_secp) <- NYO
  dim(out_sdis_sero_vacc) <- c(NYO,4)
  dim(out_sdis_sero_vacc_neg) <- c(NYO,4)
  dim(out_sdis_sero_vacc_pos) <- c(NYO,4)
  dim(out_sdis_sero_vacc_pri) <- c(NYO,4)
  dim(out_sdis_sero_vacc_secp) <- c(NYO,4)
  dim(out_disc_sdis_all_pop) <- NYO
  dim(out_disc_sdis_all_vacc) <- NYO
  
  dim(yll_sero) <- c(N_age, 3, 4)
  dim(yll_sero_vacc_pri_age) <- c(N_age, 4)
  dim(disc_yll_all) <- c(N_age, 3)
  dim(yll_sero_unvacc) <- 4
  dim(yll_sero_vacc) <- 4
  dim(yll_sero_pop) <- 4
  dim(yll_sero_vacc_neg) <- 4
  dim(yll_sero_vacc_pos) <- 4
  dim(yll_sero_vacc_pri) <- 4
  dim(yll_sero_vacc_secp) <- 4
  
  dim(out_yll_all_pop) <- NYO
  dim(out_yll_sero_pop) <- c(NYO,4)
  dim(out_yll_all_unvacc) <- NYO
  dim(out_yll_sero_unvacc) <- c(NYO,4)
  dim(out_yll_all_vacc) <- NYO
  dim(out_yll_all_vacc_neg) <- NYO
  dim(out_yll_all_vacc_pos) <- NYO
  dim(out_yll_all_vacc_pri) <- NYO
  dim(out_yll_all_vacc_secp) <- NYO
  dim(out_yll_sero_vacc) <- c(NYO,4)
  dim(out_yll_sero_vacc_neg) <- c(NYO,4)
  dim(out_yll_sero_vacc_pos) <- c(NYO,4)
  dim(out_yll_sero_vacc_pri) <- c(NYO,4)
  dim(out_yll_sero_vacc_secp) <- c(NYO,4)
  dim(out_disc_yll_all_pop) <- NYO
  dim(out_disc_yll_all_vacc) <- NYO
  
  dim(out_nvacc_all_pop) <- NYO
  dim(out_nvacc_all_neg) <- NYO
  dim(out_nvacc_all_pos) <- NYO
  dim(out_disc_nvacc_all_pop) <- NYO
  
  dim(num_child_vacc_age) <- N_age
  dim(num_cu_vacc_age) <- N_age
  dim(num_child_vacc_age_neg) <- N_age
  dim(num_cu_vacc_age_neg) <- N_age
  
  dim(vacc_ct) <- N_age
  dim(time_from_last_dose_base) <- N_age
  dim(time_from_last_dose) <- N_age
  
  dim(phi_scale) <- 4
  dim(phi_pri) <- 4
  dim(phi_sec) <- 4
  dim(phi_tert) <- 4
  dim(phi_quart) <- 4
  
  dim(rho1) <- c(N_age,3)
  dim(rho2) <- c(N_age,3)
  dim(rho3) <- c(N_age,3)
  dim(rho4) <- c(N_age,3)
  
  dim(rho21) <- c(N_age,3)
  dim(rho31) <- c(N_age,3)
  dim(rho41) <- c(N_age,3)
  dim(rho12) <- c(N_age,3)
  dim(rho32) <- c(N_age,3)
  dim(rho42) <- c(N_age,3)
  dim(rho13) <- c(N_age,3)
  dim(rho23) <- c(N_age,3)
  dim(rho43) <- c(N_age,3)
  dim(rho14) <- c(N_age,3)
  dim(rho24) <- c(N_age,3)
  dim(rho34) <- c(N_age,3)
  
  dim(rho231) <- c(N_age,3)
  dim(rho241) <- c(N_age,3)
  dim(rho341) <- c(N_age,3)
  dim(rho132) <- c(N_age,3)
  dim(rho142) <- c(N_age,3)
  dim(rho342) <- c(N_age,3)
  dim(rho123) <- c(N_age,3)
  dim(rho143) <- c(N_age,3)
  dim(rho243) <- c(N_age,3)
  dim(rho124) <- c(N_age,3)
  dim(rho134) <- c(N_age,3)
  dim(rho234) <- c(N_age,3)
  
  dim(rho2341) <- c(N_age,3)
  dim(rho1342) <- c(N_age,3)
  dim(rho1243) <- c(N_age,3)
  dim(rho1234) <- c(N_age,3)
  
  
  dim(S) <- c(N_age,3)
  
  dim(I1) <- c(N_age,3)
  dim(I2) <- c(N_age,3)
  dim(I3) <- c(N_age,3)
  dim(I4) <- c(N_age,3)
  dim(I21) <- c(N_age,3)
  dim(I31) <- c(N_age,3)
  dim(I41) <- c(N_age,3)
  dim(I12) <- c(N_age,3)
  dim(I32) <- c(N_age,3)
  dim(I42) <- c(N_age,3)
  dim(I13) <- c(N_age,3)
  dim(I23) <- c(N_age,3)
  dim(I43) <- c(N_age,3)
  dim(I14) <- c(N_age,3)
  dim(I24) <- c(N_age,3)
  dim(I34) <- c(N_age,3)
  dim(I231) <- c(N_age,3)
  dim(I241) <- c(N_age,3)
  dim(I341) <- c(N_age,3)
  dim(I132) <- c(N_age,3)
  dim(I142) <- c(N_age,3)
  dim(I342) <- c(N_age,3)
  dim(I123) <- c(N_age,3)
  dim(I143) <- c(N_age,3)
  dim(I243) <- c(N_age,3)
  dim(I124) <- c(N_age,3)
  dim(I134) <- c(N_age,3)
  dim(I234) <- c(N_age,3)
  dim(I2341) <- c(N_age,3)
  dim(I1342) <- c(N_age,3)
  dim(I1243) <- c(N_age,3)
  dim(I1234) <- c(N_age,3)
  
  dim(R1) <- c(N_age,3)
  dim(R2) <- c(N_age,3)
  dim(R3) <- c(N_age,3)
  dim(R4) <- c(N_age,3)
  dim(R12) <- c(N_age,3)
  dim(R13) <- c(N_age,3)
  dim(R14) <- c(N_age,3)
  dim(R23) <- c(N_age,3)
  dim(R24) <- c(N_age,3)
  dim(R34) <- c(N_age,3)
  dim(R123) <- c(N_age,3)
  dim(R124) <- c(N_age,3)
  dim(R134) <- c(N_age,3)
  dim(R234) <- c(N_age,3)
  dim(R1234) <- c(N_age,3)
  
  dim(Ntotal) <- c(N_age,3)
  dim(Ntotal_nv) <- N_age
  dim(Ntotal_nvS) <- N_age
  dim(Ntotal_v) <- N_age
  dim(Ntotal_vS) <- N_age
  dim(Ntotal_all) <- N_age
  
  dim(Y1) <- c(N_age,3)
  dim(Y2) <- c(N_age,3)
  dim(Y3) <- c(N_age,3)
  dim(Y4) <- c(N_age,3)
  
  dim(FOI1a) <- N_age
  dim(FOI2a) <- N_age
  dim(FOI3a) <- N_age
  dim(FOI4a) <- N_age
  
  dim(infection_pri) <- c(N_age,3)
  dim(infection_sec) <- c(N_age,3)
  dim(infection_tq) <- c(N_age,3)

  dim(cum_infection_pri) <- N_age
  dim(cum_infection_sec) <- N_age
  dim(cum_infection_tq) <- N_age
  
  dim(cum_infection_priV) <- N_age
  dim(cum_infection_secV) <- N_age
  dim(cum_infection_tqV) <- N_age
  
  
  dim(nveai) <- 3
  dim(nnveai) <- 3
  
  dim(inf_1) <- c(N_age,3)
  dim(inf_2) <- c(N_age,3)
  dim(inf_3) <- c(N_age,3)
  dim(inf_4) <- c(N_age,3)
  dim(inf_21) <- c(N_age,3)
  dim(inf_31) <- c(N_age,3)
  dim(inf_41) <- c(N_age,3)
  dim(inf_12) <- c(N_age,3)
  dim(inf_32) <- c(N_age,3)
  dim(inf_42) <- c(N_age,3)
  dim(inf_13) <- c(N_age,3)
  dim(inf_23) <- c(N_age,3)
  dim(inf_43) <- c(N_age,3)
  dim(inf_14) <- c(N_age,3)
  dim(inf_24) <- c(N_age,3)
  dim(inf_34) <- c(N_age,3)
  dim(inf_231) <- c(N_age,3)
  dim(inf_241) <- c(N_age,3)
  dim(inf_341) <- c(N_age,3)
  dim(inf_132) <- c(N_age,3)
  dim(inf_142) <- c(N_age,3)
  dim(inf_342) <- c(N_age,3)
  dim(inf_123) <- c(N_age,3)
  dim(inf_143) <- c(N_age,3)
  dim(inf_243) <- c(N_age,3)
  dim(inf_124) <- c(N_age,3)
  dim(inf_134) <- c(N_age,3)
  dim(inf_234) <- c(N_age,3)
  dim(inf_2341) <- c(N_age,3)
  dim(inf_1342) <- c(N_age,3)
  dim(inf_1243) <- c(N_age,3)
  dim(inf_1234) <- c(N_age,3)
  
  dim(infV_1) <- c(N_age,3)
  dim(infV_2) <- c(N_age,3)
  dim(infV_3) <- c(N_age,3)
  dim(infV_4) <- c(N_age,3)  
  
  dim(recov_1) <- c(N_age,3)
  dim(recov_2) <- c(N_age,3)
  dim(recov_3) <- c(N_age,3)
  dim(recov_4) <- c(N_age,3)
  dim(recov_21) <- c(N_age,3)
  dim(recov_31) <- c(N_age,3)
  dim(recov_41) <- c(N_age,3)
  dim(recov_12) <- c(N_age,3)
  dim(recov_32) <- c(N_age,3)
  dim(recov_42) <- c(N_age,3)
  dim(recov_13) <- c(N_age,3)
  dim(recov_23) <- c(N_age,3)
  dim(recov_43) <- c(N_age,3)
  dim(recov_14) <- c(N_age,3)
  dim(recov_24) <- c(N_age,3)
  dim(recov_34) <- c(N_age,3)
  dim(recov_231) <- c(N_age,3)
  dim(recov_241) <- c(N_age,3)
  dim(recov_341) <- c(N_age,3)
  dim(recov_132) <- c(N_age,3)
  dim(recov_142) <- c(N_age,3)
  dim(recov_342) <- c(N_age,3)
  dim(recov_123) <- c(N_age,3)
  dim(recov_143) <- c(N_age,3)
  dim(recov_243) <- c(N_age,3)
  dim(recov_124) <- c(N_age,3)
  dim(recov_134) <- c(N_age,3)
  dim(recov_234) <- c(N_age,3)
  dim(recov_2341) <- c(N_age,3)
  dim(recov_1342) <- c(N_age,3)
  dim(recov_1243) <- c(N_age,3)
  dim(recov_1234) <- c(N_age,3)
  
  dim(age_S) <- c(N_age_p1,3)
  
  dim(age_I1) <- c(N_age_p1,3)
  dim(age_I2) <- c(N_age_p1,3)
  dim(age_I3) <- c(N_age_p1,3)
  dim(age_I4) <- c(N_age_p1,3)
  dim(age_I21) <- c(N_age_p1,3)
  dim(age_I31) <- c(N_age_p1,3)
  dim(age_I41) <- c(N_age_p1,3)
  dim(age_I12) <- c(N_age_p1,3)
  dim(age_I32) <- c(N_age_p1,3)
  dim(age_I42) <- c(N_age_p1,3)
  dim(age_I13) <- c(N_age_p1,3)
  dim(age_I23) <- c(N_age_p1,3)
  dim(age_I43) <- c(N_age_p1,3)
  dim(age_I14) <- c(N_age_p1,3)
  dim(age_I24) <- c(N_age_p1,3)
  dim(age_I34) <- c(N_age_p1,3)
  dim(age_I231) <- c(N_age_p1,3)
  dim(age_I241) <- c(N_age_p1,3)
  dim(age_I341) <- c(N_age_p1,3)
  dim(age_I132) <- c(N_age_p1,3)
  dim(age_I142) <- c(N_age_p1,3)
  dim(age_I342) <- c(N_age_p1,3)
  dim(age_I123) <- c(N_age_p1,3)
  dim(age_I143) <- c(N_age_p1,3)
  dim(age_I243) <- c(N_age_p1,3)
  dim(age_I124) <- c(N_age_p1,3)
  dim(age_I134) <- c(N_age_p1,3)
  dim(age_I234) <- c(N_age_p1,3)
  dim(age_I2341) <- c(N_age_p1,3)
  dim(age_I1342) <- c(N_age_p1,3)
  dim(age_I1243) <- c(N_age_p1,3)
  dim(age_I1234) <- c(N_age_p1,3)
  
  dim(age_R1) <- c(N_age_p1,3)
  dim(age_R2) <- c(N_age_p1,3)
  dim(age_R3) <- c(N_age_p1,3)
  dim(age_R4) <- c(N_age_p1,3)
  dim(age_R12) <- c(N_age_p1,3)
  dim(age_R13) <- c(N_age_p1,3)
  dim(age_R14) <- c(N_age_p1,3)
  dim(age_R23) <- c(N_age_p1,3)
  dim(age_R24) <- c(N_age_p1,3)
  dim(age_R34) <- c(N_age_p1,3)
  dim(age_R123) <- c(N_age_p1,3)
  dim(age_R124) <- c(N_age_p1,3)
  dim(age_R134) <- c(N_age_p1,3)
  dim(age_R234) <- c(N_age_p1,3)
  dim(age_R1234) <- c(N_age_p1,3)
  
   
  dim(O_S) <- c(N_age,3)
  
  dim(O_I1) <- c(N_age,3)
  dim(O_I2) <- c(N_age,3)
  dim(O_I3) <- c(N_age,3)
  dim(O_I4) <- c(N_age,3)
  dim(O_I21) <- c(N_age,3)
  dim(O_I31) <- c(N_age,3)
  dim(O_I41) <- c(N_age,3)
  dim(O_I12) <- c(N_age,3)
  dim(O_I32) <- c(N_age,3)
  dim(O_I42) <- c(N_age,3)
  dim(O_I13) <- c(N_age,3)
  dim(O_I23) <- c(N_age,3)
  dim(O_I43) <- c(N_age,3)
  dim(O_I14) <- c(N_age,3)
  dim(O_I24) <- c(N_age,3)
  dim(O_I34) <- c(N_age,3)
  dim(O_I231) <- c(N_age,3)
  dim(O_I241) <- c(N_age,3)
  dim(O_I341) <- c(N_age,3)
  dim(O_I132) <- c(N_age,3)
  dim(O_I142) <- c(N_age,3)
  dim(O_I342) <- c(N_age,3)
  dim(O_I123) <- c(N_age,3)
  dim(O_I143) <- c(N_age,3)
  dim(O_I243) <- c(N_age,3)
  dim(O_I124) <- c(N_age,3)
  dim(O_I134) <- c(N_age,3)
  dim(O_I234) <- c(N_age,3)
  dim(O_I2341) <- c(N_age,3)
  dim(O_I1342) <- c(N_age,3)
  dim(O_I1243) <- c(N_age,3)
  dim(O_I1234) <- c(N_age,3)
  
  dim(O_R1) <- c(N_age,3)
  dim(O_R2) <- c(N_age,3)
  dim(O_R3) <- c(N_age,3)
  dim(O_R4) <- c(N_age,3)
  dim(O_R12) <- c(N_age,3)
  dim(O_R13) <- c(N_age,3)
  dim(O_R14) <- c(N_age,3)
  dim(O_R23) <- c(N_age,3)
  dim(O_R24) <- c(N_age,3)
  dim(O_R34) <- c(N_age,3)
  dim(O_R123) <- c(N_age,3)
  dim(O_R124) <- c(N_age,3)
  dim(O_R134) <- c(N_age,3)
  dim(O_R234) <- c(N_age,3)
  dim(O_R1234) <- c(N_age,3)
  
  
  dim(nvacc_S) <- c(N_age,3)
  
  dim(nvacc_I1) <- c(N_age,3)
  dim(nvacc_I2) <- c(N_age,3)
  dim(nvacc_I3) <- c(N_age,3)
  dim(nvacc_I4) <- c(N_age,3)
  dim(nvacc_I21) <- c(N_age,3)
  dim(nvacc_I31) <- c(N_age,3)
  dim(nvacc_I41) <- c(N_age,3)
  dim(nvacc_I12) <- c(N_age,3)
  dim(nvacc_I32) <- c(N_age,3)
  dim(nvacc_I42) <- c(N_age,3)
  dim(nvacc_I13) <- c(N_age,3)
  dim(nvacc_I23) <- c(N_age,3)
  dim(nvacc_I43) <- c(N_age,3)
  dim(nvacc_I14) <- c(N_age,3)
  dim(nvacc_I24) <- c(N_age,3)
  dim(nvacc_I34) <- c(N_age,3)
  dim(nvacc_I231) <- c(N_age,3)
  dim(nvacc_I241) <- c(N_age,3)
  dim(nvacc_I341) <- c(N_age,3)
  dim(nvacc_I132) <- c(N_age,3)
  dim(nvacc_I142) <- c(N_age,3)
  dim(nvacc_I342) <- c(N_age,3)
  dim(nvacc_I123) <- c(N_age,3)
  dim(nvacc_I143) <- c(N_age,3)
  dim(nvacc_I243) <- c(N_age,3)
  dim(nvacc_I124) <- c(N_age,3)
  dim(nvacc_I134) <- c(N_age,3)
  dim(nvacc_I234) <- c(N_age,3)
  dim(nvacc_I2341) <- c(N_age,3)
  dim(nvacc_I1342) <- c(N_age,3)
  dim(nvacc_I1243) <- c(N_age,3)
  dim(nvacc_I1234) <- c(N_age,3)
  
  dim(nvacc_R1) <- c(N_age,3)
  dim(nvacc_R2) <- c(N_age,3)
  dim(nvacc_R3) <- c(N_age,3)
  dim(nvacc_R4) <- c(N_age,3)
  dim(nvacc_R12) <- c(N_age,3)
  dim(nvacc_R13) <- c(N_age,3)
  dim(nvacc_R14) <- c(N_age,3)
  dim(nvacc_R23) <- c(N_age,3)
  dim(nvacc_R24) <- c(N_age,3)
  dim(nvacc_R34) <- c(N_age,3)
  dim(nvacc_R123) <- c(N_age,3)
  dim(nvacc_R124) <- c(N_age,3)
  dim(nvacc_R134) <- c(N_age,3)
  dim(nvacc_R234) <- c(N_age,3)
  dim(nvacc_R1234) <- c(N_age,3)
  
  dim(n_S) <- c(N_age,3)
  
  dim(n_I1) <- c(N_age,3)
  dim(n_I2) <- c(N_age,3)
  dim(n_I3) <- c(N_age,3)
  dim(n_I4) <- c(N_age,3)
  dim(n_I21) <- c(N_age,3)
  dim(n_I31) <- c(N_age,3)
  dim(n_I41) <- c(N_age,3)
  dim(n_I12) <- c(N_age,3)
  dim(n_I32) <- c(N_age,3)
  dim(n_I42) <- c(N_age,3)
  dim(n_I13) <- c(N_age,3)
  dim(n_I23) <- c(N_age,3)
  dim(n_I43) <- c(N_age,3)
  dim(n_I14) <- c(N_age,3)
  dim(n_I24) <- c(N_age,3)
  dim(n_I34) <- c(N_age,3)
  dim(n_I231) <- c(N_age,3)
  dim(n_I241) <- c(N_age,3)
  dim(n_I341) <- c(N_age,3)
  dim(n_I132) <- c(N_age,3)
  dim(n_I142) <- c(N_age,3)
  dim(n_I342) <- c(N_age,3)
  dim(n_I123) <- c(N_age,3)
  dim(n_I143) <- c(N_age,3)
  dim(n_I243) <- c(N_age,3)
  dim(n_I124) <- c(N_age,3)
  dim(n_I134) <- c(N_age,3)
  dim(n_I234) <- c(N_age,3)
  dim(n_I2341) <- c(N_age,3)
  dim(n_I1342) <- c(N_age,3)
  dim(n_I1243) <- c(N_age,3)
  dim(n_I1234) <- c(N_age,3)
  
  dim(n_R1) <- c(N_age,3)
  dim(n_R2) <- c(N_age,3)
  dim(n_R3) <- c(N_age,3)
  dim(n_R4) <- c(N_age,3)
  dim(n_R12) <- c(N_age,3)
  dim(n_R13) <- c(N_age,3)
  dim(n_R14) <- c(N_age,3)
  dim(n_R23) <- c(N_age,3)
  dim(n_R24) <- c(N_age,3)
  dim(n_R34) <- c(N_age,3)
  dim(n_R123) <- c(N_age,3)
  dim(n_R124) <- c(N_age,3)
  dim(n_R134) <- c(N_age,3)
  dim(n_R234) <- c(N_age,3)
  dim(n_R1234) <- c(N_age,3)
  
  
  dim(ncu_S) <- N_age
  
  dim(ncu_I1) <- N_age
  dim(ncu_I2) <- N_age
  dim(ncu_I3) <- N_age
  dim(ncu_I4) <- N_age
  dim(ncu_I21) <- N_age
  dim(ncu_I31) <- N_age
  dim(ncu_I41) <- N_age
  dim(ncu_I12) <- N_age
  dim(ncu_I32) <- N_age
  dim(ncu_I42) <- N_age
  dim(ncu_I13) <- N_age
  dim(ncu_I23) <- N_age
  dim(ncu_I43) <- N_age
  dim(ncu_I14) <- N_age
  dim(ncu_I24) <- N_age
  dim(ncu_I34) <- N_age
  dim(ncu_I231) <- N_age
  dim(ncu_I241) <- N_age
  dim(ncu_I341) <- N_age
  dim(ncu_I132) <- N_age
  dim(ncu_I142) <- N_age
  dim(ncu_I342) <- N_age
  dim(ncu_I123) <- N_age
  dim(ncu_I143) <- N_age
  dim(ncu_I243) <- N_age
  dim(ncu_I124) <- N_age
  dim(ncu_I134) <- N_age
  dim(ncu_I234) <- N_age
  dim(ncu_I2341) <- N_age
  dim(ncu_I1342) <- N_age
  dim(ncu_I1243) <- N_age
  dim(ncu_I1234) <- N_age
  
  dim(ncu_R1) <- N_age
  dim(ncu_R2) <- N_age
  dim(ncu_R3) <- N_age
  dim(ncu_R4) <- N_age
  dim(ncu_R12) <- N_age
  dim(ncu_R13) <- N_age
  dim(ncu_R14) <- N_age
  dim(ncu_R23) <- N_age
  dim(ncu_R24) <- N_age
  dim(ncu_R34) <- N_age
  dim(ncu_R123) <- N_age
  dim(ncu_R124) <- N_age
  dim(ncu_R134) <- N_age
  dim(ncu_R234) <- N_age
  dim(ncu_R1234) <- N_age
  
  
  dim(vw_S) <- N_age
  
  dim(vw_I1) <- N_age
  dim(vw_I2) <- N_age
  dim(vw_I3) <- N_age
  dim(vw_I4) <- N_age
  dim(vw_I21) <- N_age
  dim(vw_I31) <- N_age
  dim(vw_I41) <- N_age
  dim(vw_I12) <- N_age
  dim(vw_I32) <- N_age
  dim(vw_I42) <- N_age
  dim(vw_I13) <- N_age
  dim(vw_I23) <- N_age
  dim(vw_I43) <- N_age
  dim(vw_I14) <- N_age
  dim(vw_I24) <- N_age
  dim(vw_I34) <- N_age
  dim(vw_I231) <- N_age
  dim(vw_I241) <- N_age
  dim(vw_I341) <- N_age
  dim(vw_I132) <- N_age
  dim(vw_I142) <- N_age
  dim(vw_I342) <- N_age
  dim(vw_I123) <- N_age
  dim(vw_I143) <- N_age
  dim(vw_I243) <- N_age
  dim(vw_I124) <- N_age
  dim(vw_I134) <- N_age
  dim(vw_I234) <- N_age
  dim(vw_I2341) <- N_age
  dim(vw_I1342) <- N_age
  dim(vw_I1243) <- N_age
  dim(vw_I1234) <- N_age
  
  dim(vw_R1) <- N_age
  dim(vw_R2) <- N_age
  dim(vw_R3) <- N_age
  dim(vw_R4) <- N_age
  dim(vw_R12) <- N_age
  dim(vw_R13) <- N_age
  dim(vw_R14) <- N_age
  dim(vw_R23) <- N_age
  dim(vw_R24) <- N_age
  dim(vw_R34) <- N_age
  dim(vw_R123) <- N_age
  dim(vw_R124) <- N_age
  dim(vw_R134) <- N_age
  dim(vw_R234) <- N_age
  dim(vw_R1234) <- N_age
  
  
  
  
  