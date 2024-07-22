config_params <- function(posterior_file,ns,vacc_cov,zero_ve,vacc_age,DoVEInf,vacc_by_serostatus,max_decay_years) {


  ## set year vars
  YL <<- 365 # year length in days 
  FIRST_YEAR <<- 1950
  LAST_YEAR <<- 2045
  CALIB_YEAR <<- 2024
  VACC_YEAR <<- 2024
  VACC_END_YEAR <<- 2045
  equilib_years <<- 98  ## because there are 74 years from 1950 to 2024 (start of vacc), and want 98+74 to be divisible by 4
  
  equilib_steps <<- equilib_years*YL
  run_steps_pre_vacc <<- (VACC_YEAR-FIRST_YEAR)*YL+equilib_steps
  run_steps_post_vacc <<- (VACC_END_YEAR-FIRST_YEAR)*YL+equilib_steps
  max_ve_decay_time <- max_decay_years*YL
  
  countries <<- c("BRA","PHL")
  country_hosp <<- c(0.133,0.089) # prop of symp cases hospitalised by country
  nc <<- length(countries) # number of countries
  
  nvy <- 2100-VACC_YEAR+1
  coverage_d <- matrix(c(VACC_YEAR:2100,rep(vacc_cov,nvy)),nrow=nvy,ncol=2)
  
  # demog <- process_demog("demog/",
  #                      "demog/births_all",
  #                      "demog/pop_all",
  #                      "demog/life_ex_all",
  #                      "demog/mort_rate_all",
  #                      countries[cn],
  #                      c(0:19,20,30,40,50,60,70,80,90),
  #                      TRUE,"")
  
  demog <- process_demog("demog/",
                              "demog/births_all",
                              "demog/pop_all",
                              "demog/life_ex_all",
                              "demog/mort_rate_all",
                              countries[cn],
                              c(0:32,33,40,50,60,70,80,90),
                              TRUE,"")

  pp <- read.csv(posterior_file)
  nrp <- nrow(pp)
  if(nrp<ns) stop(paste0("Not enough posterior param samples (",ns,">",nrow(pp),")"))
  ## thin posterior samples
  st <- floor(nrp/ns)
  pp <- pp[seq(1,nrp,by=st),]
  
  # some defaults - over-written later
  dis_pri <- 0.15
  dis_sec <- 2*dis_pri
  dis_tq <- 0.25*dis_pri
  sdis_pri <- dis_pri*0.125
  sdis_sec <- dis_sec*0.125
  sdis_tq <- dis_tq*0.125

  nc0 <- matrix(c(296.1,2046.4,239.8,223.8,
                 2462.0,4749.2,1701.4,1367.9,
                 2462.0,4749.2,1701.4,1367.9),nrow = 4, ncol = 3)
  pars_nv <- list(YL=YL,
               FIRST_YEAR = 1950,
               LAST_YEAR = 2100,
               CALIB_YEAR = 2024,
               VACC_YEAR = VACC_YEAR,
               VACC_END_YEAR = VACC_END_YEAR,
               EQUILIB_YEARS = equilib_years,
               N_sim = 5000000, # set to <=0 to run with actual country pop size
        ## params for health econ calcs
               cfr = 0.005, # fraction of sdis cases who die
               daly_dis = 0.006, # 0.545*4/365
               daly_sdis = 0.014, # 0.545*10/365 - daly_dis
               disc_fact=1.03, # 3% discount
               disc_fact_vacc=1.03, # 3% discount
        ## mosquito params
               Rm = 2.69,
               omega = 1.0,
               delta = 0.1,
               epsilon = 1.0/19.0,
               Mwt = 1.5,
               sigma = 0.025,
               Kc_season = 0.3,
               kappa = 0.6,
               Beta_mh_mean = 1.0,
               Beta_mh_season = 0.0,
        ## external FOI to prevent fade-out
               extInf = 1e-8,
        ## required R0 - over-written later - Mwt also changed to match
               R0_req = 1.5, 
               Rel_R01 = 1.0, # best to break degeneracy for deterministic runs
               Rel_R02 = 1.0,
               Rel_R03 = 1.0,
               Rel_R04 = 1.0,
               phi_dis_enhance=1,
               incub = 5,
               eip = 8,
               inf_per = 4,
               dur_cross_prot = YL,
        ## disease severity params - over-written later
               dis_pri = c(dis_pri,dis_pri,dis_pri,dis_pri),
               dis_sec = c(dis_sec,dis_sec,dis_sec,dis_sec),
               dis_tert = c(dis_tq,dis_tq,dis_tq,dis_tq),
               dis_quart = c(dis_tq,dis_tq,dis_tq,dis_tq),
               sdis_pri = c(sdis_pri,sdis_pri,sdis_pri,sdis_pri),
               sdis_sec = c(sdis_sec,sdis_sec,sdis_sec,sdis_sec),
               sdis_tert = c(sdis_tq,sdis_tq,sdis_tq,sdis_tq),
               sdis_quart = c(sdis_tq,sdis_tq,sdis_tq,sdis_tq),
        ## Vacc params - over-written later
               VaccOn = 1, # leave at 1
               ZeroVE = zero_ve, # set to 1 for VE=0
               DoVEInf = DoVEInf,
               VaccRoutineAge = vacc_age,
               vacc_by_serostatus=vacc_by_serostatus,
               max_ve_decay_time = max_ve_decay_time,
               nc0 = nc0,
               hs = log(2)/c(50,100,100),
               hl = log(2)/c(400,400,400),
               ts = c(120,120,120),
               nc50_inf =  matrix(c(200,200,200,200,
                             400,400,400,400,
                             400,400,400,400),nrow = 4, ncol = 3),
               nc50_dis =  matrix(c(200,200,200,200,
                                    400,400,400,400,
                                    400,400,400,400),nrow = 4, ncol = 3),
               nc50_sdis = matrix(c(50,50,50,50,
                                    100,100,100,100,
                                    100,100,100,100),nrow = 4, ncol = 3),
               L_dis = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1),nrow = 4, ncol = 3),
               L_sdis = matrix(c(3,3,3,3,1,1,1,1,1,1,1,1),nrow = 4, ncol = 3),
               ws = c(0.75,0.75,0.75,0.75),
               nc50_age_under6 = 0,
        ## demog prarams
               age_removal_d=as.matrix(demog[[3]]),
               age_rate_d=as.matrix(demog[[4]]), 
               pop_size_d=unlist(demog[[2]]),
               births_d=as.matrix(demog[[1]]), 
               life_expec_d=as.matrix(demog[[5]]),
        ## vaccine coverage
               coverage_d=coverage_d)
  # scale ve against inf
  chi = c(rep(log(12), 4), rep(log(3), 8))
  chi = matrix(chi, nrow = 4, ncol = 3)
  ## create list of ns copies of param list
  pars <- unname(rep(list(pars_nv),ns))
  
  ## now update each param set in turn with posterior param samples
  for(i in 1:ns) {
    pars[[i]]$hl <- log(2)/c(pp$hl[i]*YL/12,pp$hl[i]*YL/12,pp$hl[i]*YL/12)
    pars[[i]]$hs <- log(2)/c(pp$hs.1.[i]*YL/12,pp$hs.2.[i]*YL/12,pp$hs.2.[i]*YL/12)
    pars[[i]]$ts <- c(pp$ts.1.[i]*YL/12,pp$ts.2.[i]*YL/12,pp$ts.2.[i]*YL/12)
    pars[[i]]$nc0[1,] <- nc0[1,]/(exp(-pars[[i]]$ts*pars[[i]]$hs)+exp(-pars[[i]]$ts*pars[[i]]$hl))
    pars[[i]]$nc0[2,] <- nc0[2,]/(exp(-pars[[i]]$ts*pars[[i]]$hs)+exp(-pars[[i]]$ts*pars[[i]]$hl))
    pars[[i]]$nc0[3,] <- nc0[3,]/(exp(-pars[[i]]$ts*pars[[i]]$hs)+exp(-pars[[i]]$ts*pars[[i]]$hl))
    pars[[i]]$nc0[4,] <- nc0[4,]/(exp(-pars[[i]]$ts*pars[[i]]$hs)+exp(-pars[[i]]$ts*pars[[i]]$hl))
    pars[[i]]$dis_sec <- c(pp$gamma[i],pp$gamma[i],pp$gamma[i],pp$gamma[i])
#    pars[[i]]$dis_pri <- c(pp$gamma[i]/pp$rho.1.[i],pp$gamma[i]/pp$rho.2.[i],pp$gamma[i]/pp$rho.3.[i],pp$gamma[i]/pp$rho.4.[i])
    pars[[i]]$dis_pri <- c(pp$gamma[i]/pp$rho.1.[i],pp$gamma[i]/pp$rho.1.[i],pp$gamma[i]/pp$rho.1.[i],pp$gamma[i]/pp$rho.1.[i])
    pars[[i]]$dis_tert <- pars[[i]]$dis_pri * pp$phi[i]
    pars[[i]]$dis_quart <- pars[[i]]$dis_tert
### code to directly use serotype specific hosp estimates (not used, as biased by country)
#    pars[[i]]$sdis_pri <-  pars[[i]]$dis_pri * c(pp$delta.1[i],pp$delta.2[i],pp$delta.3[i],pp$delta.4[i])
#    pars[[i]]$sdis_sec <-  pars[[i]]$dis_sec * c(pp$delta.1[i],pp$delta.2[i],pp$delta.3[i],pp$delta.4[i])
#    pars[[i]]$sdis_tert <-  pars[[i]]$dis_tert * c(pp$delta.1[i],pp$delta.2[i],pp$delta.3[i],pp$delta.4[i])

# code assuming no serotype differences in prop of symp cases hospitalised
# likely needs lower risk for tertiary/quaternary
    pars[[i]]$sdis_pri <-  pars[[i]]$dis_pri * country_hosp[cn]
    pars[[i]]$sdis_sec <-  pars[[i]]$dis_sec * country_hosp[cn]
    pars[[i]]$sdis_tert <-  pars[[i]]$dis_tert * country_hosp[cn]
    pars[[i]]$sdis_quart <-  pars[[i]]$sdis_tert
    pars[[i]]$ws <- c(pp$w.1.[i],pp$w.1.[i],pp$w.1.[i],pp$w.1.[i])
#    pars[[i]]$L_dis <- matrix(c(1+pp$L.1.1.[i],1+pp$L.2.1.[i],1+pp$L.3.1.[i],1+pp$L.4.1.[i],
#                                1,1,1,1,1,1,1,1),nrow = 4, ncol = 3)
#    pars[[i]]$L_sdis <- matrix(c(1+pp$L.1.1.[i]*pp$tau.1.[i],1+pp$L.2.1.[i]*pp$tau.2.[i],1+pp$L.3.1.[i]*pp$tau.3.[i],1+pp$L.4.1.[i]*pp$tau.4.[i],
#                                 1,1,1,1,1,1,1,1),nrow = 4, ncol = 3)
    pars[[i]]$L_dis <- matrix(c(1+pp$L.1.1.[i],1+pp$L.1.1.[i],1+pp$L.1.1.[i],1+pp$L.1.1.[i],
                                1,1,1,1,1,1,1,1),nrow = 4, ncol = 3)
    pars[[i]]$L_sdis <- matrix(c(1+pp$L.1.1.[i]*pp$tau.1.[i],1+pp$L.1.1.[i]*pp$tau.2.[i],1+pp$L.1.1.[i]*pp$tau.3.[i],1+pp$L.1.1.[i]*pp$tau.4.[i],
                                 1,1,1,1,1,1,1,1),nrow = 4, ncol = 3)
#    pars[[i]]$nc50_dis <- matrix(c(pp$lc.1.1.[i],pp$lc.1.2.[i],pp$lc.1.3.[i],pp$lc.1.4.[i],
#                                   pp$lc.2.1.[i],pp$lc.2.2.[i],pp$lc.2.3.[i],pp$lc.2.4.[i],
#                                   pp$lc.3.1.[i],pp$lc.3.2.[i],pp$lc.3.3.[i],pp$lc.3.4.[i]),nrow = 4, ncol = 3)
    pars[[i]]$nc50_dis <- matrix(c(pp$lc.1.1.[i],pp$lc.1.1.[i],pp$lc.1.1.[i],pp$lc.1.1.[i],
                                   pp$lc.2.1.[i],pp$lc.2.2.[i],pp$lc.2.3.[i],pp$lc.2.4.[i],
                                   pp$lc.3.1.[i],pp$lc.3.1.[i],pp$lc.3.1.[i],pp$lc.3.1.[i]),nrow = 4, ncol = 3)
    pars[[i]]$nc50_inf <- pars[[i]]$nc50_dis + chi
    pars[[i]]$nc50_sdis <- pars[[i]]$nc50_dis - pp$alpha.1.[i]
    pars[[i]]$nc50_inf <- exp(pars[[i]]$nc50_inf)
    pars[[i]]$nc50_dis <- exp(pars[[i]]$nc50_dis)
    pars[[i]]$nc50_sdis <- exp(pars[[i]]$nc50_sdis)
    ## Simple age dependence - <6 different
    pars[[i]]$nc50_age_under6 <- pp$beta.1.[i]
    }
  
  ## number of different R0 values to sample
  n_R0 <<- 9
  R0 <<- matrix(nrow=n_R0,ncol=nc) # R0 values which give 10%,20%,..,90% seroprev in 9 year olds, for each country
  R0[,1] <<- c(1.170,1.345,1.570,1.818,2.140,2.534,3.077,3.858,5.198)
  R0[,2] <<- c(1.128,1.250,1.390,1.560,1.760,2.032,2.372,2.857,3.680)
  return(pars)
}