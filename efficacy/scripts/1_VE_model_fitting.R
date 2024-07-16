# Script to run final survival models ------------------------------------------
#  - Each model variant (31 total) is the same stan model, with different
#    parameters turned on and off, using the flags.
# -  The final model (M30) is given at the top, all other models are commented 
#    out but can be run below. 
#  - n_it defines the number of iterations to model is run for, set here to 1000
#    as a demonstration, but run for 10000 for all models presented in the 
#    paper. 
# ------------------------------------------------------------------------------

# source functions
# setwd(paste0(getwd(), "/efficacy"))

file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
n_it = 1000 # change this to increase model iterations 

# create folders
dir.create("output")

# M30 - as M27 but include_pK3 = 0 ---------------------------------------------
run_model (
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M30",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "final_model.stan",
  baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/vcd_data.csv"),
  hosp = read.csv("data/hosp_data.csv") ,
  mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
)


# # M1 ---------------------------------------------------------------------------
# run_model (
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M1",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M2 - as M1 but with tau_K = 1 ------------------------------------------------
# run_model (
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M2",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M3 - as M1 but with mono_lc_MU = 1 -------------------------------------------
# run_model (
#   mono_lc_MU = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M3",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M4 - as M2 but with mono_lc_MU = 1 -------------------------------------------
# run_model (
#   mono_lc_MU = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M4",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M5 - as M1 but with mono_lc_SN = 1 -------------------------------------------
# run_model (
#   mono_lc_SN = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M5",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M6 - as M2 but with mono_lc_SN = 1 -------------------------------------------
# run_model (
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M6",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M7 - as M5 but with mono_lc_MU = 1 -------------------------------------------
# run_model (
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M7",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M8 - as M6 but with mono_lc_MU = 1 -------------------------------------------
# run_model (
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M8",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
#  # M9 - as M2 but include_beta = 2 ----------------------------------------------
# run_model (
#   include_beta = 2,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M9",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M10 - as M1 but include_beta = 2 ---------------------------------------------
# run_model (
#   include_beta = 2,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M10",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M11 - as M2 but with w_CK = 1 ------------------------------------------------
# run_model (
#   w_CK = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M11",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M12 - as M2 but with w_CK = 2 ------------------------------------------------
# run_model (
#   w_CK = 2,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M12",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M13 - as M11 but with alpha_CK = 1 -------------------------------------------
# run_model (
#   alpha_CK = 1,
#   w_CK = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M13",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M14 - as M12 but with alpha_CK = 1 -------------------------------------------
# run_model (
#   alpha_CK = 1,
#   w_CK = 2,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M14",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M15 - as M11 but with mono_lc_SN = mono_lc_MU = 1 ----------------------------
# run_model (
#   mono_lc_SN = 1,
#   mono_lc_MU = 1,
#   w_CK = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M15",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M16 - as M12 but with mono_lc_SN = mono_lc_MU = 1 ----------------------------
# run_model (
#   mono_lc_SN = 1,
#   mono_lc_MU = 1,
#   w_CK = 2,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M16",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M17 - as M13 but with mono_lc_SN = mono_lc_MU = 1-----------------------------
# run_model (
#   mono_lc_SN = 1,
#   mono_lc_MU = 1,
#   alpha_CK = 1,
#   w_CK = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M17",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M18 - as M14 but with mono_lc_SN = mono_lc_MU = 1-----------------------------
# run_model (
#   mono_lc_SN = 1,
#   mono_lc_MU = 1,
#   alpha_CK = 1,
#   w_CK = 2,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M18",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M19 - as M8 but with alpha_CK = 1 --------------------------------------------
# run_model (
#   alpha_CK = 1,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_beta = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M19",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M20 - as M19 but with include_beta = 2 ---------------------------------------
# run_model (
#   include_beta = 2,
#   alpha_CK = 1,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M20",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M21 - as M18 but with include_beta = 2 ---------------------------------------
# run_model (
#   include_beta = 2,
#   mono_lc_SN = 1,
#   mono_lc_MU = 1,
#   alpha_CK = 1,
#   w_CK = 2,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M21",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M22 - as M21 but with w_CK = 1 -----------------------------------------------
# run_model (
#   w_CK = 1,
#   include_beta = 2,
#   mono_lc_SN = 1,
#   mono_lc_MU = 1,
#   alpha_CK = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M22",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M23 - as M8 but with include_beta = 0 ----------------------------------------
# run_model (
#   include_beta = 0,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M23",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M24 - as M8 but with include_beta = 3 ----------------------------------------
# run_model (
#   include_beta = 3,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M24",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M25 - as M24 but with alpha_CK = 2 -------------------------------------------
# run_model (
#   alpha_CK = 2,
#   include_beta = 3,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M25",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M26 - as M25 but with mono_lc_SN = 0 -----------------------------------------
# run_model (
#   mono_lc_SN = 0,
#   alpha_CK = 2,
#   include_beta = 3,
#   mono_lc_MU = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   rho_K = 1,
#   MU_test_SN = 1,
#   folder = "M26",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M27 - as M24 but with rho_K = 0 ----------------------------------------------
# run_model (
#   rho_K = 0,
#   include_beta = 3,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   MU_test_SN = 1,
#   folder = "M27",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M28 - as M27 L_K = 1 ---------------------------------------------------------
# run_model (
#   L_K = 1,
#   rho_K = 0,
#   include_beta = 3,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   MU_test_SN = 1,
#   folder = "M28",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M29 - as M27 but include_eps = 1 ---------------------------------------------
# run_model (
#   include_eps = 1,
#   rho_K = 0,
#   include_beta = 3,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   include_pK3 = 1,
#   MU_test_SN = 1,
#   folder = "M29",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
# # M31 - as M30 but with L_K = 1 ------------------------------------------------
# run_model (
#   L_K = 1, 
#   include_pK3 = 0,
#   rho_K = 0,
#   include_beta = 3,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = 1,
#   MU_test_SN = 1,
#   folder = "M31",
#   n_it = n_it,
#   adapt_delta = 0.77,
#   stan_model = "final_model.stan",
#   baseline_SP = read.csv("data/seropositive_by_age_baseline.csv"),
#   VCD =  read.csv("data/vcd_data.csv"),
#   hosp = read.csv("data/hosp_data.csv") ,
#   mu =  array(read.csv("data/n0_new.csv")$mean, dim = c(2, 4))
# )
# 
