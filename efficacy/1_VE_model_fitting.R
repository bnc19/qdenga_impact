# Script to run final survival models ------------------------------------------
# Each model variant (31 total) is the same stan model, with different
# parameters turned on and off, using the flags.
# Lower bond of L is 0 for all models
# Hosp given inf
# Swap gamma and rho
# Constrain p
# Model with (1+L)/(1+(n/n50)^w)
# ------------------------------------------------------------------------------

# source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
n_it = 10000

# M1 ---------------------------------------------------------------------------
run_model (
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M1",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M2 - as M1 but with tau_K = 1 ------------------------------------------------
run_model (
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M2",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M3 - as M1 but with mono_lc_MU = 1 -------------------------------------------
run_model (
  mono_lc_MU = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M3",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M4 - as M2 but with mono_lc_MU = 1 -------------------------------------------
run_model (
  mono_lc_MU = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M4",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M5 - as M1 but with mono_lc_SN = 1 -------------------------------------------
run_model (
  mono_lc_SN = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M5",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M6 - as M2 but with mono_lc_SN = 1 -------------------------------------------
run_model (
  mono_lc_SN = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M6",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M7 - as M5 but with mono_lc_MU = 1 -------------------------------------------
run_model (
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M7",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M8 - as M6 but with mono_lc_MU = 1 -------------------------------------------
run_model (
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M8",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

 # M9 - as M2 but include_beta = 2 ----------------------------------------------
run_model (
  include_beta = 2,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M9",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M10 - as M1 but include_beta = 2 ---------------------------------------------
run_model (
  include_beta = 2,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M10",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M11 - as M2 but with w_CK = 1 ------------------------------------------------
run_model (
  w_CK = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M11",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M12 - as M2 but with w_CK = 2 ------------------------------------------------
run_model (
  w_CK = 2,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M12",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M13 - as M11 but with alpha_CK = 1 -------------------------------------------
run_model (
  alpha_CK = 1,
  w_CK = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M13",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M14 - as M12 but with alpha_CK = 1 -------------------------------------------
run_model (
  alpha_CK = 1,
  w_CK = 2,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M14",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M15 - as M11 but with mono_lc_SN = mono_lc_MU = 1 ----------------------------
run_model (
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  w_CK = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M15",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M16 - as M12 but with mono_lc_SN = mono_lc_MU = 1 ----------------------------
run_model (
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  w_CK = 2,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M16",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M17 - as M13 but with mono_lc_SN = mono_lc_MU = 1-----------------------------
run_model (
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  alpha_CK = 1,
  w_CK = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M17",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M18 - as M14 but with mono_lc_SN = mono_lc_MU = 1-----------------------------
run_model (
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  alpha_CK = 1,
  w_CK = 2,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M18",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M19 - as M8 but with alpha_CK = 1 --------------------------------------------
run_model (
  alpha_CK = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_beta = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M19",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M20 - as M19 but with include_beta = 2 ---------------------------------------
run_model (
  include_beta = 2,
  alpha_CK = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M20",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M21 - as M18 but with include_beta = 2 ---------------------------------------
run_model (
  include_beta = 2,
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  alpha_CK = 1,
  w_CK = 2,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M21",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M22 - as M21 but with w_CK = 1 -----------------------------------------------
run_model (
  w_CK = 1,
  include_beta = 2,
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  alpha_CK = 1,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M22",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M23 - as M8 but with include_beta = 0 ----------------------------------------
run_model (
  include_beta = 0,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M23",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M24 - as M8 but with include_beta = 3 ----------------------------------------
run_model (
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M24",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M25 - as M24 but with alpha_CK = 2 -------------------------------------------
run_model (
  alpha_CK = 2,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M25",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M26 - as M25 but with mono_lc_SN = 0 -----------------------------------------
run_model (
  mono_lc_SN = 0,
  alpha_CK = 2,
  include_beta = 3,
  mono_lc_MU = 1,
  tau_K = 1,
  include_pK3 = 1,
  rho_K = 1,
  MU_test_SN = 1,
  folder = "final/M26",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M27 - as M24 but with rho_K = 0 ----------------------------------------------
run_model (
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  MU_test_SN = 1,
  folder = "final/M27",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M28 - as M27 L_K = 1 ---------------------------------------------------------
run_model (
  L_K = 1,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  MU_test_SN = 1,
  folder = "final/M28",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M29 - as M27 but include_eps = 1 ---------------------------------------------
run_model (
  include_eps = 1,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 1,
  MU_test_SN = 1,
  folder = "final/M29",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M30 - as M27 but include_pK3 = 0 ---------------------------------------------
run_model (
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/M30",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M31 - as M30 but with L_K = 1 ------------------------------------------------
run_model (
  L_K = 1, 
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/M31",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "multitypic/final_model.stan",
  baseline_SP = read.csv("data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

