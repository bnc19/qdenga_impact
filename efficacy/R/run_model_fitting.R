run_model = function(n_it = 10000,
                     adapt_delta = 0.8,
                     n_chains = 4,
                     stan_model ,
                     folder,
                     init_values = function(){
                       list(
                         hs = runif(2, 2,4),
                         hl = runif(1, 108,132),
                         ts = runif(2, -2,2),
                         lambda_D = replicate(6, runif(4,0,0.02)),
                         p = runif(3,0.1,0.9), 
                         pK3 = runif(4,0.1,0.9),
                         gamma = runif(1,0.1,0.9),
                         rho = runif(4,1,5), 
                         phi = runif(1,0.1,0.5),
                         delta = runif(4,0.1,0.8),
                         L = runif(4,0,4),
                         w = runif(4,1,4),
                         alpha = runif(4,0,2),
                         beta = runif(2,0,2),
                         tau = runif(4,1,5),
                         lc = replicate(4,runif(3,2,8)),
                         sens = runif(1,.8,1),
                         spec = runif(1, 0.95,1),
                         epsilon = runif(1,1,4)
                       )},
                     VCD_years = c(12, 18, 24, 36, 48, 54) / 12,
                     time =1:54,
                     B = 2,
                     K = 4,
                     V = 2,
                     R = 2,
                     J = 3,
                     C = 3,
                     HI = 12,
                     include_pK3 = 0,
                     include_eps = 0,
                     include_beta = 0,
                     mono_lc_SN = 0,
                     mono_lc_MU = 0,
                     rho_K = 0,
                     L_K = 0,
                     w_CK = 0,
                     alpha_CK = 0,
                     tau_K = 0,
                     L_sd = 1,
                     L_mean =0,
                     lower_bound_L =0,
                     MU_test_SN = 0,
                     MU_symp = 1,
                     enhancement = 1, 
                     diagnostics = F,
                     mu,
                     baseline_SP,
                     VCD,
                     hosp,
                     titres = NULL,
                     metric = c("VCD_BVKD",
                                "VCD_BVJA",
                                "VCD_KJ2",
                                "HOSP_BVK4" ,
                                "HOSP_BVJA" ,
                                "HOSP_KJ2"),
                     serostatus = c("seronegative", "monotypic", "multitypic"),
                     BF=F) {

# set up -----------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)  
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
options(mc.cores = parallel::detectCores())

# source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

comp_model = cmdstan_model(paste0("models/", stan_model),stanc_options = list("O1"))
file_path = (paste0("output/", folder))
dir.create(file_path)
# data -----------------------------------------------------------------------
hosp = factor_VCD(hosp)
VCD = factor_VCD(VCD)
# fit Stan model -------------------------------------------------------------

list_data = format_stan_data(
  baseline_SP = baseline_SP,
  VCD = VCD,
  hosp = hosp,
  B = B,
  R = R,
  K = K,
  V = V,
  J = J,
  C = C,
  HI = HI,
  include_pK3 = include_pK3,
  include_eps = include_eps,
  include_beta = include_beta,
  mono_lc_SN = mono_lc_SN,
  mono_lc_MU = mono_lc_MU,
  rho_K = rho_K,
  L_K = L_K,
  w_CK = w_CK,
  alpha_CK = alpha_CK,
  tau_K = tau_K,
  L_sd = L_sd,
  L_mean = L_mean,
  lower_bound_L = lower_bound_L,
  MU_test_SN = MU_test_SN,
  MU_symp = MU_symp,
  enhancement = enhancement,
  VCD_years = VCD_years,
  time = time,
  mu = mu)

start = Sys.time()
stan_fit = comp_model$sample(
  data = list_data,
  chains = n_chains,
  parallel_chains = n_chains,
  iter_warmup = floor(n_it/2),
  iter_sampling = floor(n_it/2),
  thin = 1,
  init = init_values,
  seed = 14,
  refresh = 500,
  adapt_delta = adapt_delta)
end = Sys.time()
print(end - start)

# save posts -------------------------------------------------------------------

fit_ext = stan_fit$draws(format = "df")
i1 = which(names(fit_ext) == "ll")
names_select = names(fit_ext)[1:i1]
posterior_chains = fit_ext[names_select]
write.csv(posterior_chains, paste0(file_path, "/posterior_chains.csv"))

posterior = summarise_draws(posterior_chains)
write.csv(posterior, paste0(file_path, "/posterior.csv"))

# WAIC 
WAIC = waic(stan_fit$draws("log_lik"))
saveRDS(WAIC, file = paste0(file_path, "/WAIC.RDS"))

# PPC -------------------------------------------------------------------------- 
y = y_rep = out = out2 = list()
for (i in 1:length(metric)) {
  y[[i]] = c(list_data[[metric[i]]])
  y_rep[[i]] = as.matrix(fit_ext[, grepl(paste0("pred_" , metric[i]), names(fit_ext))])
  out[[i]] = ppc_dens_overlay(y[[i]], y_rep[[i]][1:50,])
  ggsave(out[[i]], file = paste0(file_path, "/PPC_", metric[i], ".jpg"))
  out2[[i]] = ppc_hist(y[[i]], y_rep[[i]][1:10,], binwidth = 5)
  ggsave(out2[[i]], file = paste0(file_path, "/PPC2_", metric[i], ".jpg"))
}

# save output ---------------------------------------------------------------------  
AR = which(grepl("AR" , names(fit_ext)))
AR_out = fit_ext[AR] %>%  as.data.frame()
saveRDS(AR_out, paste0(file_path, "/AR.RDS"))

VE = which(grepl("VE" , names(fit_ext)))
VE_out = fit_ext[VE] %>%  as.data.frame()
saveRDS(VE_out, paste0(file_path, "/VE.RDS"))

n = which(grepl("n" , names(fit_ext)))
n_out = fit_ext[n] %>%  as.data.frame()
saveRDS(n_out, paste0(file_path, "/n.RDS"))

}
