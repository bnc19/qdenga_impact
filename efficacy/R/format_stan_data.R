format_stan_data = function(baseline_SP,
                            VCD,
                            hosp,
                            mu,
                            B,
                            K,
                            R,
                            C,
                            HI,
                            V,
                            J,
                            VCD_years,
                            time,
                            include_pK3,
                            include_eps,
                            include_beta,
                            mono_lc_SN,
                            mono_lc_MU,
                            rho_K,
                            L_K,
                            w_CK,
                            alpha_CK,
                            tau_K,
                            L_sd,
                            L_mean,
                            lower_bound_L,
                            MU_test_SN,
                            MU_symp,
                            enhancement) {

T = length(time)
D = length(VCD_years)
A = sum(VCD_years < 4)

# pop by age group at trial start
pop_J = VCD %>%
  filter(age != "all", year == 1, serotype == "all") %>%
  group_by(age) %>%
  summarise(N = sum(N))
  
# calculate # VCD over time ---------------------------------------------------
VCD_D = VCD %>%
  filter(serostatus == "both",
         age == "all",
         serotype == "all",
         year %in% VCD_years) %>%
  group_by(year) %>%
  summarise(Y = sum(Y)) %>%
  select(Y)

# calculate # hosp over time ---------------------------------------------------
# year 1-3 from age data 

hosp_4 = hosp %>%
  filter(month > 36) %>%
  group_by(year) %>%
  summarise(Y = sum(Y)) %>%
  select(Y)

hosp_D =
  hosp %>%
  filter(age != "all", trial != "both", year < 48) %>%
  group_by(year) %>%
  summarise(Y = sum(Y)) %>%
  select(Y)  %>%
  bind_rows(hosp_4)
  
# calculate # pop by trial arm, serostatus, age group, over time ---------------
N_pop_BVJD = VCD %>%
  filter(serostatus != "both",
         age != "all",
         serotype == "all",
         trial != "both") %>%
  arrange(year, age, trial , serostatus)

pop_BVJD = array(N_pop_BVJD$N, dim = c(B, V, J, D))

# calculate # pop by trial arm, serostatus, over time  -------------------------
N_pop_BVD =  VCD %>%
  filter(serostatus != "both",
         age == "all",
         serotype == "all",
         trial != "both",
         year %in% VCD_years) %>%
  arrange(year, trial , serostatus)

pop_BVD = array(N_pop_BVD$N, dim = c(B, V, D))

# calculate # hosp by serostatus + trial arm + serotype, over time (36-54 months)
N_hosp_BVK4 = hosp %>%
  filter(serotype != "all", trial != "both") %>%
  arrange(year, serostatus, trial , serotype)

hosp_BVK4_m = array(N_hosp_BVK4$Y, dim = c(B * K * V, 3))  # (D=4)

# calculate # VCD by serostatus + trial arm + serotype, over time --------------
# for multinomial likelihood
N_VCD_BVKD_m = VCD %>%
  group_by(serostatus, trial, age, serotype, year) %>%
  summarise(Y = sum(Y)) %>%
  filter(serostatus != "both",
         age == "all",
         serotype != "all",
         trial != "both",
         year %in% VCD_years) %>%
  arrange(year, serostatus, trial , serotype)

VCD_BVKD_m = array(N_VCD_BVKD_m$Y, dim = c(B * K * V, D))

# calculate # hosp by serostatus + trial arm + age-group, over time -------------
N_hosp_BVJA_m =  hosp %>%
  filter(serotype == "all", trial != "both", age != "all") %>%
  arrange(year, serostatus, trial , age)

hosp_BVJA_m = array(N_hosp_BVJA_m$Y, dim = c(B * V * J, A))

# calculate # VCD by serostatus + trial arm + age-group, over time -------------
# for multinomial likelihood
N_VCD_BVJA_m = VCD %>%
  group_by(serostatus, trial, age, serotype, year) %>%
  summarise(Y = sum(Y)) %>%
  filter(serostatus != "both",
         age != "all",
         serotype == "all",
         trial != "both",
         year < 4) %>%
  arrange(year, serostatus, trial , age)

VCD_BVJA_m = array(N_VCD_BVJA_m$Y, dim = c(B * V * J, A))

# VCD age and serotype year and 2 ----------------------------------------------
N_VCD_KJ2_m = VCD %>%
  filter(age != "all", serotype != "all") %>%
  arrange(year, serotype, age)

VCD_KJ2_m = array(N_VCD_KJ2_m$Y, dim = c(K * J, 2))

# hosp age and serotype year and 2 ----------------------------------------------
N_hosp_KJ2_m = hosp %>%
  filter(age != "all", serotype != "all") %>%
  arrange(year, serotype, age)

hosp_KJ2_m = array(N_hosp_KJ2_m$Y, dim = c(K * J, 2))

# VCD by serostatus, trial arm for t_5 (to censor population t_6)
N_VCD_BV5 = VCD %>%
  filter(age == "all", year == 4, serotype == "all", serostatus != "both") %>%
  arrange(trial, serostatus)

N_VCD_BV5_m  = array(N_VCD_BV5$Y, dim = c(2, 2))

# data -------------------------------------------------------------------------

stan_data = list(
  time = time,
  T = T,
  J = J,
  K = K,
  B = B,
  V = V,
  D = D,
  A = A,
  C = C,
  HI = HI,
  R = R,
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
  SP_J = baseline_SP$Seropositive,
  pop_J = baseline_SP$Total,
  VCD_D = as.numeric(VCD_D$Y),
  pop = pop_BVJD,
  pop_BVD = pop_BVD,
  VCD_BVKD = VCD_BVKD_m,
  VCD_BVJA = VCD_BVJA_m,
  mu = mu,
  HOSP_BVJA = hosp_BVJA_m,
  HOSP_BVK4 = hosp_BVK4_m,
  HOSP_D = as.numeric(hosp_D$Y),
  HOSP_KJ2 = hosp_KJ2_m,
  VCD_KJ2 = VCD_KJ2_m,
  N_VCD_BV5 = N_VCD_BV5_m
)
  return(stan_data)
}
