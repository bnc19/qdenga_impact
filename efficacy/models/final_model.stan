data{
   int<lower = 1> C; // No. current immune status
   int<lower = 1> B; // No. baseline serostatus
   int<lower = 1> K; // No. serotypes
   int<lower = 1> D; // No. vcd data
   int<lower = 1> V; // No. trial arms 
   int<lower = 1> J; // No. age groups
   int<lower = 1> A; // No. age vcd data 
   int<lower = 1> T; // No. model time points 
   int<lower = 1> R; // No. outcomes 
   int<lower = 0> time[T]; // time 
   real<lower = 0> mu[B,K]; // titres at T0
   int<lower = 0> HI; // period of heterotypic immunity
   array[J] int<lower=0> SP_J;     // baseline seropositive
   array[J] int<lower=0> pop_J;    // baseline pop
   array[D] int<lower=0> VCD_D;    // total VCD 
   array[D] int<lower=0> HOSP_D;   // total hosp  
   array[B*V*K,D] int<lower=0> VCD_BVKD;  // VCD 
   array[B*V*J,A] int<lower=0> VCD_BVJA;  // VCD 
   array[K*J,2] int<lower=0> VCD_KJ2;     // VCD at 12 and 24 months
   array[B*V*J,A] int<lower=0> HOSP_BVJA; // hosp 
   array[B*V*K,4] int<lower=0> HOSP_BVK4; // hosp  
   array[K*J,2] int<lower=0> HOSP_KJ2;    // hosp at 12 and 24 months
   array[B,V,J,D] int<lower=0> pop;       // pop
   array[B,V,D] int<lower=0> pop_BVD;     // pop
   array[B,V] int<lower=0> N_VCD_BV5;     // cases in t_5 (not age-specific)
   // FLAGS
   int<lower = 0, upper = 1> include_pK3;      // include serotype specific p? (T/F)
   int<lower = 0, upper = 1> include_eps;      // include enhanced secondary hosp? (T/F)
   int<lower = 0, upper = 3> include_beta;     // include age-specific nc50?: 1= change age groups 1 & 2 for both outcomes / 2 = only age grp 1, sep for each outcome / 3 = only age grp 1 for both outcomes 
   int<lower = 0, upper = 1> mono_lc_SN;       // 0 for serotype specific lc SN / 1 for single 
   int<lower = 0, upper = 1> mono_lc_MU;       // 0 for serotype specific lc MO / 1 for single 
   int<lower = 0, upper = 1> rho_K;            // 0 for mono rho, 1 for serotype rho 
   int<lower = 0, upper = 1> L_K;              // 0 for mono L, 1 for serotype K 
   int<lower = 0, upper = 2> w_CK;             // 0 for mono w, 1 for serostatus w, 2 for serotype w 
   int<lower = 0, upper = 2> alpha_CK;         // 0 for mono alpha, 1 for serostatus alpha, 2 for serotype alpha 
   int<lower = 0, upper = 1> tau_K;            // 0 for mono tau, 1 for serotype tau
   real<lower = 0> L_mean;                     // mean (trunc normal) L prior
   real<lower = 0> L_sd;                       // sd (trunc normal) L prior 
   int<lower = -1, upper = 0>  lower_bound_L;  // lower bound of L  (-1, 0)
   int<lower = 0, upper = 1> MU_test_SN;       // multitypic test SN at baseline (T/F)
   int<lower = 0, upper = 1>  MU_symp;         // 0 for only 1' and 2' symp infections, 1 for post-sec symp inf   
   int<lower = 0, upper = 1> enhancement;      // 0 for no vac enhancement of SN, 1 for enhancement 
}

transformed data {
  array[B,J,2] real pop_BJD; 
  array[J,2] real pop_JD; 
  array[B,V,D] real pop_BVD_real;     
  array[B,V] real pop_BV;  
  array[B,V] real pop_BV_K;  
  array[V,D] int pop_VD;
  array[B,D] int pop_BD;
  array[D] int pop_D;
  vector<lower=0>[D] v_pop_D;
  vector[J] v_SP_J = to_vector(SP_J);     // baseline seropositive
  vector[J] v_pop_J = to_vector(pop_J);    // baseline pop
  vector[D] v_VCD_D = to_vector(VCD_D);    // total VCD 
  vector[D] v_HOSP_D = to_vector(HOSP_D);   // total hosp  
  matrix[B*V*K,D] v_VCD_BVKD=to_matrix(VCD_BVKD);  // VCD 
  matrix[B*V*J,A] v_VCD_BVJA=to_matrix(VCD_BVJA);  // VCD 
  matrix[K*J,2] v_VCD_KJ2=to_matrix(VCD_KJ2);     // VCD at 12 and 24 months
  matrix[B*V*J,A] v_HOSP_BVJA=to_matrix(HOSP_BVJA); // hosp 
  matrix[B*V*K,4] v_HOSP_BVK4=to_matrix(HOSP_BVK4); // hosp  
  matrix[K*J,2] v_HOSP_KJ2=to_matrix(HOSP_KJ2);    // hosp at 12 and 24 months

// POPULATION AGGREGATION 
for(b in 1:B)
 for(j in 1:J){
 pop_BJD[b,j,1] = sum(pop[b, ,j,1]);
 pop_BJD[b,j,2] = sum(pop[b, ,j,3]);
 }
 
for(j in 1:J)
 for(d in 1:2)
 pop_JD[j,d] = sum(pop_BJD[ ,j,d]); 

for(b in 1:B)
 for(v in 1:V)
  for(d in 1:D)
   pop_BVD_real[b,v,d] = pop_BVD[b,v,d]; 

for(b in 1:B)
 for(v in 1:V)
 pop_BV[b,v] = mean(pop_BVD_real[b,v, ]);
 
for(b in 1:B)
 for(v in 1:V)
 pop_BV_K[b,v] = mean(pop_BVD_real[b,v,4:6]); // serotype hosp data from 36months  
 
for(v in 1:V)
 for(d in 1:D)
 pop_VD[v,d] = sum(pop_BVD[,v,d]);

for(b in 1:B)
 for(d in 1:D)
  pop_BD[b,d] = sum(pop_BVD[b, ,d]); 

for(d in 1:D)
 pop_D[d] = sum(pop_BD[,d]);
 v_pop_D = to_vector(pop_D);
}

parameters{
  real<lower = 0> hs[B];                       // short term decay 
  real<lower = 12> hl;                         // long term decay 
  real<lower = -12> ts[B];                     // time switch decay rate
  real<lower = 0> lambda_D[K,D];               // FOI  
  real<lower = 0, upper = 1> p[J];             // probability of exposure before T0
  real<lower = 0, upper = 1> pK3[K];           // serotype-specific p in age group 3 
  real<lower = 1> epsilon ;                    // prob hosp 2' rel to 1'/3'/4'
  real<lower = 0, upper = 1> gamma;            // prob symp 2' 
  real<lower = 0> rho[K];                      // 1/risk symp 1' rel to 2'
  real<lower = 0, upper = 1> phi;              // prob symp 3'/4' rel to 1'
  real<lower = 0, upper = 1> delta[K];         // prob hosp
  real<lower = lower_bound_L> L[K];            // enhancement parameter 
  real<lower = 0> tau[K];                      // multiply L for hosp
  real<lower = 0> w[K];                        // shape parameter 
  real<lower = 0> lc[C,K];                     // 50% symp protection
  real<lower = 0> alpha[K];                    // reduction in titre for protection against hosp
  real beta[J-1];                              // increase in lc50 for protection in younger 
  real<lower = 0, upper = 1> sens ;            // baseline test sensitivity 
  real<lower = 0, upper = 1> spec ;            // baseline test specificity 
}

transformed parameters{
  real pSP[J]; // prob SP (baseline for likelihood)
  real ll; // log-likelihood passed to model block and then added to target
  real<lower = 0> C_KJRD[K,J,R,D];
  real<lower = 0> C_BVKJRD[B,V,K,J,R,D];
  real<lower = 0> n[B,K,T]; // titres
  real<lower = 0> RR_symp[C,V,K,J,T];             
  real<lower = 0> RR_hosp[C,V,K,J,T];
  // censor population for last time interval 
  real<lower = 0> pop_BVJD[B,V,J,D]; 
  // matrix distribution of cases
  matrix [K*J,2]   pD_KJ2 ;    
  matrix [B*V*K,D] pD_BVKD ;   
  matrix [B*V*J,A] pD_BVJA ;   
  matrix [K*J,2]   pH_KJ2 ;    
  matrix [B*V*K,4] pH_BVK4 ;   
  matrix [B*V*J,A] pH_BVJA ;
  real<lower=0> shape1; // pk3 prior
  real<lower=0> shape2; 
  //  distribution of cases
  real pC_RD[R,D];
  
  { // this { defines a block within which variables declared are local (can't have lower or upper)
  
  real lambda_mD[D];
  real lambda_m;
  real pm; 
  
  // contrain pk3 / p3 by lambda 
  for(d in 1:D) lambda_mD[d] = mean(lambda_D[ ,d]); // mean FOI at time D 
  lambda_m = (lambda_mD[1] * 12 + lambda_mD[2] * 6 + lambda_mD[3] * 6 + lambda_mD[4] * 12 + lambda_mD[5] * 12+ lambda_mD[6] * 6) / 54 ; // weighted mean across time 
  pm = 1-exp(-14*12*lambda_m); // always prob exposure to a single serotype 
  if(pm<1e-3) pm=1e-3;
  if(pm>0.99) pm=0.99;
  shape1 = 49 * pm ;
  shape2 = shape1 * (1/pm - 1);

  vector[K] rhoT; // keep code simple by always referring to vector length K
  vector[K] gammaT;
  vector[K] phiT;

  real hK3[K]; // age3 historic hazard
  real q[K];   // historic serotype proportion 
  real hK1[K]; // age1 historic hazard 
  real hK2[K]; // age2 historic hazard 
  real p_KJ[K,J]; // serotype age prob exposure 
  real TpSN[J];              // true SN
  real TpMO[K,J];            // true monotypic  
  real TpMU2[(K+2),J];       // true multitypic
  real TpMU3[K,J];           // true multitypic  
  real TpMU4[J];             // true multitypic  
  real pSN[B,V,J,(T+1)];     // prob SN time plus baseline 
  real pMO[B,V,K,J,(T+1)];   // prob monotypic time plus baseline 
  real pMU2[B,V,(K+2),J,(T+1)]; // prob multitypic time plus baseline 
  real pMU3[B,V,K,J,(T+1)];  // prob multitypic time plus baseline 
  real n_C[C,K,T];           // allow for multitypic titres 
  real L_C[C,K,R];           // MO and MU L = 1 
  real nc50[C,K,J,R];        // nc50
  real lambda[K,T];          // FOI 
  real Inc1[B,V,K,J,T];
  real Inc2[B,V,K,J,T];
  real Inc3[B,V,K,J,T];
  real Inc4[B,V,K,J,T];
  real D1[B,V,K,J,T];
  real D2[B,V,K,J,T];
  real D34[B,V,K,J,T];
  real H1[B,V,K,J,T];
  real H2[B,V,K,J,T];
  real H34[B,V,K,J,T];
  real Di[B,V,K,J,D];
  real H[B,V,K,J,D];
  real pop_VJD[V,J,D];
  real pC_BVKRD[B,V,K,R,D];
  real pC_BVJRA[B,V,J,R,A];
  real pC_BVK4[B,V,K,4];
  real pC_KJR2[K,J,R,2];
  real pJ_VCD[J]; // age dist of cases for month 37-48
  // cases 
  real C_BVJ5[B,V,J]; // serostatus trial age cases for month 37-48
  real C_BJ5[B,J]; // serostatus age cases for month 37-48
  real C_J5[J];  // age cases for month 37-48
  real C_BVKJR4[B,V,K,J,R,4];
  real C_VKJRD[V,K,J,R,D];
  real C_JRD[J,R,D];
  real C_RD[R,D];

// BIPHASIC TITRES 
  real pi_1[B];   // decay rate 1 
  real pi_2;      // decay rate 2

for(b in 1:B) pi_1[b] = -log(2) / hs[b];
pi_2 = -log(2) / hl ; 

for(b in 1:B)
 for(k in 1:K)
  for(t in 1:T)
   n[b,k,t] = mu[b,k] * (exp(pi_1[b] * time[t] + pi_2 * ts[b]) + exp(pi_2 * time[t] + pi_1[b] * ts[b])) / (exp(pi_1[b] * ts[b]) + exp(pi_2 * ts[b])) ; 

 for(k in 1:K)
  for(t in 1:T){
    n_C[1,k,t] =  n[1,k,t];
    n_C[2,k,t] =  n[2,k,t]; // MO and MU have SP titres 
    n_C[3,k,t] =  n[2,k,t];
  }
  
// MODEL SELECTION 
if(rho_K==1) {
  for(k in 1:K) {
    rhoT[k] = rho[k];
    gammaT[k] = gamma;
    phiT[k] = phi;
  }
} else {
  for(k in 1:K) {
    rhoT[k] = rho[1];
    gammaT[k] = gamma; 
    phiT[k] = phi;
  }
}
// rho is reduction in primary symp compared to sec
// gamma is prob sec symp but gammaT is prob primary (divided by rho)
    for(k in 1:K) gammaT[k] /= rhoT[k]; 
      
// INITIAL CONDITIONS
if(include_pK3 == 0) {  // p common to all serotypes
  for(k in 1:K){
    for(j in 1:J) p_KJ[k,j] = p[j]; // if include_pK3 == 0, all p are prob exposure to a single serotype 
    hK3[k] = 0;
    q[k] = 0;
    hK1[k] = 0;
    hK2[k] = 0;
    }
 } else { // serotype-specific p
  for(k in 1:K) hK3[k] = -log(1-pK3[k]) ; // pK3 prob exposure to serotype k 
  real shK3;
  shK3=sum(hK3);
  for(k in 1:K) {
    q[k] = hK3[k] / shK3 ;
    hK1[k] = q[k] * -log(1-p[1]) ; // if include_pK3 == 1, p1 and p2 is cumulative prob exposure 
    hK2[k] = q[k] * -log(1-p[2]) ;
    p_KJ[k,1] = 1 - exp(-hK1[k]) ;
    p_KJ[k,2] = 1 - exp(-hK2[k]) ;
    p_KJ[k,3] = pK3[k];
  }
}

// this defines the true serostatus populations 
 for(j in 1:J){
    TpSN[j] = (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1
    TpMO[1,j] = p_KJ[1,j]   * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 2
    TpMO[2,j] = p_KJ[2,j]   * (1-p_KJ[1,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 3
    TpMO[3,j] = p_KJ[3,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 4
    TpMO[4,j] = p_KJ[4,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 1,2
    TpMU2[1,j] = p_KJ[1,j] * p_KJ[2,j] * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1,3
    TpMU2[2,j] = p_KJ[1,j] * p_KJ[3,j] * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 1,4
    TpMU2[3,j] = p_KJ[1,j] * p_KJ[4,j] * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 2,3
    TpMU2[4,j] = p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[1,j]) * (1-p_KJ[4,j]);
    // 2,4
    TpMU2[5,j] = p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[3,j]);
    // 3,4 
    TpMU2[6,j] = p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[2,j]);
    // not 1 
    TpMU3[1,j] = p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) ;
    // not 2
    TpMU3[2,j] = p_KJ[1,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[2,j]) ;
    // not 3
    TpMU3[3,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[3,j]) ;
    // not 4 
    TpMU3[4,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[4,j]) ;
    // all 
    TpMU4[j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j];
 }

// this accounts for imperfect test performance 
for(v in 1:V)
 for(j in 1:J)
   for(k in 1:K){
    pSN[1,v,j,1] = spec * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec) * TpSN[j] ;
    pMO[1,v,k,j,1] = (1- sens) * TpMO[k,j] ;
    pMO[2,v,k,j,1] = sens * TpMO[k,j] ;
   }

if(MU_test_SN == 0){ // assume multitypic all classified SP
  for(v in 1:V)
   for(j in 1:J) {
    for(k in 1:6){  // there are 6 MU_2 combinations  
        pMU2[1,v,k,j,1] = 0;
        pMU2[2,v,k,j,1] = TpMU2[k,j] ; 
    }
    for(k in 1:K){  // there are 4 MU_3 combinations 
        pMU3[1,v,k,j,1] = 0; 
        pMU3[2,v,k,j,1] = TpMU3[k,j] ;
    } }
     } else { // allow mulitypic misclassification 
    for(v in 1:V)
     for(j in 1:J){
      for(k in 1:6){  // there are 6 MU_2 combinations  
          pMU2[1,v,k,j,1] = (1- sens) * TpMU2[k,j] ;
          pMU2[2,v,k,j,1] = sens * TpMU2[k,j] ; 
      }
      for(k in 1:K){  // there are 4 MU_3 combinations 
          pMU3[1,v,k,j,1] = (1- sens) * TpMU3[k,j] ; 
          pMU3[2,v,k,j,1] = sens * TpMU3[k,j] ;
      }}}

// probabilities of testing seropositive 
if(MU_test_SN == 0){  
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];
} else{
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * (sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j]) ;
}

// VACCINE RISK RATIO 

// outcome and age-group titre offsets 
real e_beta[2] = {exp(beta[1]), exp(beta[2])};
real e_alpha[C,K]; 

for(k in 1:K)
 for(c in 1:C){
   if(alpha_CK == 0){ 
     e_alpha[c,k] = exp(-alpha[1]);  // mono outcome offset 
     } else if(alpha_CK == 1) {
       e_alpha[c,k] = exp(-alpha[c]);  // serostatus outcome offset 
       } else {
         e_alpha[c,k] = exp(-alpha[k]);  // serotype outcome offset 
         }}
        
if(mono_lc_SN == 1){ // SN, oldest, symp 
  for(k in 1:K) nc50[1,k,3,1] = exp(lc[1,1]); 
   } else if (mono_lc_SN == 0) {
  for(k in 1:K) nc50[1,k,3,1] = exp(lc[1,k]); 
}

if(mono_lc_MU == 1){ // MU, oldest, symp
 for(k in 1:K)  nc50[3,k,3,1] = exp(lc[3,1]); 
   } else if(mono_lc_MU == 0) {
 for(k in 1:K) nc50[3,k,3,1] = exp(lc[3,k]); 
}


// MO, oldest, symp
for(k in 1:K) nc50[2,k,3,1] = exp(lc[2,k]); 

// age group offsets 
if(include_beta == 1){
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,3,1] * e_beta[1]; // youngest
   nc50[c,k,2,1] =   nc50[c,k,3,1] * e_beta[2]; // middle 
 } 
} else if(include_beta == 3) {
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,3,1] * e_beta[1]; // youngest
   nc50[c,k,2,1] =   nc50[c,k,3,1]; // middle 
   }
 } else {
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] = nc50[c,k,3,1];
   nc50[c,k,2,1] = nc50[c,k,3,1];   
  }
}

// hosp 
for(c in 1:C)
 for(k in 1:K) 
  for(j in 1:J)
    nc50[c,k,j,2] =   nc50[c,k,j,1] * e_alpha[c,k];

if(include_beta == 2){
  for(c in 1:C)
   for(k in 1:K){
     nc50[c,k,1,1] =   nc50[c,k,1,1] * e_beta[1]; 
     nc50[c,k,1,2] =   nc50[c,k,1,2] * e_beta[2]; 
   }
}

// Enhancement 
for(k in 1:K){
  if(enhancement == 1){ 
   if(L_K == 0){
   L_C[1,k,1] = 1 + L[1]; // symp 
   L_C[1,k,2] = 1 + L[1] * ((tau_K==0)?tau[1]:tau[k]);
  } else if (L_K == 1){
   L_C[1,k,1] = 1 + L[k]; // symp 
   L_C[1,k,2] = 1 + L[k] * ((tau_K==0)?tau[1]:tau[k]);
 } 
 } else if (enhancement == 0){ // no SN enhancement 
    for(r in 1:R)  L_C[1,k,r] = 1 ; 
 } 

 for(r in 1:R) { // no enhancement if seropositive 
        L_C[2,k,r] = 1; // MO
        L_C[3,k,r] = 1; // MU
      }
   }
 
// Risk ratios - w can be mono, serostatus or serotype dependent  
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J)
   for(t in 1:T){
     if(w_CK == 0){
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[1]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[1]) ;
     } else if(w_CK == 1){
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[c]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[c]) ;  
     } else{
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[k]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[k]) ;  
     }
     RR_symp[c,1,k,j,t] =  1 ;
     RR_hosp[c,1,k,j,t] =  1 ;
    }

// FOI  
for(k in 1:K){
  for(t in 1:12)  lambda[k,t] = exp(-lambda_D[k,1]);
  for(t in 13:18) lambda[k,t] = exp(-lambda_D[k,2]);
  for(t in 19:24) lambda[k,t] = exp(-lambda_D[k,3]);
  for(t in 25:36) lambda[k,t] = exp(-lambda_D[k,4]);
  for(t in 37:48) lambda[k,t] = exp(-lambda_D[k,5]);
  for(t in 49:54) lambda[k,t] = exp(-lambda_D[k,6]);
}

// SURVIVAL MODEL - prob of surviving each time point without infection 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(t in 1:T) pSN[b,v,j,(t+1)] = lambda[1,t]*lambda[2,t]*lambda[3,t]*lambda[4,t] * pSN[b,v,j,t];
   
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J) {
    for(t in 1:T) pMO[b,v,k,j,(t+1)] = lambda[1,t]*lambda[2,t]*lambda[3,t]*lambda[4,t]/lambda[k,t]* pMO[b,v,k,j,t] ; 
    for(t in (HI+1):T) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp[1,v,k,j,(t-HI)]) * (1- lambda[k,(t-HI)]) * pSN[b,v,j,(t-HI)] ;
    }  
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J) {
      for(t in 1:T){
          // MU_12
          pMU2[b,v,1,j,(t+1)] = lambda[3,t] * lambda[4,t] *  pMU2[b,v,1,j,t] ; 
           // MU_13
          pMU2[b,v,2,j,(t+1)] = lambda[2,t] * lambda[4,t] *  pMU2[b,v,2,j,t] ; 
           // MU_14
          pMU2[b,v,3,j,(t+1)] = lambda[2,t] * lambda[3,t] *  pMU2[b,v,3,j,t] ; 
          // MU_23
          pMU2[b,v,4,j,(t+1)] = lambda[1,t] * lambda[4,t] *  pMU2[b,v,4,j,t] ; 
          // MU_24
          pMU2[b,v,5,j,(t+1)] = lambda[1,t] * lambda[3,t] *  pMU2[b,v,5,j,t] ; 
          // MU_34
          pMU2[b,v,6,j,(t+1)] = lambda[1,t] * lambda[2,t] *  pMU2[b,v,6,j,t] ; 
      }
      for(t in (HI+1):T) {
          // MU_12
          pMU2[b,v,1,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_34
          pMU2[b,v,6,j,(t+1)] +=  (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
        }}

for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J){
     for(t in 1:T){ 
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] = lambda[1,t]  *  pMU3[b,v,1,j,t] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] = lambda[2,t]  *  pMU3[b,v,2,j,t] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] = lambda[3,t]  *  pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = lambda[4,t]  *  pMU3[b,v,4,j,t] ; 
     }
     for(t in (HI+1):T) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] += (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J) {
      for(t in 1:T){
         Inc3[b,v,1,j,t] = (1 - lambda[1,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,3,j,t] = (1 - lambda[3,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,5,j,t]);  
         Inc3[b,v,4,j,t] = (1 - lambda[4,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,2,j,t] + pMU2[b,v,4,j,t]);
         
         for(k in 1:K){
           Inc1[b,v,k,j,t] = (1 - lambda[k,t]) * pSN[b,v,j,t]; 
           Inc2[b,v,k,j,t] = (1 - lambda[k,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp[3,v,k,j,t];  
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * delta[k] * RR_hosp[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * delta[k]  * RR_hosp[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * delta[k] * epsilon * RR_hosp[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * delta[k] * RR_hosp[3,v,k,j,t] ;   
         }}}

// Aggregate to match published time points e.g. 1-12 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j,1] = sum(D1[b,v,k,j,1:12])  + sum(D2[b,v,k,j,1:12])   ;
   H[b,v,k,j,1]  = sum(H1[b,v,k,j,1:12])  + sum(H2[b,v,k,j,1:12])   ;
   Di[b,v,k,j,2] = sum(D1[b,v,k,j,13:18]) + sum(D2[b,v,k,j,13:18])  ;
   H[b,v,k,j,2]  = sum(H1[b,v,k,j,13:18]) + sum(H2[b,v,k,j,13:18])  ;
   Di[b,v,k,j,3] = sum(D1[b,v,k,j,19:24]) + sum(D2[b,v,k,j,19:24])  ;
   H[b,v,k,j,3]  = sum(H1[b,v,k,j,19:24]) + sum(H2[b,v,k,j,19:24])  ;
   Di[b,v,k,j,4] = sum(D1[b,v,k,j,25:36]) + sum(D2[b,v,k,j,25:36])  ;
   H[b,v,k,j,4]  = sum(H1[b,v,k,j,25:36]) + sum(H2[b,v,k,j,25:36])  ;
   Di[b,v,k,j,5] = sum(D1[b,v,k,j,37:48]) + sum(D2[b,v,k,j,37:48])  ;
   H[b,v,k,j,5]  = sum(H1[b,v,k,j,37:48]) + sum(H2[b,v,k,j,37:48])  ;
   Di[b,v,k,j,6] = sum(D1[b,v,k,j,49:54]) + sum(D2[b,v,k,j,49:54])  ;
   H[b,v,k,j,6]  = sum(H1[b,v,k,j,49:54]) + sum(H2[b,v,k,j,49:54])  ;
   }

if(MU_symp == 1){ // if post-secondary can be symptomatic 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j,1] += sum(D34[b,v,k,j,1:12])  ;
   H[b,v,k,j,1]  += sum(H34[b,v,k,j,1:12])  ;
   Di[b,v,k,j,2] += sum(D34[b,v,k,j,13:18])  ;
   H[b,v,k,j,2]  += sum(H34[b,v,k,j,13:18])  ;
   Di[b,v,k,j,3] += sum(D34[b,v,k,j,19:24])  ;
   H[b,v,k,j,3]  += sum(H34[b,v,k,j,19:24])  ;
   Di[b,v,k,j,4] += sum(D34[b,v,k,j,25:36])  ;
   H[b,v,k,j,4]  += sum(H34[b,v,k,j,25:36])  ;
   Di[b,v,k,j,5] += sum(D34[b,v,k,j,37:48])  ;
   H[b,v,k,j,5]  += sum(H34[b,v,k,j,37:48])  ;
   Di[b,v,k,j,6] += sum(D34[b,v,k,j,49:54])  ;
   H[b,v,k,j,6]  += sum(H34[b,v,k,j,49:54])  ;
   }}
   
for(v in 1:V)
 for(j in 1:J)
  for(d in 1:(D-1))
    pop_VJD[v,j,d] = sum(pop[ ,v,j,d]) ; // pop not by serostatus 

// Cases 
for(j in 1:J) {
   for(b in 1:B) {
    for(v in 1:V) {
      for(k in 1:K)
        for(d in 1:(D-1))
          C_BVKJRD[b,v,k,j,1,d] = Di[b,v,k,j,d] * pop_VJD[v,j,d]; // VCD up to month 48
          C_BVJ5[b,v,j] = sum(C_BVKJRD[b,v, ,j,1,5]) ;
    }
    C_BJ5[b,j] = sum(C_BVJ5[b, ,j]);
  }
  C_J5[j] = sum(C_BJ5[ ,j]);
}

real s_C_J5 = sum(C_J5);

for(j in 1:J) {
 pJ_VCD[j] = C_J5[j] / s_C_J5 ;
 for(v in 1:V) {
   for(b in 1:B) {
     for(d in 1:(D-1)) pop_BVJD[b,v,j,d] = pop[b,v,j,d]  ; // editable population 
     pop_BVJD[b,v,j,6] = pop[b,v,j,6] - (N_VCD_BV5[b,v] * pJ_VCD[j])  ; // censor population 
   }
   pop_VJD[v,j,6] = sum(pop_BVJD[ ,v,j,6]) ; // pop not by serostatus
   for(b in 1:B)
     for(k in 1:K)
       C_BVKJRD[b,v,k,j,1,6] = Di[b,v,k,j,6] *  pop_VJD[v,j,6]; // VCD at 54 months 
 }
}

for(j in 1:J) {
 for(k in 1:K) {
  for(v in 1:V) {
   for(b in 1:B) {
    for(d in 1:D)  
      C_BVKJRD[b,v,k,j,2,d] = H[b,v,k,j,d] * pop_VJD[v,j,d] ; // hosp
    for(r in 1:R) {
      C_BVKJR4[b,v,k,j,r,1] = sum(C_BVKJRD[b,v,k,j,r,1:3]) ;  // aggregate 1:24 months for hosp
      for(d in 2:4)
        C_BVKJR4[b,v,k,j,r,d] = C_BVKJRD[b,v,k,j,r,(d+2)] ; // 1 is 1:24, 2 is 36 (4), 3 is 48 (5), 4 is 54 (6)
    }
   }
  }
  for(r in 1:R)
   for(d in 1:D) {
     for(v in 1:V) C_VKJRD[v,k,j,r,d] = sum(C_BVKJRD[ ,v,k,j,r,d]); 
     C_KJRD[k,j,r,d] = sum(C_VKJRD[ ,k,j,r,d]) ; 
   }
 }
 for(r in 1:R)
  for(d in 1:D)
    C_JRD[j,r,d] = sum(C_KJRD[ ,j,r,d]);
}

for(r in 1:R)
 for(d in 1:D)
  C_RD[r,d] = sum(C_JRD[ ,r,d]) ; 

for(b in 1:B)
 for(v in 1:V)
  for(r in 1:R) {
   for(k in 1:K)
    for(d in 1:D)
     pC_BVKRD[b,v,k,r,d] = (sum(C_BVKJRD[b,v,k, ,r,d])) / C_RD[r,d];
   for(j in 1:J)
    for(a in 1:A)
     pC_BVJRA[b,v,j,r,a] = sum(C_BVKJRD[b,v, ,j,r,a]) / C_RD[r,a];
  }

// sum over year 2 for serotype and age cases
 for(k in 1:K)
  for(j in 1:J)
   for(r in 1:R) {
     pC_KJR2[k,j,r,1] = C_KJRD[k,j,r,1] / C_RD[r,1] ;
     pC_KJR2[k,j,r,2] = sum(C_KJRD[k,j,r,2:3]) / sum(C_RD[r,2:3]) ;
  }

// serotype hosp combined for year 1 and 2 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K) {
   pC_BVK4[b,v,k,1] = sum(C_BVKJR4[b,v,k, ,2,1]) / sum(C_RD[2,1:3]) ; // total hosp 1:24 months
   for(d in 2:4)
     pC_BVK4[b,v,k,d] = sum(C_BVKJR4[b,v,k, ,2,d]) / C_RD[2,(d+2)]; // 4, 5, 6
  }

// matrix for likelihood function 
for(j in 1:J)
  for(d in 1:2) {
   pD_KJ2[j,d]   = pC_KJR2[1,j,1,d];
   pD_KJ2[j+3,d] = pC_KJR2[2,j,1,d];
   pD_KJ2[j+6,d] = pC_KJR2[3,j,1,d];
   pD_KJ2[j+9,d] = pC_KJR2[4,j,1,d];
   pH_KJ2[j,d]   = pC_KJR2[1,j,2,d];
   pH_KJ2[j+3,d] = pC_KJR2[2,j,2,d];
   pH_KJ2[j+6,d] = pC_KJR2[3,j,2,d];
   pH_KJ2[j+9,d] = pC_KJR2[4,j,2,d];
 }

for(k in 1:K)
  for(d in 1:D){
   pD_BVKD[k,d]    =  pC_BVKRD[1,1,k,1,d];
   pD_BVKD[k+4,d]  =  pC_BVKRD[1,2,k,1,d]; 
   pD_BVKD[k+8,d]  =  pC_BVKRD[2,1,k,1,d]; 
   pD_BVKD[k+12,d] =  pC_BVKRD[2,2,k,1,d]; 
  }

// serotype hosp published for 4 time intervals
for(k in 1:K)
  for(d in 1:4){
   pH_BVK4[k,d]    =  pC_BVK4[1,1,k,d];
   pH_BVK4[k+4,d]  =  pC_BVK4[1,2,k,d]; 
   pH_BVK4[k+8,d]  =  pC_BVK4[2,1,k,d]; 
   pH_BVK4[k+12,d] =  pC_BVK4[2,2,k,d]; 
 }

for(j in 1:J)
  for(a in 1:A){
   pD_BVJA[j,a]     =  pC_BVJRA[1,1,j,1,a];
   pD_BVJA[j+3,a]   =  pC_BVJRA[1,2,j,1,a];
   pD_BVJA[j+6,a]   =  pC_BVJRA[2,1,j,1,a];
   pD_BVJA[j+9,a]   =  pC_BVJRA[2,2,j,1,a];
   pH_BVJA[j,a]     =  pC_BVJRA[1,1,j,2,a];
   pH_BVJA[j+3,a]   =  pC_BVJRA[1,2,j,2,a];
   pH_BVJA[j+6,a]   =  pC_BVJRA[2,1,j,2,a];
   pH_BVJA[j+9,a]   =  pC_BVJRA[2,2,j,2,a];
 }

for(d in 1:D)
 for(r in 1:R)
  pC_RD[r,d] = C_RD[r,d] / pop_D[d]; 

// log likelihood 
 ll=0;
 
 ll += binomial_lpmf(SP_J  | pop_J, pSP);
 ll += binomial_lpmf(VCD_D | pop_D, pC_RD[1,]);
 ll += binomial_lpmf(HOSP_D| pop_D, pC_RD[2,]);


   ll +=  multinomial_lpmf(VCD_BVKD[ ,1] | pD_BVKD[ ,1]);
  for(d in 2:D)   ll += multinomial_lpmf(VCD_BVKD[ ,d] | pD_BVKD[ ,d]);


 for(a in 1:A)    ll += multinomial_lpmf(VCD_BVJA[ ,a] | pD_BVJA[ ,a]) ;
 for(d in 1:2)    ll += multinomial_lpmf(VCD_KJ2[ ,d]  | pD_KJ2[ ,d]) ;
 for(d in 1:4)    ll += multinomial_lpmf(HOSP_BVK4[ ,d] | pH_BVK4[ ,d]);
 for(a in 1:A)    ll += multinomial_lpmf(HOSP_BVJA[ ,a] | pH_BVJA[ ,a]) ;
 for(d in 1:2)    ll += multinomial_lpmf(HOSP_KJ2[ ,d]  | pH_KJ2[ ,d]) ;

}}

model {
  target += ll;
// priors
  hs[1] ~ normal(1.65,0.5); // SN
  hs[2] ~ normal(4.20,0.5);  // SP
  hl ~ normal(84,12); // fits plus https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7557381/
  ts[1] ~ normal(-2.21,0.5) ;
  ts[2] ~ normal(0.15,0.5) ;
  p[1] ~ beta(3,5);
  p[2] ~ beta(3,5);
  
  if(include_pK3==1) {
  pK3 ~ beta(shape1,shape2) ; // prob exposure to serotype k 
  p[3] ~ beta(3,5) ;
  } else {
  pK3 ~ beta(3,5) ;
  p[3] ~ beta(shape1,shape2) ; // prob exposure to a serotype   
  }
  
  for(k in 1:K) lambda_D[k, ] ~ lognormal(-7,2);
  gamma ~ normal(0.85,0.2);  // prob secondary
  rho ~ normal(1.97,0.40);   // RR secondary
  delta ~ normal(0.25,0.10);
  phi ~ normal(0.25,0.05);
  epsilon ~ normal(1,1);
  L ~ normal(L_mean,L_sd);
  tau ~ normal(1,1);
  w ~ normal(1,2);
  lc[1,] ~ normal(4.5,1);
  lc[2,] ~ normal(6.5,1);
  lc[3,] ~ normal(6.5,1);
  alpha ~ normal(0,2);
  beta ~ normal(0,2);
  sens ~ normal(0.9,0.05) ;
  spec ~ normal(0.995, 0.01) ;
}

generated quantities{
  // AR
  real<lower = 0 > C_BVKRD[B,V,K,R,D];
  real<lower = 0 > C_BVJRD[B,V,J,R,D] ;
  real<lower = 0 > C_BVRD[B,V,R,D]    ;
  real<lower = 0 > AR_BVJRD[B,V,J,R,D] ;
  real<lower = 0 > AR_BVKRD[B,V,K,R,D] ;
  real<lower = 0 > AR_BVKHD[B,V,K,4] ;
  real<lower = 0 > AR_BVKR[B,V,K,R] ;
  real<lower = 0 > AR_BVJR[B,V,J,R] ;
  real<lower = 0 > AR_KJRD[K,J,R,2];
  real<lower = 0 > AR_BRD[B,R,D]   ;
  real<lower = 0 > AR_VRD[V,R,D]   ;
  // VE
  real<upper =1 > VE[C,K,J,R,T]   ;
  real<upper =1 > VE_BKRT[C,K,R,T]   ;
  real<upper =1 > VE_BJRT[C,J,R,T]   ;
  // LL
  vector[J+ 2*D+B*V*K*D+2*B*V*J*A+K*J*4+B*V*K*4] log_lik;
  // PPC
  array[2] int VCD_2;
  array[2] int HOSP_2;
  array[B*V*K, D] int pred_VCD_BVKD;
  array[B*V*J, A] int pred_VCD_BVJA;
  array[K*J, 2] int pred_VCD_KJ2;
  array[B*V*K,4] int pred_HOSP_BVK4;
  array[B*V*J, A] int pred_HOSP_BVJA;
  array[K*J, 2] int pred_HOSP_KJ2;

 {
// LL
vector[J] ll1;
vector[D] ll2;
vector[D] ll3;
vector[B*V*K*D] ll4;
vector[B*V*J*A] ll5;
vector[B*V*J*A] ll6;
vector[K*J*2] ll7;
vector[K*J*2] ll8;
vector[B*V*K*4] ll9;
vector[J] p_tmp1;
vector[D] p_tmp2;

p_tmp1 = to_vector(pSP);
ll1= v_SP_J .* log(p_tmp1) + (v_pop_J-v_SP_J) .* log(1-p_tmp1);
p_tmp2 = to_vector(pC_RD[1,]);
ll2= v_VCD_D .* log(p_tmp2) + (v_pop_D-v_VCD_D) .* log(1-p_tmp2);
p_tmp2 = to_vector(pC_RD[2,]);
ll3= v_HOSP_D .* log(p_tmp2) + (v_pop_D-v_HOSP_D) .* log(1-p_tmp2);

ll4 = to_vector(v_VCD_BVKD .* log(pD_BVKD));
ll5 = to_vector(v_VCD_BVJA .* log(pD_BVJA));
ll6 = to_vector(v_HOSP_BVJA .* log(pH_BVJA));
ll7 = to_vector(v_VCD_KJ2 .* log(pD_KJ2));
ll8 = to_vector(v_HOSP_KJ2 .* log(pH_KJ2));
ll9 = to_vector(v_HOSP_BVK4 .* log(pH_BVK4));

log_lik[1:J] = ll1;
log_lik[(J+1):(J+D)]=ll2;
log_lik[(J+D+1):(J+2*D)]=ll3;
log_lik[(J+2*D+1):(J+2*D+B*V*K*D)]=ll4;
log_lik[(J+2*D+B*V*K*D+1):(J+2*D+B*V*K*D+B*V*J*A)]=ll5;
log_lik[(J+2*D+B*V*K*D+B*V*J*A+1):(J+2*D+B*V*K*D+2*B*V*J*A)]=ll6;
log_lik[(J+2*D+B*V*K*D+2*B*V*J*A+1):(J+2*D+B*V*K*D+2*B*V*J*A+K*J*2)]=ll7;
log_lik[(J+2*D+B*V*K*D+2*B*V*J*A+K*J*2+1):(J+2*D+B*V*K*D+2*B*V*J*A+K*J*4)]=ll8;
log_lik[(J+2*D+B*V*K*D+2*B*V*J*A+K*J*4+1):(J+2*D+B*V*K*D+2*B*V*J*A+K*J*4+B*V*K*4)]=ll9;
}

// Attack rates
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(r in 1:R)
    for(d in 1:D)
    C_BVKRD[b,v,k,r,d] = sum(C_BVKJRD[b,v,k, ,r,d]) ;

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(r in 1:R)
    for(d in 1:D)
    C_BVJRD[b,v,j,r,d] = sum(C_BVKJRD[b,v, ,j,r,d]);

for(b in 1:B)
 for(v in 1:V)
  for(r in 1:R)
   for(d in 1:D)
    C_BVRD[b,v,r,d] = sum(C_BVKRD[b,v, ,r,d]) ;

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(r in 1:R)
    for(d in 1:D)
     AR_BVJRD[b,v,j,r,d] = C_BVJRD[b,v,j,r,d] / pop_BVJD[b,v,j,d];

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(r in 1:R)
    for(d in 1:D)
     AR_BVKRD[b,v,k,r,d] = C_BVKRD[b,v,k,r,d] / pop_BVD[b,v,d];

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(r in 1:R)
   AR_BVKR[b,v,k,r] = sum(C_BVKRD[b,v,k,r, ]) / pop_BV[b,v] ;

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(r in 1:R)
    AR_BVJR[b,v,j,r] = sum(C_BVJRD[b,v,j,r,1:4]) / mean(pop_BVJD[b,v,j,1:4]) ; // age data up to 36 months

for(v in 1:V)
 for(r in 1:R)
  for(d in 1:D)
   AR_VRD[v,r,d] = sum(C_BVRD[ ,v,r,d]) / pop_VD[v,d];

for(b in 1:B)
 for(r in 1:R)
  for(d in 1:D)
   AR_BRD[b,r,d] = sum(C_BVRD[b, ,r,d]) / pop_BD[b,d];
  
  for(k in 1:K)
   for(j in 1:J)
    for(r in 1:R){
   AR_KJRD[k,j,r,1] = C_KJRD[k,j,r,1] / pop_JD[j,1]; 
   AR_KJRD[k,j,r,2] = sum(C_KJRD[k,j,r,2:3]) / pop_JD[j,2]; 
}

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K){
    AR_BVKHD[b,v,k,1] = sum(C_BVKRD[b,v,k,2,1:3]) / pop_BVD[b,v,2];
    AR_BVKHD[b,v,k,2] = C_BVKRD[b,v,k,2,4] / pop_BVD[b,v,4];
    AR_BVKHD[b,v,k,3] = C_BVKRD[b,v,k,2,5] / pop_BVD[b,v,5];
    AR_BVKHD[b,v,k,4] = C_BVKRD[b,v,k,2,5] / pop_BVD[b,v,6];
    }
    
// VE
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J)
    for(t in 1:T) VE[c,k,j,1,t] = 1 - RR_symp[c,2,k,j,t] ;

for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J)
    for(t in 1:T) VE[c,k,j,2,t] = 1 - RR_hosp[c,2,k,j,t] ; 

for(c in 1:C)
 for(k in 1:K)
  for(r in 1:R)
   for(t in 1:T) VE_BKRT[c,k,r,t] = mean(VE[c,k, ,r,t]);

for(c in 1:C)
 for(j in 1:J)
  for(r in 1:R)
   for(t in 1:T) VE_BJRT[c,j,r,t] = mean(VE[c, ,j,r,t]);

// POSTERIOR PREDICTIVE CHECKS
VCD_2[1] = VCD_D[1];
VCD_2[2] = sum(VCD_D[2:3]);
HOSP_2[1] = HOSP_D[1];
HOSP_2[2] = sum(HOSP_D[2:3]);
for(d in 1:D)  pred_VCD_BVKD[ ,d]  = multinomial_rng(pD_BVKD[ ,d] , VCD_D[d] ) ;
for(a in 1:A)  pred_VCD_BVJA[ ,a]  = multinomial_rng(pD_BVJA[ ,a] , VCD_D[a] ) ;
for(d in 1:2)  pred_VCD_KJ2[ ,d]   = multinomial_rng(pD_KJ2 [ ,d] , VCD_2[d] ) ;
               pred_HOSP_BVK4[ ,1] = multinomial_rng(pH_BVK4[ ,1], sum(HOSP_2) ) ;
for(d in 2:4)  pred_HOSP_BVK4[ ,d] = multinomial_rng(pH_BVK4[ ,d] , HOSP_D[(d+2)]) ;
for(a in 1:A)  pred_HOSP_BVJA[ ,a] = multinomial_rng(pH_BVJA[ ,a] , HOSP_D[a]) ;
for(d in 1:2)  pred_HOSP_KJ2[ ,d]  = multinomial_rng(pH_KJ2[ ,d] , HOSP_2[d]) ;
}
