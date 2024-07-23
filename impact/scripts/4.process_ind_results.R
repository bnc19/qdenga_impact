################################################################################
# Code to calculate prop number of cases averted per million vaccinees
# in first vaccinated cohort 
################################################################################

## set up folder structure 
pathin <- "demo/" 
pathout <- "demo/processed"


### scenarios 
countries <- c("BRA")
ages <- c(6)
ve_type <- c( "VS_D15") 
sp9 = c(2,4,6,8)


## create folders
dir.create(paste0(pathout, "/" , ve_type, "/cohort_0.8vc_prop_averted"), showWarnings  = F)
dir.create(paste0(pathout, "/" , ve_type, "/cohort_0.8vc_prop_averted/csv"), showWarnings = F)


# set up matrices 
# empty matrix 
res <- matrix(nrow=9,ncol=51)
# rows are %sp9
res[,1] <- 100-10*(1:9)


res.l <- res
res.u <- res

resm <- res
resm.l <- res
resm.u <- res

resd <- res
resd.l <- res
resd.u <- res


# age structure
st <- c(1:10,21:30,41:50,61:70)

# indexes for serotype outputs
idx1 <- seq(12,48,4)
idx2 <- seq(13,49,4)
idx3 <- seq(14,50,4)
idx4 <- seq(15,51,4)

#idx for multipliers
sdis_idx <- c(2,3,4,5, 10, seq(12,27,1), seq(44,47,1))
dis_idx <- c(6,7,8,9, 11, seq(28,43,1), seq(48,51,1))


# prop sdis over all serotypes over 10 years
vars <- c("out_cohort_sdis_sero_vacc_pri",
          "out_cohort_sdis_sero_vacc_secp",
          "out_cohort_sdis_sero_vacc_pos",
          "out_cohort_sdis_sero_vacc_neg", 
          "out_cohort_dis_sero_vacc_pri",
          "out_cohort_dis_sero_vacc_secp",
          "out_cohort_dis_sero_vacc_pos",
          "out_cohort_dis_sero_vacc_neg")

# create serotype specific cols 
cols <- paste0(rep(c(vars, "out_cohort_sdis_sero_vacc", "out_cohort_dis_sero_vacc"), each = 4), "_D", rep(seq(1,4,1),11))

# overall plus serotype 
vars_out <- c(vars, "out_cohort_sdis_sero_vacc", "out_cohort_dis_sero_vacc", cols)

vars_out_mat <- c("p", gsub("sero", "all", vars_out[1:10]), vars_out[11:50])

for(h in 1:length(countries)){
  
  for(a in 1:length(ages)){
    
    for(t in 1:length(ve_type)){
      
      pout_csv <- paste0(pathout,"/",ve_type[t], "/cohort_0.8vc_prop_averted/csv/") 
      
      for(i in sp9) {
        
        # lists of 83 elements - with and without vaccine 
        test <- readRDS(paste0(pathin, countries[h],"/vacc/vc_0.8/",ve_type[t],"/vacc_res_", ages[a], "_R0",10-i,".rds"))
        ntest <- readRDS(paste0(pathin, countries[h],"/no_vacc/vc_0.8/no_vacc_res_", ages[a], "_R0",10-i,".rds"))
        
        # vaccinated in Y1 - denominator
        nd <- as.vector(ntest$out_cohort_nvacc_neg[1,,]+ntest$out_cohort_nvacc_pos[1,,])
        ndm <- ntest$out_cohort_nvacc_neg[1,,]+ntest$out_cohort_nvacc_pos[1,,]
        
        for(k in 1:length(vars)){
   
          # find element of n test equal to vars, get age structure, sum over columns serotypes
          tnn <- as.vector(colSums(ntest[[vars[k]]][st,,],dims = 1))  
          tvn <- as.vector(colSums(test[[vars[k]]][st,,],dims = 1))
          
      
          # proportion of cases averted (mean and 95% CrI)
          qn <- quantile(1-tvn/tnn, probs = c(0.025,0.975), na.rm = TRUE)
          res[i,k+1] <- mean(1-tvn/tnn, na.rm = TRUE)
          res.l[i,k+1] <- qn[1]
          res.u[i,k+1] <- qn[2]
          
          # mean variability  
          dat <- 1- apply(test[[vars[k]]][st,,], c(2,3), sum)/apply(ntest[[vars[k]]][st,,], c(2,3), sum)
          mu <- apply(dat, 2, mean)
          qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
          resm[i,k+1] <- mean(mu, na.rm = TRUE)
          resm.l[i,k+1] <- qn[1]
          resm.u[i,k+1] <- qn[2]
          
          # dynamic variability
          mu <- apply(dat, 1, mean)
          qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
          resd[i,k+1] <- mean(mu, na.rm = TRUE)
          resd.l[i,k+1] <- qn[1]
          resd.u[i,k+1] <- qn[2]
          
          if(k == 3){  # vars = out_cohort_sdis_sero_vacc_pos
            tnn1 <- (tnn + as.vector(colSums(ntest[[vars[k+1]]][st,,],dims = 1)))
            tvn1 <- (tvn + as.vector(colSums(test[[vars[k+1]]][st,,],dims = 1)))
            
            qn <- quantile(1-tvn1/tnn1, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,10] <- mean(1-tvn1/tnn1, na.rm = TRUE)
            res.l[i,10] <- qn[1]
            res.u[i,10] <- qn[2]
            
            # mean variability
            dat <- 1-apply(test[[vars[k]]][st,,]+test[[vars[k+1]]][st,,], c(2,3), sum)/apply(ntest[[vars[k]]][st,,]+ntest[[vars[k+1]]][st,,], c(2,3), sum)
            mu <- apply(dat, 2, mean, na.rm = TRUE)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,10] <- mean(mu)
            resm.l[i,10] <- qn[1]
            resm.u[i,10] <- qn[2]
            
            # dynamic variability
            mu <- apply(dat, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,10] <- mean(mu, na.rm = TRUE)
            resd.l[i,10] <- qn[1]
            resd.u[i,10] <- qn[2]
            
            # by serotype
            # overall variability
            tnn.1 <- (as.vector(colSums(ntest[[vars[k]]][1:10,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][1:10,,],dims = 1)))
            tnn.2 <- (as.vector(colSums(ntest[[vars[k]]][21:30,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][21:30,,],dims = 1)))
            tnn.3 <- (as.vector(colSums(ntest[[vars[k]]][41:50,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][41:50,,],dims = 1)))
            tnn.4 <- (as.vector(colSums(ntest[[vars[k]]][61:70,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][61:70,,],dims = 1))) 
            
            tvn.1 <- (as.vector(colSums(test[[vars[k]]][1:10,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][1:10,,],dims = 1)))
            tvn.2 <- (as.vector(colSums(test[[vars[k]]][21:30,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][21:30,,],dims = 1)))
            tvn.3 <- (as.vector(colSums(test[[vars[k]]][41:50,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][41:50,,],dims = 1)))
            tvn.4 <- (as.vector(colSums(test[[vars[k]]][61:70,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][61:70,,],dims = 1)))
            
            qn <- quantile(1-tvn.1/tnn.1, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx1[6+k]] <- mean(1-tvn.1/tnn.1, na.rm = TRUE)
            res.l[i,idx1[6+k]] <- qn[1]
            res.u[i,idx1[6+k]] <- qn[2]
            
            qn <- quantile(1-tvn.2/tnn.2, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx2[6+k]] <- mean(1-tvn.2/tnn.2, na.rm = TRUE)
            res.l[i,idx2[6+k]] <- qn[1]
            res.u[i,idx2[6+k]] <- qn[2]
            
            qn <- quantile(1-tvn.3/tnn.3, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx3[6+k]] <- mean(1-tvn.3/tnn.3, na.rm = TRUE)
            res.l[i,idx3[6+k]] <- qn[1]
            res.u[i,idx3[6+k]] <- qn[2]
            
            qn <- quantile(1-tvn.4/tnn.4, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx4[6+k]] <- mean(1-tvn.4/tnn.4, na.rm = TRUE)
            res.l[i,idx4[6+k]] <- qn[1]
            res.u[i,idx4[6+k]] <- qn[2]
            
            # mean variability
            dat1 <- 1-apply(test[[vars[k]]][1:10,,]+test[[vars[k+1]]][1:10,,], c(2,3), sum)/apply(ntest[[vars[k]]][1:10,,]+ntest[[vars[k+1]]][1:10,,], c(2,3), sum)
            dat2 <- 1-apply(test[[vars[k]]][21:30,,]+test[[vars[k+1]]][21:30,,], c(2,3), sum)/apply(ntest[[vars[k]]][21:30,,]+ntest[[vars[k+1]]][21:30,,], c(2,3), sum)
            dat3 <- 1-apply(test[[vars[k]]][41:50,,]+test[[vars[k+1]]][41:50,,], c(2,3), sum)/apply(ntest[[vars[k]]][41:50,,]+ntest[[vars[k+1]]][41:50,,], c(2,3), sum)
            dat4 <- 1-apply(test[[vars[k]]][61:70,,]+test[[vars[k+1]]][61:70,,], c(2,3), sum)/apply(ntest[[vars[k]]][61:70,,]+ntest[[vars[k+1]]][61:70,,], c(2,3), sum)
            
            mu <- apply(dat1, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx1[6+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx1[6+k]] <- qn[1]
            resm.u[i,idx1[6+k]] <- qn[2]
            
            mu <- apply(dat2, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx2[6+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx2[6+k]] <- qn[1]
            resm.u[i,idx2[6+k]] <- qn[2]
            
            mu <- apply(dat3, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx3[6+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx3[6+k]] <- qn[1]
            resm.u[i,idx3[6+k]] <- qn[2]
            
            mu <- apply(dat4, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx4[6+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx4[6+k]] <- qn[1]
            resm.u[i,idx4[6+k]] <- qn[2]
            
            # dynamic variability
            mu <- apply(dat1, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx1[6+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx1[6+k]] <- qn[1]
            resd.u[i,idx1[6+k]] <- qn[2]
            
            mu <- apply(dat2, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx2[6+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx2[6+k]] <- qn[1]
            resd.u[i,idx2[6+k]] <- qn[2]
            
            mu <- apply(dat3, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx3[6+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx3[6+k]] <- qn[1]
            resd.u[i,idx3[6+k]] <- qn[2]
            
            mu <- apply(dat4, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx4[6+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx4[6+k]] <- qn[1]
            resd.u[i,idx4[6+k]] <- qn[2]
            
          }else if(k == 7){
            # all_vacc
            tnn1 <- (tnn + as.vector(colSums(ntest[[vars[k+1]]][st,,],dims = 1)))
            tvn1 <- (tvn + as.vector(colSums(test[[vars[k+1]]][st,,],dims = 1)))
            
            qn <- quantile(1-tvn1/tnn1, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,11] <- mean(1-tvn1/tnn1, na.rm = TRUE)
            res.l[i,11] <- qn[1]
            res.u[i,11] <- qn[2]
            
            # mean variability
            dat <- 1-apply(test[[vars[k]]][st,,]+test[[vars[k+1]]][st,,], c(2,3), sum)/apply(ntest[[vars[k]]][st,,]+ntest[[vars[k+1]]][st,,], c(2,3), sum)
            mu <- apply(dat, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,11] <- mean(mu, na.rm = TRUE)
            resm.l[i,11] <- qn[1]
            resm.u[i,11] <- qn[2]
            
            # dynamic variability
            mu <- apply(dat, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,11] <- mean(mu, na.rm = TRUE)
            resd.l[i,11] <- qn[1]
            resd.u[i,11] <- qn[2]
            
            # by serotype
            tnn.1 <- (as.vector(colSums(ntest[[vars[k]]][1:10,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][1:10,,],dims = 1)))
            tnn.2 <- (as.vector(colSums(ntest[[vars[k]]][21:30,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][21:30,,],dims = 1)))
            tnn.3 <- (as.vector(colSums(ntest[[vars[k]]][41:50,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][41:50,,],dims = 1)))
            tnn.4 <- (as.vector(colSums(ntest[[vars[k]]][61:70,,],dims = 1))+
                        as.vector(colSums(ntest[[vars[k+1]]][61:70,,],dims = 1)))
            
            tvn.1 <- (as.vector(colSums(test[[vars[k]]][1:10,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][1:10,,],dims = 1)))
            tvn.2 <- (as.vector(colSums(test[[vars[k]]][21:30,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][21:30,,],dims = 1)))
            tvn.3 <- (as.vector(colSums(test[[vars[k]]][41:50,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][41:50,,],dims = 1)))
            tvn.4 <- (as.vector(colSums(test[[vars[k]]][61:70,,],dims = 1))+
                        as.vector(colSums(test[[vars[k+1]]][61:70,,],dims = 1)))
            
            qn <- quantile(1-tvn.1/tnn.1, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx1[3+k]] <- mean(1-tvn.1/tnn.1, na.rm = TRUE)
            res.l[i,idx1[3+k]] <- qn[1]
            res.u[i,idx1[3+k]] <- qn[2]
            
            qn <- quantile(1-tvn.2/tnn.2, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx2[3+k]] <- mean(1-tvn.2/tnn.2, na.rm = TRUE)
            res.l[i,idx2[3+k]] <- qn[1]
            res.u[i,idx2[3+k]] <- qn[2]
            
            qn <- quantile(1-tvn.3/tnn.3, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx3[3+k]] <- mean(1-tvn.3/tnn.3, na.rm = TRUE)
            res.l[i,idx3[3+k]] <- qn[1]
            res.u[i,idx3[3+k]] <- qn[2]
            
            qn <- quantile(1-tvn.4/tnn.4, probs = c(0.025,0.975), na.rm = TRUE)
            res[i,idx4[3+k]] <- mean(1-tvn.4/tnn.4, na.rm = TRUE)
            res.l[i,idx4[3+k]] <- qn[1]
            res.u[i,idx4[3+k]] <- qn[2]
            
            # mean variability
            dat1 <- 1-apply(test[[vars[k]]][1:10,,]+test[[vars[k+1]]][1:10,,], c(2,3), sum)/apply(ntest[[vars[k]]][1:10,,]+ntest[[vars[k+1]]][1:10,,], c(2,3), sum)
            dat2 <- 1-apply(test[[vars[k]]][21:30,,]+test[[vars[k+1]]][21:30,,], c(2,3), sum)/apply(ntest[[vars[k]]][21:30,,]+ntest[[vars[k+1]]][21:30,,], c(2,3), sum)
            dat3 <- 1-apply(test[[vars[k]]][41:50,,]+test[[vars[k+1]]][41:50,,], c(2,3), sum)/apply(ntest[[vars[k]]][41:50,,]+ntest[[vars[k+1]]][41:50,,], c(2,3), sum)
            dat4 <- 1-apply(test[[vars[k]]][61:70,,]+test[[vars[k+1]]][61:70,,], c(2,3), sum)/apply(ntest[[vars[k]]][61:70,,]+ntest[[vars[k+1]]][61:70,,], c(2,3), sum)
            
            mu <- apply(dat1, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx1[3+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx1[3+k]] <- qn[1]
            resm.u[i,idx1[3+k]] <- qn[2]
            
            mu <- apply(dat2, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx2[3+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx2[3+k]] <- qn[1]
            resm.u[i,idx2[3+k]] <- qn[2]
            
            mu <- apply(dat3, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx3[3+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx3[3+k]] <- qn[1]
            resm.u[i,idx3[3+k]] <- qn[2]
            
            mu <- apply(dat4, 2, mean)
            qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
            resm[i,idx4[3+k]] <- mean(mu, na.rm = TRUE)
            resm.l[i,idx4[3+k]] <- qn[1]
            resm.u[i,idx4[3+k]] <- qn[2]
            
            # dynamic variability
            mu <- apply(dat1, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx1[3+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx1[3+k]] <- qn[1]
            resd.u[i,idx1[3+k]] <- qn[2]
            
            mu <- apply(dat2, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx2[3+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx2[3+k]] <- qn[1]
            resd.u[i,idx2[3+k]] <- qn[2]
            
            mu <- apply(dat3, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx3[3+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx3[3+k]] <- qn[1]
            resd.u[i,idx3[3+k]] <- qn[2]
            
            mu <- apply(dat4, 1, mean)
            qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
            resd[i,idx4[3+k]] <- mean(mu, na.rm = TRUE)
            resd.l[i,idx4[3+k]] <- qn[1]
            resd.u[i,idx4[3+k]] <- qn[2]
            
          }
          
          # by serotype
          # overall variables 
          tnn.1 <- as.vector(colSums(ntest[[vars[k]]][1:10,,],dims = 1))
          tnn.2 <- as.vector(colSums(ntest[[vars[k]]][21:30,,],dims = 1))
          tnn.3 <- as.vector(colSums(ntest[[vars[k]]][41:50,,],dims = 1))
          tnn.4 <- as.vector(colSums(ntest[[vars[k]]][61:70,,],dims = 1)) 
          
          tvn.1 <- as.vector(colSums(test[[vars[k]]][1:10,,],dims = 1))
          tvn.2 <- as.vector(colSums(test[[vars[k]]][21:30,,],dims = 1))
          tvn.3 <- as.vector(colSums(test[[vars[k]]][41:50,,],dims = 1))
          tvn.4 <- as.vector(colSums(test[[vars[k]]][61:70,,],dims = 1))
          
          qn <- quantile(1-tvn.1/tnn.1, probs = c(0.025,0.975), na.rm = TRUE)
          res[i,idx1[k]] <- mean(1-tvn.1/tnn.1, na.rm = TRUE)
          res.l[i,idx1[k]] <- qn[1]
          res.u[i,idx1[k]] <- qn[2]
          
          qn <- quantile(1-tvn.2/tnn.2, probs = c(0.025,0.975), na.rm = TRUE)
          res[i,idx2[k]] <- mean(1-tvn.2/tnn.2, na.rm = TRUE)
          res.l[i,idx2[k]] <- qn[1]
          res.u[i,idx2[k]] <- qn[2]
          
          qn <- quantile(1-tvn.3/tnn.3, probs = c(0.025,0.975), na.rm = TRUE)
          res[i,idx3[k]] <- mean(1-tvn.3/tnn.3, na.rm = TRUE)
          res.l[i,idx3[k]] <- qn[1]
          res.u[i,idx3[k]] <- qn[2]
          
          qn <- quantile(1-tvn.4/tnn.4, probs = c(0.025,0.975), na.rm = TRUE)
          res[i,idx4[k]] <- mean(1-tvn.4/tnn.4, na.rm = TRUE)
          res.l[i,idx4[k]] <- qn[1]
          res.u[i,idx4[k]] <- qn[2]
          
          # mean variability
          dat1 <- 1-apply(test[[vars[k]]][1:10,,], c(2,3), sum)/apply(ntest[[vars[k]]][1:10,,], c(2,3), sum)
          dat2 <- 1-apply(test[[vars[k]]][21:30,,], c(2,3), sum)/apply(ntest[[vars[k]]][21:30,,], c(2,3), sum)
          dat3 <- 1-apply(test[[vars[k]]][41:50,,], c(2,3), sum)/apply(ntest[[vars[k]]][41:50,,], c(2,3), sum)
          dat4 <- 1-apply(test[[vars[k]]][61:70,,], c(2,3), sum)/apply(ntest[[vars[k]]][61:70,,], c(2,3), sum)
          
          mu <- apply(dat1, 2, mean)
          qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
          resm[i,idx1[k]] <- mean(mu, na.rm = TRUE)
          resm.l[i,idx1[k]] <- qn[1]
          resm.u[i,idx1[k]] <- qn[2]
          
          mu <- apply(dat2, 2, mean)
          qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
          resm[i,idx2[k]] <- mean(mu, na.rm = TRUE)
          resm.l[i,idx2[k]] <- qn[1]
          resm.u[i,idx2[k]] <- qn[2]
          
          mu <- apply(dat3, 2, mean)
          qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
          resm[i,idx3[k]] <- mean(mu, na.rm = TRUE)
          resm.l[i,idx3[k]] <- qn[1]
          resm.u[i,idx3[k]] <- qn[2]
          
          mu <- apply(dat4, 2, mean)
          qn <- quantile(mu, probs = c(0.025,0.975), na.rm = TRUE)
          resm[i,idx4[k]] <- mean(mu, na.rm = TRUE)
          resm.l[i,idx4[k]] <- qn[1]
          resm.u[i,idx4[k]] <- qn[2]
          
          # dynamic variability
          mu <- apply(dat1, 1, mean)
          qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
          resd[i,idx1[k]] <- mean(mu, na.rm = TRUE)
          resd.l[i,idx1[k]] <- qn[1]
          resd.u[i,idx1[k]] <- qn[2]
          
          mu <- apply(dat2, 1, mean)
          qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
          resd[i,idx2[k]] <- mean(mu, na.rm = TRUE)
          resd.l[i,idx2[k]] <- qn[1]
          resd.u[i,idx2[k]] <- qn[2]
          
          mu <- apply(dat3, 1, mean)
          qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
          resd[i,idx3[k]] <- mean(mu, na.rm = TRUE)
          resd.l[i,idx3[k]] <- qn[1]
          resd.u[i,idx3[k]] <- qn[2]
          
          mu <- apply(dat4, 1, mean)
          qn <- quantile(mu, probs = c(0.05,0.95), na.rm = TRUE)
          resd[i,idx4[k]] <- mean(mu, na.rm = TRUE)
          resd.l[i,idx4[k]] <- qn[1]
          resd.u[i,idx4[k]] <- qn[2]
          
        }
      }
      
      # add column names 
      colnames(res) <- vars_out_mat
      colnames(res.l) <- vars_out_mat
      colnames(res.u) <- vars_out_mat
      
      colnames(resm) <- vars_out_mat
      colnames(resm.l) <- vars_out_mat
      colnames(resm.u) <- vars_out_mat
      
      colnames(resd) <- vars_out_mat
      colnames(resd.l) <- vars_out_mat
      colnames(resd.u) <- vars_out_mat
      
      tb <- rbind(res, res.l, res.u)
      rownames(tb) <- c(rep("mean", 9), rep("lowerb", 9), rep("upperb", 9))
      write.csv(tb, paste0(pout_csv, "cohort_prop_averted_overall_uncert_", countries[h], "_", ages[a], "_95perc.csv"), quote = FALSE, row.names = TRUE)
      
      tb <- rbind(resm, resm.l, resm.u)
      rownames(tb) <- c(rep("mean", 9), rep("lowerb", 9), rep("upperb", 9))
      write.csv(tb, paste0(pout_csv, "cohort_prop_averted_mean_uncert_", countries[h], "_", ages[a],"_95perc.csv"), quote = FALSE, row.names = TRUE)

      tb <- rbind(resd, resd.l, resd.u)
      rownames(tb) <- c(rep("mean", 9), rep("lowerb", 9), rep("upperb", 9))
      write.csv(tb, paste0(pout_csv, "cohort_prop_averted_dynam_uncert_", countries[h], "_", ages[a],"_90perc.csv"), quote = FALSE, row.names = TRUE)
      
    }

  }
}

