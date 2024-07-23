################################################################################
# This code outputs the outcome in terms of prop of dis/sdis averted           #
# yearly for the first 20 years and cumulative over 5, 10 and 20 years         #
################################################################################

library(matrixStats)

# setwd(paste0(getwd(), "/impact"))

# source functions
source("R/process_simulations.R")

## set up folder structure 
pathin <- "demo/" 
pathout <- "demo/processed"
folder_type <- "/VS_D15/"  
dis.scale <- 1.2956
sdis.scale <- c(0.6833, 1.0553) # 9% hospitalisation 
new_folder <- paste0(dis.scale, "_", sdis.scale[1], "_", sdis.scale[2])

# create folders 
dir.create(pathout, showWarnings = F)
dir.create(paste0(pathout, folder_type), showWarnings = F)
dir.create(file.path(paste0(pathout, folder_type, new_folder), "csv"), showWarnings = F)

#### scenario set set up 

# set to match demo runs, but for the paper the scenarios run were:

#  - Assumed demography: Brazil or Phillipines
#  - Vaccine mechanism of protection: only disease or also against infection 
#  - Duration of vaccine decay: 15 years or 5 years
#  - Age of vaccination: 6 to 12 years
#  - Transmission intensity setting (average 9 years seroprevalence): 10% to 90%
#  - Vaccine coverage: 20%, 40%, 60%, 80%
#  - Pre-vaccination screening to identify seropositive individuals: yes/no

countries <- c("BRA")
vacc.age <- c(6)
vacc.cov <- c(0.8)
cumy <- c(5,10,20)
R0 <- c("R02","R04","R06","R08")
serop9 <- c("sp9 20%", "sp9 40%", "sp9 60%", "sp9 80%")

# populations 
v_list <- c ("dis_all_vacc_neg","dis_all_vacc_pos","dis_all_vacc","dis_all_pop","dis_all_vacc_pri", "dis_all_vacc_secp",
             "sdis_all_vacc_neg","sdis_all_vacc_pos","sdis_all_vacc","sdis_all_pop","sdis_all_vacc_pri","sdis_all_vacc_secp")

v_list_sero <- c ("dis_sero_vacc_neg","dis_sero_vacc_pos","dis_sero_vacc","dis_sero_pop","dis_sero_vacc_pri", "dis_sero_vacc_secp",
                  "sdis_sero_vacc_neg","sdis_sero_vacc_pos","sdis_sero_vacc","sdis_sero_pop","sdis_sero_vacc_pri", "sdis_sero_vacc_secp")

v_list_he <- list("out_disc_dis_all_pop", "out_disc_sdis_all_pop", "out_disc_yll_all_pop", "out_disc_nvacc_all_pop",
               "out_dis_all_pop", "out_sdis_all_pop", "out_yll_all_pop", "out_nvacc_all_pop")



################################################################################
# calculate impacts - loop over all files (countries, vacc.cov, age.vacc, R0)
################################################################################

out.list <- list()

for(h in 1:length(countries)){ 
  
  pathin.v <- paste0(pathin, countries[h], "/vacc/")
  pathin.nv <- paste0(pathin, countries[h], "/no_vacc/")
  
  for(k in 1:length(vacc.cov)){
    
    pathin.vk <- paste0(pathin.v, "vc_", vacc.cov[k], folder_type)
    pathin.nvk <- paste0(pathin.nv, "vc_", vacc.cov[k])
    
    for(j in 1:length(vacc.age)){
       for(i in 1:length(R0)){
         
         scenario <- paste0(vacc.age[j],"_", R0[i])
         
         # read in results 
         
         #  List corresponding to each outcome 
         #  Each element of the list is at minimum a matrix of 50 200
         #  Corresponding to 50 simulations of 200 draws
 
         result.nv <- readRDS(paste0(pathin.nvk,"/no_vacc_res_", scenario,".rds"))
         result.v <- readRDS(paste0(pathin.vk,"/vacc_res_", scenario,".rds"))
         
        for(v in v_list) {
          a <- paste0("out_",v)

          
          if(v %in% c("dis_all_vacc_neg",
                      "dis_all_vacc_pos",
                      "dis_all_vacc",
                      "dis_all_pop",
                      "dis_all_vacc_pri")){
            
            # yearly outcome (first 20 years)
            t.v <- dis.scale*result.v[[a]]
            t.nv <- dis.scale*result.nv[[a]]
            
          }else{
            
            # yearly outcome (first 20 years)
            t.v <- dis.scale * sdis.scale[h] * result.v[[a]]
            t.nv <- dis.scale * sdis.scale[h] * result.nv[[a]]
            
          }

          # yearly stats on prop averted (first 20 years)
          out <- prop_averted(t.nv, t.v)
          res.diff.mean <- out[[1]]
          res.diff.l95 <- out[[2]]
          res.diff.q50 <- out[[3]]
          res.diff.u95 <- out[[4]]

          # cum prop averted at 5, 10 and 20 yrs
          res.cum.nv <- apply(t.nv, c(2,3), cumsum)[cumy, ,]
          res.cum.v <- apply(t.v, c(2,3), cumsum)[cumy, ,]

          # cum stats on prop averted at 5, 10 and 20 yrs
          out <- prop_averted(res.cum.nv, res.cum.v)
          res.cum.diff.mean <- out[[1]]
          res.cum.diff.l95 <- out[[2]]
          res.cum.diff.q50 <- out[[3]]
          res.cum.diff.u95 <- out[[4]]

          res <- list()
          res[[1]] <- cbind(res.diff.mean, res.diff.l95, res.diff.q50, res.diff.u95)
          res[[2]] <- cbind(res.cum.diff.mean, res.cum.diff.l95, res.cum.diff.q50, res.cum.diff.u95)
          
          out.list[[paste0("prop_averted_", v, "_", countries[h], "_", vacc.cov[k], "_", scenario)]] <- res

          # reconstruct secp - pos from simulations 
          if(v %in% c("dis_all_vacc_secp", "sdis_all_vacc_secp")){
          
            if(v == "dis_all_vacc_secp"){
              
              # secp - pos
              secp_minus_pos.v <- t.v  - dis.scale * result.v$out_dis_all_vacc_pos
              secp_minus_pos.nv <- t.nv - dis.scale * result.nv$out_dis_all_vacc_pos
              txt <- "dis_all_vacc_"
              
            }else{
              
              secp_minus_pos.v <- t.v  - dis.scale * sdis.scale[h] * (result.v$out_sdis_all_vacc_pos)
              secp_minus_pos.nv <- t.nv - dis.scale * sdis.scale[h] * (result.nv$out_sdis_all_vacc_pos)
              txt <- "sdis_all_vacc_"
              
            }
            
            out <- prop_averted(secp_minus_pos.v, secp_minus_pos.nv)
            textstring <- paste0(txt, "secp_minus_pos_")

            res.diff.mean <- out[[1]]
            res.diff.l95 <- out[[2]]
            res.diff.q50 <- out[[3]]
            res.diff.u95 <- out[[4]]
              
            # cum prop averted at 5, 10 and 20 yrs
            # replace NAs with 0s
            idx.na <- which(is.na(secp_minus_pos.nv))
            if(length(idx.na)>0) secp_minus_pos.nv[idx.na] <- 0
            
            idx.na <- which(is.na(secp_minus_pos.v))
            if(length(idx.na)>0) secp_minus_pos.v[idx.na] <- 0
            
            res.cum.nv <- apply(secp_minus_pos.nv, c(2,3), cumsum)[c(5,10,20), ,]
            res.cum.v <- apply(secp_minus_pos.v, c(2,3), cumsum)[c(5,10,20), ,]
              
            # cum stats on prop averted at 5, 10 and 20 yrs
            out <- prop_averted(res.cum.nv, res.cum.v)
            res.cum.diff.mean <- out[[1]]
            res.cum.diff.l95 <- out[[2]]
            res.cum.diff.q50 <- out[[3]]
            res.cum.diff.u95 <- out[[4]]
              
            res <- list()
            res[[1]] <- cbind(res.diff.mean, res.diff.l95, res.diff.q50, res.diff.u95)
            res[[2]] <- cbind(res.cum.diff.mean, res.cum.diff.l95, res.cum.diff.q50, res.cum.diff.u95)
              
            out.list[[paste0("prop_averted_", textstring, countries[h], "_", vacc.cov[k], "_", scenario)]] <- res

            }
          }

        for(v in v_list_sero) {
          a <- paste0("out_",v)
          
          if(v %in% c("dis_sero_vacc_neg",
                      "dis_sero_vacc_pos",
                      "dis_sero_vacc",
                      "dis_sero_pop",
                      "dis_sero_vacc_pri",
                      "dis_sero_vacc_secp")){
            
            # yearly outcome (first 20 years)
            t.v <- dis.scale*result.v[[a]]
            t.nv <- dis.scale*result.nv[[a]]
            
          }else{
            
            # yearly outcome (first 20 years)
            t.v <- dis.scale * sdis.scale[h] * result.v[[a]]
            t.nv <- dis.scale * sdis.scale[h] *result.nv[[a]]
            
          }

          # yearly stats on prop averted by serotype (first 20 years)
          out.1 <- prop_averted(t.nv[1:20,,], t.v[1:20,,])
          out.2 <- prop_averted(t.nv[21:40,,], t.v[21:40,,])
          out.3 <- prop_averted(t.nv[41:60,,], t.v[41:60,,])
          out.4 <- prop_averted(t.nv[61:80,,], t.v[61:80,,])

          res.diff.mean.1 <- out.1[[1]]
          res.diff.mean.2 <- out.2[[1]]
          res.diff.mean.3 <- out.3[[1]]
          res.diff.mean.4 <- out.4[[1]]

          res.diff.l95.1 <- out.1[[2]]
          res.diff.l95.2 <- out.2[[2]]
          res.diff.l95.3 <- out.3[[2]]
          res.diff.l95.4 <- out.4[[2]]

          res.diff.q50.1 <- out.1[[3]]
          res.diff.q50.2 <- out.2[[3]]
          res.diff.q50.3 <- out.3[[3]]
          res.diff.q50.4 <- out.4[[3]]

          res.diff.u95.1 <- out.1[[4]]
          res.diff.u95.2 <- out.2[[4]]
          res.diff.u95.3 <- out.3[[4]]
          res.diff.u95.4 <- out.4[[4]]

          # cum prop averted at 5, 10 and 20 yrs
          res.cum.nv.1 <- apply(t.nv[1:20,,], c(2,3), cumsum)[cumy,,]
          res.cum.nv.2 <- apply(t.nv[21:40,,], c(2,3), cumsum)[cumy,,]
          res.cum.nv.3 <- apply(t.nv[41:60,,], c(2,3), cumsum)[cumy,,]
          res.cum.nv.4 <- apply(t.nv[61:80,,], c(2,3), cumsum)[cumy,,]

          res.cum.v.1 <- apply(t.v[1:20,,], c(2,3), cumsum)[cumy,,]
          res.cum.v.2 <- apply(t.v[21:40,,], c(2,3), cumsum)[cumy,,]
          res.cum.v.3 <- apply(t.v[41:60,,], c(2,3), cumsum)[cumy,,]
          res.cum.v.4 <- apply(t.v[61:80,,], c(2,3), cumsum)[cumy,,]

          # cum stats on prop averted at 5, 10 and 20 yrs
          out.1 <- prop_averted(res.cum.nv.1, res.cum.v.1)
          out.2 <- prop_averted(res.cum.nv.2, res.cum.v.2)
          out.3 <- prop_averted(res.cum.nv.3, res.cum.v.3)
          out.4 <- prop_averted(res.cum.nv.4, res.cum.v.4)

          res.cum.diff.mean.1 <- out.1[[1]]
          res.cum.diff.mean.2 <- out.2[[1]]
          res.cum.diff.mean.3 <- out.3[[1]]
          res.cum.diff.mean.4 <- out.4[[1]]

          res.cum.diff.l95.1 <- out.1[[2]]
          res.cum.diff.l95.2 <- out.2[[2]]
          res.cum.diff.l95.3 <- out.3[[2]]
          res.cum.diff.l95.4 <- out.4[[2]]

          res.cum.diff.q50.1 <- out.1[[3]]
          res.cum.diff.q50.2 <- out.2[[3]]
          res.cum.diff.q50.3 <- out.3[[3]]
          res.cum.diff.q50.4 <- out.4[[3]]

          res.cum.diff.u95.1 <- out.1[[4]]
          res.cum.diff.u95.2 <- out.2[[4]]
          res.cum.diff.u95.3 <- out.3[[4]]
          res.cum.diff.u95.4 <- out.4[[4]]

          res.s <- list()
          res.s[[1]] <- cbind(res.diff.mean.1, res.diff.l95.1, res.diff.q50.1, res.diff.u95.1,
                              res.diff.mean.2, res.diff.l95.2, res.diff.q50.2, res.diff.u95.2,
                              res.diff.mean.3, res.diff.l95.3, res.diff.q50.3, res.diff.u95.3,
                              res.diff.mean.4, res.diff.l95.4, res.diff.q50.4, res.diff.u95.4)

          res.s[[2]] <- cbind(res.cum.diff.mean.1, res.cum.diff.l95.1, res.cum.diff.q50.1, res.cum.diff.u95.1,
                              res.cum.diff.mean.2, res.cum.diff.l95.2, res.cum.diff.q50.2, res.cum.diff.u95.2,
                              res.cum.diff.mean.3, res.cum.diff.l95.3, res.cum.diff.q50.3, res.cum.diff.u95.3,
                              res.cum.diff.mean.4, res.cum.diff.l95.4, res.cum.diff.q50.4, res.cum.diff.u95.4)

          out.list[[paste0("prop_averted_", v, "_", countries[h], "_", vacc.cov[k], "_", scenario)]] <- res.s

          # estimate secp - pos from simulations 
          if(v %in% c("dis_sero_vacc_secp", "sdis_sero_vacc_secp")){
            
            if(v == "dis_sero_vacc_secp"){
              
              secp_minus_pos.v  <- dis.scale * (result.v$out_dis_sero_vacc_secp - result.v$out_dis_sero_vacc_pos)
              secp_minus_pos.nv <- dis.scale * (result.nv$out_dis_sero_vacc_secp - result.nv$out_dis_sero_vacc_pos)
              txt <- "dis_sero_vacc_"
              
            }else{
              
              secp_minus_pos.v  <- dis.scale * sdis.scale[h] * (result.v$out_dis_sero_vacc_secp - result.v$out_sdis_sero_vacc_pos)
              secp_minus_pos.nv <- dis.scale * sdis.scale[h] * (result.nv$out_dis_sero_vacc_secp - result.nv$out_sdis_sero_vacc_pos)
              txt <- "sdis_sero_vacc_"
              
            }

            # yearly stats on prop averted by serotype (first 20 years)
            out.1 <- prop_averted(secp_minus_pos.nv[1:20,,], secp_minus_pos.v[1:20,,])
            out.2 <- prop_averted(secp_minus_pos.nv[21:40,,], secp_minus_pos.v[21:40,,])
            out.3 <- prop_averted(secp_minus_pos.nv[41:60,,], secp_minus_pos.v[41:60,,])
            out.4 <- prop_averted(secp_minus_pos.nv[61:80,,], secp_minus_pos.v[61:80,,])
                
            textstring <- paste0(txt, "secp_minus_pos_")
      
            res.diff.mean.1 <- out.1[[1]]
            res.diff.mean.2 <- out.2[[1]]
            res.diff.mean.3 <- out.3[[1]]
            res.diff.mean.4 <- out.4[[1]]
            
            res.diff.l95.1 <- out.1[[2]]
            res.diff.l95.2 <- out.2[[2]]
            res.diff.l95.3 <- out.3[[2]]
            res.diff.l95.4 <- out.4[[2]]
            
            res.diff.q50.1 <- out.1[[3]]
            res.diff.q50.2 <- out.2[[3]]
            res.diff.q50.3 <- out.3[[3]]
            res.diff.q50.4 <- out.4[[3]]
            
            res.diff.u95.1 <- out.1[[4]]
            res.diff.u95.2 <- out.2[[4]]
            res.diff.u95.3 <- out.3[[4]]
            res.diff.u95.4 <- out.4[[4]]
            
            # cum prop averted at 5, 10 and 20 yrs
            idx.na <- which(is.na(secp_minus_pos.nv))
            if(length(idx.na)>0) secp_minus_pos.nv[idx.na] <- 0
            
            idx.na <- which(is.na(secp_minus_pos.v))
            if(length(idx.na)>0) secp_minus_pos.v[idx.na] <- 0
            
            res.cum.nv.1 <- apply(secp_minus_pos.nv[1:20,,], c(2,3), cumsum)[cumy,,]
            res.cum.nv.2 <- apply(secp_minus_pos.nv[21:40,,], c(2,3), cumsum)[cumy,,]
            res.cum.nv.3 <- apply(secp_minus_pos.nv[41:60,,], c(2,3), cumsum)[cumy,,]
            res.cum.nv.4 <- apply(secp_minus_pos.nv[61:80,,], c(2,3), cumsum)[cumy,,]
            
            res.cum.v.1 <- apply(secp_minus_pos.v[1:20,,], c(2,3), cumsum)[cumy,,]
            res.cum.v.2 <- apply(secp_minus_pos.v[21:40,,], c(2,3), cumsum)[cumy,,]
            res.cum.v.3 <- apply(secp_minus_pos.v[41:60,,], c(2,3), cumsum)[cumy,,]
            res.cum.v.4 <- apply(secp_minus_pos.v[61:80,,], c(2,3), cumsum)[cumy,,]
            
            # cum stats on prop averted at 5, 10 and 20 yrs
            out.1 <- prop_averted(res.cum.nv.1, res.cum.v.1)
            out.2 <- prop_averted(res.cum.nv.2, res.cum.v.2)
            out.3 <- prop_averted(res.cum.nv.3, res.cum.v.3)
            out.4 <- prop_averted(res.cum.nv.4, res.cum.v.4)
            
            res.cum.diff.mean.1 <- out.1[[1]]
            res.cum.diff.mean.2 <- out.2[[1]]
            res.cum.diff.mean.3 <- out.3[[1]]
            res.cum.diff.mean.4 <- out.4[[1]]
            
            res.cum.diff.l95.1 <- out.1[[2]]
            res.cum.diff.l95.2 <- out.2[[2]]
            res.cum.diff.l95.3 <- out.3[[2]]
            res.cum.diff.l95.4 <- out.4[[2]]
            
            res.cum.diff.q50.1 <- out.1[[3]]
            res.cum.diff.q50.2 <- out.2[[3]]
            res.cum.diff.q50.3 <- out.3[[3]]
            res.cum.diff.q50.4 <- out.4[[3]]
            
            res.cum.diff.u95.1 <- out.1[[4]]
            res.cum.diff.u95.2 <- out.2[[4]]
            res.cum.diff.u95.3 <- out.3[[4]]
            res.cum.diff.u95.4 <- out.4[[4]]
            
            res.s <- list()
            res.s[[1]] <- cbind(res.diff.mean.1, res.diff.l95.1, res.diff.q50.1, res.diff.u95.1,
                                res.diff.mean.2, res.diff.l95.2, res.diff.q50.2, res.diff.u95.2,
                                res.diff.mean.3, res.diff.l95.3, res.diff.q50.3, res.diff.u95.3,
                                res.diff.mean.4, res.diff.l95.4, res.diff.q50.4, res.diff.u95.4)
            
            res.s[[2]] <- cbind(res.cum.diff.mean.1, res.cum.diff.l95.1, res.cum.diff.q50.1, res.cum.diff.u95.1,
                                res.cum.diff.mean.2, res.cum.diff.l95.2, res.cum.diff.q50.2, res.cum.diff.u95.2,
                                res.cum.diff.mean.3, res.cum.diff.l95.3, res.cum.diff.q50.3, res.cum.diff.u95.3,
                                res.cum.diff.mean.4, res.cum.diff.l95.4, res.cum.diff.q50.4, res.cum.diff.u95.4)
            
            out.list[[paste0("prop_averted_", textstring, countries[h], "_", vacc.cov[k], "_", scenario)]] <- res.s
            }
          }
        
     
      }
    }
  }
}

saveRDS(out.list, paste0(pathout, folder_type, dis.scale, "_", sdis.scale[1], "_", sdis.scale[2], "_all_results.rds"))

################################################################################
# create data frames for plotting impacts 
################################################################################

v_list <- c ("dis_all_vacc_neg","dis_all_vacc_pos","dis_all_vacc","dis_all_pop","dis_all_vacc_pri", "dis_all_vacc_secp", "dis_all_vacc_secp_minus_pos",
             "sdis_all_vacc_neg","sdis_all_vacc_pos","sdis_all_vacc","sdis_all_pop","sdis_all_vacc_pri", "sdis_all_vacc_secp", "sdis_all_vacc_secp_minus_pos")

v_list_sero <- c ("dis_sero_vacc_neg","dis_sero_vacc_pos","dis_sero_vacc","dis_sero_pop","dis_sero_vacc_pri", "dis_sero_vacc_secp", "dis_sero_vacc_secp_minus_pos",
                  "sdis_sero_vacc_neg","sdis_sero_vacc_pos","sdis_sero_vacc","sdis_sero_pop","sdis_sero_vacc_pri", "sdis_sero_vacc_secp", "sdis_sero_vacc_secp_minus_pos")


mat <- matrix(ncol = 0, nrow = 0)
df.out.year <- data.frame(mat)
df.out.cum <- data.frame(mat)
df.out.year.sero <- data.frame(mat)
df.out.cum.sero <- data.frame(mat)

for(h in 1:length(countries)){
    for(j in 1:length(vacc.age)){
      for(k in 1:length(vacc.cov)){
        for(i in 1:length(R0)){
          for(v in 1:length(v_list)){
            
          scenario <- paste0(countries[h], "_", vacc.cov[k], "_", vacc.age[j], "_", R0[i])
          
          # read in results - all serotypes
          #res <- readRDS(paste0(pathout, folder_type, new_folder, "prop_averted_", v_list[v], "_", scenario,".rds"))
          res <- out.list[[paste0("prop_averted_", v_list[v], "_", scenario)]]
          res.year <- res[[1]]
          res.cum <- res[[2]]
          
          # read in results - by  serotypes
          #res.sero <- readRDS(paste0(pathout, folder_type, new_folder, "prop_averted_", v_list_sero[v], "_", scenario,".rds"))
          res.sero <- out.list[[paste0("prop_averted_", v_list_sero[v], "_", scenario)]]
          res.sero.year <- res.sero[[1]]
          res.sero.cum <- res.sero[[2]]
          
          if(v <= length(v_list)/2){
            type_dis <- "symptomatic disease"
          }else{
            type_dis <- "hospitalisation"
          }
          
          if(v == 4 || v == 11){
            pop <- "pop"
          }else{
            pop <- "vacc"
          }
          
          if(v == 1 || v == 8){
            base <- "neg"
          }else if(v == 2 || v == 9){
            base <- "pos"
          }else if(v == 5 || v == 12){
            base <- "pri"
          }else if(v == 6 || v == 13){
            base <- "secp"
          }else if(v == 7 || v == 14){
            base <- "secp_minus_pos"
          }else{
            base <- "all"
          }
          
          # output yearly only for variables that are not secp_minus_pos
          if(v %in% c(1,2,3,4,5,6,8,9,10,11,12,13)){

            df.year <- data.frame(year = seq(1,20,1), 
                                  mean = res.year[,1],
                                  low = res.year[,2], 
                                  median = res.year[,3], 
                                  upper = res.year[,4], 
                                  serop9 = serop9[i],
                                  coverage = vacc.cov[k], 
                                  age = vacc.age[j],
                                  type = type_dis,
                                  pop = pop, 
                                  base = base, 
                                  country = countries[h])
            
            df.year.sero <- data.frame(year = rep(seq(1,20,1), 4), 
                                       mean = c(res.sero.year[,c(1, 5, 9, 13)]),
                                       low = c(res.sero.year[,c(2, 6, 10, 14)]), 
                                       median = c(res.sero.year[,c(3, 7, 11, 15)]), 
                                       upper = c(res.sero.year[,c(4, 8, 12, 16)]), 
                                       serotype = rep(c("D1", "D2", "D3", "D4"), each = 20), 
                                       serop9 = serop9[i],
                                       coverage = vacc.cov[k], 
                                       age = vacc.age[j],
                                       type = type_dis,
                                       pop = pop, 
                                       base = base, 
                                       country = countries[h])
          }

        
          df.cum <- data.frame(year = c(5,10,20), 
                               mean = res.cum[,1],
                               low = res.cum[,2], 
                               median = res.cum[,3], 
                               upper = res.cum[,4], 
                               serop9 = serop9[i],
                               coverage = vacc.cov[k], 
                               age = vacc.age[j],
                               type = type_dis,
                               pop = pop, 
                               base = base, 
                               country = countries[h])
          
          df.cum.sero <- data.frame(year = rep(c(5,10,20), 4), 
                                    mean = c(res.sero.cum[,c(1, 5, 9, 13)]),
                                    low = c(res.sero.cum[,c(2, 6, 10, 14)]), 
                                    median = c(res.sero.cum[,c(3, 7, 11, 15)]), 
                                    upper = c(res.sero.cum[,c(4, 8, 12, 16)]), 
                                    serotype = rep(c("D1", "D2", "D3", "D4"), each = 3), 
                                    serop9 = serop9[i],
                                    coverage = vacc.cov[k], 
                                    age = vacc.age[j],
                                    type = type_dis,
                                    pop = pop, 
                                    base = base, 
                                    country = countries[h])
          
          df.out.year <- rbind(df.out.year, df.year)
          df.out.cum <- rbind(df.out.cum, df.cum)
          
          df.out.year.sero <- rbind(df.out.year.sero, df.year.sero)
          df.out.cum.sero <- rbind(df.out.cum.sero, df.cum.sero)
          
        }
      }
    }
  }
}

# output csv files with impact estimates 
write.csv(df.out.year, paste0(pathout, folder_type, new_folder, "/csv/Summary_yearly_impacts_all_serotypes.csv"), row.names = FALSE, quote = FALSE)
write.csv(df.out.cum, paste0(pathout, folder_type, new_folder, "/csv/Summary_cumulative_impacts_all_serotypes.csv"), row.names = FALSE, quote = FALSE)
write.csv(df.out.year.sero, paste0(pathout, folder_type, new_folder, "/csv/Summary_yearly_impacts_by_serotype.csv"), row.names = FALSE, quote = FALSE)
write.csv(df.out.cum.sero, paste0(pathout, folder_type, new_folder, "/csv/Summary_cumulative_impacts_by_serotype.csv"), row.names = FALSE, quote = FALSE)

