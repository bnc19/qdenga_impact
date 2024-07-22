library(hipercow)
hipercow_configure(driver = "windows")

packages = c("dust", "odin", "odin.dust", "denguetak")
my_sources = c("R/run_who_cluster.r","R/process_demog.r","R/config_params.r")

hipercow_provision()

hipercow_environment_create(sources = my_sources,
                            packages = packages)
                          
ns <- 200
np <- 50
nt <- 50
pf <- "post.samp.csv"
outf <- "Z:/outputsBCD/impact_BCD/outputs"
outf2 <- "Z:/outputsBCD/impact_BCD/outputs"
outf <- outf2

## vacc runs

job <- list()
k <- 0
vbs <- 0 

# age for VS D15 0.8
for(vc in c(0.8)) {   
  for(va in c(6:12)) { # make 6 - 12 if want age differences 
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(0)) {
        for(mdy in c(15)) {
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}


# VS D5 

job <- list()
k <- 0
vbs <- 0 

for(vc in c(0.8)) {   
  for(va in c(6)) {
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(0)) {
        for(mdy in c(5)) {
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}


#  VI coverage for age 6, decay 15
job <- list()
k <- 0
vbs <- 0 

for(vc in c(0.2,0.4, 0.6, 0.8)) {   
  for(va in c(6)) { 
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(1)) {
        for(mdy in c(15)) {
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}

# SA on VI 3

outf <- "Z:/outputsBCD/impact_BCD/outputs_3VI"
outf2 <- "Z:/outputsBCD/impact_BCD/outputs_3VI"


job <- list()
k <- 0
vbs <- 0 

for(vc in c(0.8)) {   
  for(va in c(6)) {  
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(1)) {
        for(mdy in c(15)) {
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}




# SA on VI 12

outf <- "Z:/outputsBCD/impact_BCD/outputs_12VI"
outf2 <- "Z:/outputsBCD/impact_BCD/outputs_12VI"


job <- list()
k <- 0
vbs <- 0 

for(vc in c(0.8)) {   
  for(va in c(6)) {  
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(1)) {
        for(mdy in c(15)) {
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}


## RUN SCREENING 
outf <- "Z:/outputsBCD/impact_BCD/outputs_vbs"
outf2 <- "Z:/outputsBCD/impact_BCD/outputs_vbs"



job <- list()
k <- 0
vbs <- 1

# age for VS D15 0.8
for(vc in c(0.8)) {   
  for(va in c(6)) { # make 6 - 12 if want age differences 
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(0)) {
        for(mdy in c(15)) {
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}
