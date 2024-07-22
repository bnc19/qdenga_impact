
########### Set up to run on the high performance computing cluster ############
# library(hipercow)
# hipercow_configure(driver = "windows")
# packages = c("dust", "odin", "odin.dust", "denguetak")
# my_sources = c("R/run_who_cluster.r","R/process_demog.r","R/config_params.r")
# hipercow_provision()
# hipercow_environment_create(sources = my_sources,
#                             packages = packages)
################################################################################ 

# install packages

install.packages(
  "odin.dust",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))

install.packages(
  "odin",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))

# install.packages("drat") # -- if you don't have drat installed
drat:::add("ncov-ic")
install.packages("dust")
conan::conan_sources(c("mrc-ide/dust", "mrc-ide/odin", "mrc-ide/odin.dust", "mrc-ide/denguetak"))

install.packages("dust")
install.packages("denguetak")
install.packages("odin.dust")

library(odin.dust)
# setwd(paste0(getwd(), "/impact"))

# source files  
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)


ns <- 200
np <- 50
nt <- 50
pf <- "ps_new24.csv"
outf <- "//wpia-hn-app/denguetak/outputs24"
outf2 <- "//wpia-hn.hpc.dide.ic.ac.uk/denguetak/outputs24"
outf <- outf2

### equilib runs

job <- list()
k<-0

for(i in 1:9) {
  k <- k+1
  job[[k]] <- hpcq$enqueue(run_WHO_eq(1,i,ns,np,pf,nt,outf))
  k <- k+1
  job[[k]] <- hpcq$enqueue(run_WHO_eq(2,i,ns,np,pf,nt,outf))
}


## vacc runs

job <- list()
k <- 0
vbs <- 0

for(vc in c(0.8)) {   
  for(va in c(4:18)) {
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- hpcq$enqueue(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(0,1)) {
        for(mdy in c(5,15)) {
          k <- k+1
          job[[k]] <- hpcq$enqueue(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- hpcq$enqueue(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}


job <- list()
k <- 0
vbs <- 0

for(vc in c(0.2,0.4,0.6)) {   
  for(va in c(6)) {
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- hpcq$enqueue(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(1)) {
        for(mdy in c(5,15)) {
          k <- k+1
          job[[k]] <- hpcq$enqueue(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- hpcq$enqueue(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}

## with sero-screening

outf <- "//wpia-hn-app/denguetak/outputs24vbs"
outf2 <- "//wpia-hn.hpc.dide.ic.ac.uk/denguetak/outputs24vbs"
#outf <- outf2

job <- list()
k <- 0
vbs <- 1

for(vc in c(0.8)) {   
  for(va in c(4:18)) {
    ## no vacc runs
    for(i in 1:9) {
      k <- k+1
      job[[k]] <- hpcq$enqueue(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
      k <- k+1
      job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
      ## vacc runs
      for(vei in c(0,1)) {
        for(mdy in c(5,15)) {
          k <- k+1
          job[[k]] <- hpcq$enqueue(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
          k <- k+1
          job[[k]] <- hpcq$enqueue(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
        }
      }
    }
  }
}


k <- k+1
job[[k]] <- hpcq$enqueue(run_WHO_vacc(1,1,ns,np,pf,nt,15,0.8,1,vbs,5,outf))
k <- k+1
job[[k]] <- hpcq$enqueue(run_WHO_vacc(2,9,ns,np,pf,nt,14,0.8,1,vbs,5,outf))



job[[6]]$status()
job[[8]]$result()
rm(hpcq)



k <- k+1
job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,4,ns,np,pf,nt,9,0.8,outf))
k <- k+1
job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,5,ns,np,pf,nt,9,0.8,outf))
k <- k+1
job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,7,ns,np,pf,nt,9,0.8,outf))
k <- k+1
job[[k]] <- hpcq$enqueue(run_WHO_novacc(2,8,ns,np,pf,nt,9,0.8,outf))



###############

source("R/config_params.r")
source("R/process_demog.r")
ns <- 200
np <- 50
nt <- 50
pf <- "ps_md57.csv"
outf <- "//wpia-hn-app/denguetak/outputs3"
outf2 <- "//wpia-hn.hpc.dide.ic.ac.uk/denguetak/outputs3"
cn <-1
pars_nv <- config_params(pf,200,0.8,1,9,0,0,10)
gen_dengue <- denguetak::model
mod <- gen_dengue$new(pars = pars_nv, 
                      time = 0,
                      n_particles = 200,
                      n_threads = 50,
                      seed=1L,
                      pars_multi = TRUE,
                      deterministic = FALSE)
mod.var.list <- mod$info()[[1]]
saveRDS(mod.var.list,"U:/outputs3/var_list.rds")

has context menu