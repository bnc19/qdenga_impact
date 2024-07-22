process_demog <- function(fpath,fbirths,fpop,flife,fmortr,country.code,age_l,save.out=FALSE,fdesc="") {
  
  if(file.exists(paste0(fpath,"init_pop_",country.code,fdesc,".rds")) &&
     file.exists(paste0(fpath,"age_rate_",country.code,fdesc,".rds")) &&
     file.exists(paste0(fpath,"mort_rate_",country.code,fdesc,".rds")) &&
     file.exists(paste0(fpath,"life_ex_",country.code,fdesc,".rds")) &&
     file.exists(paste0(fpath,"births_",country.code,fdesc,".rds"))) {
    
    init.pop <- readRDS(paste0(fpath,"init_pop_",country.code,fdesc,".rds"))
    age.r <- readRDS(paste0(fpath,"age_rate_",country.code,fdesc,".rds"))
    mort.r <- readRDS(paste0(fpath,"mort_rate_",country.code,fdesc,".rds"))
    life.r <- readRDS(paste0(fpath,"life_ex_",country.code,fdesc,".rds"))
    births.r <- readRDS(paste0(fpath,"births_",country.code,fdesc,".rds"))
    
  } else {
    
    if(!file.exists(paste0(fbirths,".RData"))) {
      births.lf <- read.csv(paste0(fbirths,".csv"))
      tmp <- read.csv("demog/births.csv")
      births.lf <- rbind(births.lf,tmp)
      births.lf <- births.lf[,c(2,6,8)]
      save(births.lf,file=paste0(fbirths,".RData"),compress=TRUE)
    } else {
      load(file=paste0(fbirths,".RData"))
    }
    births.wf <- births.lf[births.lf$country_code==country.code,c(2,3)]
    colnames(births.wf) <- c("year","births")
    births.wf <- births.wf[c(1,1:nrow(births.wf),nrow(births.wf)),]
    births.wf[1,1] <- births.wf[1,1]-1
    births.wf[nrow(births.wf),1] <- births.wf[nrow(births.wf),1]+1
    
    if(!file.exists(paste0(fpop,".RData"))) {
      pop.lf <- read.csv(paste0(fpop,".csv"))
      pop.lf <- pop.lf[pop.lf$year>=1950,]
      tmp <- read.csv("demog/pop.csv")
      pop.lf <- rbind(pop.lf,tmp)
      pop.lf <- pop.lf[,c(2,6,4,8)]
      save(pop.lf,file=paste0(fpop,".RData"),compress=TRUE)
    } else {
      load(file=paste0(fpop,".RData"))
    }
    pop.lf <- pop.lf[pop.lf$country_code==country.code,c(2,3,4)]
    pop.wf <- aggregate(value ~ year, pop.lf, I,simplify=TRUE)
    pop.wf <- as.data.frame(cbind(pop.wf$year,pop.wf$value))
    colnames(pop.wf) <- c("year",0:100)
    pop.wf <- pop.wf[c(1,1:nrow(pop.wf),nrow(pop.wf)),]
    pop.wf[1,1] <- pop.wf[1,1]-1
    pop.wf[nrow(pop.wf),1] <- pop.wf[nrow(pop.wf),1]+1  
    
    if(!file.exists(paste0(flife,".RData"))) {
      life.lf <- read.csv(paste0(flife,".csv"))
      tmp <- read.csv("demog/life_ex.csv")
      life.lf <- rbind(life.lf,tmp)
      life.lf <- life.lf[,c(2,6,4,8)]
      save(life.lf,file=paste0(flife,".RData"),compress=TRUE)
    } else {
      load(file=paste0(flife,".RData"))
    }
    life.lf <- life.lf[life.lf$country_code==country.code,c(2,3,4)]
    life.wf <- aggregate(value ~ year, life.lf, I,simplify=TRUE)
    life.wf <- as.data.frame(cbind(life.wf$year,life.wf$value))
    colnames(life.wf) <- c("year",0:100)
    life.wf <- life.wf[c(1,1:nrow(life.wf),nrow(life.wf)),]
    life.wf[1,1] <- life.wf[1,1]-1
    life.wf[nrow(life.wf),1] <- life.wf[nrow(life.wf),1]+1  
    
    if(!file.exists(paste0(fmortr,".RData"))) {
      mortr.lf <- read.csv(paste0(fmortr,".csv"))
      tmp <- read.csv("demog/mort_rate.csv")
      mortr.lf <- rbind(mortr.lf,tmp)
      mortr.lf <- mortr.lf[,c(2,6,4,8)]
      save(mortr.lf,file=paste0(fmortr,".RData"),compress=TRUE)
    } else {
      load(file=paste0(fmortr,".RData"))
    }
    mortr.lf <- mortr.lf[mortr.lf$country_code==country.code,c(2,3,4)]
    mortr.wf <- aggregate(value ~ year, mortr.lf, I,simplify=TRUE)
    mortr.wf <- as.data.frame(cbind(mortr.wf$year,mortr.wf$value))
    colnames(mortr.wf) <- c("year",0:100)
    mortr.wf <- mortr.wf[c(1,1:nrow(mortr.wf),nrow(mortr.wf)),]
    mortr.wf[1,1] <- mortr.wf[1,1]-1
    mortr.wf[nrow(mortr.wf),1] <- mortr.wf[nrow(mortr.wf),1]+1
    
    pop.mort.wf <- pop.wf * mortr.wf
    pop.mort.life.wf <- pop.mort.wf * life.wf
    
    n.r <- nrow(pop.wf)
    n.a <- length(age_l)
    age_u <-c(age_l[-1],101)
    age_c <- age_u-age_l
    
    pop.r <- pop.wf[,1:(n.a+1)]
    colnames(pop.r) <- c("year", age_l)
    age.r <- pop.r
    life.r <- pop.r
    for(i in 1:n.a) {
      pop.r[,i+1] <- rowSums(cbind(pop.wf[,(2+age_l[i]):(1+age_u[i])]),rep(0,n.r))
      age.r[,i+1] <- pop.wf[,1+age_u[i]]/pop.r[,i+1]
      life.r[,i+1] <- rowSums(cbind(pop.mort.life.wf[,(2+age_l[i]):(1+age_u[i])]),rep(0,n.r))/rowSums(cbind(pop.mort.wf[,(2+age_l[i]):(1+age_u[i])]),rep(0,n.r))
    }
    age.up <- pop.r
    age.up[,3:(1+n.a)] <- pop.r[,2:n.a]*age.r[,2:n.a]
    age.up[,2] <- c(births.wf[c(2:n.r,n.r),2])
    age.stay <- age.up
    age.stay[,2:(1+n.a)] <- pop.r[,2:(1+n.a)]*(1-age.r[,2:(1+n.a)])
    age.up <- age.up[c(1,1:(n.r-1)),]
    age.stay <- age.stay[c(1,1:(n.r-1)),]
    mort.r <- pop.r
    mort.r[,2:(1+n.a)] <- -log(pop.r[,2:(1+n.a)]/(age.stay[,2:(1+n.a)]+age.up[,2:(1+n.a)]))
  #  mort.r[,2] <- mort.r[,2]*2
    zero.r <- pop.r
    zero.r[,2:(1+n.a)] <- 0
    mort.r <- pmax(mort.r,zero.r)
    init.pop <- pop.r[2,]
    age.r <- age.r[2:(n.r-1),]
    mort.r <- mort.r[2:(n.r-1),]
    life.r <- life.r[2:(n.r-1),]
    births.r <- births.wf[2:(n.r-1),]
    if(save.out) {
      saveRDS(init.pop,paste0(fpath,"init_pop_",country.code,fdesc,".rds"))
      saveRDS(age.r,paste0(fpath,"age_rate_",country.code,fdesc,".rds"))
      saveRDS(mort.r,paste0(fpath,"mort_rate_",country.code,fdesc,".rds"))
      saveRDS(life.r,paste0(fpath,"life_ex_",country.code,fdesc,".rds"))
      saveRDS(births.r,paste0(fpath,"births_",country.code,fdesc,".rds"))
    }
  }
  return (list(births.r,init.pop,mort.r,age.r,life.r))
}

process_demog_gavi <- function(fpath,fbirths,fpop,flife,fmortr,country.code,age_l,save.out=FALSE,fdesc="") {
  
  if(file.exists(paste0(fpath,"/",country.code,"/init_pop",fdesc,".rds")) &&
     file.exists(paste0(fpath,"/",country.code,"/age_rate",fdesc,".rds")) &&
     file.exists(paste0(fpath,"/",country.code,"/mort_rate",fdesc,".rds")) &&
     file.exists(paste0(fpath,"/",country.code,"/life_ex",fdesc,".rds")) &&
     file.exists(paste0(fpath,"/",country.code,"/births",fdesc,".rds"))) {
    
    init.pop <- readRDS(paste0(fpath,"/",country.code,"/init_pop",fdesc,".rds"))
    age.r <- readRDS(paste0(fpath,"/",country.code,"/age_rate",fdesc,".rds"))
    mort.r <- readRDS(paste0(fpath,"/",country.code,"/mort_rate",fdesc,".rds"))
    life.r <- readRDS(paste0(fpath,"/",country.code,"/life_ex",fdesc,".rds"))
    births.r <- readRDS(paste0(fpath,"/",country.code,"/births",fdesc,".rds"))
    
  } else {
    
    if(!file.exists(paste0(fbirths,".RData"))) {
      births.lf <- read.csv(paste0(fbirths,".csv"))
      births.lf <- births.lf[,c(2,6,8)]
      save(births.lf,file=paste0(fbirths,".RData"),compress=TRUE)
    } else {
      load(file=paste0(fbirths,".RData"))
    }
    births.wf <- births.lf[births.lf$country_code==country.code,c(2,3)]
    colnames(births.wf) <- c("year","births")
    births.wf <- births.wf[c(1,1:nrow(births.wf),nrow(births.wf)),]
    births.wf[1,1] <- births.wf[1,1]-1
    births.wf[nrow(births.wf),1] <- births.wf[nrow(births.wf),1]+1
    
    if(!file.exists(paste0(fpop,".RData"))) {
      pop.lf <- read.csv(paste0(fpop,".csv"))
      pop.lf <- pop.lf[pop.lf$year>=1950,]
      pop.lf <- pop.lf[,c(2,6,4,8)]
      save(pop.lf,file=paste0(fpop,".RData"),compress=TRUE)
    } else {
      load(file=paste0(fpop,".RData"))
    }
    pop.lf <- pop.lf[pop.lf$country_code==country.code,c(2,3,4)]
    pop.wf <- aggregate(value ~ year, pop.lf, I,simplify=TRUE)
    pop.wf <- as.data.frame(cbind(pop.wf$year,pop.wf$value))
    colnames(pop.wf) <- c("year",0:100)
    pop.wf <- pop.wf[c(1,1:nrow(pop.wf),nrow(pop.wf)),]
    pop.wf[1,1] <- pop.wf[1,1]-1
    pop.wf[nrow(pop.wf),1] <- pop.wf[nrow(pop.wf),1]+1  
    
    if(!file.exists(paste0(flife,".RData"))) {
      life.lf <- read.csv(paste0(flife,".csv"))
      life.lf <- life.lf[,c(2,6,4,8)]
      save(life.lf,file=paste0(flife,".RData"),compress=TRUE)
    } else {
      load(file=paste0(flife,".RData"))
    }
    life.lf <- life.lf[life.lf$country_code==country.code,c(2,3,4)]
    life.wf <- aggregate(value ~ year, life.lf, I,simplify=TRUE)
    life.wf <- as.data.frame(cbind(life.wf$year,life.wf$value))
    colnames(life.wf) <- c("year",0:100)
    life.wf <- life.wf[c(1,1:nrow(life.wf),nrow(life.wf)),]
    life.wf[1,1] <- life.wf[1,1]-1
    life.wf[nrow(life.wf),1] <- life.wf[nrow(life.wf),1]+1  
    
    if(!file.exists(paste0(fmortr,".RData"))) {
      mortr.lf <- read.csv(paste0(fmortr,".csv"))
      mortr.lf <- mortr.lf[,c(2,6,4,8)]
      save(mortr.lf,file=paste0(fmortr,".RData"),compress=TRUE)
    } else {
      load(file=paste0(fmortr,".RData"))
    }
    mortr.lf <- mortr.lf[mortr.lf$country_code==country.code,c(2,3,4)]
    mortr.wf <- aggregate(value ~ year, mortr.lf, I,simplify=TRUE)
    mortr.wf <- as.data.frame(cbind(mortr.wf$year,mortr.wf$value))
    colnames(mortr.wf) <- c("year",0:100)
    mortr.wf <- mortr.wf[c(1,1:nrow(mortr.wf),nrow(mortr.wf)),]
    mortr.wf[1,1] <- mortr.wf[1,1]-1
    mortr.wf[nrow(mortr.wf),1] <- mortr.wf[nrow(mortr.wf),1]+1
    
    pop.mort.wf <- pop.wf * mortr.wf
    pop.mort.life.wf <- pop.mort.wf * life.wf
    
    n.r <- nrow(pop.wf)
    n.a <- length(age_l)
    age_u <-c(age_l[-1],101)
    age_c <- age_u-age_l
    
    pop.r <- pop.wf[,1:(n.a+1)]
    colnames(pop.r) <- c("year", age_l)
    age.r <- pop.r
    life.r <- pop.r
    for(i in 1:n.a) {
      pop.r[,i+1] <- rowSums(cbind(pop.wf[,(2+age_l[i]):(1+age_u[i])]),rep(0,n.r))
      age.r[,i+1] <- pop.wf[,1+age_u[i]]/pop.r[,i+1]
      life.r[,i+1] <- rowSums(cbind(pop.mort.life.wf[,(2+age_l[i]):(1+age_u[i])]),rep(0,n.r))/rowSums(cbind(pop.mort.wf[,(2+age_l[i]):(1+age_u[i])]),rep(0,n.r))
    }
    age.up <- pop.r
    age.up[,3:(1+n.a)] <- pop.r[,2:n.a]*age.r[,2:n.a]
    age.up[,2] <- c(births.wf[c(2:n.r,n.r),2])
    age.stay <- age.up
    age.stay[,2:(1+n.a)] <- pop.r[,2:(1+n.a)]*(1-age.r[,2:(1+n.a)])
    age.up <- age.up[c(1,1:(n.r-1)),]
    age.stay <- age.stay[c(1,1:(n.r-1)),]
    mort.r <- pop.r
    mort.r[,2:(1+n.a)] <- -log(pop.r[,2:(1+n.a)]/(age.stay[,2:(1+n.a)]+age.up[,2:(1+n.a)]))
    #  mort.r[,2] <- mort.r[,2]*2
    zero.r <- pop.r
    zero.r[,2:(1+n.a)] <- 0
    mort.r <- pmax(mort.r,zero.r)
    init.pop <- pop.r[2,]
    age.r <- age.r[2:(n.r-1),]
    mort.r <- mort.r[2:(n.r-1),]
    life.r <- life.r[2:(n.r-1),]
    births.r <- births.wf[2:(n.r-1),]
    if(save.out) {
      if(!dir.exists(paste0(fpath,"/",country.code))) dir.create(paste0(fpath,"/",country.code))
      saveRDS(init.pop,paste0(fpath,"/",country.code,"/init_pop",fdesc,".rds"))
      saveRDS(age.r,paste0(fpath,"/",country.code,"/age_rate",fdesc,".rds"))
      saveRDS(mort.r,paste0(fpath,"/",country.code,"/mort_rate",fdesc,".rds"))
      saveRDS(life.r,paste0(fpath,"/",country.code,"/life_ex",fdesc,".rds"))
      saveRDS(births.r,paste0(fpath,"/",country.code,"/births",fdesc,".rds"))
    }
  }
  return (list(births.r,init.pop,mort.r,age.r,life.r))
}
