process_demog <- function(fpath,fbirths,fpop,flife,fmortr,country.code,age_l,save.out=FALSE,fdesc="") {
  
    init.pop <- readRDS(paste0(fpath,"init_pop_",country.code,fdesc,".rds"))
    age.r <- readRDS(paste0(fpath,"age_rate_",country.code,fdesc,".rds"))
    mort.r <- readRDS(paste0(fpath,"mort_rate_",country.code,fdesc,".rds"))
    life.r <- readRDS(paste0(fpath,"life_ex_",country.code,fdesc,".rds"))
    births.r <- readRDS(paste0(fpath,"births_",country.code,fdesc,".rds"))
    

  return (list(births.r,init.pop,mort.r,age.r,life.r))
}
