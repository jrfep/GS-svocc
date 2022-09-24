#! R --vanilla

args <- commandArgs(TRUE)
sp <- args[1]

input.file <- sprintf("Rdata/svocc/explore/%s.rda",sp)

require(dplyr)
require(tidyr)
require(raster)
require(chron)
require(detect)
GIS.data <- "Rdata/GIS.rda"
load(GIS.data)
mi.rda <- sprintf("Rdata/padata/%s.rda",sp)
if (file.exists(mi.rda))
   load(file=mi.rda)
if (file.exists(input.file))
   load(file=input.file)

for (kk in 1:4) {
  out.file <- sprintf("Rdata/svocc/best-%s/%s.rda",kk,sp)
  if (!file.exists(out.file)) {
    ss <- switch(kk,`1`=rep(T,nrow(pa.data)), `2`=pa.data$muestreo, `3`=pa.data$metodo, `4`=pa.data$metodo & pa.data$muestreo)
    lkf <- "cloglog"

    slc_params <-  params %>% filter(test1,test2==5,AIC<AICnull,k==kk,linkfuns==lkf) %>% arrange(AIC-AICnull)
    if (nrow(slc_params)>0) {
       fml <- slc_params %>% slice(1) %>% pull(fml)
       nfml <- slc_params %>% slice(1) %>% pull(nfml)

       fit.full <- svocc(formula(fml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
       fit.null <- svocc(formula(nfml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
       save(file=out.file,fit.null,fit.full)
  } else {
    slc_params <-  nulls %>% filter(k==kk,linkfuns==lkf) %>% arrange(AICnull) %>% slice(1)
    nfml <- slc_params %>% pull(nfml)
     fit.null <- svocc(formula(nfml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
     save(file=out.file,fit.null)
  }
}
