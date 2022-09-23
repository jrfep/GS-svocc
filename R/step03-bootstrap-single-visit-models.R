#! R --vanilla
require(dplyr)
require(tidyr)
require(raster)
require(chron)
require(detect)

GIS.data <- "Rdata/GIS.rda"
load(GIS.data)

input.dir <- sprintf("Rdata/svocc/step/",sp)

for (sprda in dir(input.dir)) {
   print(sprda)
   mi.rda <- sprintf("Rdata/padata/%s",sprda)
   input.file <- sprintf("%s/%s",input.dir,sprda)
   out.file <- sprintf("Rdata/svocc/bootstrap/%s",sprda)

   if (file.exists(mi.rda))
      load(file=mi.rda)
   if (file.exists(input.file))
      load(file=input.file)

   if (!file.exists(out.file)) {
      if (nrow(sparams)>0) {

        if ( with(sparams,any(AIC<AICnull))) {
           slc_params <-  sparams %>% filter(test1,test2==5,AIC<AICnull,k==min(k)) %>% arrange(AIC-AICnull)
        } else if (with(sparams,any(AICs<AICnull))){
           slc_params <-  sparams %>% filter(test1,test2==5,AICs<AICnull,k==min(k)) %>% arrange(AIC-AICnull)
        } else {
           slc_params <-  sparams %>% filter(test1,test2==3,AIC<AICnull | AICs<AICnull,k==min(k)) %>% arrange(AIC-AICnull)
        }
      if (nrow(slc_params)>0) {
         j<-1
         fml <- slc_params %>% slice(j) %>% pull(fml)
         nfml <- slc_params %>% slice(j) %>% pull(nfml)
         lkf <- slc_params %>% slice(j) %>% pull(linkfuns)
         ss <- slc_params %>% slice(j) %>% pull(k) %>% switch(`1`=rep(T,nrow(pa.data)), `2`=pa.data$muestreo, `3`=pa.data$metodo, `4`=pa.data$metodo & pa.data$muestreo)
         fit <- svocc(formula(fml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
         fit.step <- svocc.step(fit, model="sta")
         fit.null <- svocc(formula(nfml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
         assign(sprintf("%s.full",sub(".rda","",sprda)),fit)
         bfit <- bootstrap(fit.null, B=50)
         assign(sprintf("%s.null",sp),bfit)
         bfit <- bootstrap(fit, B=50)
         assign(sprintf("%s.boot",sub(".rda","",sprda)),bfit)
         bfit <- bootstrap(fit.step, B=50)
         assign(sprintf("%s.step",sub(".rda","",sprda)),bfit)

         save(file=out.file,list=ls(pattern=".boot|.full|.null|.step"))
         rm(list=ls(pattern=".boot|.full|.null|.step"))
      }
   }
   }
}
