#! R --vanilla

args <- commandArgs(TRUE)
sp <- args[1]

input.file <- sprintf("Rdata/svocc/explore/%s.rda",sp)
out.file <- sprintf("Rdata/svocc/step/%s.rda",sp)

if (!file.exists(out.file)) {
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

     
   testfit <- function(x) {
      tsts <- c(
        all(abs(x$coefficients$sta)<10),
        all(!is.na(x$std.error$sta)),
        all(x$std.error$sta<6),
        all(abs(x$coefficients$sta)<6),
        all(x$std.error$sta<4)
      )
      return(tsts)
   }

slc_params <-  params %>% filter(test1,test2==5,AIC<AICnull,k==min(k)) %>% arrange(AIC-AICnull)
if (nrow(slc_params)>0) {
   j<-1
   fml <- slc_params %>% slice(j) %>% pull(fml)
   nfml <- slc_params %>% slice(j) %>% pull(nfml)
   lkf <- slc_params %>% slice(j) %>% pull(linkfuns)
   ss <- slc_params %>% slice(j) %>% pull(k) %>% switch(`1`=rep(T,nrow(pa.data)), `2`=pa.data$muestreo, `3`=pa.data$metodo, `4`=pa.data$metodo & pa.data$muestreo)
   sparams <- slc_params %>% slice(j)

   fit <- svocc(formula(fml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
   fit.step <- svocc.step(fit, model="sta")

   sparams[1,"AICs"]=AIC(fit.step)
   sparams[1,"sfml"]=paste(as.character(fit$formula$full)[c(2,1,3)],collapse=" ")

} else {
   sparams <-  params %>% filter(test1,test2==max(test2),k==min(k)) %>% arrange(AIC-AICnull) %>% mutate(test1=F,test2=0)
   if (nrow(sparams)>0) {
      for (j in 1:nrow(sparams)) {
        fml <- sparams %>% slice(j) %>% pull(fml)
        nfml <- sparams %>% slice(j) %>% pull(nfml)
        lkf <- sparams %>% slice(j) %>% pull(linkfuns)
        ss <- sparams %>% slice(j) %>% pull(k) %>% switch(`1`=rep(T,nrow(pa.data)), `2`=pa.data$muestreo, `3`=pa.data$metodo, `4`=pa.data$metodo & pa.data$muestreo)


        prefit <- svocc(formula(fml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
        fit <- svocc.step(prefit, model="sta")
        tst <- testfit(fit)

        if (sum(tst[1:3],na.rm=T)==3) {
          sparams[j,"test1"]=T
          sparams[j,"test2"]=sum(tst[1:5],na.rm=T)
          sparams[j,"AICs"]=AIC(fit)
          sparams[j,"sfml"]=paste(as.character(fit$formula$full)[c(2,1,3)],collapse=" ")
        }
      }
   }
}



   #sparams %>% filter(test1) %>% dplyr::select(sfml,k,linkfuns,test2:AICs) %>% arrange(AICs)

 sparams <- {sparams %>% filter(test1)}
 save(file=out.file,sparams)
}
