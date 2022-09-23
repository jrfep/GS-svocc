#! R --vanilla

args <- commandArgs(TRUE)
sp <- args[1]
k <- args[2]

out.file <- sprintf("Rdata/svocc-%s-%s.rda",sp,k)

if (!file.exists(out.file)) {
  require(dplyr)
  require(tidyr)
  require(raster)
  require(chron)
  require(detect)

  GIS.data <- "Rdata/GIS.rda"
  load(GIS.data)
  mi.rda <- sprintf("Rdata/padata-%s.rda",sp)
  if (file.exists(mi.rda))
     load(file=mi.rda)


  exito <- "NOT"

  linkfuns <- c("cloglog","probit")
  detvars <- c("muestreo + (walk * cam)","muestreo + walk + cam","walk * cam","walk + cam")
  occvars <- c("bsq + dcon + dcom + frs + dbsq",
  "bsq + dcon + frs + dbsq",
  "bsq + dcon + frs + dbsq",
  "bsq + dcom + frs + dbsq",
  "bsq + dcon + dcom + frs",
  "bsq + dcon + dbsq",
  "bsq + dcon + frs",
  "bsq + frs + dbsq",
  "bsq + dcon + dcom",
  "bsq + dcom + dbsq",
  "bsq + dcom + frs",
  "bsq + dbsq",
  "bsq + frs",
  "bsq + dcon",
  "bsq + dcom"
  )
  params <- expand_grid(occvars,detvars,linkfuns,k=1:4) %>%
  mutate(fml=sprintf("pa ~ %s |  %s",occvars,detvars),
    nfml=sprintf("pa ~ bloque |  %s",detvars),test1=F,test2=0L)


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

  for (j in 1:nrow(params)) {
    fml <- params %>% slice(j) %>% pull(fml)
    nfml <- params %>% slice(j) %>% pull(nfml)
    lkf <- params %>% slice(j) %>% pull(linkfuns)
    ss <- params %>% slice(j) %>% pull(k) %>% switch(`1`=rep(T,nrow(pa.data)), `2`=pa.data$muestreo, `3`=pa.data$metodo, `4`=pa.data$metodo & pa.data$muestreo)

    fit <- svocc(formula(fml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
    tst <- testfit(fit)

    if (sum(tst[1:3],na.rm=T)==3) {
      fit.null <- svocc(formula(nfml), data=pa.data[ss,],link.sta = lkf, link.det = "logit", penalized = FALSE, method = c( "optim"))
      params[j,"test1"]=T
      params[j,"test2"]=sum(tst[1:5],na.rm=T)
      params[j,"AIC"]=AIC(fit)
      params[j,"AICnull"]=AIC(fit.null)
      if (sum(tst[1:5],na.rm=T)==5 & AIC(fit)<AIC(fit.null)) stop("we made it!")

    }
  }

  params %>% filter(test1) %>% dplyr::select(fml,k,linkfuns,test2:AICnull)

sparams <- params %>% filter(test1) %>% mutate(test1=F,test2=0)

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
    if (sum(tst[1:5],na.rm=T)==5 & AIC(fit)<sparams[j,"AICnull"]) stop("we made it!")
  }
}

sparams %>% filter(test1) %>% dplyr::select(sfml,k,linkfuns,test2:AICs) %>% arrange(AICs)



  assign(sprintf("%s.full",sp),fit)
  bfit <- bootstrap(fit, B=50)
  assign(sprintf("%s.boot",sp),bfit)
  bfit <- bootstrap(fit.null, B=50)
  assign(sprintf("%s.null",sp),bfit)
  exito <- "OK"

  function(x) {
    }
    for (occvar in occvars) {
       for(mi.link in linkfuns) {
          for (detvar in detvars) {
             if (exito != "OK") {
                fit <- svocc(formula(sprintf("pa ~ %s |  %s",occvar,detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                if (all(abs(fit$coefficients$sta)<6) & all(!is.na(fit$std.error$sta)) & all(fit$std.error$sta<4) ) {
                   fit.null <- svocc(formula(sprintf("pa ~ bloque |  %s",detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                   if (AIC(fit)<AIC(fit.null)) {
                     assign(sprintf("%s.full",sp),fit)
                     bfit <- bootstrap(fit, B=50)
                     assign(sprintf("%s.boot",sp),bfit)
                     bfit <- bootstrap(fit.null, B=50)
                     assign(sprintf("%s.null",sp),bfit)
                     exito <- "OK"
                   }
                }
             }
          }
       }
    }


    for (occvar in occvars) {
       for(mi.link in linkfuns) {
          for (detvar in detvars) {
             if (exito != "OK") {
                fit <- svocc(formula(sprintf("pa ~ %s |  %s",occvar,detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                if (all(abs(fit$coefficients$sta)<10) & all(!is.na(fit$std.error$sta)) & all(fit$std.error$sta<6) ) {
                   fit.null <- svocc(formula(sprintf("pa ~ bloque |  %s",detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                   if (AIC(fit)<AIC(fit.null)) {
                     assign(sprintf("%s.alt",sp),fit)
                     bfit <- bootstrap(fit, B=50)
                     assign(sprintf("%s.boot",sp),bfit)
                     bfit <- bootstrap(fit.null, B=50)
                     assign(sprintf("%s.null",sp),bfit)
                     exito <- "OK"
                   }
                }
             }
          }
       }
    }

    for (occvar in occvars) {
       for(mi.link in linkfuns) {
          for (detvar in detvars) {
             if (exito != "OK") {
                prefit <- svocc(formula(sprintf("pa ~ %s |  %s",occvar,detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                fit <- svocc.step(prefit, model="sta")
                if (all(abs(fit$coefficients$sta)<6) & all(!is.na(fit$std.error$sta)) & all(fit$std.error$sta<4) ) {
                   fit.null <- svocc(formula(sprintf("pa ~ bloque |  %s",detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                   if (AIC(fit)<AIC(fit.null)) {
                     assign(sprintf("%s.step",sp),fit)
                     bfit <- bootstrap(fit, B=50)
                     assign(sprintf("%s.boot",sp),bfit)
                     bfit <- bootstrap(fit.null, B=50)
                     assign(sprintf("%s.null",sp),bfit)
                     exito <- "OK"
                   }
                }
             }
          }
       }
    }

    for (occvar in occvars) {
       for(mi.link in linkfuns) {
          for (detvar in detvars) {
             if (exito != "OK") {
                prefit <- svocc(formula(sprintf("pa ~ %s |  %s",occvar,detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                fit <- svocc.step(prefit, model="sta")
                if (all(abs(fit$coefficients$sta)<10) & all(!is.na(fit$std.error$sta)) & all(fit$std.error$sta<6) ) {
                   fit.null <- svocc(formula(sprintf("pa ~ bloque |  %s",detvar)), data=pa.data[ss,],link.sta = mi.link, link.det = "logit", penalized = FALSE, method = c( "optim"))
                   if (AIC(fit)<AIC(fit.null)) {
                     assign(sprintf("%s.step",sp),fit)
                     bfit <- bootstrap(fit, B=50)
                     assign(sprintf("%s.boot",sp),bfit)
                     bfit <- bootstrap(fit.null, B=50)
                     assign(sprintf("%s.null",sp),bfit)
                     exito <- "OK"
                   }
                }
             }
          }
       }
    }

    if (exito != "OK") {
       bfit <- bootstrap(fit.null, B=50)
       assign(sprintf("%s.null",sp),bfit)
    }
  save(file=out.file,list=ls(pattern=".boot|.full|.null|.alt|.step"))

}
