#!R --vanilla

require(dplyr)
require(tidyr)
require(raster)
require(chron)
require(detect)
require(vegan)
require(stringr)
require(magrittr)
require(ggplot2)
require(forcats)
require(tidyr)



GIS.data <- sprintf("Rdata/GIS.rda")
load(GIS.data)

input.dir <- sprintf("Rdata/svocc/explore")

allmodels <- tibble()

for (k in dir(input.dir,pattern="*rda",full.names=T)) {

      (load(k))
      spname <- basename(k) %>% str_replace("\\.rda","")
      allmodels %<>% bind_rows(params %>% mutate(sp=spname))
      rm(params)
}

## Nagelkerke R2 does not work for this example (they need to be NESTED models)
allmodels %>% filter(test1,test2==5,k==1,AIC<AICnull) %>% dplyr::select(sp,fml,test2,linkfuns,k,AIC,AICnull) %>% mutate(deltaAIC=AIC-AICnull) %>% print(n=100)

allmodels %>% filter(test1,test2==5,k==1,AIC<AICnull) %>% mutate(deltaAIC=AIC-AICnull) %>% group_by(sp,linkfuns) %>% summarise(n=n(),deltaAIC=min(deltaAIC)) %>% print(n=100)

# better number of models per species with cloglog
allmodels %>% filter(test1,test2==5,AIC<AICnull,linkfuns=="cloglog") %>% mutate(deltaAIC=AIC-AICnull) %>% group_by(sp,k) %>% summarise(n=n(),deltaAIC=min(deltaAIC)) %>% pivot_wider(id_cols="sp",names_from=c(k),values_from=n) %>% print(n=100)



## prepare data frame for prediction
dts <- data.frame(bsq=values(vbsq),dcon=values(dist.conucos),dcom=values(dist.comunidades),frs=values(dist.frs),dbsq=values(dist.dbsq),dcaz1=values(dist.caza1),dcaz2=values(dist.caza2),grid=values(rgrd)) %>% group_by(grid) %>% summarise(bsq=mean(bsq),dbsq=mean(dbsq),dcon=mean(dcon),dcom=mean(dcom),dpob=mean(min(dcon,dcom)),dhum=mean(min(dcon,dcom,dcaz1,dcaz2)),dcaz1=mean(dcaz1),dcaz2=mean(dcaz2),dcaz=mean(min(dcaz1,dcaz2)),frs=mean(frs))

dts %>% dplyr::select(bsq:frs) %>% decostand(.,"stand") -> mi.data

## fit best model for each method/data combination, or use null model for prediction

input.dir <- sprintf("Rdata/svocc/bootstrap")
mcoefs <- tibble()
preds <- tibble()
for (j in pull(allmodels,sp)) {
   input.file <- sprintf("%s/%s.rda",input.dir,j)
   if (file.exists(input.file)) {
      mdls <- (load(input.file))
      mdl <- get(grep("\\.boot",mdls,value=T))
      cfs <- coefficients(mdl)
      cis <- confint(mdl,type="boot") # very wide intervals
      cis <- confint(mdl)
      nms <- str_split_fixed(names(cfs), "_", n=2)
      mcoefs %<>% bind_rows(tibble(species=j, component=nms[,1],variable=nms[,2],mu=cfs,lower=cis[,1],upper=cis[,2]))
     
      # do not back transform, rather use error propagation
      # sqrt(sum(prd0$se.fit^2))
      prd0 <- predict(mdl,mi.data,type="response",se.fit=T)
      #A <- 1.85 * boot::logit(prd0$se.fit)
      #A[is.infinite(A)] <- -60
      #B <- boot::logit(prd0$fit)
      #ci.min <- boot::inv.logit(B+A)
      #ci.max <- boot::inv.logit(B-A)
      
      preds %<>% bind_rows(tibble(species=j, best=sum(prd0$fit),lower=sum(ci.min),upper=sum(ci.max)))
      
      rm(mdls,mdl,cfs,cis)
      
   }
}


mcoefs %>% filter(component=="sta",variable=="bsq") %>% print(n=100)

cols <- c("sig. pos." = "slateblue4", "sig. neg." = "orange", "not significant" = "grey66")

plotModCoef <- function(varcode,varname) {
  dat <- {mcoefs %>% filter(component=="sta",variable==varcode) %>% mutate(
    species = fct_reorder(str_replace(species,"\\.",". "), mu),
    response=case_when(
      lower>0 ~ "sig. pos.",
      upper<0 ~ "sig. neg.",
      TRUE ~ "not significant",
    )) %>% arrange(mu)}
  p <- ggplot(data=dat) +
    geom_point(aes(x=mu,y=species,colour=response),cex=3) +
    geom_errorbar(aes(y=species,xmin=lower,xmax=upper,colour=response), width = 0.2,lwd=1.05) +
    scale_colour_manual(values = cols) +
    theme_minimal() + theme(legend.position="none", axis.text.y=element_text(face = "italic")) + geom_vline(xintercept=0,lty=2,col="black") + ylab("") + xlab("Estimated model coefficient") + labs(title=varname)
  return(p)
}

plotModCoef("bsq","Forest")
plotModCoef("dbsq","Dist. to deforestation")
plotModCoef("frs","Dist. to fire")
plotModCoef("dcom","Dist. to communities")
plotModCoef("dcon","Dist. to conucos")

# + coord_cartesian(xlim=c(-6,6))


dat <- {preds %>% mutate(
  species = fct_reorder(str_replace(species,"\\.",". "), best)) %>% arrange(best)}
p <- ggplot(data=dat) +
  geom_point(aes(x=best,y=species),cex=3) +
  geom_errorbar(aes(y=species,xmin=lower,xmax=upper), width = 0.2,lwd=1.05) +
  theme_minimal() + theme(legend.position="none", axis.text.y=element_text(face = "italic")) + geom_vline(xintercept=0,lty=2,col="black") + ylab("") + xlab("Nr. of cells occupied") 

#   dev.copy(png,file=sprintf("%s-confint-spps.png",vv))
#   dev.off()
