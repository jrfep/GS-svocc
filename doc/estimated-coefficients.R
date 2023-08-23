require(dplyr)
#require(raster)
#require(chron)
require(detect)
#require(vegan)
require(stringr)
require(magrittr)
require(ggplot2)
require(forcats)
require(tidyr)
library(DescTools)
library(doParallel)
here::i_am("doc/svocc-analysis-GS.Rmd")

mcoefs <- tibble()
for (k in 1:4) {
  input.dir <- here::here("Rdata","svocc",sprintf("best-%s",k))
  for (j in dir(input.dir,full.names=T)) {
    mdls <- (load(j))
    sp = basename(j) %>% str_replace("\\.rda","") %>%
      str_replace("\\.",". ")
    if ("fit.boot" %in% mdls) {
      cfs <- coefficients(fit.boot)
      cis <- confint(fit.boot)
      nms <- str_split_fixed(names(cfs), "_", n=2)
      mcoefs %<>% bind_rows(tibble(species=sp, ss=k,component=nms[,1],variable=nms[,2],mu=cfs,lower=cis[,1],upper=cis[,2]))
    } 
  }
}
mcoefs %<>% mutate(method=if_else(ss %in% 1:2,"camera + observations","camera only"),
                   sampling=if_else(ss %% 2 == 0,"Warapata","Warapata + Kavanayen"))

#We will use colours to highlight significant results:
cols <- c("sig. pos." = "slateblue4", "sig. neg." = "orange", "not significant" = "grey66")

#And we create this function to layout the figures for each variable:
  
plotModCoef <- function(varcode,varname) {
  dat <- {mcoefs %>% filter(component=="sta",variable==varcode) %>% mutate(
    species = fct_reorder(species, mu,.fun=first),
    response=case_when(
      lower>0 ~ "sig. pos.",
      upper<0 ~ "sig. neg.",
      TRUE ~ "not significant",
    )) %>% arrange(mu)}
  p <- ggplot(data=dat) +
    geom_point(aes(x=mu,y=species,colour=response),cex=3) +
    geom_errorbar(aes(y=species,xmin=lower,xmax=upper,colour=response), width = 0.2,lwd=1.05) + facet_grid(method~sampling) +
    scale_colour_manual(values = cols) +
    theme_minimal() + theme(legend.position="none", axis.text.y=element_text(face = "italic")) + geom_vline(xintercept=0,lty=2,col="black") + ylab("") + xlab("Estimated model coefficient") + labs(title=varname)
  return(p)
}

plotBestModCoef <- function(varcode,varname) {
  dat <- {mcoefs %>% filter(component=="sta",variable==varcode,ss==1) %>% mutate(
    species = fct_reorder(species, mu,.fun=first),
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

## Lin's CCC
foreach(v=c("bsq", "dbsq", "frs", "dcom", "dcon"), .combine=bind_rows) %do% {
  xys <- mcoefs %>% filter(variable %in% v) %>% 
    pivot_wider(id_cols="species", values_from=mu,names_from=ss)  
  xys %<>% filter(!if_all(`1`:`4`,~is.na(.)))
  x <- xys %>% pull(`1`) %>% coalesce(0L) # camera + observations Warapata + Kavanayen
  y2 <- xys %>% pull(`2`) %>% coalesce(0L) # camera + observations Warapata 
  y3 <- xys %>% pull(`3`) %>% coalesce(0L) # camera only           Warapata + Kavanayen
  y4 <- xys %>% pull(`4`) %>% coalesce(0L) # camera only           Warapata 
  
  Ltab <- bind_rows(CCC(x,y2)$rho.c, CCC(x,y3)$rho.c, CCC(x,y4)$rho.c) %>% mutate(vars=v, groups = 2:4)
} %>% 
  pivot_wider(id_cols = vars, values_from=c("est","lwr.ci","upr.ci"), names_from = groups, names_vary="slowest")


#We start with the covariate forest cover, which is include in all models:
  
plotModCoef("bsq","Forest")

plotModCoef("dbsq","Dist. to deforestation")

plotModCoef("frs","Dist. to fire")
plotModCoef("dcom","Dist. to communities")
plotModCoef("dcon","Dist. to conucos")




require(cowplot)
left_plot <- plot_grid(plotBestModCoef("bsq","Tree cover"),
                        plotBestModCoef("frs","Dist. to fire"),
                        labels = c("A","D"), ncol=1, rel_heights = c(2,1))
right_plot <- plot_grid(plotBestModCoef("dbsq","Dist. to deforestation"),
                        plotBestModCoef("dcon","Dist. to conucos"),
                        plotBestModCoef("dcom","Dist. to communities"),
           labels = c('B',"C","E"), ncol=1)
plot_grid(left_plot,right_plot)

ggsave(here::here("figs","Estimated-parameters-best-model.png"), width = 7, height = 8)
