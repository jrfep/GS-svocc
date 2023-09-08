require(dplyr)
require(raster)
#require(chron)
require(detect)
require(vegan)
require(stringr)
require(magrittr)
require(ggplot2)
require(forcats)
#require(tidyr)
here::i_am("doc/svocc-analysis-GS.Rmd")

GIS.data <- sprintf(here::here("Rdata","GIS.rda"))
load(GIS.data)

dts <- data.frame(bsq=values(vbsq),
                  dcon=values(dist.conucos),
                  dcom=values(dist.comunidades),
                  frs=values(dist.frs),
                  dbsq=values(dist.dbsq),
                  dcaz1=values(dist.caza1),
                  dcaz2=values(dist.caza2),
                  grid=values(rgrd)) %>% 
  group_by(grid) %>% 
  summarise(bsq=mean(bsq),
            dbsq=mean(dbsq),
            dcon=mean(dcon),
            dcom=mean(dcom),
            dpob=mean(min(dcon,dcom)),
            dhum=mean(min(dcon,dcom,dcaz1,dcaz2)),
            dcaz1=mean(dcaz1),
            dcaz2=mean(dcaz2),
            dcaz=mean(min(dcaz1,dcaz2)),
            frs=mean(frs))

mi.data <- dts %>% dplyr::select(bsq:frs) %>% decostand(.,"stand")

preds <- tibble()
for (k in 1:4) {
  input.dir <- here::here("Rdata","svocc",sprintf("best-%s",k))
  for (j in dir(input.dir,full.names=T)) {
    mdls <- (load(j))
    sp = basename(j) %>% str_replace("\\.rda","") %>%
      str_replace("\\.",". ")
    if ("fit.boot" %in% mdls) {
      # 
      prd0 <- predict(fit.boot,mi.data,type="response",se.fit=T)
      preds %<>% bind_rows(tibble(
        species=sp, ss=k, model="full",
        best=sum(prd0$fit),
        se=sqrt(sum(prd0$se.fit^2))))
    } else {
      prd0 <- predict(fit.null,
                      data.frame(bloque=rep(fit.null$levels$bloque,25)),
                      type="response",se.fit=T)
      preds %<>% bind_rows(tibble(
        species=sp, ss=k, model="null",
        best=sum(prd0$fit),
        se=sqrt(sum(prd0$se.fit^2))))   # propagate errors
    }
    
  }
}

preds %<>% mutate(
  method=if_else(ss %in% 1:2,"camera + off-camera","camera only"),
  sampling=if_else(ss %% 2 == 0,"Warapata","Warapata + Kavanayen"),
  ssplot = case_when(
    ss == 1 ~ "(b)",
    ss == 2 ~ "(a)",
    ss == 3 ~ "(d)",
    ss == 4 ~ "(c)"
  ))


dat <- {preds %>%
    mutate(species = fct_reorder(species, best,.fun=first)) %>% arrange(best)}
grps = tibble(xmin=c(0),xmax=c(50),ymin=c(1),ymax=c(13))
ggplot(data=dat) +
  geom_point(aes(x=best,y=species,colour=ss),cex=3) +
  geom_text(aes(label=ssplot), x=1,y=27) +
  geom_errorbar(aes(y=species,xmin=best-(2*se),xmax=best+(2*se),colour=ss), width = 0.2,lwd=1.05) + 
  facet_grid(method~sampling) +
  theme_minimal() + theme(legend.position="none", axis.text.y=element_text(face = "italic")) +
  annotate("rect", xmin = 0, xmax = 50, ymin = .5, ymax = 15.5,
           alpha = .2, lty=2, col="grey62") + 
  annotate("rect", xmin = 50, xmax = 125, ymin = 15.5, ymax = 23.5,
           alpha = .2, lty=2, col="grey62") + 
  annotate("rect", xmin = 125, xmax = 225, ymin = 23.5, ymax = 28.5,
           alpha = .2, lty=2, col="grey62") + 
  ylab("") + xlab("Predicted Nr. of cells")
ggsave(here::here("figs","Predicted-cells-data-selections.png"), width = 6, height = 7)

# Lin's concordance correlation coefficient 
# L.I.-K. Lin "A concordance correlation coefficient to evaluate reproducibility" 
# Biometrics, 45 (1989), pp. 255-268 https://doi.org/10.2307/2532051
# recommended by Watson and Petrie https://doi.org/10.1016/j.theriogenology.2010.01.003
blandAltmanPlot <- function(x,y) {
  ds <- x-y
  ms <- (x+y)/2
  mdifs <- mean(ds)
  sdifs <- sd(ds)
  plot(ds~ms)
  abline(h=0, col="olivedrab", lwd=3, lty=2)
  abline(h=mdifs - c(-2*sdifs,0,2*sdifs), col="maroon")
}
corLin <- function(x,y) {
  n <- length(x)
  r <- cor(x,y, method="pearson")
  sx2 <- var(x) * (n-1)/n 
  sy2 <- var(y) * (n-1)/n
  xm <- mean(x)
  ym <- mean(y)
  (2*r*sqrt(sx2)*sqrt(sy2)) / (sx2 + sy2 + (xm-ym)^2)
}
preds %>% distinct(ss,method,sampling)
x <- preds %>% filter(ss==1) %>% pull(best) # camera + observations Warapata + Kavanayen
y2 <- preds %>% filter(ss==2) %>% pull(best) # camera + observations Warapata 
y3 <- preds %>% filter(ss==3) %>% pull(best) # camera only           Warapata + Kavanayen
y4 <- preds %>% filter(ss==4) %>% pull(best) # camera only           Warapata 
blandAltmanPlot(x,y2)
corLin(x,y2)

blandAltmanPlot(x,y4)
corLin(x,y4)

## More information available if we use function from DescTools package
library(DescTools)
CCC(x,y2)$rho.c
CCC(x,y3)$rho.c
CCC(x,y4)$rho.c
