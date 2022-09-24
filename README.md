# GS-svocc
Fit single visit occupancy models to data from camera trap and walking transects in the Gran Sabana, Venezuela



```sh
SPECIE=P.maximus
mkdir -p Rdata/padata/
mkdir -p Rdata/svocc/explore/
#mkdir -p Rdata/svocc/step/
#mkdir -p Rdata/svocc/bootstrap/
Rscript --vanilla R/step00-presence-absence-dataframe.R $SPECIE
Rscript --vanilla R/step01-explore-single-visit-models.R $SPECIE
#Rscript --vanilla R/step02-reduce-single-visit-models.R $SPECIE
```


We can run these scripts in sequential order:

```sh
mkdir -p Rdata/padata/
mkdir -p Rdata/svocc/explore/

for SPECIE in C.alector C.olivaceus C.paca C.thous C.unicinctus D.imperfecta D.kappleri D.leporina D.marsupialis D.novemcinctus E.barbara H.hydrochaeris L.pardalis L.rufaxilla L.tigrinus L.wiedii M.americana M.gouazoubira M.pratti M.tridactyla N.nasua O.virginianus P.concolor P.maximus P.onca P.tajacu S.venaticus T.major T.pecari T.terrestris T.tetradactyla
do
   Rscript --vanilla R/step00-presence-absence-dataframe.R $SPECIE
   nohup Rscript --vanilla R/step01-explore-single-visit-models.R $SPECIE > nohup-${SPECIE}.out &
#   nohup Rscript --vanilla R/step02-reduce-single-visit-models.R $SPECIE > nohup-${SPECIE}.out &
done

```

Synchronize local copy with remote

```sh
rsync -gloptrunv terra.ad.unsw.edu.au:/home/jferrer/proyectos/IVIC/GS-svocc/Rdata/ Rdata

```
