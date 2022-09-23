# GS-svocc
Fit single visit occupancy models to data from camera trap and walking transects in the Gran Sabana, Venezuela

We can run these scripts in sequential order:

```sh
for SPECIE in C.alector C.olivaceus #C.paca C.thous C.unicinctus D.imperfecta D.kappleri D.leporina D.marsupialis D.novemcinctus E.barbara H.hydrochaeris L.pardalis L.rufaxilla L.tigrinus L.wiedii M.americana M.gouazoubira M.pratti M.tridactyla N.nasua O.virginianus P.concolor P.maximus P.onca P.tajacu S.venaticus T.major T.pecari T.terrestris T.tetradactyla
do
    Rscript --vanilla R/presence-absence-dataframe.R $SPECIE
    for k in $(seq 1 4)
    do
      nohup Rscript --vanilla R/single_visit_models.R $SPECIE $k &
    done
done
```
