# Running all the scripts

First prepare all folder:

```sh
mkdir -p Rdata/padata/
mkdir -p Rdata/svocc/explore/

for j in $(seq 1 4)
do
  mkdir -p Rdata/svocc/best-$j
done

```

There are three R-scripts to be run on each species, for example, this will run all scripts in sequential order for one species:

```sh
SPECIE=P.maximus
Rscript --vanilla R/step00-presence-absence-dataframe.R $SPECIE
Rscript --vanilla R/step01-explore-single-visit-models.R $SPECIE
Rscript --vanilla R/step02-best-single-visit-models.R $SPECIE
```

The scripts will skip the species if the output files already exist, thus it is necessary to delete any old files in order to update the results.

If we want to run all the scripts for all the species, we will first declare our list of species:

```sh
ALLSPECIES="C.alector C.olivaceus C.paca C.thous C.unicinctus D.imperfecta D.kappleri D.leporina D.marsupialis D.novemcinctus E.barbara H.hydrochaeris L.pardalis L.rufaxilla L.tigrinus L.wiedii M.americana M.gouazoubira M.pratti M.tridactyla N.nasua O.virginianus P.concolor P.maximus P.onca P.tajacu S.venaticus T.major T.pecari T.terrestris T.tetradactyla"

```

We can now run step 0 and step 1 for each species, the step 0 is necessary for step 1 so we have to wait for it to complete before running the script of step 1, but this one can be run in the background so that multiple species can be run in parallel.

```sh
for SPECIE in $ALLSPECIES
do
   Rscript --vanilla R/step00-presence-absence-dataframe.R $SPECIE
   nohup Rscript --vanilla R/step01-explore-single-visit-models.R $SPECIE > nohup-${SPECIE}.out &
done

```

We have to wait for step 1 to complete before we launch step 2, but this can also run in the background for each species:

```sh
for SPECIE in $ALLSPECIES
do
   nohup Rscript --vanilla R/step02-best-single-visit-models.R $SPECIE > nohup-${SPECIE}.out &
done

```

# Synchronizing copies 

We ran all the scripts in one computer (terra), but want to analyse the results in another one (roraima)?
No problem, let's synchronize our local copy with the remote:

```sh
rsync -gloptrunv terra.ad.unsw.edu.au:/home/jferrer/proyectos/IVIC/GS-svocc/Rdata/ Rdata

```
