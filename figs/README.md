Previsualizar archivo `.dot` en RStudio:

```r
DiagrammeR::grViz("figs/GranSabana-analisis.dot")
```

Convertir archivo `.dot` a `.png` o `.svg` :

```sh
dot -Tpng GranSabana-analisis.dot > GranSabana-analisis.png
dot -Tsvg GranSabana-analisis.dot > GranSabana-analisis.svg
```