# Resources selection model application for wildlife monitoring in Gran Sabana, Venezuela

Fit a Zero Inflated Binomial model to accumulated evidence of species detection from on camera and off camera records. 

## Repository structure

This repository includes scripts that can be run within ***R*** or using ***Rscript***. Use the Rproj file to open the project in RStudio or VScode. 

The [R](/R) folder includes three scripts and a markdown file that describes the workflow to run the scripts using ***Rscript***.

Data is downloaded and results saved under the [Rdata](/Rdata) and [output](/output) folders. These folders are not included in version control. Final results will be uploaded to an OSF repository.

Folders [doc](/doc) and [fig](/fig) include documentation of data and steps of the analysis:

- A [dashboard](/doc/dashboard.html) sumamrising sampling effort per region, blocks and cells
- A [document](/doc/svocc-analysis-GS.pdf) with the main results of the analysis for all species combined

## Reproducible research

We aim to provide all data and code necessary to reproduce our analysis and replicate our results. 

Documentation and comments in code have been written in English and/or Spanish. 

The data used in this analysis can be downloaded from:

- Ferrer-Paris, José R.; Stachowicz, Izabela (2022). Camera trap data and spatial covariates for monitoring vertebrates in the Gran Sabana, Venezuela. figshare. Dataset. https://doi.org/10.6084/m9.figshare.21174589.v1

This code repository is linked to the following component project:

- Ferrer-Paris, J. R., Stachowicz, I., & Sánchez-Mercado, A. (2023, August 20). Resources selection model application for wildlife monitoring in Gran Sabana, Venezuela. https://doi.org/10.17605/OSF.IO/6QTEP

Related repositories and resources are available under the following overarching project:

- Stachowicz, I., Ferrer-Paris, J. R., & Sánchez-Mercado, A. (2023, August 21). Monitoring vertebrates with camera traps in the Gran Sabana, Venezuela. https://doi.org/10.17605/OSF.IO/EY8FT


## References 

Our approach is based on:

> Solymos, P., Lele, S. R. 2016. Revisiting resource selection probability functions and single-visit methods: clarification and extensions. Methods in Ecology and Evolution, 7, 196–205. <doi:10.1111/2041-210X.12432>

