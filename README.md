# Integrating histopathology and transcriptomics for spatial profiling of the tumor microenvironment: a melanoma case study

This repository holds the code for all analyses related to

> t.b.a

Our pipeline (SPoTLIghT) to derive spatial graph-based interpretable features from H&E tissue slides is packaged as a docker container, which is available here: 
https://github.com/SysBioOncology/spotlight_docker.

## Spatial features for used cohorts
Spatial features for the TCGA and CPTAC melanoma cohorts are provided here:

- `data/spatial_features_matrix_TCGA.csv`: TCGA-SKCM spatial features
- `data/spatial_features_matrix_CPTAC.csv`: CPTAC-CM spatial features 

## Notebook to reproduce the figures from the manuscript

### Required .xlsx files

Download the following files and place them in the `data` directory:

- `bagaev_file`: bagaev microenvironment subtypes for TCGA patients (available from [Bagaev et al. 2021, Table S6,  sheet 1](https://ars.els-cdn.com/content/image/1-s2.0-S1535610821002221-mmc6.xlsx))
- `tcga_survival_file`: survival data for TCGA patients (available from [Liu et al. 2018, Table S1, sheet 1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx)).

```bash

# First, activate conda environment

conda env create -f spatial_features.yml
conda activate spatial_features

# Second, render the notebook

Rscript -e "rmarkdown::render('spatial_features.Rmd, params=list(tcga_bagaev_subtypes = bagaev_file, tcga_survival = tcga_survival_file)')"

```

## Notebook to reproduce figures for validation with CPTAC using logistic modeling


```bash

# First, activate conda environment

conda env create -f spatial_features_manuscript.yml
conda activate spatial_features_manuscript

# Second, render the notebook

Rscript -e "rmarkdown::render('CPTAC_validation_plots.rmd')"

```