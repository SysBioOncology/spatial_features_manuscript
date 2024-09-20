# Integrating histopathology and transcriptomics for spatial profiling of the tumor microenvironment: a melanoma case study

This repository holds the code for all analyses related to

> t.b.a

Our pipeline (SPoTLIghT) to derive spatial graph-based interpretable features from H&E tissue slides is packaged as a docker container, which is available here:
https://github.com/SysBioOncology/spotlight_docker.

## Spatial features for used cohorts

Spatial features for the TCGA and CPTAC melanoma cohorts are provided here:

* `data/spatial_features_matrix_TCGA.csv`: TCGA-SKCM spatial features
* `data/spatial_features_matrix_CPTAC.csv`: CPTAC-CM spatial features

## Notebook to reproduce the figures from the manuscript

### Required .xlsx files

Download the following files and place them in the `data` directory:

* `bagaev_file`: bagaev microenvironment subtypes for TCGA patients (available from [Bagaev et al. 2021, Table S6, sheet 1](https://ars.els-cdn.com/content/image/1-s2.0-S1535610821002221-mmc6.xlsx))
* `tcga_survival_file`: survival data for TCGA patients (available from [Liu et al. 2018, Table S1, sheet 1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx)).

```bash

# First, activate conda environment

conda env create -f spatial_features.yml
conda activate spatial_features

```

Enter R and install the following packages prior to rendering the notebook.

```r
# These packages aren't included in the conda environment

install.packages("pak")
pak::pkg_install(c("tidyverse", "GaitiLab/GaitiLabUtils", "patchwork", "Seurat", "patchwork", "ggplot2", "data.table", "glue", "stringr", "rstatix", "prabhakarlab/Banksy", "satijalab/seurat-wrappers", "ggrastr", "gridExtra" ) )

```

Alternatively, these packages are also loaded through `pacman::p_load` so if package is missing it will automatically be installed.

After installation you can render a notebook as follows in the terminal or using a bash script.

```bash

Rscript -e "rmarkdown::render('analysis/analysis/Fig1de_SuppFig3c_spatial_validation_maps_xenium_skin_panel_w_add_on.Rmd')"

# If you want to run all notebooks, you can run: 
Rscript -e "analysis/run_notebooks.R"
```
