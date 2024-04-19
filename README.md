# Integrating histopathology and transcriptomics for spatial profiling of the tumor microenvironment: a melanoma case study

This repository holds the code for all analyses related to

> t.b.a

Our pipeline (SPoTLIghT) to derive spatial graph-based interpretable features from H&E tissue slides is packaged as a docker container, which is available here: 
https://github.com/SysBioOncology/spotlight_docker.

### Notebook to reproduce the figures from the manuscript

```bash

# First, activate conda environment

conda env create -f spatial_features.yml
conda activate spatial_features

# Second, render the notebook

Rscript -e "rmarkdown::render('inputfile.Rmd')"

```
