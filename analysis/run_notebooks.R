# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)


notebooks <- list.files(file.path(here::here(), "analysis"), pattern = ".Rmd", full.names = TRUE)


for(notebook in notebooks) {
    rmarkdown::render(notebook, envir = new.env())
}
