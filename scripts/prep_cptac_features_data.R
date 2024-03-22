# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Validation with cptac",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
cptac_cm_clinical <- data.table::fread(file = "data/TCIA_CPTAC_Pathology_Portal_v2.csv") %>%
    data.frame() %>%
    filter(Vital_status_at_12months_follow_up != "Not Reported" & Last_Known_vital_status != "Not Reported" &
        Specimen_Type == "tumor_tissue" &
        Tumor_Segment_Acceptable == "Yes")
nrow(cptac_cm_clinical)

features <- readRDS("data/features_matrix_validation_cptac_cm.rds") %>%
    data.frame() %>%
    select(-pat_id) %>%
    select(-grep("sol_clust|round_clust", colnames(.))) %>%
    mutate_at(vars(grep("prox_clust|n_shortest_paths|LCC|ND|ND_effsize|Coloc", colnames(.))), ~ replace(., is.na(.), 0)) %>%
    rownames_to_column("Slide_ID")
head(features)

survival <- cptac_cm_clinical %>%
    select(Slide_ID, Vital_status_at_12months_follow_up)

prepped_spatial_features_data <- features %>%
    left_join(survival, by = "Slide_ID") %>%
    column_to_rownames("Slide_ID")

# ---- Clinical features ---- #
numeric_features <- c("Weight", "Percent_Necrosis", "Percent_Tumor_Nuclei", "Percent_Total_Cellularity", "Percent_Necrosis")
character_features <- c("Slide_ID", "Tumor_Site", "Tumor_Histological_Type", "Gender", "Race", "Vital_status_at_12months_follow_up")
clinical_feature_names <- c(numeric_features, character_features)
# clinical_feature_names
clinical_features <- cptac_cm_clinical %>%
    select(all_of(clinical_feature_names))

clinical_character <- do.call(cbind, lapply(character_features, function(colname) {
    return(cptac_cm_clinical %>%
        pull(!!sym(colname)) %>%
        replace(., "N/A", NA))
})) %>%
    as.data.frame() %>%
    remove_rownames()
colnames(clinical_character) <- character_features
nrow(clinical_character)
clinical_numeric <- cptac_cm_clinical %>%
    select(all_of(numeric_features), Slide_ID) %>%
    data.frame()

clinical_features <- clinical_numeric %>%
    left_join(clinical_character) %>%
    select_if(~ !any(is.na(.))) %>%
    column_to_rownames("Slide_ID")


saveRDS(prepped_spatial_features_data, file = glue("{args$output_dir}/CPTAC_prepped_spatial_features_data.rds"))
saveRDS(clinical_features, file = glue("{args$output_dir}/CPTAC_prepped_clinical_features_data.rds"))
