# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require("GaitiLabUtils")

set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Prepare TCGA spatial features",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/data")
    args$spatial_features <- glue("{here::here()}/data/features_matrix_validation_TCGA_SKCM.rds")
    args$survival_data <- "data/survival.xlsx"
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
pacman::p_load(readxl)

log_info("Load tcga features...")
tcga_features <- readRDS(args$spatial_features) %>%
    rownames_to_column("bcr_patient_barcode") %>%
    mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
    tidyr::pivot_longer(
        cols = !bcr_patient_barcode,
        names_to = "feature",
        values_to = "feature_score"
    )

log_info("Load survival data, and filter...")
survival <- read_excel(args$survival_data) %>%
    filter(
        !ajcc_pathologic_tumor_stage %in% c(
            "[Not Available]",
            "[Not Applicable]",
            "Unknown",
            "[Discrepancy]"
        ), !str_detect(ajcc_pathologic_tumor_stage, "Stage IV"),
        type == "SKCM"
    ) %>%
    mutate(Vital_Status = factor(case_when(
        (death_days_to >= 365) ~ "Living",
        (last_contact_days_to >= 365) ~ "Living",
        death_days_to < 365 ~ "Deceased",
        .default = NA
    ), levels = c("Living", "Deceased"))) %>%
    # Only select columns of interest
    select(bcr_patient_barcode, Vital_Status)

print(table(survival %>% pull(Vital_Status)))
#    Alive Deceased
#      350       24

# Combine spatial featurs and survival data
table_with_survival_df <- inner_join(tcga_features, survival) %>%
    tidyr::pivot_wider(
        id_cols = c(bcr_patient_barcode, Vital_Status),
        names_from = feature,
        values_from = feature_score
    ) %>%
    tibble::column_to_rownames(var = "bcr_patient_barcode") %>%
    select(-grep("sol_clust|round_clust", colnames(.))) %>%
    mutate_at(vars(grep("prox_clust|n_shortest_paths|LCC|ND|ND_effsize|Coloc", colnames(.))), ~ replace(., is.na(.), 0)) %>%
    drop_na() %>%
    filter(complete.cases(.))

print(table_with_survival_df %>%
    pull(Vital_Status) %>%
    table())
# .
# Deceased   Living
#       15      265

log_info("Save data...")
saveRDS(table_with_survival_df, file = glue("{args$output_dir}/tcga_spatial_features.rds"))



log_info("COMPLETED!")
