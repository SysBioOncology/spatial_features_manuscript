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
        description = "Prepare TCGA spatial features",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/testing")
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
pacman::p_load(ggplot2, ggpubr, scales, ComplexHeatmap, ggpubr)

# Load tcga features
tcga_features <- readRDS(glue("{here::here()}/data/features_matrix_validation_TCGA_SKCM.rds")) %>%
    rownames_to_column("pat_id") %>%
    mutate(pat_id = substr(pat_id, 1, 12)) %>%
    tidyr::pivot_longer(
        cols = !pat_id,
        names_to = "feature",
        values_to = "feature_score"
    )

# load Bagaev data
bagaev_mfp <- readxl::read_xlsx(path = glue("{here::here()}/data/MFP.xlsx"), skip = 1) %>%
    filter(Cohort == "TCGA-SKCM") %>%
    rename(pat_id = Sample)
head(bagaev_mfp)

# load survival data
survival_tcga_liu <- data.table::fread(file = glue("{here::here()}/data/tcga_survival_data.csv"))
head(survival_tcga_liu)


### survival info from bagaev
survival_bagaev <- bagaev_mfp %>%
    filter(pat_id %in% unique(tcga_features$pat_id)) %>%
    select(pat_id, OS, OS_FLAG, DFS, DFS_FLAG)

### survival info from liu (both PFI and OS are recommended for SKCM)
survival_liu <- survival_tcga_liu %>%
    filter(type == "SKCM", bcr_patient_barcode %in% unique(tcga_features$pat_id)) %>%
    select(bcr_patient_barcode, OS, OS.time, PFI, PFI.time) %>%
    dplyr::rename("pat_id" = "bcr_patient_barcode") %>%
    mutate(
        OS = ifelse(OS == "#N/D", NA, OS),
        OS.time = ifelse(OS.time == "#N/D", NA, OS.time),
        PFI = ifelse(PFI == "#N/D", NA, PFI),
        PFI.time = ifelse(PFI.time == "#N/D", NA, PFI.time)
    ) %>%
    mutate(
        OS = as.numeric(OS),
        OS.time = as.numeric(OS.time),
        PFI = as.numeric(PFI),
        PFI.time = as.numeric(PFI.time)
    )

# Combine spatial featurs and survival data
table_with_survival_df <- inner_join(tcga_features, survival_liu) %>%
    tidyr::pivot_wider(
        id_cols = c(pat_id, OS, OS.time, PFI, PFI.time),
        names_from = feature,
        values_from = feature_score
    ) %>%
    tibble::column_to_rownames(var = "pat_id") %>%
    filter(!(is.na(OS) & is.na(OS.time) & is.na(PFI) & is.na(PFI.time))) %>%
    mutate(
        OS = as.numeric(OS),
        OS.time = as.numeric(OS.time),
        PFI = as.numeric(PFI),
        PFI.time = as.numeric(PFI.time)
    ) %>%
    mutate_at(vars(contains("prox_clust")), ~ replace(., is.na(.), 0))
table_with_survival_df %>% nrow()

# table_with_survival_df[, grep("prox_clust", colnames(table_with_survival_df))]

# table_with_survival_df %>% head()
# print(ncol(table_with_survival_df))
# # TODO Remove shape-related features (-8 features), remove prox_ based features (-36 features)
table_with_survival_df <- table_with_survival_df[, -grep("sol_clust|round_clust", colnames(table_with_survival_df))] %>%
    drop_na() %>%
    filter(complete.cases(.))

table_with_survival_df <- table_with_survival_df[is.finite(rowSums(table_with_survival_df)), ]


saveRDS(table_with_survival_df, file = glue("{args$output_dir}/tcga_spatial_features.rds"))
# # print(ncol(table_with_survival_df))

# sum(is.na(table_with_survival_df))

# # TODO keep only features that have at least 80% valid values -> complete.cases()
only_features <- table_with_survival_df %>% select(-c(OS, OS.time, PFI, PFI.time))

head(only_features)
# Grab proximity features
# mat <- apply(data.matrix(only_features), 2, function(x) {
#     as.numeric(x)
# })
# # mat[is.na(mat), ] <- NA
# Hmisc::rcorr(as.matrix(mat))
# corr_mat <- cor(mat, use = "pairwise.complete.obs")
# corr_df <- corr_mat %>%
#     data.frame(check.names = FALSE) %>%
#     rownames_to_column("variable_x") %>%
#     reshape2::melt(id.vars = "variable_x")

# highly_corr_var <-  corr_df %>% filter(value > 0.8) %>% distinct() %>% filter(variable_x != variable)
# nrow(highly_corr_var)

# #  %>% reshape2::melt()


# corrplot::corrplot(corr_mat)



frac_valid <- apply(only_features, 2, function(col) {
    sum(!is.na(col)) / length(col)
})
only_features %>%
    complete.cases() %>%
    sum()

# minfrac_mask <- frac_valid >= 0.85
# survival_cols <- c("OS", "OS.time", "PFI", "PFI.time")
# cols_to_select <- c(names(minfrac_mask)[minfrac_mask], survival_cols)
# table_with_survival_df <- table_with_survival_df %>% select(all_of(cols_to_select))

# sum(minfrac_mask)
# n_patients <- only_features %>%
#     select(all_of(names(minfrac_mask)[minfrac_mask])) %>%
#     filter(complete.cases(.)) %>%
#     nrow()
# print(n_patients)
# n_patients_retained <- sapply(seq(0.80, 1, 0.05), function(alpha) {
#     minfrac_mask <- frac_valid >= alpha
#     n_patients <- only_features %>%
#         select(all_of(names(minfrac_mask)[minfrac_mask])) %>%
#         filter(complete.cases(.)) %>%
#         nrow()

#     return(c(alpha = alpha, n_patients = n_patients, n_features = sum(minfrac_mask)))
# })
# df <- as.data.frame(t(n_patients_retained)) %>% reshape2::melt(id.vars = "alpha")

# features_vs_patients <- ggplot(data = df) +
#     geom_bar(aes(x = alpha, y = value, fill = variable), stat = "identity", position = "dodge") +
#     custom_theme() +
#     labs(y = "Counts", x = "threshold")
# ggsave(plot = features_vs_patients, filename = glue("{args$output_dir}/features_vs_patients_retained.pdf"))


# plt_minfrac <- ggplot() +
#     geom_histogram(data = data.frame(frac_valid), aes(x = frac_valid)) +
#     custom_theme() +
#     labs(title = "Fraction of valid feature values", x = "Fraction", y = "Number of features")
# plt_minfrac
# ggsave(plot = plt_minfrac, filename = glue("{args$output_dir}/frac_not_NA_hist.pdf"))


# plt_cumul <- ggplot(data.frame(frac_valid), aes(frac_valid)) +
#     stat_bin(aes(y = length(only_features) - cumsum(after_stat(count)))) +
#     custom_theme() +
#     labs(title = "Inverse cumulative sum of features w/ fraction of valid feature values", x = "Fraction of valid feature values", y = "Number of features")
# plt_cumul

# plt_combi <- ggarrange(plt_minfrac, plt_cumul)
# ggsave(plot = plt_combi, filename = glue("{args$output_dir}/valid_values_distribution.pdf"), width = 12)


# Only keep patients who have all features
table_with_survival_os <- table_with_survival_df %>%
    filter(!is.na(OS), !is.na(OS.time))
table_with_survival_pfi <- table_with_survival_df %>%
    filter(!is.na(PFI), !is.na(PFI.time))
table_with_survival_os %>%
    nrow()

data_os <- list(
    x = table_with_survival_os %>% select(-c(OS, OS.time, PFI, PFI.time)),
    y = table_with_survival_os %>%
        select(OS.time, OS) %>%
        rename(time = OS.time, status = OS)
)

data_pfi <- list(
    x = table_with_survival_pfi %>% filter() %>% select(-c(OS, OS.time, PFI, PFI.time)),
    y = table_with_survival_pfi %>%
        select(PFI.time, PFI) %>%
        rename(time = PFI.time, status = PFI)
)

saveRDS(data_os, glue("{args$output_dir}/tcga_data_os.rds"))
saveRDS(data_pfi, glue("{args$output_dir}/tcga_data_pfi.rds"))

# ---- Check missing values ---- #
spatial_features_only <- table_with_survival_os %>% select(-c(OS, OS.time, PFI, PFI.time))
spatial_features_only[!is.na(spatial_features_only)] <- TRUE
spatial_features_only[is.na(spatial_features_only)] <- FALSE
mat <- data.matrix(spatial_features_only)
head(mat)
lgd <- Legend(labels = c("Available", "Missing"))

colors <- c("blue", "red")
names(colors) <- c(1, 0)
hm <- Heatmap(mat, col = colors, show_row_names = FALSE, column_title = "Spatial features", row_title = "Patients", name = "Presence")

pdf(glue("{args$output_dir}/tcga_data_missing.pdf"), width = 20, height = 15)
draw(hm)
dev.off()

# ---- Visualization ---- #
survival_only <- table_with_survival_df %>%
    select(OS, OS.time, PFI, PFI.time) %>%
    mutate(
        OS = factor(ifelse(OS == 0, "Alive", "Dead")),
        PFI = factor(ifelse(PFI == 0, "Alive", "Dead")),
    )

# Create long-format combined OS and PFI
survival_only_os <- survival_only %>%
    select(OS, OS.time) %>%
    rename(state = OS, time = OS.time) %>%
    mutate(type = "OS")

survival_only_pfi <- survival_only %>%
    select(PFI, PFI.time) %>%
    rename(state = PFI, time = PFI.time) %>%
    mutate(type = "PFI")
survival_only <- rbind(survival_only_os, survival_only_pfi)


plt_hist <- ggplot(data = survival_only) +
    geom_histogram(aes(x = time, fill = state), alpha = 0.8, position = "identity", bins = 50) +
    scale_fill_manual(values = c("Alive" = "royalblue", "Dead" = "darkred")) +
    custom_theme() +
    facet_grid(. ~ type) +
    scale_y_continuous(breaks = pretty_breaks()) +
    labs(y = "Number of patients")

# Number of patients per outcome and type
plt_bar <- ggplot(data = survival_only) +
    geom_bar(aes(x = state, fill = state)) +
    scale_fill_manual(values = c("Alive" = "royalblue", "Dead" = "darkred")) +
    custom_theme() +
    labs(y = "Number of patients") +
    scale_y_continuous(breaks = pretty_breaks()) +
    facet_grid(. ~ type) +
    theme(legend.position = "none")

# Save plots
ggsave(plot = plt_hist, filename = glue("{args$output_dir}/clinical_hist.pdf"))
ggsave(plot = plt_bar, filename = glue("{args$output_dir}/clinical_bar.pdf"))
