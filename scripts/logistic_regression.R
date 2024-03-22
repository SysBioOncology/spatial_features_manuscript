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
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$min_lambda <- -6
    args$min_alpha <- -5
    args$n_alpha <- 10
    args$n_lambda <- 10
    args$top_n <- 20
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
pacman::p_load(caret)
data <- readRDS("output/CPTAC_prepped_data.rds") %>%
    column_to_rownames("Slide_ID") %>%
    select(-Vital_status_at_24months_follow_up) %>%
    filter(Vital_status_at_12months_follow_up != "Not Reported")

# Standardize features
preProcValues <- preProcess(data, method = c("center", "scale"))
data <- predict(preProcValues, data)

set.seed(998)
fitControl <- trainControl(
    method = "repeatedcv",
    # Number of folds
    number = 5,
    ## repeated ten times
    repeats = 5,
    # search = "random",
    savePredictions = TRUE,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
)

hyper_param_grid <- data.frame(expand.grid(
    lambda = rev(10^seq(args$min_lambda, 0, length = args$n_lambda)),
    alpha = 10^seq(args$min_alpha, 0, length = args$n_alpha)
))

model <- train(Vital_status_at_12months_follow_up ~ .,
    data = data,
    method = "glmnet",
    trControl = fitControl,
    verbose = FALSE,
    tuneGrid = hyper_param_grid,
    allowParallel = TRUE,
    preProc = c("center", "scale")
)

ggplot(model) +
    theme(legend.position = "top") +
    custom_theme()

# Assess Performance
res <- evalm(model)
roc_plt <- res$roc + custom_theme()
roc_plt
# Feature Importance Plot
feature_importance <- varImp(model)
feature_importance_df <- data.frame(feature_importance$importance) %>%
    arrange(desc(Overall)) %>%
    rename(feature_importance = Overall) %>%
    rownames_to_column("Feature") %>%
    head(args$top_n)

ggplot(feature_importance_df, aes(x = feature_importance, y = fct_reorder(Feature, feature_importance))) +
    geom_point() +
    geom_segment(aes(
        y = Feature,
        yend = Feature,
        x = 0,
        xend = feature_importance
    )) +
    custom_theme() +
    labs(
        x = "Feature Importance",
        y = "Feature",
        title = glue("Top-{args$top_n} features")
    )


# cv_res_df <- data.frame(model$pred)
# lift_df <- data.frame()

# for (fold in unique(cv_res_df$resample)) {
#     fold_df <- dplyr::filter(for_lift, resample == fold)
#     lift_obj_data <- lift(Class ~ rf, data = fold_df, class = "R")$data
#     lift_obj_data$fold <- fold
#     lift_df <- rbind(lift_df, lift_obj_data)
# }
