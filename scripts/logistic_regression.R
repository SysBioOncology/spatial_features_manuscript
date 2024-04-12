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
        description = "Fit a logistic regression model",
    )
    # parser$add_argument("--min_lambda", type = "numeric", help = "Choose smallest lambda to test, 10^min_lambda (integer < 0)", default = -6)
    parser$add_argument("--min_alpha", type = "numeric", help = "Choose smallest alpha to test, 10^min_alpha (integer < 0)", default = 0.1)
    parser$add_argument("--n_alpha", type = "numeric", help = "Number of alphas to test", default = 10)
    parser$add_argument("--n_lambda", type = "numeric", help = "Number of lambdas to test", default = 100)
    parser$add_argument("--k_folds", type = "numeric", help = "Number of folds", default = 5)
    parser$add_argument("--n_repeats", type = "numeric", help = "Number of repeats", default = 5)

    parser$add_argument("--class_weight", type = "numeric", help = "Class weight of minority class", default = 0.1)
    parser$add_argument("--input_file", type = "character", help = "Input file with data for modeling", default = "")
    parser$add_argument("--optimization_metric", type = "character", help = "Optimization metric: 'ROC' or 'AUC'", default = "")
    parser$add_argument("--n_cores", type = "numeric", help = "Numer of cores to use", default = 1)
    parser$add_argument("--model_name", type = "character", help = "Type of model used for naming output file")
    parser$add_argument("--searchtype", type = "character", help = "Hyperparmameter tuning approach: 'random' or 'grid', if random choose also tuneLength (default='random')", default = "random")
    parser$add_argument("--tunelength", type = "numeric", help = "If hyperparameter tuning approach is random choose number of combinations to test (default=1)", default = 1)

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$n_cores <- 4
    args$output_dir <- glue("{here::here()}/output/TESTING_PIPELINE")
    args$input_file <- "data/tcga_spatial_features.rds"
    args$model_name <- "spatial_features"

    # Hyper parameter tunign
    args$searchtype <- "random"
    args$min_lambda <- -6
    args$min_alpha <- -4
    args$n_alpha <- 3
    args$n_lambda <- 3

    # Repeated cross-validation
    args$k_folds <- 3
    args$n_repeats <- 3
    args$optimization_metric <- "AUC"

    # Class weights
    # args$min_weight <- 0.1
    # args$max_weight <- 0.9
    # args$n_weights <- 2
    args$class_weight <- 0.2
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
pacman::p_load(caret, MLeval, MLmetrics)

# For reproducibility
# set.seed(998)

log_info("Load data...")
data <- readRDS(args$input_file)

log_info("Standardize features...")
# preProcValues <- preProcess(data, method = c("center", "scale"))
# data <- predict(preProcValues, data)

log_info(glue("Create folds with kfolds={args$k_folds} and repeats={args$n_repeats}..."))
cvIndex <- createMultiFolds(factor(data$Vital_Status), k = args$k_folds, times = args$n_repeats)


fitControl <- trainControl(
    # Pre-set folds
    index = cvIndex,
    method = "repeatedcv",
    # Number of folds
    number = args$k_folds,
    ## repeated N times
    repeats = args$n_repeats,
    savePredictions = TRUE,
    # Type of evaluation metrics:
    # 'twoClassSummary': ROC, Sens, Spec
    # 'prSummary': AUC, Precision, Recall & F
    summaryFunction = ifelse(
        # Condition
        args$optimization_metric == "ROC",
        list(twoClassSummary)[[1]],
        # Else: AUC
        list(prSummary)[[1]]
    ),
    classProbs = TRUE,
    search = args$searchtype
)

# -1 because of response variable
n_features <- ncol(data) - 1
n_obs <- nrow(data)
min_lambda <- ifelse(n_obs < n_features, -2, -4)

log_info("Set up hyperparameter grid...")
hyper_param_grid <- data.frame(expand.grid(
    lambda = rev(10^seq(min_lambda, 0, length = args$n_lambda)),
    alpha = seq(0, 1, length = args$n_alpha)
))

log_info("Setup class weights...")
# Increasingly make 'Deceased' more important
# class_weight_grid <- seq(args$min_weight, args$max_weight, length = args$n_weights)

log_info("Fitting models...")
fit_model <- function(class_weight, data, trControl, tuneGrid, metric, tunelength = 1) {
    weights <- ifelse(data$Vital_Status == "Deceased",
        (1 / table(data$Vital_Status)["Deceased"]) * class_weight,
        (1 / table(data$Vital_Status)["Living"]) * (1 - class_weight)
    )
    model <- train(Vital_Status ~ .,
        data = data,
        method = "glmnet",
        trControl = trControl,
        verbose = FALSE,
        tuneGrid = tuneGrid,
        allowParallel = TRUE,
        # preProc = c("center", "scale"),
        metric = metric,
        weights = weights,
        tunelength = tunelength
    )
    # Storing used class weight
    model$custom_class_weight <- class_weight
    return(model)
}

fitted_model <- fit_model(class_weight = args$class_weight, data = data, trControl = fitControl, tuneGrid = hyper_param_grid, metric = args$optimization_metric, tunelength = args$tunelength)

log_info("Save models...")
saveRDS(
    fitted_model,
    file = glue("{args$output_dir}/models_{args$model_name}__{args$optimization_metric}__cw_{round(args$class_weight, 4)}.rds")
)

log_info("COMPLETED!")
