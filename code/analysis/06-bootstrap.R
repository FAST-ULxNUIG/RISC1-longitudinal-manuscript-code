# ------------------------------------------------------------------------#
# Do bootstrap of subjects.
# ------------------------------------------------------------------------#
# Load packages: ----------------------------------------------------------
library(lme4)       # CRAN v1.1-30
library(data.table) # CRAN v1.14.2
library(splines)
# Load functions -----------------------------------------------------------
functions_path <- here::here("code", "functions")
source(file.path(functions_path, "add_natural_splines_to_df.R"))
source(file.path(functions_path, "fit_spline_subject_ri_side.R"))
source(file.path(functions_path, "bootstrap_of_subjects.R"))
source(file.path(functions_path, "extract_fixef_coef.R"))

# Load data: --------------------------------------------------------------
model_fit_results <- readRDS(here::here("outputs", "model-fit-results.rds"))
basis_transformation_results <- readRDS(here::here("outputs", 
                                                   "basis-transformation-results.rds"))

# Unpack results: ---------------------------------------------------------
# only extract what is neccesary.
covariates_and_scores_dt_train <- model_fit_results$covariates_and_scores_dt_train
k_retain <- basis_transformation_results$k_retain

# Set some parameters: ----------------------------------------------------
n_cores <- parallel::detectCores() - 1
B <- 250
fixef_formula_full <- "ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent"

# -------------------------------------------------------------------------
bootstrap_time <- system.time(
  bootstrap_results <- bootstrap_of_subjects(
    df_for_bootstrap = covariates_and_scores_dt_train,
    k_retain = k_retain,
    model = "spline_subject_ri_side",
    df = 4, 
    fixef_formula = fixef_formula_full,
    B = B, 
    diagonal_covariance = FALSE,
    par_mc = TRUE, 
    n_cores = n_cores)
)

saveRDS(object = list(bootstrap_results = bootstrap_results,
                      bootstrap_time = bootstrap_time,
                      session_info = sessionInfo(),
                      sys_time = Sys.time()),
        here::here("outputs", "bootstrap-results.rds"))

       