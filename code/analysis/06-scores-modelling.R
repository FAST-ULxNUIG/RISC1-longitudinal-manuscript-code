# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(splines)
source(here::here("code/functions/source_all_analysis_functions.R"))

# Path to save the outputs of analysis: ------------------------------------
outputs_path <- here::here("outputs")
# Read in results and unpack: ---------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
mfpca <- basis_transformation_results$mfpca
k_retain <- basis_transformation_results$k_retain
covariates_dt_train <- basis_transformation_results$covariates_dt_train
N_train <- basis_transformation_results$N_train


# Create scores: ----------------------------------------------------------
# scores for mv-fpca are in a 
# N_{Total} \times K \times P matrix
# Sum over the hip knee and ankle to get overall scores:
scores_train <- apply(mfpca$scores, c(1, 2), sum)
colnames(scores_train) <- paste0("score_", seq_len(k_retain))

# Some basic checks:
stopifnot(dim(scores_train) == c(N_train, k_retain))
stopifnot(nrow(scores_train) == nrow(covariates_dt_train))

# And join back up into data.table for modelling:
covariates_and_scores_dt_train <- cbind(covariates_dt_train, scores_train)
covariates_and_scores_dt_train[, time := long_time] # our functions require longitudinal time variable be named time

# Modelling: --------------------------------------------------------------
# And join back up into data.table for modelling:

fixef_formula_full <- "ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent"

# Fit Spline Model: -------------------------------------------------------
# Fit Models: -------------------------------------------------------------
# #Spline: ----------------------------------------------------------------
spline_time <- system.time(spline_model <- fit_spline(df_scores = covariates_and_scores_dt_train,
                                                      K_retain = k_retain,
                                                      df = 4,
                                                      fixef_formula = fixef_formula_full,
                                                      diagonal_covariance = FALSE))

spline_ri_time <- system.time(
  spline_ri_model <- fit_spline_subject_ri_side(df_scores = covariates_and_scores_dt_train,
                                                K_retain = k_retain, 
                                                df = 4, 
                                                fixef_formula = fixef_formula_full, 
                                                diagonal_covariance = FALSE))

# ml-FPCA: ----------------------------------------------------------------
fpca_time <- system.time(fpca_model <- fit_fpca(df_scores = covariates_and_scores_dt_train, 
                                                K_retain = k_retain, 
                                                pve = 0.995, 
                                                df_spline = 4, 
                                                fixef_formula = fixef_formula_full))
## Simple naive model -----------------------------------------------------
naive_time <- system.time(naive_model <- fit_naive_spline_intercept(df_scores = covariates_and_scores_dt_train, 
                                                   K_retain = k_retain,
                                                   df = 4, 
                                                   fixef_formula = fixef_formula_full))




# Save Results: -----------------------------------------------------------

saveRDS(
  object = list(fpca_model = fpca_model,
                      naive_model = naive_model,
                      spline_model = spline_model,
                      spline_ri_model = spline_ri_model,
                      times = list(
                        spline_time = spline_time,
                        spline_ri_time = spline_ri_time,
                        fpca_time = fpca_time,
                        naive_time = naive_time),
                covariates_and_scores_dt_train = covariates_and_scores_dt_train,
                session_info = sessionInfo()),
        here::here("outputs","model-fit-results.rds"))


