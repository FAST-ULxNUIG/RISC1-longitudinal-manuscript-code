# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
source(here::here("code/functions/source_all_analysis_functions.R"))

# Load Results: -----------------------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
model_fit_results <- readRDS(
  here::here("outputs", "model-fit-results.rds"))

# Path to save the outputs of analysis: ------------------------------------
outputs_path <- here::here("outputs")

# Extract Objects from Results: -------------------------------------------
# From basis transformation results:
mfpca <- basis_transformation_results$mfpca
k_retain <- basis_transformation_results$k_retain
covariates_dt_test <- basis_transformation_results$covariates_dt_test
covariates_dt_test[, time := long_time] # needs this name for prediction.
mfd_obj_test <- basis_transformation_results$mfd_obj_test
# From model fitting results:
spline_ri_model <- model_fit_results$spline_ri_model

scores_point_est_mat <- extract_fixef_coef(lme_list_object = spline_ri_model$lme_fit_list, K_retain = k_retain)
scores_diag_var_mat <- extract_diag_var(lme_fit_list = spline_ri_model$lme_fit_list)
Psi_array <- eval.fd(evalarg = 0:100, fdobj = mfpca$harmonics)
Psi_mat <- rbind(Psi_array[,,1], Psi_array[,,2], Psi_array[,,3])
effects_names <- names(fixef(spline_ri_model$lme_fit_list[[1]]))
names(effects_names) <- effects_names

parameter_results_df <- purrr::map_dfr(.x = effects_names, 
               .f = function(x) {
                 coef_point_est_x <- scores_point_est_mat[x,, drop = TRUE]
                 coef_covar_mat_wald_x <- diag(scores_diag_var_mat[x,, drop = TRUE])
                 ci_list <- mvn_sim(coef_point_est = coef_point_est_x, 
                               coef_covar_mat = coef_covar_mat_wald_x, 
                               Psi_basis = Psi_mat, 
                               N_simulation_mvn = 10000,
                               coverage_level = 0.95)
                 # put into data frame:
                 data.frame(
                   t = rep(0:100, times = 3),
                   dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
                   point_est = ci_list[["point_est"]],
                   se_wald = ci_list[["fun_se"]],
                   sim_wald_lower = ci_list[["sim"]][["lower"]],
                   sim_wald_upper = ci_list[["sim"]][["upper"]],
                   pw_wald_lower = ci_list[["pw"]][["lower"]],
                   pw_wald_upper = ci_list[["pw"]][["upper"]]
                 )
               },
               .id = "parameter")

parameter_results_dt <- data.table(parameter_results_df)

# Save results: -----------------------------------------------------------
saveRDS(object = parameter_results_dt, file.path(outputs_path, "parameter_results_dt"))
