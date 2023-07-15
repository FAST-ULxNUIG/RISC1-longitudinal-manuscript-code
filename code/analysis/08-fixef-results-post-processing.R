# ------------------------------------------------------------------------#
# Post processing of fixed effects from model fitting results.
# ------------------------------------------------------------------------#

# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
source(here::here("code", "functions",
                  "source_all_analysis_functions.R"))

# Load Results: -----------------------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
model_fit_results <- readRDS(
  here::here("outputs", "model-fit-results.rds"))
bootstrap_results_list <- readRDS(
  file = here::here("outputs", "bootstrap-results.rds"))

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
# And from Bootstrap:
bootstrap_results <- bootstrap_results_list$bootstrap_results
B <- 250 # no. of bootstrap replicates. ideally should have been saved!

# Compute some necessary parameters: --------------------------------------
# point estimates:
scores_point_est_mat <- extract_fixef_coef(
  lme_list_object = spline_ri_model$lme_fit_list,
  K_retain = k_retain)
# basis function evaluations:
Psi_array <- eval.fd(evalarg = 0:100, fdobj = mfpca$harmonics)
Psi_mat <- rbind(Psi_array[,,1], Psi_array[,,2], Psi_array[,,3])
# names of the parameters:
effects_names <- names(fixef(spline_ri_model$lme_fit_list[[1]]))
names(effects_names) <- effects_names



# 1) Wald intervals: ------------------------------------------------------
# Just combine model-based varaince estimates across scores.

# extract variances
scores_diag_var_mat <- extract_diag_var(lme_fit_list = spline_ri_model$lme_fit_list)
# compute pointwise and joint intervals based off this:
parameter_results_wald_df <- purrr::map_dfr(.x = effects_names, 
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

parameter_results_wald_dt <- data.table(parameter_results_wald_df)



# 2) Bootstrap Intervals: -------------------------------------------------
# Re-organise bootstrap samples into a list of matrices
# (1 list element = 1 fixed effect)
fixef_coef_samples_boot <- lapply(
  effects_names,
  FUN = function(y) {
    t(sapply(bootstrap_results, function(x) {
      x$fixef[y, ]
    }))
  }
)

# Sense check dimensions of resultant matrix:
stopifnot(all(sapply(fixef_coef_samples_boot, function(x) {
  all(dim(x) == c(B, k_retain))
})))


# Sense check dimensions of resultant matrix:
stopifnot(all(sapply(fixef_coef_samples_boot, function(x) {
  all(dim(x) == c(B, k_retain))
})))


# Compute empirical Bootstrap covariance matrices for coefficients:
fixef_coef_covar_boot <- lapply(fixef_coef_samples_boot, var) 

parameter_results_boot_df <- purrr::map_dfr(.x = effects_names, 
                                            .f = function(x) {
                                              coef_point_est_x <- scores_point_est_mat[x,, drop = TRUE]
                                              coef_covar_mat_boot_x <- fixef_coef_covar_boot[[x]]
                                              ci_list <- mvn_sim(coef_point_est = coef_point_est_x, 
                                                                 coef_covar_mat = coef_covar_mat_boot_x, 
                                                                 Psi_basis = Psi_mat, 
                                                                 N_simulation_mvn = 10000,
                                                                 coverage_level = 0.95)
                                              # put into data frame:
                                              data.frame(
                                                t = rep(0:100, times = 3),
                                                dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
                                                point_est = ci_list[["point_est"]],
                                                se_boot = ci_list[["fun_se"]],
                                                sim_boot_lower = ci_list[["sim"]][["lower"]],
                                                sim_boot_upper = ci_list[["sim"]][["upper"]],
                                                pw_boot_lower = ci_list[["pw"]][["lower"]],
                                                pw_boot_upper = ci_list[["pw"]][["upper"]]
                                              )
                                            },
                                            .id = "parameter")

parameter_results_boot_dt <- data.table(parameter_results_boot_df)




# Join up estimates: ------------------------------------------------------
parameter_results_dt <- merge.data.table(x = parameter_results_wald_dt, 
                                         y = parameter_results_boot_dt, 
                                         by = c("t", "dimension", "parameter", "point_est"), 
                                         all = TRUE)
# check join didn't mis-join or create duplicates
stopifnot(nrow(parameter_results_dt) == nrow(parameter_results_wald_dt))
stopifnot(nrow(parameter_results_dt) == nrow(parameter_results_boot_dt))
# and check again:
stopifnot(fsetequal(
  unique(parameter_results_boot_dt[, c("t", "dimension", "parameter", "point_est")]),
  unique(parameter_results_dt[, c("t", "dimension", "parameter", "point_est")])))

# Save results: -----------------------------------------------------------
saveRDS(object = parameter_results_dt, file.path(outputs_path, "parameter_results_dt"))
