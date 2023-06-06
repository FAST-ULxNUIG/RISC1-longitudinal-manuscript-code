# Script to source all necceary functions for the simulation: -------------
library(fda) # CRAN v5.5.1
source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "add_natural_splines_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code", "functions", "generate-basis-coefficient-matrix.R"))
source(here::here("code", "functions", "generate-mvmllfd.R"))
source(here::here("code", "functions", "construct_fd_from_scores.R"))
source(here::here("code", "functions", "decenter_fd_around_new_mean.R"))
source(here::here("code", "functions", "function-generate-smooth-noise.R"))
source(here::here("code", "functions", "pca.fd_pve_cutoff.R"))
source(here::here("code", "functions", "add_pca.fd_scores_to_df.R"))
source(here::here("code", "functions", "fit_poly.R"))
source(here::here("code", "functions", "fit_naive.R"))
source(here::here("code", "functions", "split_train_test.R"))
source(here::here("code", "functions", "calculate_prediction_error.R"))

