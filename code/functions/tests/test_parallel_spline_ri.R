library(splines)
source("code/functions/generate-basis-coefficient-matrix.R")
source("code/functions/add_natural_splines_to_df.R")
source("code/functions/add_poly_to_df.R")
source("code/functions/fit_spline_subject_ri_side.R")
source("code/functions/generate_design.R")
source("code/functions/generate_polynomial_model_basis_coefficient.R")
N_test <- 280
design_df <- generate_design_multiple_subjects(N = N_test, n_i = 80, speed_sd = 1.5)
design_df <- add_poly_to_df(df = design_df, poly_object = poly(x = seq(0, 1, length.out = 101),
                                                               degree = 2,
                                                               raw = FALSE))
sim_param_final <- readRDS(file = "outputs/simulation/simulation-parameters-final.rds")
N_sim <- 20
times_mat <- matrix(NA, nrow = N_sim, ncol = 2)
fixef_array <- array(NA, dim = c(2, 10, N_sim))
seeds_list <- vector("list", length = N_sim)
set.seed(1)
for(i in seq_len(N_sim)) {
  print(paste("iteration", i))
  seeds_list[[i]] <- .Random.seed
  # generate data:
  Y_star <- generate_basis_coefficient_matrix(design_df = design_df,
                                              N = N_test, 
                                              K = sim_param_final$K_true, 
                                              beta_poly_1_vec = sim_param_final$beta_poly_1_true, 
                                              beta_poly_2_vec = sim_param_final$beta_poly_2_true,
                                              beta_poly_3_vec = sim_param_final$beta_poly_3_true,
                                              beta_sex_vec = sim_param_final$beta_sex_true,
                                              beta_speed_cent_vec = sim_param_final$beta_speed_cent_true,
                                              Q_star_array = sim_param_final$Q_star_true,
                                              R_star_array =  sim_param_final$R_star_true, 
                                              s_vec = sim_param_final$s_k_true)
  # Fit normally:
  times_mat[i, 1] <- system.time(
    test1 <- fit_spline_subject_ri_side(df_scores = Y_star, K_retain = 10, df = 3, diagonal_covariance = FALSE))["elapsed"]
  # Fit in parallel:
  times_mat[i, 2] <- system.time(
    test2 <- fit_spline_subject_ri_side(df_scores = Y_star,
                                        K_retain = 10, 
                                        df = 3, 
                                        diagonal_covariance = FALSE,
                                        parallel = TRUE,
                                        ncores = parallel::detectCores() - 1))["elapsed"]
  # and make sure it passess all checks:
  testthat::expect_equal(sapply(test1$lme_fit_list, fitted),
                         sapply(test2$lme_fit_list, fitted))
  testthat::expect_equal(sapply(test1$lme_fit_list, resid),
                         sapply(test2$lme_fit_list, resid))
  testthat::expect_equal(sapply(test1$lme_fit_list, ranef),
                         sapply(test2$lme_fit_list, ranef))
  testthat::expect_equal(sapply(test1$lme_fit_list, fixef),
                         sapply(test2$lme_fit_list, fixef))
  testthat::expect_equal(sapply(test1$lme_fit_list, fixef),
                         sapply(test2$lme_fit_list, fixef))
  testthat::expect_equal(sapply(test1$lme_fit_list, REMLcrit),
                         sapply(test2$lme_fit_list, REMLcrit))
  testthat::expect_equal(sapply(test1$lme_fit_list, AIC),
                         sapply(test2$lme_fit_list, AIC))
  # and save fixed effects:
  fixef_array[,,i] <- sapply(test1$lme_fit_list, function(x) {fixef(x)[c("sexfemale", "speed_cent")]})
  }

# passed all tests + about 3x speed up!
boxplot(times_mat[,1]/times_mat[,2])






                       
