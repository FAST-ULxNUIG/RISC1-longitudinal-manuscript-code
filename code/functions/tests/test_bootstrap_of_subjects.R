library(splines)
functions_path <- here::here("code", "functions")
source(file.path(functions_path, "generate-basis-coefficient-matrix.R"))
source(file.path(functions_path, "add_natural_splines_to_df.R"))
source(file.path(functions_path, "add_poly_to_df.R"))
source(file.path(functions_path, "fit_spline_subject_ri_side.R"))
source(file.path(functions_path, "fit_naive_spline_intercept.R"))
source(file.path(functions_path, "generate_design.R"))
source(file.path(functions_path, "generate_polynomial_model_basis_coefficient.R"))
source(file.path(functions_path, "bootstrap_of_subjects.R"))
source(file.path(functions_path, "extract_fixef_coef.R"))

N_test <- 284
K <- 1
n_i_test <- 80
n_cores <- parallel::detectCores() - 1
design_df <- generate_design_multiple_subjects(N = N_test, n_i = n_i_test, speed_sd = 1.5)
design_df <- add_poly_to_df(df = design_df, poly_object = poly(x = seq(0, 1, length.out = 101),
                                                               degree = 2,
                                                               raw = FALSE))
sim_param_final <- readRDS(file = "outputs/simulation/simulation-parameters-final.rds")

set.seed(1996)
Y_star_df <- generate_basis_coefficient_matrix(design_df = design_df,
                                            N = N_test,
                                            K = K,
                                            beta_poly_1_vec = sim_param_final$beta_poly_1_true[seq_len(K)],
                                            beta_poly_2_vec = sim_param_final$beta_poly_2_true[seq_len(K)],
                                            beta_poly_3_vec = sim_param_final$beta_poly_3_true[seq_len(K)],
                                            beta_sex_vec = sim_param_final$beta_sex_true[seq_len(K)],
                                            beta_speed_cent_vec = sim_param_final$beta_speed_cent_true[seq_len(K)],
                                            Q_star_array = sim_param_final$Q_star_true[,,seq_len(K), drop = FALSE],
                                            R_star_array =  sim_param_final$R_star_true[,,seq_len(K), drop = FALSE],
                                            s_vec = sim_param_final$s_k_true[seq_len(K)])

test_fit <- fit_spline_subject_ri_side(df_scores = Y_star_df, K_retain = K, df = 3, diagonal_covariance = FALSE)  


boot_time <- system.time(boot_results <- bootstrap_of_subjects(df_for_bootstrap = Y_star_df,
                      k_retain = K,
                      model = "spline_subject_ri_side",
                      df = 3,
                      diagonal_covariance = FALSE, 
                      B = 250, 
                      par_mc = TRUE,
                      n_cores = n_cores))
print(paste0(100 *  mean(sapply(boot_results, function(x) {x[["singular"]]})),
             "% of fits were singular!"))
print(paste("bootstrap took", boot_time["elapsed"]/60, "mins"))


score_1_results <- sapply(boot_results, function(x) {x$fixef[,1]})
par(mfrow = c(1, 2))
plot(x = sqrt(diag(vcov(test_fit$lme_fit_list[[1]]))), 
     y = apply(score_1_results, 1, sd), 
     col = 1:6,
     pch = paste0(1:6),
     xlab = "Model SEs",
     ylab = "B.o.B SEs")
abline(0, 1, col = "grey")
title("bootstrap vs. model SEs")

boxplot(t(score_1_results), col = 1:6)
abline(h = fixef(test_fit$lme_fit_list[[1]]), col = 1:6, lty = 2)
title("Bootstrap distributions of estimators")
