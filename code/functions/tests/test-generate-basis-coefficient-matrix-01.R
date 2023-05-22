# Simple test to see if function to generate matrix of mv-FPC scores ------
# actually runs. No checking the result.
source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code", "functions", "generate-basis-coefficient-matrix.R"))
library(data.table) # CRAN v1.14.2

# Some preliminary parameters: --------------------------------------------
N_test <- 300
n_i_test <- 50
sim_param <- readRDS(file = "outputs/simulation/simulation-parameters.rds")
K_test <- sim_param$k_sim
R_star_test <- Q_star_test <- array(NA, c(3, 3, K_test))
for(k in seq_len(K_test)) {
  Q_star_test[,,k] <- diag(sim_param$Q_star[,k])
  R_star_test[,,k] <- diag(sim_param$R_star[,k])
}


s_k_test <- sim_param$s_k_vec
Beta_k_test <- sim_param$Beta
beta_poly_1_test <- Beta_k_test["$\\beta_{0, 1, k}^*$", ] 
beta_poly_2_test <- Beta_k_test["$\\beta_{0, 2, k}^*$", ]
beta_poly_3_test <- Beta_k_test["$\\beta_{0, 3, k}^*$", ]
beta_sex_test <- Beta_k_test["$\\beta_{1, k}^* $", ]
beta_speed_cent_test <- Beta_k_test["$\\beta_{2, k}^* $", ]
speed_sd_test <- 1.5
longitudinal_grid <- seq(0, 1, by = 0.01)
poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)
test <- generate_design_multiple_subjects(N = N_test, n_i = n_i_test, speed_sd = speed_sd_test)
# Add polynomials to it:
test_add_poly <- data.table(add_poly_to_df(df = test, poly_object = poly_basis))
test <- generate_basis_coefficient_matrix(design_df = test_add_poly, 
                                          N = N_test,
                                          K = K_test,
                                          beta_poly_1_vec = beta_poly_1_test, 
                                          beta_poly_2_vec = beta_poly_2_test,
                                          beta_poly_3_vec = beta_poly_3_test,
                                          beta_sex_vec = beta_sex_test,
                                          beta_speed_cent_vec = beta_speed_cent_test,
                                          Q_star_array = Q_star_test,
                                          R_star_array = R_star_test, 
                                          s_vec = s_k_test)









