source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code", "functions", "generate-basis-coefficient-matrix.R"))
source(here::here("code", "functions", "generate-mvmllfd.R"))
source(here::here("code", "functions", "construct_fd_from_scores.R"))
source(here::here("code", "functions", "decenter_fd_around_new_mean.R"))
source(here::here("code", "functions", "function-generate-smooth-noise.R"))
library(data.table) # CRAN v1.14.2
library(fda)

# Some preliminary parameters: --------------------------------------------
N_test <- 280
n_i_test <- 80
N_total_test <- N_test * n_i_test * 2
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
speed_sd_test <- 1.54

pca_fd_test <- sim_param$mfpca
pca_fd_test$harmonics <- pca_fd_test$harmonics[seq_len(K_test),]
pca_fd_test$scores <- pca_fd_test$scores[, seq_len(K_test),]
pca_fd_test$varprop <- pca_fd_test$varprop[seq_len(K_test)]

N_sim <- 20
cum_varpop <- array(NA, dim = c(20, 15, 3))
for(i in seq_len(20)) {
  print(paste0("Iteration ", i))
  test_mvmlfd <- generate_mvmllfd(N = N_test, 
                                  pca_fd_obj = pca_fd_test,
                                  n_i = n_i_test, 
                                  speed_sd = speed_sd_test, 
                                  K = K_test, 
                                  beta_poly_1_vec = beta_poly_1_test, 
                                  beta_poly_2_vec = beta_poly_2_test,
                                  beta_poly_3_vec = beta_poly_3_test,
                                  beta_sex_vec = beta_sex_test,
                                  beta_speed_cent_vec = beta_speed_cent_test,
                                  Q_star_array = Q_star_test,
                                  R_star_array = R_star_test, 
                                  s_vec = s_k_test
  )
  
  test_noisy_0.85 <- add_smooth_noise_to_fd_obj(fd_obj = test_mvmlfd$fd_obj, sigma = 0.85, l = 0.25)
  test_noisy_pca_0.85 <- pca.fd(test_noisy_0.85, nharm = 15)
  cum_varpop[i,,1] <- cumsum(test_noisy_pca_0.85$varprop)
  
  test_noisy_0.9 <- add_smooth_noise_to_fd_obj(fd_obj = test_mvmlfd$fd_obj, sigma = 0.9, l = 0.25)
  test_noisy_pca_0.9 <- pca.fd(test_noisy_0.9, nharm = 15)
  cum_varpop[i,,2] <- cumsum(test_noisy_pca_0.9$varprop)
  
  test_noisy_0.95 <- add_smooth_noise_to_fd_obj(fd_obj = test_mvmlfd$fd_obj, sigma = 0.95, l = 0.25)
  test_noisy_pca_0.95 <- pca.fd(test_noisy_0.95, nharm = 15)
  cum_varpop[i,,3] <- cumsum(test_noisy_pca_0.95$varprop)
}


boxplot(cum_varpop[,11,])
abline(h=0.995)



# -------------------------------------------------------------------------
N_sim <- 100
cum_varpop_02 <- matrix(NA, nrow = N_sim, ncol = 15)
for(i in seq_len(N_sim)) {
  print(paste0("Iteration ", i))
  test_mvmlfd <- generate_mvmllfd(N = N_test, 
                                  pca_fd_obj = pca_fd_test,
                                  n_i = n_i_test, 
                                  speed_sd = speed_sd_test, 
                                  K = K_test, 
                                  beta_poly_1_vec = beta_poly_1_test, 
                                  beta_poly_2_vec = beta_poly_2_test,
                                  beta_poly_3_vec = beta_poly_3_test,
                                  beta_sex_vec = beta_sex_test,
                                  beta_speed_cent_vec = beta_speed_cent_test,
                                  Q_star_array = Q_star_test,
                                  R_star_array = R_star_test, 
                                  s_vec = s_k_test
  )
  

  test_noisy_0.9 <- add_smooth_noise_to_fd_obj(fd_obj = test_mvmlfd$fd_obj, sigma = 0.9, l = 0.25)
  remove_ind <- sample(seq_len(N_total_test), size = 0.5 * N_total_test, replace = FALSE)
  test_noisy_0.9_sparse <- test_noisy_0.9[-remove_ind, ]
  test_noisy_pca_0.9 <- pca.fd(test_noisy_0.9, nharm = 15)
  cum_varpop_02[i,] <- cumsum(test_noisy_pca_0.9$varprop)
}

boxplot(cum_varpop_02[,10])
abline(h = 0.995)




