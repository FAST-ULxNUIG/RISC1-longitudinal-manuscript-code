source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code", "functions", "generate-basis-coefficient-matrix.R"))
source(here::here("code", "functions", "generate-mvmllfd.R"))
source(here::here("code", "functions", "construct_fd_from_scores.R"))
source(here::here("code", "functions", "decenter_fd_around_new_mean.R"))
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



test_mvmlfd <- generate_mvmllfd(N = N_test, 
                 pca_fd_obj = pca_fd_test,
                 n_i = n_i_test, 
                 speed_sd = speed_sd_test, 
                 K = K_test, 
                 beta_poly_1_vec = beta_poly_1_test, 
                 beta_poly_2_vec = beta_poly_2_test,
                 beta_poly_3_vec = beta_poly_3_test,
                 beta_sex_vec = c(0, beta_sex_test[2], rep(0, 8)),
                 beta_speed_cent_vec = c(beta_speed_cent_test[1], rep(0, 9)),
                 Q_star_array = Q_star_test,
                 R_star_array = R_star_test, 
                 s_vec = s_k_test
                 )


test_pca <- pca.fd(fdobj = test_mvmlfd$fd_obj, nharm = 10)






