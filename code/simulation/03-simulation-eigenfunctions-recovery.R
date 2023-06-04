source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code", "functions", "generate-basis-coefficient-matrix.R"))
source(here::here("code", "functions", "generate-mvmllfd.R"))
source(here::here("code", "functions", "construct_fd_from_scores.R"))
source(here::here("code", "functions", "decenter_fd_around_new_mean.R"))
source(here::here("code", "functions", "theme_gunning.R"))
library(data.table) # CRAN v1.14.2
library(fda)        # CRAN v5.5.1


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
speed_sd_test <- 1.4

pca_fd_test <- sim_param$mfpca
pca_fd_test$harmonics <- pca_fd_test$harmonics[seq_len(K_test),]
pca_fd_test$scores <- pca_fd_test$scores[, seq_len(K_test),]
pca_fd_test$varprop <- pca_fd_test$varprop[seq_len(K_test)]

N_sim <- 500
test_pca_fd_list_zero_fixef <- vector(mode = "list", length = N_sim)
for(i in seq_len(N_sim)) {
  print(paste0("Iteration ", i))
  test_mvmlfd <- generate_mvmllfd(N = N_test, 
                                  pca_fd_obj = pca_fd_test,
                                  n_i = n_i_test, 
                                  speed_sd = speed_sd_test, 
                                  K = K_test, 
                                  beta_poly_1_vec = rep(0, 10), 
                                  beta_poly_2_vec = rep(0, 10),
                                  beta_poly_3_vec = rep(0, 10),
                                  # beta_sex_vec = c(0, beta_sex_test[2], rep(0, 8)),
                                  # beta_speed_cent_vec = c(beta_speed_cent_test[1], rep(0, 9)),
                                  beta_sex_vec = rep(0, 10),
                                  beta_speed_cent_vec = rep(0, 10),
                                  Q_star_array = Q_star_test,
                                  R_star_array = R_star_test, 
                                  s_vec = s_k_test
  )
  test_pca_fd_list_zero_fixef[[i]] <- pca.fd(fdobj = test_mvmlfd$fd_obj, nharm = 10)
}


test_pca_fd_list_true_fixef <- vector(mode = "list", length = N_sim)
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
  test_pca_fd_list_true_fixef[[i]] <- pca.fd(fdobj = test_mvmlfd$fd_obj, nharm = 10)
}



longitudinal_grid <- seq(0, 1, by = 0.01)
poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)

generated_scores_list <- vector(mode = "list", length = N_sim)
for(i in seq_len(N_sim)) {
  print(paste0("Iteration ", i))
  design_df <- generate_design_multiple_subjects(N = N_test,
                                                  n_i = n_i_test, 
                                                  speed_sd = speed_sd_test)
  design_df_poly <- add_poly_to_df(df = design_df, poly_object = poly_basis)
  generated_scores_list[[i]] <- generate_basis_coefficient_matrix(N = N_test,  
                                    design_df = design_df_poly,
                                    K = K_test, 
                                    beta_poly_1_vec = beta_poly_1_test, 
                                    beta_poly_2_vec = beta_poly_2_test,
                                    beta_poly_3_vec = beta_poly_3_test,
                                    beta_sex_vec = beta_sex_test,
                                    beta_speed_cent_vec = beta_speed_cent_test,
                                    Q_star_array = Q_star_test,
                                    R_star_array = R_star_test, 
                                    s_vec = s_k_test)
  }


cov_array <- cor_array <- array(NA, dim = c(10, 10, N_sim))


for(i in seq_len(N_sim)) {
  print(paste0("Iteration ", i))
  cov_array[,,i] <- cov(generated_scores_list[[i]][, paste0("score_", 1:10)])
  cor_array[,,i] <- cor(generated_scores_list[[i]][, paste0("score_", 1:10)])
}

apply(cor_array, c(1, 2), mean)
average_marginal_cov <- apply(cov_array, c(1, 2), mean)
eigen_marginal_cov <- eigen(average_marginal_cov)


pca_fd_rotated <- pca_fd_test
j<-1
for(j in 1:3) {
  pca_fd_rotated$harmonics$coefs[,,j] <- pca_fd_rotated$harmonics$coefs[,,j] %*% eigen_marginal_cov$vectors
}


for(k in 1:10) {
  harmonics_k <- pca_fd_rotated$harmonics[k, ]
  true_k <- pca_fd_test$harmonics[k, ]
  ip1 <- sum(diag(inprod(harmonics_k- true_k, harmonics_k- true_k)))
  ip2 <- sum(diag(inprod(harmonics_k + true_k, harmonics_k + true_k)))
  if(ip2 < ip1) {
    pca_fd_rotated$harmonics$coefs[,k,] <- -1* pca_fd_rotated$harmonics$coefs[,k,]
  }
}



mvFPC_df <- purrr::map_dfr(.x = test_pca_fd_list_zero_fixef, .f = function(x) {
  for(k in 1:6) {
    harmonics_k <- x$harmonics[k, ]
    true_k <- pca_fd_test$harmonics[k, ]
    ip1 <- sum(diag(inprod(harmonics_k- true_k, harmonics_k - true_k)))
    ip2 <- sum(diag(inprod(harmonics_k + true_k, harmonics_k + true_k)))
    if(ip2 < ip1) {
      x$harmonics$coefs[,k,] <- -1* x$harmonics$coefs[,k,]
    }
  }
  Psi_hat_array <- eval.fd(0:100, x$harmonics[1:6, ])
  Psi_hat <- rbind(Psi_hat_array[,,1], Psi_hat_array[,,2], Psi_hat_array[,,3])
  colnames(Psi_hat) <- paste0("mvFPC", 1:6)
  data.frame(t = rep(0:100, times = 3), dimension = rep(c("Hip", "Knee", "Ankle"), each = 101), Psi_hat)
}, .id = "sim_rep")


mvFPC_dt <- as.data.table(mvFPC_df)
mvFPC_dt_lng <- melt.data.table(mvFPC_dt,
                                id.vars = c("sim_rep", "t", "dimension"),
                                value.name = "estimate",
                                variable.name = "num",
                                measure.vars = paste0("mvFPC", 1:6),
                                variable.factor = FALSE,
                                value.factor = FALSE,
                                verbose = TRUE)

Psi_true_array <- eval.fd(0:100, pca_fd_test$harmonics[1:6, ])
Psi_true <- rbind(Psi_true_array[,,1], Psi_true_array[,,2], Psi_true_array[,,3])
colnames(Psi_true) <- paste0("mvFPC", 1:6)
true_mvFPC_dt <- data.table(t = rep(0:100, times = 3), dimension = rep(c("Hip", "Knee", "Ankle"), each = 101), Psi_true)
true_mvFPC_dt_lng <- melt.data.table(true_mvFPC_dt,
                                     id.vars = c("t", "dimension"),
                                     value.name = "estimate",
                                     variable.name = "num",
                                     measure.vars = paste0("mvFPC", 1:6),
                                     variable.factor = FALSE,
                                     value.factor = FALSE,
                                     verbose = TRUE)

Psi_true_rotated_array <- eval.fd(0:100, pca_fd_rotated$harmonics[1:6, ])
Psi_true_rotated <- rbind(Psi_true_rotated_array[,,1], Psi_true_rotated_array[,,2], Psi_true_rotated_array[,,3])
colnames(Psi_true_rotated) <- paste0("mvFPC", 1:6)
true_rotated_mvFPC_dt <- data.table(t = rep(0:100, times = 3), dimension = rep(c("Hip", "Knee", "Ankle"), each = 101), Psi_true_rotated)
true_rotated_mvFPC_dt_lng <- melt.data.table(true_rotated_mvFPC_dt,
                                     id.vars = c("t", "dimension"),
                                     value.name = "estimate",
                                     variable.name = "num",
                                     measure.vars = paste0("mvFPC", 1:6),
                                     variable.factor = FALSE,
                                     value.factor = FALSE,
                                     verbose = TRUE)











# -------------------------------------------------------------------------

mvFPC_covariates_df <- purrr::map_dfr(.x = test_pca_fd_list_true_fixef, .f = function(x) {
  for(k in 1:6) {
    harmonics_k <- x$harmonics[k, ]
    true_k <- pca_fd_test$harmonics[k, ]
    ip1 <- sum(diag(inprod(harmonics_k - true_k, harmonics_k - true_k)))
    ip2 <- sum(diag(inprod(harmonics_k + true_k, harmonics_k + true_k)))
    if(ip2 < ip1) {
      x$harmonics$coefs[,k,] <- - 1 * x$harmonics$coefs[,k,]
    }
  }
  Psi_hat_array <- eval.fd(0:100, x$harmonics[1:6, ])
  Psi_hat <- rbind(Psi_hat_array[,,1], Psi_hat_array[,,2], Psi_hat_array[,,3])
  colnames(Psi_hat) <- paste0("mvFPC", 1:6)
  data.frame(t = rep(0:100, times = 3), dimension = rep(c("Hip", "Knee", "Ankle"), each = 101), Psi_hat)
}, .id = "sim_rep")


mvFPC_covariates_dt <- as.data.table(mvFPC_covariates_df)
mvFPC_covariates_dt_lng <- melt.data.table(mvFPC_covariates_dt,
                                id.vars = c("sim_rep", "t", "dimension"),
                                value.name = "estimate",
                                variable.name = "num",
                                measure.vars = paste0("mvFPC", 1:6),
                                variable.factor = FALSE,
                                value.factor = FALSE,
                                verbose = TRUE)

# Compare: ----------------------------------------------------------------
mvFPC_dt_lng_mean <- mvFPC_dt_lng[, .(estimate = mean(estimate)), by = .(t, dimension, num)]
mvFPC_covariates_dt_lng_mean <- mvFPC_covariates_dt_lng[, .(estimate = mean(estimate)), by = .(t, dimension, num)]

mvFPC_dt_lng[, num_label := factor(x = num,
                                   levels = paste0("mvFPC", 1:6),
                                   labels = paste0("$\\boldsymbol{\\Psi}_", 1:6, "(t)$"),
                                   ordered = FALSE)]
mvFPC_covariates_dt_lng[, num_label := factor(x = num,
                                   levels = paste0("mvFPC", 1:6),
                                   labels = paste0("$\\boldsymbol{\\Psi}_", 1:6, "(t)$"),
                                   ordered = FALSE)]
mvFPC_dt_lng_mean[, num_label := factor(x = num,
                                   levels = paste0("mvFPC", 1:6),
                                   labels = paste0("$\\boldsymbol{\\Psi}_", 1:6, "(t)$"),
                                   ordered = FALSE)]
mvFPC_covariates_dt_lng_mean[, num_label := factor(x = num,
                                              levels = paste0("mvFPC", 1:6),
                                              labels = paste0("$\\boldsymbol{\\Psi}_", 1:6, "(t)$"),
                                              ordered = FALSE)]
true_mvFPC_dt_lng[, num_label := factor(x = num,
                                   levels = paste0("mvFPC", 1:6),
                                   labels = paste0("$\\boldsymbol{\\Psi}_", 1:6, "(t)$"),
                                   ordered = FALSE)]
true_rotated_mvFPC_dt_lng[, num_label := factor(x = num,
                                              levels = paste0("mvFPC", 1:6),
                                              labels = paste0("$\\boldsymbol{\\Psi}_", 1:6, "(t)$"),
                                              ordered = FALSE)]


saveRDS(object = list(
  mvFPC_dt_lng = mvFPC_dt_lng,
  mvFPC_dt_lng_mean = mvFPC_dt_lng_mean,
  mvFPC_covariates_dt_lng = mvFPC_covariates_dt_lng,
  mvFPC_covariates_dt_lng_mean = mvFPC_covariates_dt_lng_mean,
  true_mvFPC_dt_lng = true_mvFPC_dt_lng,
  true_rotated_mvFPC_dt_lng = true_rotated_mvFPC_dt_lng
), file = here::here("outputs", "simulation", "eigenfunction-estimation.rds"))





