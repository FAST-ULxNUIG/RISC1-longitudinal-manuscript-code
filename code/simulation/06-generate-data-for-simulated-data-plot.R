# ------------------------------------------------------------------------#
# Generate data for plot of simulated vs test data: -----------------------
# ------------------------------------------------------------------------#

# Packages and script: ----------------------------------------------------
source(here::here("code", "functions", "source_all_simulation_functions.R"))
library(data.table) # CRAN v1.14.2

# Some preliminary parameters: --------------------------------------------
N_demo <- 280
n_i_demo <- 80
N_total_demo <- N_demo * n_i_demo * 2
sim_param <- readRDS(file = "outputs/simulation/simulation-parameters.rds")
mfd_obj_test <- readRDS(file = here::here("outputs", "basis-transformation-results.rds"))$mfd_obj_test
K_demo <- sim_param$k_sim
R_star_demo <- Q_star_demo <- array(NA, c(3, 3, K_demo))
for(k in seq_len(K_demo)) {
  Q_star_demo[,,k] <- diag(sim_param$Q_star[,k] * c(1, 1, 1))
  R_star_demo[,,k] <- diag(sim_param$R_star[,k] * c(1, 1, 1))
}

s_k_demo <- sim_param$s_k_vec
Beta_k_demo <- sim_param$Beta
beta_poly_1_demo <- Beta_k_demo["$\\beta_{0, 1, k}^*$", ] 
beta_poly_2_demo <- Beta_k_demo["$\\beta_{0, 2, k}^*$", ]
beta_poly_3_demo <- Beta_k_demo["$\\beta_{0, 3, k}^*$", ]
beta_sex_demo <- Beta_k_demo["$\\beta_{1, k}^* $", ]
beta_speed_cent_demo <- Beta_k_demo["$\\beta_{2, k}^* $", ]
speed_sd_demo <- 1.54
pca_fd_demo <- sim_param$mfpca
pca_fd_demo$harmonics <- pca_fd_demo$harmonics[seq_len(K_demo),]
pca_fd_demo$scores <- pca_fd_demo$scores[, seq_len(K_demo),]
pca_fd_demo$varprop <- pca_fd_demo$varprop[seq_len(K_demo)]


# Generate Data: ----------------------------------------------------------
set.seed(1)
demo_mvmllfd <- generate_mvmllfd(N = N_demo, 
                                 pca_fd_obj = pca_fd_demo,
                                 n_i = n_i_demo, 
                                 speed_sd = speed_sd_demo, 
                                 K = K_demo, 
                                 beta_poly_1_vec = beta_poly_1_demo, 
                                 beta_poly_2_vec = beta_poly_2_demo,
                                 beta_poly_3_vec = beta_poly_3_demo,
                                 beta_sex_vec = beta_sex_demo,
                                 beta_speed_cent_vec = beta_speed_cent_demo,
                                 Q_star_array = Q_star_demo,
                                 R_star_array = R_star_demo, 
                                 s_vec = s_k_demo)
demo_noisy_mvmllfd <- add_smooth_noise_to_fd_obj(fd_obj = demo_mvmllfd$fd_obj,
                                                 sigma = 0.9,
                                                 l = 0.25)

# Do test train splitting as usual (although we could just sample 200 curves):
test_train_split <- split_train_test(N_total = N_total_demo, test_prop = 0.1)
demo_fd_obj_train <- demo_mvmllfd$fd_obj[test_train_split$train_inds,]
demo_fd_obj_test <- demo_mvmllfd$fd_obj[test_train_split$test_inds,]

# Nopw sample 200 obs randomly from true and simulated test sets:
sample_inds_true <- sample(x = seq_len(ncol(mfd_obj_test$coefs)), size = 200)
sample_inds_simulated <- sample(x = seq_len(ncol(demo_fd_obj_test$coefs)), size = 200)

par(mfrow = c(3, 2))
p<-1
plot(mfd_obj_test[sample_inds_true, p], ylim = c(-30, 80))
plot(demo_fd_obj_test[sample_inds_simulated,p], ylim = c(-30, 80))

p<-2
plot(mfd_obj_test[sample_inds_true,p])
abline(h = 50)
plot(demo_fd_obj_test[sample_inds_simulated, p])
abline(h = 50)

p<-3
plot(mfd_obj_test[sample_inds_true,p])
plot(demo_fd_obj_test[sample_inds_simulated,p])

saveRDS(object = list(mfd_obj_test = mfd_obj_test[sample_inds_true, ],
                      demo_fd_obj_test = demo_fd_obj_test[sample_inds_simulated, ]),
        file = here::here("outputs", "simulation", "simulated-data-plot.rds"))
