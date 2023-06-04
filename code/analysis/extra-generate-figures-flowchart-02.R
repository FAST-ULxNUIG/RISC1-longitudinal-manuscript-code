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
library(tikzDevice) # CRAN v0.12.3.1
theme_gunning()


# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
plots_path <- here::here("outputs", "figures")

# Some preliminary parameters: --------------------------------------------
N_test <- 3
n_i_test <- 10
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


set.seed(1996)
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



y1 <- eval.fd(0:100, test_mvmlfd$fd_obj[1,])

pdf(file = file.path(plots_path, "y1-hip.pdf"), width = 3, height = 3)
plot(y1[,1], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "y1-knee.pdf"), width = 3, height = 3)
plot(y1[,2], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "y1-ankle.pdf"), width = 3, height = 3)
plot(y1[,3], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()


yN <- eval.fd(0:100, test_mvmlfd$fd_obj[60,])

pdf(file = file.path(plots_path, "yN-hip.pdf"), width = 3, height = 3)
plot(yN[,1], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "yN-knee.pdf"), width = 3, height = 3)
plot(yN[,2], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "yN-ankle.pdf"), width = 3, height = 3)
plot(yN[,3], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()





sub_01_right_inds <- which(test_mvmlfd$df$subject_id == 1 & test_mvmlfd$df$side == "right" & test_mvmlfd$df$stride_ind %in% 2:11)


sub_01_right_hip <- eval.fd(0:100, fdobj = test_mvmlfd$fd_obj)[, sub_01_right_inds,1]
sub_01_right_knee <- eval.fd(0:100, fdobj = test_mvmlfd$fd_obj)[, sub_01_right_inds,2]
sub_01_right_ankle <- eval.fd(0:100, fdobj = test_mvmlfd$fd_obj)[, sub_01_right_inds,3]


pdf(file = file.path(plots_path, "long-trad-pred.pdf"), width = 7, height = 3)
par(mfrow = c(3, 1))
par(mar = rep(0, 4))
not_missing_hip <- c(sub_01_right_hip)
not_missing_hip[303:404] <- NA
plot(not_missing_hip, type = "l", lty = 1, lwd = 2, col = scales::alpha(colour = "red", alpha = 0.75), yaxt = "n",  xaxt = "n", xlab = NA, ylab = NA, axes = F)
lines(c(sub_01_right_hip), lwd = 2, lty =3, col = scales::alpha(colour = "red", alpha = 1))
abline(v = seq(0, 1000, by = 101), col = "grey")

not_missing_knee <- c(sub_01_right_knee)
not_missing_knee[303:404] <- NA
plot(not_missing_knee, type = "l", lty = 1, lwd = 2, col = scales::alpha(colour = "orange", alpha = 0.75), yaxt = "n",  xaxt = "n", xlab = NA, ylab = NA, axes = F)
lines(c(sub_01_right_knee), lwd = 2, lty =3, col = scales::alpha(colour = "orange", alpha = 1))
abline(v = seq(0, 1000, by = 101), col = "grey")

not_missing_ankle <- c(sub_01_right_ankle)
not_missing_ankle[303:404] <- NA
plot(not_missing_ankle, type = "l", lty = 1, lwd = 2, col = scales::alpha(colour = "blue", alpha = 0.75), yaxt = "n",  xaxt = "n", xlab = NA, ylab = NA, axes = F)
lines(c(sub_01_right_ankle), lwd = 2, lty =3, col = scales::alpha(colour = "blue", alpha = 1))
abline(v = seq(0, 1000, by = 101), col = "grey")
dev.off()







# And fixed effects: ------------------------------------------------------
Psi_eval <- eval.fd(0:100, sim_param$mfpca$harmonics[1:10,])
Psi_mat <- rbind(Psi_eval[,,1], Psi_eval[,,2], Psi_eval[,,3])
Beta_sex <- (sim_param$Beta[4,, drop = FALSE] %*% t(Psi_mat))[1, ]



pdf(file = file.path(plots_path, "beta-sex-hip.pdf"), width = 3, height = 3)
plot(Beta_sex[1:101], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "beta-sex-knee.pdf"), width = 3, height = 3)
plot(Beta_sex[102:202], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "beta-sex-ankle.pdf"), width = 3, height = 3)
plot(Beta_sex[203:303], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()

Beta_speed <- (sim_param$Beta[5,, drop = FALSE] %*% t(Psi_mat))[1, ]

pdf(file = file.path(plots_path, "beta-speed-hip.pdf"), width = 3, height = 3)
plot(Beta_speed[1:101], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "beta-speed-knee.pdf"), width = 3, height = 3)
plot(Beta_speed[102:202], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "beta-speed-ankle.pdf"), width = 3, height = 3)
plot(Beta_speed[203:303], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()

