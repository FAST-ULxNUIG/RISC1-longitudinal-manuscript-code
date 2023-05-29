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


test_mvmlfd_array <- eval.fd(evalarg = 0:100, fdobj = test_mvmlfd$fd_obj)
test_mvmlfd_mat <- rbind(test_mvmlfd_array[,,1],
                         test_mvmlfd_array[,,2],
                         test_mvmlfd_array[,,3])
test_mvmlfd_dt <- data.table(test_mvmlfd$df)
test_mvmlfd_dt[, id := paste(subject_id, side, stride_ind, sep = "_")]
colnames(test_mvmlfd_mat) <- test_mvmlfd_dt$id


# Create dataset for plotting: --------------------------------------------
plot_dt <- data.table(t_grid = rep(0:100, times = 3),
                      dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
                      test_mvmlfd_mat)
plot_dt_lng <- melt.data.table(plot_dt, id.vars = c("t_grid", "dimension"))
plot_dt_lng[, subject_id := stringr::str_extract(variable, "\\d{1,2}")]
plot_dt_lng[, side := stringr::str_extract(variable, "(right|left)")]
plot_dt_lng[, stride_ind := stringr::str_extract(variable, pattern = "(?<=\\d{1,2}_(right|left)_)\\d+")]
plot_dt_lng[, long_time := (as.numeric(stride_ind) - 1) / (max(as.numeric(stride_ind)) - 1)]
plot_dt_lng[, participant_id :=  paste("Participant", subject_id)]

p <- ggplot(data = plot_dt_lng) +
  aes(x = t_grid, y = value, group = interaction(stride_ind, side), color = long_time) +
  facet_grid(dimension ~ participant_id, scales = "free_y") +
  geom_line() +
  scale_color_viridis_c(option = "D",
                        name = "Longitudinal Time:",
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        label = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust = 0.8)) +
  labs(x = "$t$", y = "$y^{(p)}_{ijl}(t)$")

p

tikz(file.path(plots_path, "simulated-mvmlfd.tex"),
     width = 1 * doc_width_inches, 
     height = 0.85 * doc_width_inches)
print(p)
dev.off()
