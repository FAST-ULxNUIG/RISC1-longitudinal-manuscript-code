library(fda) # CRAN v5.5.1
library(data.table)
library(tikzDevice)

plots_path <- here::here("outputs", "figures")
file.exists(plots_path)
sim_param <- readRDS(file = "outputs/simulation/simulation-parameters.rds")

Psi_1 <- eval.fd(evalarg = 0:100, fdobj = sim_param$mfpca$harmonics[1, ])

pdf(file = file.path(plots_path, "psi-1-hip.pdf"), width = 3, height = 3)
plot(Psi_1[,1], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "psi-1-knee.pdf"), width = 3, height = 3)
plot(Psi_1[,2], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "psi-1-ankle.pdf"), width = 3, height = 3)
plot(Psi_1[,3], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()

Psi10 <- eval.fd(evalarg = 0:100, fdobj = sim_param$mfpca$harmonics[10, ])


pdf(file = file.path(plots_path, "psi-10-hip.pdf"), width = 3, height = 3)
plot(Psi10[,1], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "psi-10-knee.pdf"), width = 3, height = 3)
plot(Psi10[,2], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "psi-10-ankle.pdf"), width = 3, height = 3)
plot(Psi10[,3], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()


mu <- eval.fd(0:100, fdobj = sim_param$mfpca$meanfd)

pdf(file = file.path(plots_path, "mu-hip.pdf"), width = 3, height = 3)
plot(mu[,1,1], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "red")
dev.off()

pdf(file = file.path(plots_path, "mu-knee.pdf"), width = 3, height = 3)
plot(mu[,1,2], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "orange")
dev.off()

pdf(file = file.path(plots_path, "mu-ankle.pdf"), width = 3, height = 3)
plot(mu[,1,3], type = "l", xaxt = "n", lwd = 5, yaxt = "n", xlab = NA, axes = FALSE, ylab = NA, col = "blue")
dev.off()




# -------------------------------------------------------------------------

source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code/functions/theme_gunning.R"))

# Some settings: ----------------------------------------------------------
# Some preliminary parameters: --------------------------------------------
N_test <- 3
n_i_test <- 10
sim_param <- readRDS(file = "outputs/simulation/simulation-parameters.rds")
Q_star_k_test <- diag(sim_param$Q_star[,1] * c(1, 1, 1))
R_star_k_test <- diag(sim_param$R_star[,1] * c(1, 1, 1))

s_k_test <- sim_param$s_k_vec[1]
Beta_k_test <- sim_param$Beta[, 1]
speed_sd_test <- 1.5

longitudinal_grid <- seq(0, 1, by = 0.01)
poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)

set.seed(1996)
test <- generate_design_multiple_subjects(N = N_test,
                                          n_i = n_i_test, 
                                          speed_sd = speed_sd_test)

# Add polynomials to it:
test_add_poly <- data.table(add_poly_to_df(df = test, poly_object = poly_basis))

# Generate random score:
test_add_poly$score_k_test <- generate_polynomial_model_basis_coefficient(
  design_df = test_add_poly, 
  N = N_test, 
  beta_poly_1_k = Beta_k_test["$\\beta_{0, 1, k}^*$"], 
  beta_poly_2_k = Beta_k_test["$\\beta_{0, 2, k}^*$"],
  beta_poly_3_k = Beta_k_test["$\\beta_{0, 3, k}^*$"],
  beta_sex_k = Beta_k_test["$\\beta_{1, k}^* $"],
  beta_speed_cent_k = Beta_k_test["$\\beta_{2, k}^* $"],
  Q_star_k = Q_star_k_test,
  R_star_k = R_star_k_test,
  s_k = s_k_test)

library(ggplot2)

p <- ggplot(data = test_add_poly) +
  aes(x = time, y = score_k_test, group = interaction(subject_id, side), color = factor(subject_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  theme_light() +
  xlab("Longitudinal Time") +
  ylab("") +
  theme(legend.position = "none") +
  theme(#axis.line = element_line(colour = "black", size=1),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 4,colour = 1),
        panel.background = element_blank(),
        axis.title = element_text(size = 20)) 

p

tikz(file.path(plots_path, "long-model-01.tex"),
     width = 5, 
     height = 3,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(file = file.path(plots_path, "long-model-01.tex"))



# Generate random score:
test_add_poly$score_k_test<-NULL
test_add_poly$score_2_test <- generate_polynomial_model_basis_coefficient(
  design_df = test_add_poly, 
  N = N_test, 
  beta_poly_1_k = Beta_k_test["$\\beta_{0, 1, k}^*$"], 
  beta_poly_2_k = Beta_k_test["$\\beta_{0, 2, k}^*$"],
  beta_poly_3_k = Beta_k_test["$\\beta_{0, 3, k}^*$"],
  beta_sex_k = Beta_k_test["$\\beta_{1, k}^* $"],
  beta_speed_cent_k = Beta_k_test["$\\beta_{2, k}^* $"],
  Q_star_k = Q_star_k_test,
  R_star_k = R_star_k_test,
  s_k = s_k_test)


p <- ggplot(data = test_add_poly) +
  aes(x = time, y = score_2_test, group = interaction(subject_id, side), color = factor(subject_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  theme_light() +
  xlab("Longitudinal Time") +
  ylab("") +
  theme(legend.position = "none") +
  theme(#axis.line = element_line(colour = "black", size=1),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 4,colour = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 20)) 

p

tikz(file.path(plots_path, "long-model-02.tex"),
     width = 5, 
     height = 3,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(file = file.path(plots_path, "long-model-02.tex"))


# Generate random score:
test_add_poly$score_2_test<-NULL
test_add_poly$score_3_test <- generate_polynomial_model_basis_coefficient(
  design_df = test_add_poly, 
  N = N_test, 
  beta_poly_1_k = Beta_k_test["$\\beta_{0, 1, k}^*$"], 
  beta_poly_2_k = Beta_k_test["$\\beta_{0, 2, k}^*$"],
  beta_poly_3_k = Beta_k_test["$\\beta_{0, 3, k}^*$"],
  beta_sex_k = Beta_k_test["$\\beta_{1, k}^* $"],
  beta_speed_cent_k = Beta_k_test["$\\beta_{2, k}^* $"],
  Q_star_k = Q_star_k_test,
  R_star_k = R_star_k_test,
  s_k = s_k_test)


p <- ggplot(data = test_add_poly) +
  aes(x = time, y = score_3_test, group = interaction(subject_id, side), color = factor(subject_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  theme_light() +
  xlab("Longitudinal Time") +
  ylab("") +
  theme(legend.position = "none") +
  theme(#axis.line = element_line(colour = "black", size=1),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 4,colour = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 20)) 

p

tikz(file.path(plots_path, "long-model-03.tex"),
     width = 5, 
     height = 3,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(file = file.path(plots_path, "long-model-03.tex"))
