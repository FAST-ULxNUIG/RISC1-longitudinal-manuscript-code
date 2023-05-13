source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))
source(here::here("code", "functions", "generate_polynomial_model_basis_coefficient.R"))
source(here::here("code/functions/theme_gunning.R"))
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice)

# Path to save the outputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")

# Some settings for the plots: --------------------------------------------
theme_gunning() # use default theme
theme_update(
  axis.text = element_text(size = 10.5),
  axis.title = element_text(size = 10.5),
  strip.text = element_text(size = 10.5),
  plot.title = element_text(size = 11.5),
  legend.text = element_text(size = 10.5),
  legend.title = element_text(size = 11)
)
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Some settings: ----------------------------------------------------------
# Some preliminary parameters: --------------------------------------------
N_test <- 6
n_i_test <- 80
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

test_add_poly[, participant_id :=  paste("Participant", subject_id)]
p <- ggplot(data = test_add_poly) +
    aes(x = time, y = score_k_test, group = side, color = side) +
    facet_wrap(~ participant_id, nrow = 2, ncol = 3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.75) +
    scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
    xlab("Longitudinal Time ($T$)") +
    ylab("Simulated mv-FPC Score ($y_{ijl, 1}^*$)") +
    theme(legend.position = "bottom")

tikz(file.path(plots_path, "simulated-trajectory-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.725 * doc_width_inches)
print(p)
dev.off()






