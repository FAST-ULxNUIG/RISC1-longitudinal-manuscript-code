# ------------------------------------------------------------------------#
# Obtain parameters for simulation using an initial simple model fit -----#
# ------------------------------------------------------------------------#

# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(xtable)     # CRAN v1.8-4 
library(tikzDevice) # CRAN v0.12.3.1
library(ggplot2)    # CRAN v3.4.0
source(here::here("code/functions/theme_gunning.R"))

# Some settings for the plots: --------------------------------------------
theme_gunning() # use default theme
theme_update(
  axis.text = element_text(size = 10.5),
  axis.title = element_text(size = 10.5),
  strip.text = element_text(size = 10.5),
  plot.title = element_text(size = 11.5),
  legend.text = element_text(size = 10.5)
)
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))

# Path to save the outputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
outputs_path <- here::here("outputs")

# Read in results and unpack: ---------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
mfpca <- basis_transformation_results$mfpca
k_retain <- basis_transformation_results$k_retain
covariates_dt_train <- basis_transformation_results$covariates_dt_train
N_train <- basis_transformation_results$N_train


# Set some parameters for the simulation: ---------------------------------
k_sim <- 10 # Number of basis functions used in the simulation:
grid_points <- 0:100


# Extract the mv-FPCA scores: ---------------------------------------------
# scores for mv-fpca are in a 
# N_{Total} \times K \times P matrix
# Sum over the hip knee and ankle to get overall scores:
scores_train <- apply(mfpca$scores[,seq_len(k_sim),], c(1, 2), sum)
colnames(scores_train) <- paste0("score_", seq_len(k_sim))

# Some basic checks:
stopifnot(dim(scores_train) == c(N_train, k_sim))
stopifnot(nrow(scores_train) == nrow(covariates_dt_train))

# And join back up into data.table for modelling:
covariates_and_scores_dt_train <- cbind(covariates_dt_train, scores_train)


# Extract the mv-FPC basis functions: -------------------------------------
Psi_array <- eval.fd(evalarg = grid_points, fdobj = mfpca$harmonics)[,seq_len(k_sim),]
mu_array <- eval.fd(evalarg = grid_points,
                    fdobj = mfpca$meanfd)



# Construct the longitudinally varying basis functions: -------------------
longitudinal_grid <- seq(0, 1, by = 0.01)
poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)
matplot(longitudinal_grid, poly_basis, type = "l")
dev.off()
poly_coefs <- attr(poly_basis, "coefs")


# Double check construction of longitudinally-varying basis functi --------
xi_1 <- function(T) {
  1 / sqrt(poly_coefs$norm2[2])
}
xi_2 <- function(T) {
  (T - poly_coefs$alpha[1]) / sqrt(poly_coefs$norm2[3])
}
xi_3 <- function(T) {
  ((T - poly_coefs$alpha[2]) * sqrt(poly_coefs$norm2[3]) * xi_2(T) - poly_coefs$norm2[3] / sqrt(poly_coefs$norm2[2]) * xi_1(T)) / sqrt(poly_coefs$norm2[4])
}
  
par(mfrow = c(1, 2))
plot(longitudinal_grid, poly_basis[, 1],type = "l", lwd = 2)
lines(longitudinal_grid, xi_2(longitudinal_grid), lty = 2, col = "pink", lwd = 2)

plot(longitudinal_grid, poly_basis[, 2], type = "l", lwd = 2)
lines(longitudinal_grid, xi_3(longitudinal_grid), lty = 2, col = "pink", lwd = 2)
dev.off()



# Set up data.table for modelling fitting: --------------------------------
poly_longtime <- poly(covariates_and_scores_dt_train$long_time, degree = 2, coefs = poly_coefs)
colnames(poly_longtime) <- paste0("poly_", 1:2)
covariates_and_scores_dt_train <- cbind(covariates_and_scores_dt_train, poly_longtime)



# Fit series of scalar mixed models: --------------------------------------
lme_simulation_list <- vector(mode = "list", length = k_sim)
for(k in seq_len(k_sim)) {
  print(paste("Score", k))
  formula_k <- paste0("score_", k, " ~ poly_1 + poly_2 + sex + speed_cent + (poly_1 + poly_2||subject_id) + (poly_1 + poly_2||subject_id:side)")
  formula_k <- formula(formula_k)
  lme_simulation_list[[k]] <- lmer(
    formula = formula_k,
    data = covariates_and_scores_dt_train,
    control= lmerControl(optCtrl = list(xtol_rel=0,
                                                   xtol_abs=1e-10,
                                                   ftol_rel=0, 
                                                   ftol_abs=1e-10)))
}


# Extract Fixed-Effects Coefficients: -------------------------------------
Beta <- round(sapply(lme_simulation_list, fixef))
colnames(Beta) <- paste(seq_len(k_sim))
rownames(Beta) <- c("$\\beta_{0, 1, k}^*$",
                    "$\\beta_{0, 2, k}^*$",
                    "$\\beta_{0, 3, k}^*$",
                    "$\\beta_{1, k}^* $",
                    "$\\beta_{2, k}^* $")

## Format into table for paper: -------------------------------------------
bold <- function(x) {
  paste0("{\\bfseries ", x, "}") 
}
Beta_table <- xtable(Beta, 
                     digits = 0, 
                     label = "tab:fixed-parameters",
                     caption = "Values of the fixed-effects basis coefficients used in the simulation.")
align(Beta_table)[1] <- "l"
print(Beta_table, 
      file = file.path(outputs_path, "tables", "Beta-simulation.tex"),
      sanitize.text.function = function(x){x},
      sanitize.colnames.function = bold,
      booktabs = TRUE)
fwrite(Beta, file = file.path(outputs_path, "tables", "Beta-simulation.csv"))


# Extract random-effects parameters: --------------------------------------
## Extract Parameters of Q^*: ---------------------------------------------
subject_id_names <- c("subject_id", "subject_id.1", "subject_id.2")
Q_star <- sapply(lme_simulation_list, function(x) {
  sapply(lme4::VarCorr(x)[subject_id_names], function(y){y[1,1]})[c("subject_id.2", "subject_id.1", "subject_id")]
})
Q_star <- round(Q_star)
colnames(Q_star) <- paste(seq_len(k_sim))
rownames(Q_star) <- c("$q_{11, k}",
                       "$q_{22, k}",
                       "$q_{33, k}")
fwrite(Q_star, file = file.path(outputs_path, "tables", "Q-star-simulation.csv"))


## Extract Parameters of R^*: ---------------------------------------------

subject_id_names <- c("subject_id.side", "subject_id.side.1", "subject_id.side.2")
R_star <- sapply(lme_simulation_list, function(x) {
  sapply(lme4::VarCorr(x)[subject_id_names], function(y){y[1,1]})[c("subject_id.side.2", "subject_id.side.1", "subject_id.side")]
})
R_star <- round(R_star)
colnames(R_star) <- paste(seq_len(k_sim))
rownames(R_star) <- c("$r_{11, k}",
                      "$r_{22, k}",
                      "$r_{33, k}")
fwrite(R_star, file = file.path(outputs_path, "tables", "R-star-simulation.csv"))

## Extract s_k: -----------------------------------------------------------
s_k_vec <- sapply(lme_simulation_list, function(x) {
  attr(lme4::VarCorr(x), "sc")^2
})
stopifnot(length(s_k_vec) == k_sim)
names(s_k_vec) <- paste(seq_len(k_sim))
sk_vec <- round(s_k_vec)

fwrite(data.table(
  k = seq_len(k_sim),
  s_l = sk_vec),
  file = file.path(outputs_path, "tables", "s-simulation.csv"))


## Create combined table for the paper: -----------------------------------
random_parameters_mat <- rbind(
  Q_star,
  R_star,
  s_k_vec
)
rownames(random_parameters_mat)[7] <- c("$s_k$")

random_parameters_table <- xtable(
  random_parameters_mat,
  digits = 0,
  label = "tab:random-parameters",
  caption = "Values of the random-effects and random-error parameters used in the simulation.")
align(random_parameters_table)[1] <- "l"
print(random_parameters_table,
      file = file.path(outputs_path, "tables", "random_parameters-simulation.tex"),
      sanitize.text.function = function(x){x},
      sanitize.colnames.function = bold,
      booktabs = TRUE)
















# Finally -- Make plot of the ten mv-FPCs ---------------------------------
plot_psi_dt <- data.table(
  t_grid = rep(grid_points, times = 3),
  dimension = factor(rep(c("Hip", "Knee", "Ankle"), each = length(grid_points)), levels = c("Hip", "Knee", "Ankle")),
  rbind(Psi_array[,,1], Psi_array[,,2], Psi_array[,,3])
)
names(plot_psi_dt)[- c(1:2)] <- paste0("$\\boldsymbol{\\Psi}_{", seq_len(k_sim), "}(t)$")

plot_psi_dt_lng <- melt.data.table(data = plot_psi_dt,
                                   id.vars = c("t_grid", "dimension"),
                                    measure.vars = paste0("$\\boldsymbol{\\Psi}_{", seq_len(k_sim), "}(t)$"),
                                    variable.name = "k",
                                    value.name = "psi_t",
                                    variable.factor = FALSE,
                                    value.factor = FALSE, 
                                    verbose = TRUE)
plot_psi_dt_lng[, k := factor(k, levels =  paste0("$\\boldsymbol{\\Psi}_{", seq_len(k_sim), "}(t)$"))]


palette <- colorRampPalette(RColorBrewer::brewer.pal(8,name = 'Set1'))(10)
Psi_plot <- ggplot(data = plot_psi_dt_lng) +
  aes(x = t_grid, y = psi_t, color = k, group = k) +
  facet_wrap(~ dimension) +
  geom_line() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "$\\boldsymbol{\\Psi}(t)$") +
  scale_color_manual(values = palette)


tikz(file.path(plots_path, "Psi-simulation-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 0.5 * doc_width_inches)
print(Psi_plot)
dev.off()

