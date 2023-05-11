# ------------------------------------------------------------------------#
# Obtain parameters for simulation using an initial simple model fit -----#
# ------------------------------------------------------------------------#

# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30

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


# Create scores: ----------------------------------------------------------
# scores for mv-fpca are in a 
# N_{Total} \times K \times P matrix
# Sum over the hip knee and ankle to get overall scores:
scores_train <- apply(mfpca$scores, c(1, 2), sum)
colnames(scores_train) <- paste0("score_", seq_len(k_retain))

# Some basic checks:
stopifnot(dim(scores_train) == c(N_train, k_retain))
stopifnot(nrow(scores_train) == nrow(covariates_dt_train))

# And join back up into data.table for modelling:
covariates_and_scores_dt_train <- cbind(covariates_dt_train, scores_train)

k_sim <- 10 # Number of basis functions used in the simulation:
grid_points <- 0:100
Psi_array <- eval.fd(evalarg = grid_points, fdobj = mfpca$harmonics)[,seq_len(k_sim),]
mu_array <- eval.fd(evalarg = grid_points,
                    fdobj = mfpca$meanfd)



longitudinal_grid <- seq(0, 1, by = 0.01)
poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)
matplot(poly_basis, type = "l")
poly_longtime <- poly(covariates_and_scores_dt_train$long_time, degree = 2, coefs = attr(poly_basis, "coefs"))
attr(poly_longtime, "coefs")
colnames(poly_longtime) <- paste0("poly_", 1:2)
covariates_and_scores_dt_train <- cbind(covariates_and_scores_dt_train, poly_longtime)



lme_simulation_list <- vector(mode = "list", length = k_sim)
for(k in seq_len(k_sim)) {
  print(paste("Score", k))
  formula_k <- paste0("score_", k, " ~ sex + speed_cent + poly_1 + poly_2 + (poly_1 + poly_2||subject_id) + (poly_1 + poly_2||subject_id:side)")
  formula_k <- formula(formula_k)
  lme_simulation_list[[k]] <- lmer(
    formula = formula_k,
    data = covariates_and_scores_dt_train,
    control= lmerControl(optCtrl = list(xtol_rel=0,
                                                   xtol_abs=1e-10,
                                                   ftol_rel=0, 
                                                   ftol_abs=1e-10)))
}



fixef <- sapply(lme_simulation_list, fixef)
matplot(seq_len(k_sim), t(fixef), type = "b", lty = 1)

plot((fixef[2, ] %*% t(Psi_array[,,1]))[1,], type = "l", ylim = c(-5, 5))
abline(h = 0)



