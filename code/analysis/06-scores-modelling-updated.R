# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
source(here::here("code/functions/source_all_analysis_functions.R"))

# Path to save the outputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
outputs_path <- here::here("outputs")

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


# Modelling ---------------------------------------------------------------
set.seed(22081996)
(sample_subjects <- sample(unique(covariates_and_scores_dt_train$subject_id), size = 6))
# [1] P_4197 P_4158 P_4287 P_4288 P_4212 P_4100
sample_plot_dt <- covariates_and_scores_dt_train[subject_id %in% sample_subjects]

sample_plot_dt[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
sample_plot_dt[, participant_id := paste("Participant", participant_id)]
sample_plot_dt[, side := factor(side,
                                levels = c("left", "right"),
                                labels = c("Left", "Right"))]


(p <- ggplot(data = sample_plot_dt) +
    aes(x = long_time, y = score_1, group = side, color = side) +
    facet_wrap(~ participant_id, nrow = 2, ncol = 3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.75) +
    scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
    xlab("Longitudinal Time ($T$)") +
    ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
    theme(legend.position = "bottom"))

tikz(file.path(plots_path, "mv-FPC1-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.725 * doc_width_inches)
print(p)
dev.off()




# Modelling: --------------------------------------------------------------
# And join back up into data.table for modelling:
covariates_and_scores_dt_train <- cbind(covariates_dt_train, scores_train)
splines <- ns(covariates_and_scores_dt_train$long_time, df = 5)
spline_names <- colnames(splines) <- paste0("spline_", 1:5)
covariates_and_scores_dt_train <- cbind(covariates_and_scores_dt_train, splines)

covariates_and_scores_dt_train[, time := long_time]

fixef_formula_full <- "ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent"



system.time(poly_2_model <- fit_poly(df_scores = covariates_and_scores_dt_train,
                                     degree = 2,
                                     K_retain = k_retain,
                                     fixef_formula = fixef_formula_full,
                                     diagonal_covariance = FALSE))

system.time(fpca_model <- fit_fpca(df_scores = covariates_and_scores_dt_train, 
                       K_retain = k_retain, 
                       pve = 0.995, 
                       df_spline = 4, 
                       fixef_formula = fixef_formula_full))



system.time(naive_model <- fit_naive(df_scores = covariates_and_scores_dt_train, 
                                      K_retain = k_retain, degree = 2, fixef_formula = fixef_formula_full))



par(mfrow = c(3, 3))
for(j in seq_len(9)) {
  plot(fixef(fpca_model$lme_fit_list[[j]])[-c(1:5)], fixef(poly_2_model$lme_fit_list[[j]])[-c(1:3)])
  abline(0, 1)
}

par(mfrow = c(3, 3))
for(j in seq_len(9)) {
  plot(fixef(fpca_model$lme_fit_list[[j]])[-c(1:5)], fixef(naive_model$lme_fit_list[[j]])[-c(1:3)])
  abline(0, 1)
}


fpca_fixef <- compute_fixef_functions(fit_object = fpca_model, pca_fd_obj = mfpca)
poly_2_fixef <- compute_fixef_functions(fit_object = poly_2_model, pca_fd_obj = mfpca)
naive_fixef <- compute_fixef_functions(fit_object = naive_model, pca_fd_obj = mfpca)


# -------------------------------------------------------------------------

mfd_obj_test <- basis_transformation_results$mfd_obj_test
covariates_dt_test <- basis_transformation_results$covariates_dt_test
covariates_dt_test[, time := long_time]


predicted_fd_obj_fpca <- predict_fd_fpca(fit_fpca_object  = fpca_model,
                                              newdata = covariates_dt_test,
                                              pca_fd_obj = mfpca)
calculate_prediction_error(pred_fd_obj = predicted_fd_obj_fpca, 
                           true_fd_obj = mfd_obj_test)

predicted_fd_obj_naive <- predict_fd_naive(fit_naive_object  = naive_model,
                                         newdata = covariates_dt_test,
                                         pca_fd_obj = mfpca)
calculate_prediction_error(pred_fd_obj = predicted_fd_obj_naive, 
                           true_fd_obj = mfd_obj_test)


predicted_fd_obj_poly_2 <- predict_fd_poly(fit_poly_object = poly_2_model,
                                           newdata = covariates_dt_test,
                                           pca_fd_obj = mfpca)

calculate_prediction_error(pred_fd_obj = predicted_fd_obj_poly_2, 
                           true_fd_obj = mfd_obj_test)



# -------------------------------------------------------------------------

predicted_scores_poly <- predict_scores_poly(fit_poly_object = poly_2_model, 
                    newdata = covariates_and_scores_dt_train)

plot_dt_fit <- data.table(cbind(covariates_and_scores_dt_train, predicted_scores_poly))
sample_plot_dt_02 <- plot_dt_fit[subject_id %in% sample_subjects]
sample_plot_dt_02[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
sample_plot_dt_02[, participant_id := paste("Participant", participant_id)]
sample_plot_dt_02[, side := factor(side,
                                levels = c("left", "right"),
                                labels = c("Left", "Right"))]

(p2 <- ggplot(data = sample_plot_dt_02) +
    aes(x = long_time, y = score_1, group = side, color = side) +
    facet_wrap(~ participant_id, nrow = 2, ncol = 3) +
    geom_line(linewidth = 0.5, alpha = 0.3) +
    geom_point(size = 0.75, alpha = 0.3) +
    scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
    xlab("Longitudinal Time ($T$)") +
    ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
    theme(legend.position = "bottom")) +
    geom_line(aes(y = pred_score_1))

predicted_scores_poly <- predict_scores_poly(fit_poly_object = poly_2_model, 
                                             newdata = covariates_and_scores_dt_train)

plot_dt_fit <- data.table(cbind(covariates_and_scores_dt_train, predicted_scores_poly))
sample_plot_dt_02 <- plot_dt_fit[subject_id %in% sample_subjects]
sample_plot_dt_02[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
sample_plot_dt_02[, participant_id := paste("Participant", participant_id)]
sample_plot_dt_02[, side := factor(side,
                                   levels = c("left", "right"),
                                   labels = c("Left", "Right"))]


predicted_scores_fpca <- predict_scores_fpca(fit_fpca_object  = fpca_model, 
                                             newdata = covariates_and_scores_dt_train)

plot_dt_fit_02 <- data.table(cbind(covariates_and_scores_dt_train, predicted_scores_fpca))
sample_plot_dt_03 <- plot_dt_fit_02[subject_id %in% sample_subjects]
sample_plot_dt_03[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
sample_plot_dt_03[, participant_id := paste("Participant", participant_id)]
sample_plot_dt_03[, side := factor(side,
                                   levels = c("left", "right"),
                                   labels = c("Left", "Right"))]




(p3 <- ggplot(data = sample_plot_dt_03) +
    aes(x = long_time, y = score_1, group = side, color = side) +
    facet_wrap(~ participant_id, nrow = 2, ncol = 3) +
    geom_line(linewidth = 0.5, alpha = 0.25) +
    geom_point(size = 0.75, alpha = 0.25) +
    scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
    xlab("Longitudinal Time ($T$)") +
    ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
    theme(legend.position = "bottom") +
  geom_line(aes(y = pred_score_1, linetype = "ml-FPCA")) +
  geom_line(aes(y = pred_score_1, linetype = "Polynomial"), data = sample_plot_dt_02) +
  scale_linetype_discrete(name = NULL))
  

tikz(file.path(plots_path, "fitted-mv-FPC1-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.725 * doc_width_inches,
     standAlone = TRUE)
print(p3)
dev.off()

tinytex::lualatex(file.path(plots_path, "fitted-mv-FPC1-plot.tex"))






