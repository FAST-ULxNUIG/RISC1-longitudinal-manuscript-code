# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
library(data.table) # CRAN v1.14.2 
library(ggpubr)     # CRAN v0.4.0

# Load all helper functions: ----------------------------------------------
source(here::here("code", "functions", "source_all_analysis_functions.R"))



# Graphics: ---------------------------------------------------------------
# Set graphics settings; --------------------------------------------------
theme_gunning()
plots_path <- here::here("outputs", "figures")
theme_update(
  plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"),
  axis.text = element_text(size = 9))
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937


# Load Results: -----------------------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
model_fit_results <- readRDS(
  here::here("outputs", "model-fit-results.rds"))

# Extract Objects from Results: -------------------------------------------
# From basis transformation results:
mfpca <- basis_transformation_results$mfpca
k_retain <- basis_transformation_results$k_retain
covariates_dt_test <- basis_transformation_results$covariates_dt_test
covariates_dt_test[, time := long_time] # needs this name for prediction.
mfd_obj_test <- basis_transformation_results$mfd_obj_test
# From model fitting results:
fpca_model <- model_fit_results$fpca_model
naive_model <- model_fit_results$naive_model
spline_ri_model <- model_fit_results$spline_ri_model
covariates_and_scores_dt_train <- model_fit_results$covariates_and_scores_dt_train
# Create data to look at change: ------------------------------------------
# Want to create data to evaluate every subject at on a longitudinal grid.
# So create new dataset to predict at long. times 0, 0.01, ..., 0.99, 1.
covariates_dt_train <- basis_transformation_results$covariates_dt_train
newdata_for_change <- unique(covariates_dt_train[, .(subject_id,
                                                     side, 
                                                     ris, 
                                                     sex, 
                                                     speed_cent,
                                                     age_cent,
                                                     height_cm_cent,
                                                     weight_kg_cent)])
newdata_for_change <- newdata_for_change[,
                                         .(time = seq(0, 1, length.out = 101)),
                                         by = .(subject_id,
                                                side, 
                                                ris, 
                                                sex, 
                                                speed_cent,
                                                age_cent,
                                                height_cm_cent,
                                                weight_kg_cent)]

# Calculate subject-and-side specific predictions from each modek
newdata_for_change[, paste0("fpca_prediction_for_change_", seq_len(k_retain)) := {
  get_list_of_columns(predict_scores_fpca(fit_fpca_object = fpca_model, newdata = newdata_for_change))
  }]

newdata_for_change[, paste0("spline_ri_prediction_for_change_", seq_len(k_retain)) := {
  get_list_of_columns(predict_scores_spline(fit_spline_object = spline_ri_model, newdata = newdata_for_change))
  }]

newdata_for_change[, paste0("naive_prediction_for_change_", seq_len(k_retain)) := {
  get_list_of_columns(predict_scores_naive_spline_intercept(fit_naive_spline_intercept_object = naive_model, newdata = newdata_for_change))
  }]


# Comute Rates of Change: -------------------------------------------------
# uses custiom function to compute derivative of fitted trajectory 
# this could, somehow, be done in parallell ()
spline_ri_rate_of_change <- newdata_for_change[,
                   {
                     predictions_mat <- as.matrix(.SD)
                     as.list(apply(predictions_mat, MARGIN = 2, calculate_rate_of_change, arg_vals = seq(0, 1, length.out = 101)))
                     },
                   .SDcols = paste0("spline_ri_prediction_for_change_", seq_len(k_retain)),
                   by = .(subject_id, side)]

fpca_rate_of_change <- newdata_for_change[,
                   {
                     predictions_mat <- as.matrix(.SD)
                     as.list(apply(predictions_mat, MARGIN = 2, calculate_rate_of_change, arg_vals = seq(0, 1, length.out = 101)))
                   },
                   .SDcols = paste0("fpca_prediction_for_change_", seq_len(k_retain)),
                   by = .(subject_id, side)]

rate_of_change_dt <- merge.data.table(x = spline_ri_rate_of_change,
                                      y = fpca_rate_of_change, 
                                      by = c("subject_id", "side"),
                                      all = TRUE)
stopifnot(nrow(rate_of_change_dt) == nrow(fpca_rate_of_change) & nrow(rate_of_change_dt) == nrow(spline_ri_rate_of_change))


# Compute Overall Totals:
rate_of_change_dt[, total_rate_of_change_fpca := {
  apply(.SD, 1, sum)
  },
  .(subject_id, side),
  .SDcols = paste0("fpca_prediction_for_change_", seq_len(k_retain))]
rate_of_change_dt[, total_rate_of_change_spline_ri := {
  apply(.SD, 1, sum)
},
.(subject_id, side),
.SDcols = paste0("spline_ri_prediction_for_change_", seq_len(k_retain))]


# Now look at overall change: ---------------------------------------------

overall_change_dt <- newdata_for_change[, {
  ind_0 <- which(time == 0)
  ind_1 <- which(time == 1)
  (.SD[ind_1] - .SD[ind_0])^2},
  .SDcols = c(paste0("fpca_prediction_for_change_", seq_len(k_retain)),
              paste0("spline_ri_prediction_for_change_", seq_len(k_retain))),
  by = .(subject_id, side)]


overall_change_dt[, total_overall_change_spline_ri := {
  apply(.SD, 1, sum)
  },
  by = .(subject_id, side),
  .SDcols = paste0("spline_ri_prediction_for_change_", seq_len(k_retain))]

overall_change_dt[, total_overall_change_fpca := {
  apply(.SD, 1, sum)
},
by = .(subject_id, side),
.SDcols = paste0("fpca_prediction_for_change_", seq_len(k_retain))]






# Visualisation: ----------------------------------------------------------
## Rate of Change: --------------------------------------------------------
# Visualisation: ----------------------------------------------------------
rate_of_change_subjects <- rate_of_change_dt[rev(order(total_rate_of_change_spline_ri)),
                                             as.character(unique(subject_id)[1:4])]
rate_of_change_scores_plot_dt <- covariates_and_scores_dt_train[subject_id %in% rate_of_change_subjects & side == "left"]
rate_of_change_fitted_plot_dt <- newdata_for_change[subject_id %in% rate_of_change_subjects & side == "left"]

rate_of_change_scores_plot_dt[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
rate_of_change_scores_plot_dt[, participant_id := paste("Participant", participant_id)]
rate_of_change_fitted_plot_dt[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
rate_of_change_fitted_plot_dt[, participant_id := paste("Participant", participant_id)]

# re-labael the subjects to order them by rate of change:
rate_of_change_partcipants <- paste("Participant", as.numeric(stringr::str_extract(rate_of_change_subjects, "4\\d+")) - 4000)
# and re-order:
rate_of_change_scores_plot_dt[, participant_id := factor(participant_id, 
                                                         levels = rate_of_change_partcipants)]
rate_of_change_fitted_plot_dt[, participant_id := factor(participant_id, 
                                                         levels = rate_of_change_partcipants)]


p1 <- ggplot(data = rate_of_change_scores_plot_dt) +
  aes(x = time) +
  facet_wrap(~ participant_id, ncol = 4) +
  geom_point(aes(y = score_1), alpha = 0.2, size = 0.75) +
  geom_line(aes(y = score_1), alpha = 0.2)  +
  geom_line(data = rate_of_change_fitted_plot_dt, aes(y = fpca_prediction_for_change_1, colour = "ml-FPCA")) +
  geom_line(data = rate_of_change_fitted_plot_dt, aes(y = naive_prediction_for_change_1, colour = "Naive")) +
  geom_line(data = rate_of_change_fitted_plot_dt, aes(y = spline_ri_prediction_for_change_1, colour = "Spline")) +
  xlab("Longitudinal Time ($T$)") +
  ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
  scale_color_discrete(name = "Model:") +
  ggtitle("(a) Rate of Change") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste0(breaks = c(0, 0.25, 0.5, 0.75, 1))) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))
p1



overall_change_subjects <- overall_change_dt[rev(order(total_overall_change_spline_ri)), as.character(unique(subject_id)[1:6])]
overall_change_subjects <- overall_change_subjects[!(overall_change_subjects %in% rate_of_change_subjects)]
overall_change_scores_plot_dt <- covariates_and_scores_dt_train[subject_id %in% overall_change_subjects & side == "left"]

overall_change_fitted_plot_dt <- newdata_for_change[subject_id %in% overall_change_subjects & side == "left"]
overall_change_scores_plot_dt[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
overall_change_scores_plot_dt[, participant_id := paste("Participant", participant_id)]
overall_change_fitted_plot_dt[, participant_id := as.numeric(stringr::str_extract(subject_id, "4\\d+")) - 4000]
overall_change_fitted_plot_dt[, participant_id := paste("Participant", participant_id)]

# re-labael the subjects to order them by rate of change:
overall_change_partcipants <- paste("Participant", as.numeric(stringr::str_extract(overall_change_subjects, "4\\d+")) - 4000)
# and re-order:
overall_change_scores_plot_dt[, participant_id := factor(participant_id, 
                                                         levels = overall_change_partcipants)]
overall_change_fitted_plot_dt[, participant_id := factor(participant_id, 
                                                         levels = overall_change_partcipants)]

p2 <- ggplot(data = overall_change_scores_plot_dt) +
  aes(x = time) +
  facet_wrap(~ participant_id, ncol = 4) +
  geom_point(aes(y = score_1), alpha = 0.2, size = 0.75) +
  geom_line(aes(y = score_1), alpha = 0.2)  +
  geom_line(data = overall_change_fitted_plot_dt, aes(y = fpca_prediction_for_change_1, colour = "ml-FPCA")) +
  geom_line(data = overall_change_fitted_plot_dt, aes(y = naive_prediction_for_change_1, colour = "Naive")) +
  geom_line(data = overall_change_fitted_plot_dt, aes(y = spline_ri_prediction_for_change_1, colour = "Spline")) +
  xlab("Longitudinal Time ($T$)") +
  ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
  scale_color_discrete(name = "Model:") +
  ggtitle("(b) Overall Change") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste0(breaks = c(0, 0.25, 0.5, 0.75, 1))) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))
p2

combined_plot <- ggarrange(p1,
                           p2,
                           nrow = 2, 
                           common.legend = TRUE, legend = "bottom"
                           )
combined_plot

tikz(file.path(plots_path, "change-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.75 * doc_width_inches, 
     standAlone = TRUE)
print(combined_plot)
dev.off()

tinytex::lualatex(file.path(plots_path, "change-plot.tex"))
