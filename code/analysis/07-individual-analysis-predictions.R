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
spline_ri_model <- model_fit_results$spline_ri_model



P_4237_inds <- covariates_dt_test[, which(subject_id=="P_4237" & side == "right")]

P_4237_eval_array <- eval.fd(evalarg = 0:100, fdobj = mfd_obj_test[P_4237_inds])
P_4237_eval_mat <- rbind(
  P_4237_eval_array[,,1],
  P_4237_eval_array[,,2],
  P_4237_eval_array[,,3]
)

P_4237_eval_dt <- data.table(
  t = rep(0:100, times = 3),
  dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
  P_4237_eval_mat
)
names(P_4237_eval_dt)[-c(1:2)] <- covariates_dt_test[P_4237_inds, stride_num]


P_4237_eval_dt_lng <- melt.data.table(data = P_4237_eval_dt,
                                      id.vars = c("t", "dimension"),
                                      measure.vars = paste(covariates_dt_test[P_4237_inds, stride_num]),
                                      variable.name = "stride_num", value.name = "angle",
                                      variable.factor = FALSE, value.factor = FALSE)
P_4237_eval_dt_lng[, `:=`(
  dimension = factor(dimension, levels = c("Hip", "Knee", "Ankle")),
  stride_num = as.numeric(stride_num)
)]




P_4237_predicted_fd <- predict_fd_spline(fit_spline_object = spline_ri_model,
                                       newdata = covariates_dt_test[P_4237_inds], 
                                       pca_fd_obj = mfpca)

P_4237_predicted_eval_array <- eval.fd(evalarg = 0:100,
                             fdobj = P_4237_predicted_fd)
P_4237_predicted_eval_mat <- rbind(
  P_4237_predicted_eval_array[,,1],
  P_4237_predicted_eval_array[,,2],
  P_4237_predicted_eval_array[,,3]
)

P_4237_predicted_eval_dt <- data.table(
  t = rep(0:100, times = 3),
  dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
  P_4237_predicted_eval_mat
)
names(P_4237_predicted_eval_dt)[-c(1:2)] <- covariates_dt_test[P_4237_inds, stride_num]


P_4237_predicted_eval_dt_lng <- melt.data.table(data = P_4237_predicted_eval_dt,
                                      id.vars = c("t", "dimension"),
                                      measure.vars = paste(covariates_dt_test[P_4237_inds, stride_num]),
                                      variable.name = "stride_num", value.name = "angle",
                                      variable.factor = FALSE, value.factor = FALSE)
P_4237_predicted_eval_dt_lng[, `:=`(
  dimension = factor(dimension, levels = c("Hip", "Knee", "Ankle")),
  stride_num = as.numeric(stride_num)
)]


p1 <- ggplot(data = P_4237_eval_dt_lng) +
  aes(x = t, y = angle, group = stride_num, colour= factor(stride_num)) +
  facet_wrap(~ dimension, scales = "free_y") +
  geom_line() +
  ggtitle("(a) Held-Out Strides") +
  labs(y = "Angle ($^{\\circ}$)",
       x = "Normalised Time ($\\%$ of Stride)",
       colour = "Stride Number:") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

p2 <- ggplot(data = P_4237_predicted_eval_dt_lng) +
  aes(x = t, y = angle, group = stride_num, colour= factor(stride_num)) +
  facet_wrap(~ dimension, scales = "free_y") +
  geom_line() +
  ggtitle("(b) Predictions") +
  scale_y_continuous() +
  labs(y = "Angle ($^{\\circ}$)", 
       x = "Normalised Time ($\\%$ of Stride)", 
       colour = "Stride Number:") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

combined_plot <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "bottom") 

tikz(file.path(plots_path, "p4237-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.75 * doc_width_inches, 
     standAlone = TRUE)
print(combined_plot)
dev.off()

tinytex::pdflatex(file.path(plots_path, "p4237-plot.tex"))


