# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
library(data.table) # CRAN v1.14.2 
library(ggpubr)

# Load all helper functions: ----------------------------------------------
source(here::here("code", "functions", "source_all_analysis_functions.R"))


# Set graphics settings; --------------------------------------------------
theme_gunning()
plots_path <- here::here("outputs", "figures")
theme_gunning() # use default theme
theme_update(plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"))
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

# Calculate Prediction Errors: --------------------------------------------
## Spline RI model: --------------------------------------------------------
predicted_fd_obj_spline_ri <- predict_fd_spline(
  fit_spline_object = spline_ri_model, 
  newdata = covariates_dt_test,
  pca_fd_obj = mfpca
)
calculate_prediction_error(pred_fd_obj = predicted_fd_obj_spline_ri, 
                           true_fd_obj = mfd_obj_test)
individual_test_pe_spline_ri <- calculate_individual_prediction_errors(
  pred_fd_obj = predicted_fd_obj_spline_ri, true_fd_obj = mfd_obj_test)

## FPCA model: ------------------------------------------------------------
predicted_fd_obj_fpca <- predict_fd_fpca(fit_fpca_object = fpca_model,
                                         newdata = covariates_dt_test,
                                         pca_fd_obj = mfpca)
calculate_prediction_error(pred_fd_obj = predicted_fd_obj_fpca, 
                           true_fd_obj = mfd_obj_test)
individual_test_pe_fpca <- calculate_individual_prediction_errors(
  pred_fd_obj = predicted_fd_obj_fpca, true_fd_obj = mfd_obj_test)

## Naive Model: -----------------------------------------------------------
predicted_fd_obj_naive <- predict_fd_naive(fit_naive_object  = naive_model,
                                           newdata = covariates_dt_test,
                                           pca_fd_obj = mfpca)
calculate_prediction_error(pred_fd_obj = predicted_fd_obj_naive, 
                           true_fd_obj = mfd_obj_test)
individual_test_pe_naive <- calculate_individual_prediction_errors(
  pred_fd_obj = predicted_fd_obj_naive, true_fd_obj = mfd_obj_test)




# Inspect test errors: ----------------------------------------------------

boxplot(individual_test_pe_spline_ri / individual_test_pe_naive,
        individual_test_pe_fpca / individual_test_pe_naive)


# Create Publishable Boxplot: ---------------------------------------------
## By Individual Stride: ---------------------------------------------------
individual_pe_dt <- copy(covariates_dt_test)
individual_pe_dt[, `:=`(
  pe_spline_ri = individual_test_pe_spline_ri,
  pe_fpca = individual_test_pe_fpca,
  pe_naive = individual_test_pe_naive)]

individual_pe_dt[,  `:=`(
  pe_spline_ri_ratio = pe_spline_ri / pe_naive,
  pe_fpca_ratio = pe_fpca / pe_naive)]

individual_pe_dt_lng <- melt.data.table(individual_pe_dt[, .(subject_id, stride_num, side, pe_fpca_ratio, pe_spline_ri_ratio)],
                                        id.vars = c("subject_id", "stride_num", "side"), verbose = TRUE,
                                        measure.vars = c("pe_fpca_ratio", "pe_spline_ri_ratio"),
                                        variable.factor = FALSE,
                                        value.factor = FALSE,
                                        variable.name = "method",
                                        value.name = "pe_ratio")

individual_pe_dt_lng[, method := factor(method, 
                                        levels = c("pe_spline_ri_ratio", "pe_fpca_ratio"),
                                        c("Spline Model", "ml-FPCA model"))]

p1 <- ggplot(data = individual_pe_dt_lng) +
  aes(y = pe_ratio, x = method, colour = method) +
  geom_hline(yintercept = c(0.95,1), lty = c(2, 1), col = c("darkgrey", 1)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = "ISPE$_{method}$ / ISPE$_{naive}$", 
       title = "(a) Test-Set Prediction Error",
       subtitle = "Individual Strides",
       x = "Method") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#CC79A7", "#009E73")) 


## by Subject -------------------------------------------------------------
individual_pe_dt_subject <- individual_pe_dt[, .(pe_spline_ri = mean(pe_spline_ri),
                                                 pe_fpca = mean(pe_fpca),
                                                 pe_naive = mean(pe_naive)),
                                             by = subject_id]
individual_pe_dt_subject[,  `:=`(
  pe_spline_ri_ratio = pe_spline_ri / pe_naive,
  pe_fpca_ratio = pe_fpca / pe_naive)]

individual_pe_dt_subject_lng <- melt.data.table(individual_pe_dt_subject[, .(subject_id, pe_fpca_ratio, pe_spline_ri_ratio)],
                                        id.vars = c("subject_id"), verbose = TRUE,
                                        measure.vars = c("pe_fpca_ratio", "pe_spline_ri_ratio"),
                                        variable.factor = FALSE,
                                        value.factor = FALSE,
                                        variable.name = "method",
                                        value.name = "pe_ratio")

individual_pe_dt_subject_lng[, method := factor(method, 
                                        levels = c("pe_spline_ri_ratio", "pe_fpca_ratio"),
                                        c("Spline Model", "ml-FPCA model"))]


p2 <- ggplot(data = individual_pe_dt_subject_lng) +
  aes(y = pe_ratio, x = method, colour = method) +
  geom_hline(yintercept = c(0.95,1), lty = c(2, 1), col = c("darkgrey", 1)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = "Av(ISPE$_{method}$) / Av(ISPE$_{naive}$)", 
       title = "(b) Test-Set Prediction Error",
       subtitle = "Subject Averages",
       x = "Method")  +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#CC79A7", "#009E73")) 



## Combine: ----------------------------------------------------------------
ggarrange(p1, p2)
tikz(file.path(plots_path, "test-set-PE.tex"),
     width = (1) * doc_width_inches, 
     height = 0.5 * doc_width_inches,
     standAlone = TRUE)
ggarrange(p1, p2)
dev.off()

tinytex::lualatex(file.path(plots_path, "test-set-PE.tex"))



# Save Output: ------------------------------------------------------------
saveRDS(object = individual_pe_dt, 
        file = here::here("outputs", "individual_prediction_error_dt.rds"))

