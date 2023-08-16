# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice) # CRAN v0.12.3.1

# Helper functions for looking at results: --------------------------------
source(here::here("code", "functions", "simulation_results_functions.R"))

# Graphics Settings: ------------------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning()
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))
theme_update(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"))

# Path: -------------------------------------------------------------------
simulation_results_path <- "/Users/edwardgunning/Dropbox/simulation"
plots_path <- here::here("outputs", "figures")
outputs_path <- here::here("outputs")


# Varying Long-strength ---------------------------------------------------
N_sub_280_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_1.rds"
  )
)

N_sub_280_prop_missing_0.1_long_strength_2 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_2.rds"
  )
)

N_sub_280_prop_missing_0.1_long_strength_3 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_3.rds"
  )
)

# -------------------------------------------------------------------------
stopifnot(any(!sapply(N_sub_280_prop_missing_0.1_long_strength_1[1:500], is.null)))
stopifnot(any(!sapply(N_sub_280_prop_missing_0.1_long_strength_2[1:500], is.null)))
stopifnot(any(!sapply(N_sub_280_prop_missing_0.1_long_strength_3[1:500], is.null)))


N_sub_280_prop_missing_0.1_long_strength_1$n_cores
N_sub_280_prop_missing_0.1_long_strength_2$n_cores
N_sub_280_prop_missing_0.1_long_strength_3$n_cores

results_1 <- N_sub_280_prop_missing_0.1_long_strength_1[1:500]
results_2 <- N_sub_280_prop_missing_0.1_long_strength_2[1:500]
results_3 <- N_sub_280_prop_missing_0.1_long_strength_3[1:500]
# -------------------------------------------------------------------------

# Check for warnings: -----------------------------------------------------
sum(unlist(lapply(results_1, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))

sum(unlist(lapply(results_2, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))

sum(unlist(lapply(results_3, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))


# Check FPCA Times: -------------------------------------------------------
boxplot(sapply(results_1, function(x) {
  x$pca$time
}),
sapply(results_2, function(x) {
  x$pca$time
}),
sapply(results_3, function(x) {
  x$pca$time
}))

pca_times_dt <- data.table(rep = c(rep(1:500, times = 2), 1:500), 
                           long_strength = rep(c(1:3), each = 500),
                           time_secs = c(
                             sapply(results_1, function(x) {x$pca$time}),
                             sapply(results_2, function(x) {x$pca$time}),
                             sapply(results_3, function(x) {x$pca$time})
                           ))

pca_time_plot <- ggplot(data = pca_times_dt) +
  aes(x = long_strength, group = long_strength, y = time_secs) + 
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "",
       colour = "Model:",
       y = "Time (seconds)",
       title = "(a) Computation Time (mv-FPCA)") +
  scale_x_continuous(breaks = 1:3)

pca_time_plot

# Check Modelling Times: --------------------------------------------------
model_times_1 <- extract_simulation_element(results_1, extractor_function = extract_time)
model_times_1$rep <- 1:500
model_times_1$long_strength <- 1

model_times_2 <- extract_simulation_element(results_2, extractor_function = extract_time)
model_times_2$rep <- 1:500
model_times_2$long_strength <- 2

model_times_3 <- extract_simulation_element(results_3, extractor_function = extract_time)
model_times_3$rep <- 1:500
model_times_3$long_strength <- 3


model_times_dt <- data.table(rbind(model_times_1, model_times_2, model_times_3))
model_times_dt_lng <- melt.data.table(model_times_dt, 
                                      id.vars = c("long_strength", "rep"),
                                      variable.factor = FALSE,
                                      value.factor = FALSE,
                                      variable.name = "model", 
                                      value.name = "time_sec")
model_times_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]

model_fit_time_plot <- ggplot(data = model_times_dt_lng) +
  aes(x = long_strength, group = interaction(long_strength, model), colour = model, y = time_sec) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x ="Longitudinal Variation Strength",
       colour = "Model:",
       y = "Time (seconds)",
       title = "(b) Computation Time (Model Fit)") +
  scale_x_continuous(breaks = 1:3)

model_fit_time_plot


# -------------------------------------------------------------------------

## Prediction Error: ------------------------------------------------------
par(mfrow = c(1, 3))
boxplot(extract_simulation_element(results_1, extractor_function = extract_prediction_error))
boxplot(extract_simulation_element(results_2, extractor_function = extract_prediction_error))
boxplot(extract_simulation_element(results_3, extractor_function = extract_prediction_error))


pred_error_1 <- extract_simulation_element(results_1, extractor_function = extract_prediction_error)
pred_error_1$rep <- 1:500
pred_error_1$long_strength <- 1

pred_error_2 <- extract_simulation_element(results_2, extractor_function = extract_prediction_error)
pred_error_2$rep <- 1:500
pred_error_2$long_strength <- 2

pred_error_3 <- extract_simulation_element(results_3, extractor_function = extract_prediction_error)
pred_error_3$rep <- 1:500
pred_error_3$long_strength <- 3

pred_error_dt <- data.table(rbind(pred_error_1, pred_error_2, pred_error_3))
pred_error_dt_lng <- melt.data.table(pred_error_dt, 
                                     id.vars = c("long_strength", "rep"),
                                     variable.factor = FALSE,
                                     value.factor = FALSE,
                                     variable.name = "model", 
                                     value.name = "ISPE")
pred_error_dt_lng[, model := factor(model, 
                                    levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                    labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]

pred_error_ISPE_plot <- ggplot(data = pred_error_dt_lng) +
  aes(x = long_strength, group = interaction(long_strength, model), colour = model, y = ISPE) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "",
       colour = "Model:",
       y = "ISPE",
       title = "(c) Individual Prediction Error") +
  scale_x_continuous(breaks = 1:3)

pred_error_ISPE_plot 

## Intercept: --------------------------------------------------------------

par(mfrow = c(3, 3))
boxplot(extract_simulation_element(results_1, extractor_function = extract_intercept_error))
boxplot(extract_simulation_element(results_2, extractor_function = extract_intercept_error))
boxplot(extract_simulation_element(results_3, extractor_function = extract_intercept_error))

fixef_error_intercept_1 <- extract_simulation_element(results_1, extractor_function = extract_intercept_error)
fixef_error_intercept_1$rep <- 1:500
fixef_error_intercept_1$long_strength <- 1

fixef_error_intercept_2 <- extract_simulation_element(results_2, extractor_function = extract_intercept_error)
fixef_error_intercept_2$rep <- 1:500
fixef_error_intercept_2$long_strength <- 2

fixef_error_intercept_3 <- extract_simulation_element(results_3, extractor_function = extract_intercept_error)
fixef_error_intercept_3$rep <- 1:500
fixef_error_intercept_3$long_strength <- 3

fixef_error_intercept_dt <- data.table(rbind(fixef_error_intercept_1,
                                             fixef_error_intercept_2,
                                             fixef_error_intercept_3))
fixef_error_intercept_dt_lng <- melt.data.table(fixef_error_intercept_dt, 
                                                id.vars = c("long_strength", "rep"),
                                                variable.factor = FALSE,
                                                value.factor = FALSE,
                                                variable.name = "model", 
                                                value.name = "intercept")

## Speed: ----------------------------------------------------------------

boxplot(extract_simulation_element(results_1, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_2, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_3, extractor_function = extract_speed_error))


fixef_error_speed_1 <- extract_simulation_element(results_1, extractor_function = extract_speed_error)
fixef_error_speed_1$rep <- 1:500
fixef_error_speed_1$long_strength <- 1

fixef_error_speed_2 <- extract_simulation_element(results_2, extractor_function = extract_speed_error)
fixef_error_speed_2$rep <- 1:500
fixef_error_speed_2$long_strength <- 2

fixef_error_speed_3 <- extract_simulation_element(results_3, extractor_function = extract_speed_error)
fixef_error_speed_3$rep <- 1:500
fixef_error_speed_3$long_strength <- 3

fixef_error_speed_dt <- data.table(rbind(fixef_error_speed_1, fixef_error_speed_2, fixef_error_speed_3))
fixef_error_speed_dt_lng <- melt.data.table(fixef_error_speed_dt, 
                                            id.vars = c("long_strength", "rep"),
                                            variable.factor = FALSE,
                                            value.factor = FALSE,
                                            variable.name = "model", 
                                            value.name = "speed")


## Sex: --------------------------------------------------------------------

boxplot(extract_simulation_element(results_1, extractor_function = extract_sex_error))
boxplot(extract_simulation_element(results_2, extractor_function = extract_sex_error))
boxplot(extract_simulation_element(results_3, extractor_function = extract_sex_error))

fixef_error_sex_1 <- extract_simulation_element(results_1, extractor_function = extract_sex_error)
fixef_error_sex_1$rep <- 1:500
fixef_error_sex_1$long_strength <- 1

fixef_error_sex_2 <- extract_simulation_element(results_2, extractor_function = extract_sex_error)
fixef_error_sex_2$rep <- 1:500
fixef_error_sex_2$long_strength <- 2

fixef_error_sex_3 <- extract_simulation_element(results_3, extractor_function = extract_sex_error)
fixef_error_sex_3$rep <- 1:500
fixef_error_sex_3$long_strength <- 3

fixef_error_sex_dt <- data.table(rbind(fixef_error_sex_1, fixef_error_sex_2, fixef_error_sex_3))
fixef_error_sex_dt_lng <- melt.data.table(fixef_error_sex_dt, 
                                          id.vars = c("long_strength", "rep"),
                                          variable.factor = FALSE,
                                          value.factor = FALSE,
                                          variable.name = "model", 
                                          value.name = "sex")


fixef_error_dt_wide <- merge.data.table(merge.data.table(fixef_error_intercept_dt_lng,
                                                         fixef_error_speed_dt_lng,
                                                         by = c("long_strength", "rep", "model"),
                                                         all = TRUE), 
                                        y = fixef_error_sex_dt_lng,
                                        by = c("long_strength", "rep", "model"),
                                        all = TRUE)

fixef_error_dt_lng <- melt.data.table(fixef_error_dt_wide, 
                                      id.vars = c("long_strength", "rep", "model"),
                                      variable.name = "parameter",
                                      measure.vars = c("intercept", "speed", "sex"), 
                                      value.name = "MISE", 
                                      variable.factor = FALSE, 
                                      value.factor = FALSE)

fixef_error_dt_lng[, parameter := factor(parameter,
                                         levels = c("intercept", "speed", "sex"),
                                         labels = c("Intercept $\\boldsymbol{\\beta}_0 (t, T)$",
                                                    "Speed $\\boldsymbol{\\beta}_1 (t)$",
                                                    "Sex $\\boldsymbol{\\beta}_2 (t)$"))]

fixef_error_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]


fixef_error_plot <- ggplot(data = fixef_error_dt_lng) +
  facet_wrap(~ parameter, nrow = 1, ncol = 3, scales = "free_y") +
  aes(x = long_strength, group = interaction(long_strength, model),
      colour = model,
      y = sqrt(MISE)) +
  geom_boxplot(outlier.size = 0.5)  +
  labs(x = "Longitudinal Variation Strength",
       colour = "Model:",
       y = "$\\sqrt{\\text{MISE}}$",
       title = "(d) Fixed Effects Estimation Error") +
  scale_x_continuous(breaks = 1:3)


# Combine: ----------------------------------------------------------------
top_panel <- ggarrange(pca_time_plot,
                       model_fit_time_plot,
                       pred_error_ISPE_plot, 
                       nrow = 1, 
                       legend = "none")

tikz(file.path(plots_path, "simulation-long-strength-plot.tex"),
     width = 1.2 * doc_width_inches, 
     height = 1.3 * (1/1.618) * doc_width_inches, # near golden ratio...
     standAlone = TRUE)
ggarrange(top_panel, fixef_error_plot,
          nrow = 2, 
          common.legend = TRUE,
          legend = "bottom",
          heights = c(0.48, 0.52))
dev.off()

tinytex::lualatex(file.path(plots_path, "simulation-long-strength-plot.tex"))
