# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice) # CRAN v0.12.3.1

# Graphics Settings: ------------------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
source(here::here("code", "functions", "simulation_results_functions.R"))

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


# Varying Prop. Missing ---------------------------------------------------
N_sub_280_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_1.rds"
  )
)

N_sub_280_prop_missing_0.2_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.2_long_strength_1.rds"
  )
)

N_sub_280_prop_missing_0.5_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.5_long_strength_1.rds"
  )
)

# -------------------------------------------------------------------------
stopifnot(any(!sapply(N_sub_280_prop_missing_0.1_long_strength_1[1:500], is.null)))
stopifnot(any(!sapply(N_sub_280_prop_missing_0.2_long_strength_1[1:500], is.null)))
stopifnot(any(!sapply(N_sub_280_prop_missing_0.5_long_strength_1[1:500], is.null)))


N_sub_280_prop_missing_0.1_long_strength_1$n_cores
N_sub_280_prop_missing_0.2_long_strength_1$n_cores
N_sub_280_prop_missing_0.5_long_strength_1$n_cores

results_0.1 <- N_sub_280_prop_missing_0.1_long_strength_1[1:500]
results_0.2 <- N_sub_280_prop_missing_0.2_long_strength_1[1:500]
results_0.5 <- N_sub_280_prop_missing_0.5_long_strength_1[1:500]
# -------------------------------------------------------------------------

# Check for warnings: -----------------------------------------------------
sum(unlist(lapply(results_0.1, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))

sum(unlist(lapply(results_0.2, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))

sum(unlist(lapply(results_0.5, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))


# Check FPCA Times: -------------------------------------------------------
boxplot(sapply(results_0.1, function(x) {
  x$pca$time
}),
sapply(results_0.2, function(x) {
  x$pca$time
}),
sapply(results_0.5, function(x) {
  x$pca$time
}))

pca_times_dt <- data.table(rep = c(rep(1:500, times = 2), 1:500), 
                           prop_missing = rep(c(0.1, 0.2, 0.5), each = 500),
                           time_secs = c(
                             sapply(results_0.1, function(x) {x$pca$time}),
                             sapply(results_0.2, function(x) {x$pca$time}),
                             sapply(results_0.5, function(x) {x$pca$time})
                           ))

pca_time_plot <- ggplot(data = pca_times_dt) +
  aes(x = prop_missing, group = prop_missing, y = time_secs) + 
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "",
       colour = "Model:",
       y = "Time (seconds)",
       title = "(a) Computation Time (mv-FPCA)") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.5))

pca_time_plot

# Check Modelling Times: --------------------------------------------------

# Had an error here: ------------------------------------------------------
 which(sapply(results_0.2, function(x) {
  is.null(x$models$fpca$time)
}))

which(sapply(results_0.5, function(x) {
  is.null(x$models$fpca$time)
}))


model_times_0.1 <- extract_simulation_element(results_0.1, extractor_function = extract_time)
model_times_0.1$rep <- 1:500
model_times_0.1$prop_missing <- 0.1

model_times_0.2 <- extract_simulation_element(results_0.2, extractor_function = extract_time) # Note: Added because of 1 failing fpca replicate.
model_times_0.2$rep <- 1:500
model_times_0.2$prop_missing <- 0.2

model_times_0.5 <- extract_simulation_element(results_0.5, extractor_function = extract_time)
model_times_0.5$rep <- 1:500
model_times_0.5$prop_missing <- 0.5

model_times_dt <- data.table(rbind(model_times_0.1, model_times_0.2, model_times_0.5))
model_times_dt_lng <- melt.data.table(model_times_dt, 
                                      id.vars = c("prop_missing", "rep"),
                                      variable.factor = FALSE,
                                      value.factor = FALSE,
                                      variable.name = "model", 
                                      value.name = "time_sec")
model_times_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]

model_fit_time_plot <- ggplot(data = model_times_dt_lng) +
  aes(x = prop_missing, group = interaction(prop_missing, model), colour = model, y = time_sec) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x ="Proportion of Missing Strides",
       colour = "Model:",
       y = "Time (seconds)",
       title = "(b) Computation Time (Model Fit)") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.5))

model_fit_time_plot



# -------------------------------------------------------------------------

## Prediction Error: ------------------------------------------------------
par(mfrow = c(1, 3))
boxplot(extract_simulation_element(results_0.1, extractor_function = extract_prediction_error))
boxplot(extract_simulation_element(results_0.2, extractor_function = extract_prediction_error))
boxplot(extract_simulation_element(results_0.5, extractor_function = extract_prediction_error))


pred_error_0.1 <- extract_simulation_element(results_0.1, extractor_function = extract_prediction_error)
pred_error_0.1$rep <- 1:500
pred_error_0.1$prop_missing <- 0.1

pred_error_0.2 <- extract_simulation_element(results_0.2, extractor_function = extract_prediction_error)
pred_error_0.2$rep <- 1:500
pred_error_0.2$prop_missing <- 0.2

pred_error_0.5 <- extract_simulation_element(results_0.5, extractor_function = extract_prediction_error)
pred_error_0.5$rep <- 1:500
pred_error_0.5$prop_missing <- 0.5

pred_error_dt <- data.table(rbind(pred_error_0.1, pred_error_0.2, pred_error_0.5))
pred_error_dt_lng <- melt.data.table(pred_error_dt, 
                                     id.vars = c("prop_missing", "rep"),
                                     variable.factor = FALSE,
                                     value.factor = FALSE,
                                     variable.name = "model", 
                                     value.name = "ISPE")
pred_error_dt_lng[, model := factor(model, 
                                    levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                    labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]

pred_error_ISPE_plot <- ggplot(data = pred_error_dt_lng) +
  aes(x = prop_missing, group = interaction(prop_missing, model), colour = model, y = ISPE) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "",
       colour = "Model:",
       y = "ISPE",
       title = "(c) Individual Prediction Error") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.5))

pred_error_ISPE_plot 

## Intercept: --------------------------------------------------------------

par(mfrow = c(3, 3))
boxplot(extract_simulation_element(results_0.1, extractor_function = extract_intercept_error))
boxplot(extract_simulation_element(results_0.2, extractor_function = extract_intercept_error))
boxplot(extract_simulation_element(results_0.5, extractor_function = extract_intercept_error))

fixef_error_intercept_0.1 <- extract_simulation_element(results_0.1, extractor_function = extract_intercept_error)
fixef_error_intercept_0.1$rep <- 1:500
fixef_error_intercept_0.1$prop_missing <- 0.1

fixef_error_intercept_0.2 <- extract_simulation_element(results_0.2, extractor_function = extract_intercept_error)
fixef_error_intercept_0.2$rep <- 1:500
fixef_error_intercept_0.2$prop_missing <- 0.2

fixef_error_intercept_0.5 <- extract_simulation_element(results_0.5, extractor_function = extract_intercept_error)
fixef_error_intercept_0.5$rep <- 1:500
fixef_error_intercept_0.5$prop_missing <- 0.5

fixef_error_intercept_dt <- data.table(rbind(fixef_error_intercept_0.1, fixef_error_intercept_0.2, fixef_error_intercept_0.5))
fixef_error_intercept_dt_lng <- melt.data.table(fixef_error_intercept_dt, 
                                                id.vars = c("prop_missing", "rep"),
                                                variable.factor = FALSE,
                                                value.factor = FALSE,
                                                variable.name = "model", 
                                                value.name = "intercept")

## Speed: ----------------------------------------------------------------

boxplot(extract_simulation_element(results_0.1, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_0.2, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_0.5, extractor_function = extract_speed_error))

fixef_error_speed_0.1 <- extract_simulation_element(results_0.1, extractor_function = extract_speed_error)
fixef_error_speed_0.1$rep <- 1:500
fixef_error_speed_0.1$prop_missing <- 0.1

fixef_error_speed_0.2 <- extract_simulation_element(results_0.2, extractor_function = extract_speed_error)
fixef_error_speed_0.2$rep <- 1:500
fixef_error_speed_0.2$prop_missing <- 0.2

fixef_error_speed_0.5 <- extract_simulation_element(results_0.5, extractor_function = extract_speed_error)
fixef_error_speed_0.5$rep <- 1:500
fixef_error_speed_0.5$prop_missing <- 0.5

fixef_error_speed_dt <- data.table(rbind(fixef_error_speed_0.1, fixef_error_speed_0.2, fixef_error_speed_0.5))
fixef_error_speed_dt_lng <- melt.data.table(fixef_error_speed_dt, 
                                            id.vars = c("prop_missing", "rep"),
                                            variable.factor = FALSE,
                                            value.factor = FALSE,
                                            variable.name = "model", 
                                            value.name = "speed")


## Sex: --------------------------------------------------------------------
boxplot(extract_simulation_element(results_0.1, extractor_function = extract_sex_error))
boxplot(extract_simulation_element(results_0.2, extractor_function = extract_sex_error))
boxplot(extract_simulation_element(results_0.5, extractor_function = extract_sex_error))

fixef_error_sex_0.1 <- extract_simulation_element(results_0.1, extractor_function = extract_sex_error)
fixef_error_sex_0.1$rep <- 1:500
fixef_error_sex_0.1$prop_missing <- 0.1

fixef_error_sex_0.2 <- extract_simulation_element(results_0.2, extractor_function = extract_sex_error)
fixef_error_sex_0.2$rep <- 1:500
fixef_error_sex_0.2$prop_missing <- 0.2

fixef_error_sex_0.5 <- extract_simulation_element(results_0.5, extractor_function = extract_sex_error)
fixef_error_sex_0.5$rep <- 1:500
fixef_error_sex_0.5$prop_missing <- 0.5

fixef_error_sex_dt <- data.table(rbind(fixef_error_sex_0.1, fixef_error_sex_0.2, fixef_error_sex_0.5))
fixef_error_sex_dt_lng <- melt.data.table(fixef_error_sex_dt, 
                                          id.vars = c("prop_missing", "rep"),
                                          variable.factor = FALSE,
                                          value.factor = FALSE,
                                          variable.name = "model", 
                                          value.name = "sex")


fixef_error_dt_wide <- merge.data.table(merge.data.table(fixef_error_intercept_dt_lng,
                                                         fixef_error_speed_dt_lng,
                                                         by = c("prop_missing", "rep", "model"),
                                                         all = TRUE), 
                                        y = fixef_error_sex_dt_lng,
                                        by = c("prop_missing", "rep", "model"),
                                        all = TRUE)

fixef_error_dt_lng <- melt.data.table(fixef_error_dt_wide, 
                                      id.vars = c("prop_missing", "rep", "model"),
                                      variable.name = "parameter",
                                      measure.vars = c("intercept", "speed", "sex"), 
                                      value.name = "MISE", 
                                      variable.factor = FALSE, 
                                      value.factor = FALSE)

fixef_error_dt_lng[, parameter := factor(parameter,
                                         levels = c("intercept",  "sex", "speed"),
                                         labels = c("Intercept $\\boldsymbol{\\beta}_0 (t, T)$",
                                                    "Sex $\\boldsymbol{\\beta}_1 (t)$",
                                                    "Speed $\\boldsymbol{\\beta}_2 (t)$"))]

fixef_error_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]


fixef_error_plot <- ggplot(data = fixef_error_dt_lng) +
  facet_wrap(~ parameter, nrow = 1, ncol = 3, scales = "free_y") +
  aes(x = prop_missing, group = interaction(prop_missing, model), colour = model, y = sqrt(MISE)) +
  geom_boxplot(outlier.size = 0.5)  +
  labs(x = "Proportion of Missing Strides",
       colour = "Model:",
       y = "$\\sqrt{\\text{ISE}}$",
       title = "(d) Fixed Effects Estimation Error") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.5))

fixef_error_plot



# Combine: ----------------------------------------------------------------

top_panel <- ggarrange(pca_time_plot, model_fit_time_plot, pred_error_ISPE_plot, nrow = 1, legend = "none")

tikz(file.path(plots_path, "simulation-prop-missing-plot.tex"),
     width = 1.2 * doc_width_inches, 
     height = 1.3 * (1/1.618) * doc_width_inches,
     standAlone = TRUE)
ggarrange(top_panel, fixef_error_plot, nrow = 2, common.legend = TRUE, legend = "bottom", heights = c(0.48, 0.52))
dev.off()


tinytex::lualatex(file.path(plots_path, "simulation-prop-missing-plot.tex"))
