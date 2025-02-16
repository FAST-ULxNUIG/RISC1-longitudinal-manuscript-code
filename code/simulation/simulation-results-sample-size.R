# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice) # CRAN v0.12.3.1

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

N_sub_500_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_500_prop_missing_0.1_long_strength_1.rds"
  )
)

N_sub_1000_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_1000_prop_missing_0.1_long_strength_1.rds"
  )
)

N_sub_1000_prop_missing_0.1_long_strength_1_part_2 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_1000_prop_missing_0.1_long_strength_part_21.rds"
  )
)

# -------------------------------------------------------------------------
stopifnot(any(!sapply(N_sub_280_prop_missing_0.1_long_strength_1[1:500], is.null)))
stopifnot(any(!sapply(N_sub_500_prop_missing_0.1_long_strength_1[1:500], is.null)))
stopifnot(any(!sapply(N_sub_1000_prop_missing_0.1_long_strength_1[1:250], is.null)))
stopifnot(any(!sapply(N_sub_1000_prop_missing_0.1_long_strength_1_part_2[1:250], is.null)))

N_sub_280_prop_missing_0.1_long_strength_1$n_cores
N_sub_500_prop_missing_0.1_long_strength_1$n_cores
N_sub_1000_prop_missing_0.1_long_strength_1$n_cores
N_sub_1000_prop_missing_0.1_long_strength_1_part_2$n_cores

results_280 <- N_sub_280_prop_missing_0.1_long_strength_1[1:500]
results_500 <- N_sub_500_prop_missing_0.1_long_strength_1[1:500]
results_1000 <- N_sub_1000_prop_missing_0.1_long_strength_1[1:250]
results_1000[251:500] <- N_sub_1000_prop_missing_0.1_long_strength_1[1:250]

# -------------------------------------------------------------------------

# Check for warnings: -----------------------------------------------------
sum(unlist(lapply(results_280, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))

sum(unlist(lapply(results_500, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))

sum(unlist(lapply(results_1000, function(x) {
  sapply(x$models$spline_subject_ri_side$model_checks$warnings, length)
})))


# Check FPCA Times: -------------------------------------------------------
boxplot(sapply(results_280, function(x) {
  x$pca$time
}),
sapply(results_500, function(x) {
  x$pca$time
}),
sapply(results_1000, function(x) {
  x$pca$time
}))

pca_times_dt <- data.table(rep = c(rep(1:500, times = 2), 1:500), 
           N = rep(c(280, 500, 1000), each = 500),
           time_secs = c(
             sapply(results_280, function(x) {x$pca$time}),
             sapply(results_500, function(x) {x$pca$time}),
             sapply(results_1000, function(x) {x$pca$time})
           ))

pca_time_plot <- ggplot(data = pca_times_dt) +
  aes(x = N, group = N, y = time_secs) + 
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "",
      colour = "Model:",
      y = "Time (seconds)",
      title = "(a) Computation Time (mv-FPCA)") +
  scale_x_continuous(breaks = c(280, 500, 1000), labels = paste0("$N=", c(280, 500, 1000),"$"))

pca_time_plot

# Check Modelling Times: --------------------------------------------------


model_times_280 <- extract_simulation_element(results_280, extractor_function = extract_time)
model_times_280$rep <- 1:500
model_times_280$N <- 280

model_times_500 <- extract_simulation_element(results_500, extractor_function = extract_time)
model_times_500$rep <- 1:500
model_times_500$N <- 500

model_times_1000 <- extract_simulation_element(results_1000, extractor_function = extract_time)
model_times_1000$rep <- 1:500
model_times_1000$N <- 1000


model_times_dt <- data.table(rbind(model_times_280, model_times_500, model_times_1000))
model_times_dt_lng <- melt.data.table(model_times_dt, 
                                      id.vars = c("N", "rep"),
                                      variable.factor = FALSE,
                                      value.factor = FALSE,
                                      variable.name = "model", 
                                      value.name = "time_sec")
model_times_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]

model_fit_time_plot <- ggplot(data = model_times_dt_lng) +
  aes(x = N, group = interaction(N, model), colour = model, y = time_sec) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "Number of Subjects",
       colour = "Model:",
       y = "Time (seconds)",
       title = "(b) Computation Time (Model Fit)") +
  scale_x_continuous(breaks = c(280, 500, 1000), labels = paste0("$N=", c(280, 500, 1000),"$"))
  
model_fit_time_plot

## Prediction Error: ------------------------------------------------------
par(mfrow = c(1, 3))
boxplot(extract_simulation_element(results_280, extractor_function = extract_prediction_error))
boxplot(extract_simulation_element(results_500, extractor_function = extract_prediction_error))
boxplot(extract_simulation_element(results_1000, extractor_function = extract_prediction_error))


pred_error_280 <- extract_simulation_element(results_280, extractor_function = extract_prediction_error)
pred_error_280$rep <- 1:500
pred_error_280$N <- 280

pred_error_500 <- extract_simulation_element(results_500, extractor_function = extract_prediction_error)
pred_error_500$rep <- 1:500
pred_error_500$N <- 500

pred_error_1000 <- extract_simulation_element(results_1000, extractor_function = extract_prediction_error)
pred_error_1000$rep <- 1:500
pred_error_1000$N <- 1000

pred_error_dt <- data.table(rbind(pred_error_280, pred_error_500, pred_error_1000))
pred_error_dt_lng <- melt.data.table(pred_error_dt, 
                                      id.vars = c("N", "rep"),
                                      variable.factor = FALSE,
                                      value.factor = FALSE,
                                      variable.name = "model", 
                                      value.name = "ISPE")
pred_error_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]

pred_error_ISPE_plot <- ggplot(data = pred_error_dt_lng) +
  aes(x = N, group = interaction(N, model), colour = model, y = ISPE) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "",
       colour = "Model:",
       y = "ISPE",
       title = "(c) Individual Prediction Error") +
  scale_x_continuous(breaks = c(280, 500, 1000), labels = paste0("$N=", c(280, 500, 1000),"$"))

pred_error_ISPE_plot 

## Intercept: --------------------------------------------------------------

par(mfrow = c(3, 3))
boxplot(extract_simulation_element(results_280, extractor_function = extract_intercept_error))
boxplot(extract_simulation_element(results_500, extractor_function = extract_intercept_error))
boxplot(extract_simulation_element(results_1000, extractor_function = extract_intercept_error))

fixef_error_intercept_280 <- extract_simulation_element(results_280, extractor_function = extract_intercept_error)
fixef_error_intercept_280$rep <- 1:500
fixef_error_intercept_280$N <- 280

fixef_error_intercept_500 <- extract_simulation_element(results_500, extractor_function = extract_intercept_error)
fixef_error_intercept_500$rep <- 1:500
fixef_error_intercept_500$N <- 500

fixef_error_intercept_1000 <- extract_simulation_element(results_1000, extractor_function = extract_intercept_error)
fixef_error_intercept_1000$rep <- 1:500
fixef_error_intercept_1000$N <- 1000

fixef_error_intercept_dt <- data.table(rbind(fixef_error_intercept_280, 
                                             fixef_error_intercept_500,
                                             fixef_error_intercept_1000))
fixef_error_intercept_dt_lng <- melt.data.table(fixef_error_intercept_dt, 
                                     id.vars = c("N", "rep"),
                                     variable.factor = FALSE,
                                     value.factor = FALSE,
                                     variable.name = "model", 
                                     value.name = "intercept")

## Speed: ----------------------------------------------------------------
boxplot(extract_simulation_element(results_280, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_500, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_1000, extractor_function = extract_speed_error))

fixef_error_speed_280 <- extract_simulation_element(results_280, extractor_function = extract_speed_error)
fixef_error_speed_280$rep <- 1:500
fixef_error_speed_280$N <- 280

fixef_error_speed_500 <- extract_simulation_element(results_500, extractor_function = extract_speed_error)
fixef_error_speed_500$rep <- 1:500
fixef_error_speed_500$N <- 500

fixef_error_speed_1000 <- extract_simulation_element(results_1000, extractor_function = extract_speed_error)
fixef_error_speed_1000$rep <- 1:500
fixef_error_speed_1000$N <- 1000

fixef_error_speed_dt <- data.table(rbind(fixef_error_speed_280, fixef_error_speed_500, fixef_error_speed_1000))
fixef_error_speed_dt_lng <- melt.data.table(fixef_error_speed_dt, 
                                                id.vars = c("N", "rep"),
                                                variable.factor = FALSE,
                                                value.factor = FALSE,
                                                variable.name = "model", 
                                                value.name = "speed")


## Sex: --------------------------------------------------------------------
boxplot(extract_simulation_element(results_280, extractor_function = extract_sex_error))
boxplot(extract_simulation_element(results_500, extractor_function = extract_speed_error))
boxplot(extract_simulation_element(results_1000, extractor_function = extract_speed_error))


fixef_error_sex_280 <- extract_simulation_element(results_280, extractor_function = extract_sex_error)
fixef_error_sex_280$rep <- 1:500
fixef_error_sex_280$N <- 280

fixef_error_sex_500 <- extract_simulation_element(results_500, extractor_function = extract_sex_error)
fixef_error_sex_500$rep <- 1:500
fixef_error_sex_500$N <- 500

fixef_error_sex_1000 <- extract_simulation_element(results_1000, extractor_function = extract_sex_error)
fixef_error_sex_1000$rep <- 1:500
fixef_error_sex_1000$N <- 1000

fixef_error_sex_dt <- data.table(rbind(fixef_error_sex_280, fixef_error_sex_500, fixef_error_sex_1000))
fixef_error_sex_dt_lng <- melt.data.table(fixef_error_sex_dt, 
                                            id.vars = c("N", "rep"),
                                            variable.factor = FALSE,
                                            value.factor = FALSE,
                                            variable.name = "model", 
                                            value.name = "sex")


fixef_error_dt_wide <- merge.data.table(merge.data.table(fixef_error_intercept_dt_lng,
                                        fixef_error_speed_dt_lng,
                                        by = c("N", "rep", "model"),
                                        all = TRUE), 
                                        y = fixef_error_sex_dt_lng,
                                        by = c("N", "rep", "model"),
                                        all = TRUE)

fixef_error_dt_lng <- melt.data.table(fixef_error_dt_wide, 
                                      id.vars = c("N", "rep", "model"),
                                      variable.name = "parameter",
                                      measure.vars = c("intercept", "speed", "sex"), 
                                      value.name = "MISE", 
                                      variable.factor = FALSE, 
                                      value.factor = FALSE)

fixef_error_dt_lng[, parameter := factor(parameter,
                                         levels = c("intercept", "sex", "speed"),
                                         labels = c("Intercept $\\boldsymbol{\\beta}_0 (t, T)$",
                                                    "Sex $\\boldsymbol{\\beta}_1 (t)$",
                                                    "Speed $\\boldsymbol{\\beta}_2 (t)$"))]

fixef_error_dt_lng[, model := factor(model, 
                                     levels = c("poly", "naive", "spline_subject_ri_side", "fpca"),
                                     labels = c("Polynomial", "Naive", "Spline", "ml-FPCA"))]


fixef_error_plot <- ggplot(data = fixef_error_dt_lng) +
  facet_wrap(~ parameter, nrow = 1, ncol = 3, scales = "free_y") +
  aes(x = N, group = interaction(N, model), colour = model, y = sqrt(MISE)) +
  geom_boxplot(outlier.size = 0.5)  +
  labs(x = "Number of Subjects",
       colour = "Model:",
       y = "$\\sqrt{\\text{ISE}}$",
       title = "(d) Fixed Effects Estimation Error") +
  scale_x_continuous(breaks = c(280, 500, 1000), labels = paste0("$N=", c(280, 500, 1000),"$"))

# Combine: ----------------------------------------------------------------
top_panel <- ggarrange(pca_time_plot,
                       model_fit_time_plot, 
                       pred_error_ISPE_plot,
                       nrow = 1,
                       legend = "none")
tikz(file.path(plots_path, "simulation-sample-size-plot.tex"),
     width = 1.2 * doc_width_inches, 
     height = 1.3 * (1/1.618) * doc_width_inches,
     standAlone = TRUE)
ggarrange(top_panel, fixef_error_plot,
          nrow = 2, 
          common.legend = TRUE,
          legend = "bottom",
          heights = c(0.48, 0.52))
dev.off()
tinytex::lualatex(file.path(plots_path, "simulation-sample-size-plot.tex"))
