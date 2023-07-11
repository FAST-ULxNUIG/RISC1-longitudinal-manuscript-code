# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0

# Functions Used: ---------------------------------------------------------
source(here::here("code/functions/theme_gunning.R"))
source("code/functions/fit_spline_subject_ri_side.R")
source("code/functions/extract_fixef_coef.R")
source("code/functions/construct_fd_from_scores.R")
source("code/functions/center_fd_around_new_mean.R")
source("code/functions/decenter_fd_around_new_mean.R")

# Load Results: -----------------------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
model_fit_results <- readRDS(
  here::here("outputs", "model-fit-results.rds"))

# Path to save the outputs of analysis: ------------------------------------
outputs_path <- here::here("outputs")
plots_path <- file.path(outputs_path, "figures")


# Graphics Settings: ------------------------------------------------------
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
theme_gunning()
theme_update(
  axis.text = element_text(size = 10.5),
  axis.title = element_text(size = 10.5),
  strip.text = element_text(size = 10.5),
  plot.title = element_text(size = 11.5),
  legend.text = element_text(size = 10.5),
  legend.title = element_text(size = 11)
)

# Extract Objects from Results: -------------------------------------------
# From basis transformation results:
mfpca <- basis_transformation_results$mfpca
k_retain <- basis_transformation_results$k_retain
# From model fitting results:
spline_ri_model <- model_fit_results$spline_ri_model



# Do it for overall intercept function: -----------------------------------
intercept_fd_on_grid <- extract_intercept_fd_spline(fit_spline_object = spline_ri_model,
                                                    pca_fd_obj = mfpca)
intercept_eval_on_grid_array <- eval.fd(0:100, intercept_fd_on_grid)
intercept_eval_on_grid_mat <- rbind(intercept_eval_on_grid_array[,,1],
                                    intercept_eval_on_grid_array[,,2],
                                    intercept_eval_on_grid_array[,,3])

intercept_eval_on_grid_dt <- data.table(t = rep(0:100, times = 3),
                                        dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
                                        intercept_eval_on_grid_mat)
names(intercept_eval_on_grid_dt)[-c(1:2)] <- seq(0, 1, length.out = 101)


intercept_eval_on_grid_dt_lng <- melt.data.table(intercept_eval_on_grid_dt, 
                                                 id.vars = c("t", "dimension"), 
                                                 measure.vars = paste0(seq(0, 1, length.out = 101)),
                                                 variable.name = "long_time",
                                                 value.name = "angle", 
                                                 variable.factor = FALSE, 
                                                 value.factor = FALSE)
intercept_eval_on_grid_dt_lng[, long_time := as.numeric(long_time)]
intercept_eval_on_grid_dt_lng[, dimension := factor(dimension, levels = c("Hip", "Knee", "Ankle"))]

p1 <- ggplot(data = intercept_eval_on_grid_dt_lng[long_time %in% seq(0, 1, by = 0.05)]) +
  aes(x = t, y = angle, group = long_time, colour = long_time) +
  facet_wrap(~ dimension, scales = "free_y") +
  geom_line(linewidth = 0.4) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Angle ($^{\\circ}$)",
       colour = "Longitudinal Time ($T$):",
       title = "(a) Intercept Function") +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust = 0.75)) +
  scale_color_gradientn(colours = rainbow(10),
                        breaks = c(0, 0.5, 1),
                        labels = paste(c(0, 0.5, 1)))

p1



# Centered Plot: ----------------------------------------------------------
intercept_fd_on_grid_centered <- center_fd_around_new_mean(fdobj = intercept_fd_on_grid, mean.fd.obj = mfpca$meanfd)

intercept_eval_on_grid_centered_array <- eval.fd(0:100, intercept_fd_on_grid_centered)
intercept_eval_on_grid_centered_mat <- rbind(intercept_eval_on_grid_centered_array[,,1],
                                    intercept_eval_on_grid_centered_array[,,2],
                                    intercept_eval_on_grid_centered_array[,,3])

intercept_eval_on_grid_centered_dt <- data.table(t = rep(0:100, times = 3),
                                        dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
                                        intercept_eval_on_grid_centered_mat)
names(intercept_eval_on_grid_centered_dt)[-c(1:2)] <- seq(0, 1, length.out = 101)


intercept_eval_on_grid_centered_dt_lng <- melt.data.table(intercept_eval_on_grid_centered_dt, 
                                                 id.vars = c("t", "dimension"), 
                                                 measure.vars = paste0(seq(0, 1, length.out = 101)),
                                                 variable.name = "long_time",
                                                 value.name = "angle", 
                                                 variable.factor = FALSE, 
                                                 value.factor = FALSE)
intercept_eval_on_grid_centered_dt_lng[, long_time := as.numeric(long_time)]
intercept_eval_on_grid_centered_dt_lng[, dimension := factor(dimension, levels = c("Hip", "Knee", "Ankle"))]

p2 <- ggplot(data = intercept_eval_on_grid_centered_dt_lng[long_time %in% seq(0, 1, by = 0.05)]) +
  aes(x = t, y = angle, group = long_time, colour = long_time) +
  facet_wrap(~ dimension, scales = "free_y") +
  geom_line() +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Angle ($^{\\circ}$)",
       colour = "Longitudinal Time ($T$):",
       title = "(b) Intercept Function (Centred around Overall Mean)") +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust = 0.75)) +
  scale_color_gradientn(colours = rainbow(10),
                        breaks = c(0, 0.5, 1),
                        labels = paste(c(0, 0.5, 1)))

p2


(p <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "bottom"))

tikz(file.path(plots_path, "time-varying-mean.tex"),
     width = (4/3) * doc_width_inches, 
     height = 1 * doc_width_inches)
print(p)
dev.off()

