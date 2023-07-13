# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0

# Helper functions: -------------------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
source(here::here("code", "functions", "fit_spline_subject_ri_side.R"))

# Paths to files: ---------------------------------------------------------
plots_path <- here::here("outputs", "figures")
outputs_path <- here::here("outputs")

# Graphics Settings: ------------------------------------------------------
theme_gunning()
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Read in Results: --------------------------------------------------------
parameter_results_dt <- readRDS(file.path(outputs_path, "parameter_results_dt"))

parameter_results_dt[,
                     beta_label_part_1 := fcase(
                       parameter == "(Intercept)", "Intercept",
                       parameter == "spline_1", "Spline Coefficient 1",
                       parameter == "spline_2", "Spline Coefficient 2",
                       parameter == "spline_3", "Spline Coefficient 3",
                       parameter == "spline_4", "Spline Coefficient 4",
                       parameter == "risinjured_greater_than_2_yr", "RIS: Injured $>2$ yr.",
                       parameter == "risinjured_1_to_2_yr", "RIS: Injured $1-2$ yr.",
                       parameter == "risinjured_less_than_1_yr", "RIS: Injured $<1$ yr.",
                       parameter == "sexfemale", "Sex",
                       parameter == "speed_cent", "Speed",
                       parameter == "weight_kg_cent", "Weight",
                       parameter == "height_cm_cent", "Height",
                       parameter == "age_cent", "Age"
                     )]
parameter_results_dt[, beta_label_part_1 := factor(
  beta_label_part_1, 
  levels = c(
    "Intercept",
    "Spline Coefficient 1",
    "Spline Coefficient 2",
    "Spline Coefficient 3",
    "Spline Coefficient 4",
    "Speed",
    "Sex",
    "RIS: Injured $>2$ yr.",
    "RIS: Injured $1-2$ yr.",
    "RIS: Injured $<1$ yr.",
    "Age",
    "Weight",
    "Height"
  ))]


parameter_results_dt[, dimension := factor(dimension, levels = c("Hip", "Knee", "Ankle"))]
spline_names <- paste0("spline_", 1:4)
p <- ggplot(data = parameter_results_dt[!(parameter %in% c(spline_names, "(Intercept)"))]) +
  aes(x = t) +
  facet_wrap(dimension ~ beta_label_part_1, scales = "free_y", ncol = 4) +
  geom_line(aes(y = point_est)) +
  geom_ribbon(mapping = aes(ymin = sim_wald_lower, ymax = sim_wald_upper), alpha = 0.25) +
geom_hline(yintercept = 0, col = "darkgrey") +
  geom_line(aes(y = pw_wald_lower), lty = 3) +
  geom_line(aes(y = pw_wald_upper), lty = 3) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "\\leftarrow $\\beta_{a}^{(p)} (t)$ \\rightarrow")

p_hip <- ggplot(data = parameter_results_dt[dimension == "Hip" &!(parameter %in% c(spline_names, "(Intercept)"))]) +
  aes(x = t) +
  geom_hline(yintercept = 0, col = "darkgrey") +
  facet_wrap(~beta_label_part_1, scales = "free_y", ncol = 4) +
  geom_line(aes(y = point_est)) +
  geom_ribbon(mapping = aes(ymin = sim_wald_lower, ymax = sim_wald_upper), alpha = 0.25, fill = "cornflowerblue") +
  geom_line(aes(y = pw_wald_lower), lty = 3) +
  geom_line(aes(y = pw_wald_upper), lty = 3) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "$\\leftarrow$ Extension \\hspace{0.25cm} $\\beta_{a}^{(hip)} (t)$ \\hspace{0.25cm}  Flexion $\\rightarrow$",
       title = "Hip")

p_knee <- ggplot(data = parameter_results_dt[dimension == "Knee" &!(parameter %in% c(spline_names, "(Intercept)"))]) +
  aes(x = t) +
  geom_hline(yintercept = 0, col = "darkgrey") +
  facet_wrap(~beta_label_part_1, scales = "free_y", ncol = 4) +
  geom_line(aes(y = point_est)) +
  geom_ribbon(mapping = aes(ymin = sim_wald_lower, ymax = sim_wald_upper), alpha = 0.25, fill = "cornflowerblue") +
  geom_line(aes(y = pw_wald_lower), lty = 3) +
  geom_line(aes(y = pw_wald_upper), lty = 3) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "$\\leftarrow$ Extension \\hspace{0.25cm} $\\beta_{a}^{(hip)} (t)$ \\hspace{0.25cm}  Flexion $\\rightarrow$",
       title = "Knee")

p_ankle <- ggplot(data = parameter_results_dt[dimension == "Ankle" &!(parameter %in% c(spline_names, "(Intercept)"))]) +
  aes(x = t) +
  geom_hline(yintercept = 0, col = "darkgrey") +
  facet_wrap(~beta_label_part_1, scales = "free_y", ncol = 4) +
  geom_line(aes(y = point_est)) +
  geom_ribbon(mapping = aes(ymin = sim_wald_lower, ymax = sim_wald_upper), alpha = 0.25, fill = "cornflowerblue") +
  geom_line(aes(y = pw_wald_lower), lty = 3) +
  geom_line(aes(y = pw_wald_upper), lty = 3) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "$\\leftarrow$ Plantarflexion \\hspace{0.05cm} $\\beta_{a}^{(ankle)} (t)$ \\hspace{0.05cm}  Dorsiflexion $\\rightarrow$",
       title = "Ankle")

combined_plot <- ggarrange(p_hip, p_knee, p_ankle, nrow = 3)



tikz(file.path(plots_path, "fixef_coef_plot.tex"),
     width = 1 * doc_width_inches, 
     height = (6/4) * doc_width_inches,
     standAlone = TRUE)
print(combined_plot)
dev.off()

tinytex::lualatex(file.path(plots_path, "fixef_coef_plot.tex"))
