library(tikzDevice) # CRAN v0.12.3.1
library(ggplot2)    # CRAN v3.4.0
library(data.table) # CRAN v1.14.2
source(here::here("code/functions/theme_gunning.R"))

# Some settings for the plots: --------------------------------------------
theme_gunning() # use default theme
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))


# Path to save the outputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")

outputs_path <- here::here("outputs")
eigenfunction_results <- readRDS(file = file.path(outputs_path, "simulation", "eigenfunction-estimation.rds"))


# Unpack: -----------------------------------------------------------------
mvFPC_dt_lng <- eigenfunction_results$mvFPC_dt_lng
mvFPC_dt_lng_mean <- eigenfunction_results$mvFPC_dt_lng_mean
mvFPC_covariates_dt_lng <- eigenfunction_results$mvFPC_covariates_dt_lng
mvFPC_covariates_dt_lng_mean <- eigenfunction_results$mvFPC_covariates_dt_lng_mean
true_mvFPC_dt_lng <- eigenfunction_results$true_mvFPC_dt_lng
true_rotated_mvFPC_dt_lng <- eigenfunction_results$true_rotated_mvFPC_dt_lng


covariates_zero_plot <- ggplot(data = mvFPC_dt_lng[num %in% paste0("mvFPC", 1:3)]) +
  aes(x = t, y = estimate) +
  facet_grid(num_label ~ factor(dimension, levels = c("Hip", "Knee", "Ankle"))) +
  geom_line(aes(group = sim_rep),col = "darkgrey", alpha = 0.25) +
  geom_line(data = true_mvFPC_dt_lng[num %in% paste0("mvFPC", 1:3)], col = "black") +
  geom_line(data = mvFPC_dt_lng_mean[num %in% paste0("mvFPC", 1:3)], col = "black", lty = 2) +
  labs(x = "$t$", y = "$\\psi_k^{(p)}(t)$")

covariates_zero_plot

tikz(file.path(plots_path, "Eigenfunction-estimation-01.tex"),
     width = 1 * doc_width_inches, 
     height = 0.8 * doc_width_inches)
print(covariates_zero_plot)
dev.off()

covariates_true_plot <-  ggplot(data = mvFPC_covariates_dt_lng[num %in% paste0("mvFPC", 1:3)]) +
  aes(x = t, y = estimate) +
  facet_grid(num_label ~ factor(dimension, levels = c("Hip", "Knee", "Ankle"))) +
  geom_line(aes(group = sim_rep),col = "darkgrey", alpha = 0.25) +
  geom_line(data = true_mvFPC_dt_lng[num %in% paste0("mvFPC", 1:3)], col = "black") + 
  geom_line(data = true_rotated_mvFPC_dt_lng[num %in% paste0("mvFPC", 1:3)], col = "red") +
  geom_line(data = mvFPC_covariates_dt_lng_mean[num %in% paste0("mvFPC", 1:3)], col = "black", lty = 2) +
  labs(x = "$t$", y = "$\\psi^{(p)}_k (t)$")

covariates_true_plot 

tikz(file.path(plots_path, "Eigenfunction-estimation-02.tex"),
     width = 1 * doc_width_inches, 
     height = 0.8 * doc_width_inches)
print(covariates_true_plot)
dev.off()
  
  
  
