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



# Plot of longitudinal trajectories: --------------------------------------
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
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
    xlab("Longitudinal Time ($T$)") +
    ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
    theme(legend.position = "bottom"))

tikz(file.path(plots_path, "mv-FPC1-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.725 * doc_width_inches)
print(p)
dev.off()