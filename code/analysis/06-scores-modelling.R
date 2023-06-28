# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
source(here::here("code/functions/theme_gunning.R"))

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


# Modelling ---------------------------------------------------------------
set.seed(22081996)
(sample_subjects <- sample(unique(covariates_and_scores_dt_train$subject_id), size = 6))
                    # [1] P_4035 P_4206 P_4226 P_4283 P_4240 P_4046
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
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  xlab("Longitudinal Time ($T$)") +
  ylab("mv-FPC1 Score ($y_{ijl, 1}^*$)") +
  theme(legend.position = "bottom"))

tikz(file.path(plots_path, "mv-FPC1-plot.tex"),
     width = (1) * doc_width_inches, 
     height = 0.725 * doc_width_inches)
print(p)
dev.off()

dev.off()


# Modelling: --------------------------------------------------------------
# And join back up into data.table for modelling:
covariates_and_scores_dt_train <- cbind(covariates_dt_train, scores_train)
splines <- ns(covariates_and_scores_dt_train$long_time, df = 5)
spline_names <- colnames(splines) <- paste0("spline_", 1:5)
covariates_and_scores_dt_train <- cbind(covariates_and_scores_dt_train, splines)

paste(spline_names, collapse = " + ")
fixef_structure <- "spline_1 + spline_2 + spline_3 + spline_4 + spline_4 + spline_5 + ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent"
ranef_structure <- "+ (spline_1 + spline_2 + spline_3 + spline_4 + spline_5 | subject_id) + (1|subject_id:side)"
lme_time <- lme_list <- vector(mode = "list", length = k_retain)









k<-1
for(k in seq_len(k_retain)) {
  print(paste("Iteration", k))
  lmer_formula <- formula(paste0("score_", k, "~", fixef_structure, ranef_structure))
  lme_time[[k]] <- system.time({lme_list[[k]] <- lmer(formula = lmer_formula, data = covariates_and_scores_dt_train,
               control= lmerControl(optCtrl = list(xtol_rel=0,
                                                   xtol_abs=1e-10,
                                                   ftol_rel=0, 
                                                   ftol_abs=1e-10)))})
}


summary(lme_list[[21]])

sum(sapply(lme_time, FUN = function(x) {x["elapsed"]})) / 60

lme_list[[2]]@optinfo



set.seed(22081996)
(sample_subjects <- sample(unique(covariates_and_scores_dt_train$subject_id), size = 6))
covariates_and_scores_dt_train$score_1_fitted <- fitted(lme1)
covariates_and_scores_dt_train$score_2_fitted <- fitted(lme_list[[2]])
covariates_and_scores_dt_train$score_3_fitted <- fitted(lme_list[[3]])
covariates_and_scores_dt_train$score_4_fitted <- fitted(lme_list[[4]])

(sample_subjects <- sample(unique(covariates_and_scores_dt_train$subject_id), size = 6))
# [1] P_4199 P_4048 P_4131 P_4248 P_4294 P_4164
ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id) +
  geom_point(aes(y = score_1), alpha = 0.5) +
  geom_line(aes(y = score_1_fitted))

ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id, scales = "free_y") +
  geom_point(aes(y = score_1), alpha = 0.5) +
  geom_line(aes(y = score_1_fitted))


ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id) +
  geom_point(aes(y = score_2), alpha = 0.5) +
  geom_line(aes(y = score_2_fitted))

ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id, scales = "free_y") +
  geom_point(aes(y = score_2), alpha = 0.5) +
  geom_line(aes(y = score_2_fitted))


ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id) +
  geom_point(aes(y = score_3), alpha = 0.5) +
  geom_line(aes(y = score_3_fitted))

ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id, scales = "free_y") +
  geom_point(aes(y = score_3), alpha = 0.5) +
  geom_line(aes(y = score_3_fitted))


# -------------------------------------------------------------------------

ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id) +
  geom_point(aes(y = score_4), alpha = 0.5) +
  geom_line(aes(y = score_4_fitted))

ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id, scales = "free_y") +
  geom_point(aes(y = score_4), alpha = 0.5) +
  geom_line(aes(y = score_4_fitted))


# -------------------------------------------------------------------------


sample_subjects <- "P_4219"
sample_subjects <- c("P_4145", "P_4289", "P_4248", "P_4138", "P_4130", "P_4196")

change<- apply(ranef(lme_list[[4]])[[2]][, 2:6], 1, FUN = function(x) {sum(x^2)})
sample_subjects <- names(change[rev(order(change))[1:6]])




# -------------------------------------------------------------------------

library(refund)
covariates_and_scores_dt_train$score_1_r
covariates_and_scores_dt_train_score_1 <- covariates_and_scores_dt_train[, .(subject_id, side, stride_num, long_time, score_1=score_3_resid)]
covariates_and_scores_dt_train_score_1[, long_time_bin := round(long_time, 2)]

covariates_and_scores_dt_train_score_1 <- covariates_and_scores_dt_train_score_1[, .(score_1_bin = mean(score_1)),
                                                                                 by = .(subject_id,side, long_time_bin)]

covariates_and_scores_dt_train_score_1[, .N, by = .(long_time_bin)][, plot(long_time_bin, N)]


wide_dt <- dcast.data.table(data = covariates_and_scores_dt_train_score_1, subject_id + side ~ long_time_bin, value.var = "score_1_bin")
names(wide_dt)

Y <- as.matrix(wide_dt[, paste0(seq(0, 1, by = 0.01))])




wide_dt[, subject_id_int := as.integer(subject_id)]
wide_dt[, side_int := as.integer(side)]

mfpca_face <- mfpca.face(Y = Y, id = wide_dt$subject_id_int,
                         visit = wide_dt$side_int,
                         weight = "subj", 
                         pve = 0.995,
                         knots = 10,
                         argvals = seq(0, 1, by = 0.01))

mfpca_face$evalues

mfpca_face$sigm

matplot(mfpca_face$Xhat.subject)

covariates_and_scores_dt_train$phi_u_1 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level1[, 1],
                                                 xout = covariates_and_scores_dt_train$long_time)$y
covariates_and_scores_dt_train$phi_u_2 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level1[, 2],
                                                 xout = covariates_and_scores_dt_train$long_time)$y
covariates_and_scores_dt_train$phi_u_3 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level1[, 3],
                                                 xout = covariates_and_scores_dt_train$long_time)$y
covariates_and_scores_dt_train$phi_u_4 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level1[, 4],
                                                 xout = covariates_and_scores_dt_train$long_time)$y
covariates_and_scores_dt_train$phi_u_5 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level1[, 5],
                                                 xout = covariates_and_scores_dt_train$long_time)$y

covariates_and_scores_dt_train$phi_v_1 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level2[, 1],
                                                 xout = covariates_and_scores_dt_train$long_time)$y
covariates_and_scores_dt_train$phi_v_2 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level2[, 2],
                                                 xout = covariates_and_scores_dt_train$long_time)$y
covariates_and_scores_dt_train$phi_v_3 <- approx(x = seq(0, 1, by = 0.01),
                                                 y = mfpca_face$efunctions$level2[, 3],
                                                 xout = covariates_and_scores_dt_train$long_time)$y


test_lme_fpca <- lmer(score_1 ~ spline_1 + spline_2 + spline_3 + spline_4 + spline_4 + spline_5 + ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent +
                        (0 + phi_u_1 + phi_u_2 + phi_u_3 + phi_u_4 + phi_u_5  || subject_id) + (0 + phi_v_1 + phi_v_2 + phi_v_3 || subject_id:side), 
                      data = covariates_and_scores_dt_train)
test_lme_bspline <- lme_list[[1]]

test_lme_very_simple <- lmer(score_1 ~ spline_1 + spline_2 + spline_3 + spline_4 + spline_4 + spline_5 + ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent +
                               (1| subject_id) + (1| subject_id:side), 
                             data = covariates_and_scores_dt_train)

covariates_and_scores_dt_train$score_1_fitted_bspline <- fitted(test_lme_bspline)
covariates_and_scores_dt_train$score_1_fitted_fpca <- fitted(test_lme_fpca)

plot(fitted(test_lme_bspline), fitted(test_lme_fpca))
plot(fitted(), fitted(test_lme_fpca))

plot(fixef(test_lme_bspline), fixef(test_lme_fpca))
abline(a = 0, b = 1)


change<- apply(ranef(lme_list[[1]])[[2]][, 2:6], 1, FUN = function(x) {sum(x^2)})
sample_subjects <- names(change[rev(order(change))[12:18]])


ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id) +
  geom_point(aes(y = score_1), alpha = 0.3) +
  geom_line(aes(y = score_1_fitted_bspline)) +
  geom_line(aes(y = score_1_fitted_fpca), linetype = 2)

ggplot(data = covariates_and_scores_dt_train[subject_id %in% sample_subjects]) +
  aes(x = long_time, group = side, colour = side) +
  scale_color_manual(values = c("red3", "#619CFF"), name = "Side:") +
  facet_wrap(~ subject_id, scales = "free_y") +
  geom_point(aes(y = score_1), alpha = 0.3) +
  geom_line(aes(y = score_1_fitted_bspline)) +
  geom_line(aes(y = score_1_fitted_fpca), linetype = 2)








