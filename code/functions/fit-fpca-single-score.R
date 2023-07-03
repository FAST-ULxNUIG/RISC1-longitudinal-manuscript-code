# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(splines)
library(refund)     # CRAN v0.1-26
source(here::here("code/functions/theme_gunning.R"))
source(here::here("code/functions/fit_spline.R"))
source(here::here("code/functions/fit_spline_subject_ri_side.R"))
source(here::here("code/functions/fit_poly.R"))
source(here::here("code/functions/add_natural_splines_to_df.R"))



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

covariates_and_scores_dt_train[, time := long_time]

score_ind <- 8; df <- covariates_and_scores_dt_train; fixef_formula = "ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent"

fit_fpca_single_score <- function(score_ind = 1, df = covariates_and_scores_dt_train, 
                                  fixef_formula = "ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent") {
  
  time_arg_vals <- seq(0, 1, length.out = 101)
  bin_digits <- 2 # for now, round time values to 2 digits.
  
  require(data.table)
  if(!is.data.table(df)) df <- data.table(df)
  
  # Step 1 -- Construct Spline Basis:
  print("Fitting Working Independence Model:")
  spline_basis <- ns(x = time_arg_vals, df = 3)
  df <- add_natural_splines_to_df(df = df, spline_object = spline_basis)
  spline_seq <- paste0("spline_", seq_len(ncol(spline_basis)))
  fixef_formula_full <- paste(paste0(spline_seq, collapse = " + "), fixef_formula, sep = " + ")
  ind_model_formula  <- paste0("score_", score_ind, " ~ ", fixef_formula_full)
  lm_ind <- lm(formula = ind_model_formula, data = df)
  score_resid_name <- paste0(c("score", score_ind, "resid"), collapse = "_")
  df[, score_resid_name] <- resid(lm_ind)
  

  # Bin data and reshape to apply FACE: -------------------------------------
  df$time_binned <- round(df$time, bin_digits)
  stopifnot(setequal(df$time_binned, round(time_arg_vals, bin_digits))) # need to make this more general!!
  score_resid_binned_name <- paste0(score_resid_name, "_binned")
  
  df_lng <- df[, .(mean(get(score_resid_name))), by = .(subject_id, side, time_binned)]
  setnames(df_lng, old = "V1", new = score_resid_binned_name)
 
  df_wide <- dcast.data.table(data = df_lng,
                              formula = subject_id + side ~ time_binned,
                              value.var = paste0(score_resid_binned_name))
  time_arg_vals_char <- paste(time_arg_vals)
  Y_face <- as.matrix(df_wide[, ..time_arg_vals_char])
  mfpca_face <- mfpca.face(Y = Y_face, 
                           id = df_wide$subject_id,
                           visit = df_wide$side,
                           twoway = FALSE, 
                           knots = 10,
                           weight = "subj",
                           pve = 0.995)

  (npc_u <- mfpca_face$npc$level1)
  (npc_v <- mfpca_face$npc$level2)
  if(npc_u == 1 & max(abs(mfpca_face$efunctions$level1[, 1]) - 1) < 0.05) {mfpca_face$efunctions$level1[, 1] <- 1}
  if(npc_v == 1 & max(abs(mfpca_face$efunctions$level2[, 1]) - 1) < 0.05) {mfpca_face$efunctions$level2[, 1] <- 1}
  # if(max(abs(mfpca_face$efunctions$level1[, 1]) - 1) < 0.03) {mfpca_face$efunctions$level1[, 1] <- 1}
  # if(max(abs(mfpca_face$efunctions$level2[, 1]) - 1) < 0.03) {mfpca_face$efunctions$level2[, 1] <- 1}
  
  # Extract Eigenfunctions:
  efuns_list <- vector(mode = "list", length = 2)
  names(efuns_list) <- c("level_1", "level_2")
  efuns_list[["level_1"]] <- vector(mode = "list", length = npc_u)
  efuns_list[["level_2"]] <- vector(mode = "list", length = npc_v)
  for(k in seq_len(npc_u)) {
    efuns_list[["level_1"]][[k]] <- approxfun(time_arg_vals, mfpca_face$efunctions$level1[, k])
  }
  for(k in seq_len(npc_v)) {
    efuns_list[["level_2"]][[k]] <- approxfun(time_arg_vals, mfpca_face$efunctions$level2[, k])
  }
  

  # -------------------------------------------------------------------------
  for(k in seq_len(npc_u)) {
    df[, paste0("phi_u_", k)] <- efuns_list[["level_1"]][[k]](v = df$time)
  }
  for(k in seq_len(npc_v)) {
    df[, paste0("phi_v_", k)] <- efuns_list[["level_2"]][[k]](v = df$time)
  }
  
  ranef_formula_subject <- paste("(0 + ",
                                 paste(paste0("phi_u_", seq_len(npc_u)), collapse = " + "),
                                 " || subject_id)")
  ranef_formula_side <- paste("(0 + ",
                                 paste(paste0("phi_v_", seq_len(npc_v)), collapse = " + "),
                                 " || subject_id:side)")
  
  ranef_formula <- paste(ranef_formula_subject, ranef_formula_side,sep = " + ")
  
  full_lme_formula <- paste(ind_model_formula, ranef_formula, sep = " + ")
  
  lme_fit <- lmer(formula = full_lme_formula, data = df)
  
}



poly_model <- lmer(score_1 ~ spline_1 + spline_2 + spline_3 + ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent + (spline_1 + spline_2 + spline_3 ||subject_id) + (1 + time | subject_id:side), data = df)
ri_model <- lmer(score_1 ~ spline_1 + spline_2 + spline_3 + ris + sex + speed_cent + age_cent + height_cm_cent + weight_kg_cent + (1 |subject_id) + (1|subject_id:side), data = df)
AIC(poly_model); AIC(lme_fit); AIC(ri_model)


plot(fixef(ri_model), fixef(lme_fit))
abline(0, 1)
plot(fixef(poly_model), fixef(lme_fit))
abline(0, 1)
plot(coef(lm_ind), fixef(lme_fit))
abline(0, 1)
plot(coef(lm_ind), fixef(ri_model))
abline(0, 1)
plot(fixef(poly_model), fixef(ri_model))
abline(0, 1)


matplot(mfpca_face$efunctions$level1, type = "l", lty = 1)
matplot(mfpca_face$efunctions$level2, type = "l", lty = 1)
