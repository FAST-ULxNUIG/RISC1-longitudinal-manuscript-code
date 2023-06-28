fixInNamespace("face.Cov.mfpca", "refund")
require(refund)
fit_fpca_single_score <- function(df_scores,
                                  score_ind = 1, 
                                  pve = 0.99,
                                  df_spline = 3,
                                  fixef_formula = "sex + speed_cent",
                                  silent = FALSE,
                                  control = lmerControl(optCtrl = list(xtol_rel=0,
                                                                                 xtol_abs=1e-10,
                                                                                 ftol_rel=0, 
                                                                                 ftol_abs=1e-10)),
                                  REML = TRUE) {
  
  time_arg_vals <- seq(0, 1, length.out = 101)
  bin_digits <- 2 # for now, round time values to 2 digits.
  
  require(data.table)
  if(!is.data.table(df_scores)) df_scores <- data.table(df_scores)
  

  # Fit working independence model: -----------------------------------------
  if(!silent) print("Fitting Working Independence Model:")
  spline_basis <- ns(x = time_arg_vals, df = df_spline)
  df_scores <- add_natural_splines_to_df(df = df_scores, spline_object = spline_basis)
  spline_seq <- paste0("spline_", seq_len(ncol(spline_basis)))
  fixef_formula_full <- paste(paste0(spline_seq, collapse = " + "), fixef_formula, sep = " + ")
  ind_model_formula  <- paste0("score_", score_ind, " ~ ", fixef_formula_full)
  lm_ind <- lm(formula = ind_model_formula, data = df_scores)
  score_resid_name <- paste0(c("score", score_ind, "resid"), collapse = "_")
  df_scores[, score_resid_name] <- resid(lm_ind)
  
  
  # Bin data and reshape to apply FACE: -------------------------------------
  if(!silent) print("Reshaping Data and Applying mfpca.face()")
  df_scores$time_binned <- round(df_scores$time, bin_digits)
  stopifnot(setequal(df_scores$time_binned, round(time_arg_vals, bin_digits))) # need to make this more general!!
  score_resid_binned_name <- paste0(score_resid_name, "_binned")
  
  df_scores_lng <- df_scores[, .(mean(get(score_resid_name))), by = .(subject_id, side, time_binned)]
  setnames(df_scores_lng, old = "V1", new = score_resid_binned_name)
  
  df_scores_wide <- dcast.data.table(data = df_scores_lng,
                              formula = subject_id + side ~ time_binned,
                              value.var = paste0(score_resid_binned_name))
  time_arg_vals_char <- paste(time_arg_vals)
  Y_face <- as.matrix(df_scores_wide[, ..time_arg_vals_char])
  mfpca_face <- mfpca.face(Y = Y_face, 
                           id = df_scores_wide$subject_id,
                           visit = df_scores_wide$side,
                           twoway = FALSE, 
                           knots = 10,
                           weight = "subj",
                           pve = pve)
  
  (npc_u <- mfpca_face$npc$level1)
  (npc_v <- mfpca_face$npc$level2)
  if(npc_u == 1 & max(abs(mfpca_face$efunctions$level1[, 1]) - 1) < 0.05) {mfpca_face$efunctions$level1[, 1] <- 1}
  if(npc_v == 1 & max(abs(mfpca_face$efunctions$level2[, 1]) - 1) < 0.05) {mfpca_face$efunctions$level2[, 1] <- 1}
  
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
  for(k in seq_len(npc_u)) {
    df_scores[, paste0("phi_u_", k)] <- efuns_list[["level_1"]][[k]](v = df_scores$time)
  }
  for(k in seq_len(npc_v)) {
    df_scores[, paste0("phi_v_", k)] <- efuns_list[["level_2"]][[k]](v = df_scores$time)
  }
  
  if(!silent) print("Doing final mixed model fit")
  ranef_formula_subject <- paste("(0 + ",
                                 paste(paste0("phi_u_", seq_len(npc_u)), collapse = " + "),
                                 " || subject_id)")
  ranef_formula_side <- paste("(0 + ",
                              paste(paste0("phi_v_", seq_len(npc_v)), collapse = " + "),
                              " || subject_id:side)")
  ranef_formula <- paste(ranef_formula_subject, ranef_formula_side,sep = " + ")
  full_lme_formula <- paste(ind_model_formula, ranef_formula, sep = " + ")
  lme_fit <- lmer(formula = full_lme_formula, data = df_scores, control = control, REML = REML)
  
  # Return Results:
  list(lme_fit = lme_fit, efuns_list = efuns_list, spline_basis = spline_basis)
}


fit_fpca <- function(df_scores,
                     K_retain,
                     pve = 0.99,
                     df_spline = 3,
                     fixef_formula = "sex + speed_cent",
                     control = lmerControl(optCtrl = list(xtol_rel=0,
                                                          xtol_abs=1e-10,
                                                          ftol_rel=0, 
                                                          ftol_abs=1e-10)),
                     REML = TRUE) {
  
  # Step 0 - some routine checks before we start:
  if(!(is.data.frame(df_scores))) stop("df_scores must be a data frame.")
  if(!all(paste0("score_", seq_len(K_retain)) %in% names(df_scores))) {
    stop("df_scores must contain columns score_1, ..., score_{K_retain}")
  }
  
  
  # Step 3 - Fit series of scalar longitudinal models:
  pb <- progress_bar$new(total = K_retain)
  efuns_list <- lme_fit_list <- vector(mode = "list", length = K_retain)
  for(k in seq_len(K_retain)) {
    pb$tick()
    fpca_fit_k <- fit_fpca_single_score(score_ind = k,
                          pve = pve,
                          df_scores = df_scores, 
                          df_spline = df_spline,
                          fixef_formula = fixef_formula, 
                          control = control,
                          REML = REML, 
                          silent = TRUE)
    lme_fit_list[[k]] <- fpca_fit_k$lme_fit
    efuns_list[[k]] <- fpca_fit_k$efuns_list
  }
  spline_basis <- fpca_fit_k$spline_basis
  
  list(lme_fit_list = lme_fit_list,
       efuns_list = efuns_list,
       spline_basis = spline_basis, 
       K_retain = K_retain)
}


predict_scores_fpca <- function(fit_fpca_object, newdata, pca_fd_obj) {
  
  # Extract K from spline fit object:
  # Step 1 -- Unpack objects
  spline_basis <- fit_fpca_object$spline_basis
  lme_fit_list <- fit_fpca_object$lme_fit_list
  K_retain <- fit_fpca_object$K_retain
  efuns_list <- fit_fpca_object$efuns_list
  # routine checks:
  stopifnot(sapply(lme_fit_list, function(x) {
    class(x) == "lmerMod"
  }))
  stopifnot(length(lme_fit_list) == K_retain)
  
  # Step 2 - Do predictions
  newdata_add_spline <- add_natural_splines_to_df(df = newdata, spline_object = spline_basis)
  
  predictions_matrix <- sapply(seq_len(K_retain),
                               function(k) {
                                 newdata_add_spline_k <- newdata_add_spline
                                 for(ku in seq_along(efuns_list[[k]][["level_1"]])) {
                                   newdata_add_spline_k[, paste0("phi_u_", ku)] <- efuns_list[[k]][["level_1"]][[ku]](v = newdata_add_spline_k$time)
                                 }
                                 for(kv in seq_along(efuns_list[[k]][["level_2"]])) {
                                   newdata_add_spline_k[, paste0("phi_v_", kv)] <- efuns_list[[k]][["level_2"]][[kv]](v = newdata_add_spline_k$time)
                                 }
                                 predict(object = lme_fit_list[[k]], newdata = newdata_add_spline_k)
                                 }
                               )
  
  # Step 3 - check and name output
  stopifnot(dim(predictions_matrix) == c(nrow(newdata), K_retain))
  colnames(predictions_matrix) <- paste0("pred_score_", seq_len(K_retain))
  predictions_matrix
}


predict_fd_fpca <- function(fit_fpca_object, newdata, pca_fd_obj) {
  
  # Extract K from spline fit object:
  K_retain <- fit_spline_object$K_retain
  
  # Predict scores:
  scores_predictions_matrix <- predict_scores_fpca(fit_fpca_object, newdata)
  
  # Use to predict functional data object:
  construct_fd_from_scores(
    pca_fd_obj = pca_fd_obj,
    scores_matrix = scores_predictions_matrix,
    K = K_retain)
}



predict_fd_fpca <- function(fit_fpca_object, newdata, pca_fd_obj) {
  
  # Extract K from fpca fit object:
  K_retain <- fit_fpca_object$K_retain
  
  # Predict scores:
  scores_predictions_matrix <- predict_scores_fpca(fit_fpca_object, newdata)
  
  # Use to predict functional data object:
  construct_fd_from_scores(
    pca_fd_obj = pca_fd_obj,
    scores_matrix = scores_predictions_matrix,
    K = K_retain)
}

extract_intercept_fd_spline <- function(fit_spline_object, pca_fd_obj) {
  # Extracts multivariate intercept function at a grid of 101 equally-spaced points
  # on the longitudinal domain [0,1].
  spline_basis <- fit_spline_object$spline_basis
  spline_design_matrix <- cbind(1, spline_basis)
  spline_basis_coefficients <- extract_fixef_coef(
    lme_list_object = fit_spline_object$lme_fit_list,
    K_retain = fit_spline_object$K_retain)
  spline_basis_coefficients <- spline_basis_coefficients[c("(Intercept)", paste0("spline_", seq_len(ncol(spline_basis)))), ]
  stopifnot(ncol(spline_design_matrix) == nrow(spline_basis_coefficients))
  spline_scores <- spline_design_matrix %*% spline_basis_coefficients
  intercept_on_longitudinal_grid <- construct_fd_from_scores(pca_fd_obj = pca_fd_obj,
                                                             scores_matrix = spline_scores,
                                                             K = fit_spline_object$K_retain, 
                                                             add_back_mean = TRUE)
  intercept_on_longitudinal_grid
}


extract_intercept_fd_fpca <- function(fit_fpca_object, pca_fd_obj) {
  # Extracts multivariate intercept function at a grid of 101 equally-spaced points
  # on the longitudinal domain [0,1].
  spline_basis <- fit_fpca_object$spline_basis
  spline_design_matrix <- cbind(1, spline_basis)
  spline_basis_coefficients <- extract_fixef_coef(
    lme_list_object = fit_fpca_object$lme_fit_list,
    K_retain = fit_fpca_object$K_retain)
  spline_basis_coefficients <- spline_basis_coefficients[c("(Intercept)", paste0("spline_", seq_len(ncol(spline_basis)))), ]
  stopifnot(ncol(spline_design_matrix) == nrow(spline_basis_coefficients))
  spline_scores <- spline_design_matrix %*% spline_basis_coefficients
  
  intercept_on_longitudinal_grid <- construct_fd_from_scores(pca_fd_obj = pca_fd_obj,
                                                             scores_matrix = spline_scores,
                                                             K = fit_fpca_object$K_retain, 
                                                             add_back_mean = TRUE)
  intercept_on_longitudinal_grid
}


