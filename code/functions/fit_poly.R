require(progress) # CRAN v1.2.2
require(lme4) # CRAN v1.1-30

fit_poly <- function(df_scores,
                     K_retain,
                     degree, 
                     fixef_formula = "sex + speed_cent",
                     control = lmerControl(optCtrl = list(xtol_rel=0,
                                                         xtol_abs=1e-10,
                                                         ftol_rel=0, 
                                                         ftol_abs=1e-10)),
                     REML = TRUE,
                     diagonal_covariance = FALSE) {
  # Step 0 - some routine checks before we start:
  if(!(is.data.frame(df_scores))) stop("df_scores must be a data frame.")
  if(!all(paste0("score_", seq_len(K_retain)) %in% names(df_scores))) {
    stop("df_scores must contain columns score_1, ..., score_{K_retain}")
  }
  if(!(degree >= 2)) {
    stop("degree must be >= 2. If degree = 1, just use slope model.")
  }
  
  # Step 1 -- Construct Polynomial Basis:
  poly_basis <- poly(seq(0, 1, length.out = 101), degree = degree, raw = FALSE)
  df_scores <- add_poly_to_df(df = df_scores, poly_object = poly_basis)
  poly_seq <- paste0("poly_", seq_len(degree))
  
  # Step 2 - Set up formula for model fitting:
  poly_formula <- paste(poly_seq, collapse = " + ")
  if(!diagonal_covariance) {
    ranef_formula <- paste0("(", poly_formula, " | subject_id/side)")
  } else if(diagonal_covariance) {
    ranef_formula <- paste0("(", poly_formula, " || subject_id/side)")
  }
  rhs_formula <- paste(poly_formula, fixef_formula, ranef_formula, sep = " + ")
  
  # Step 3 - Fit series of scalar longitudinal models:
  pb <- progress_bar$new(total = K_retain)
  lme_fit_list <- vector(mode = "list", length = K_retain)
  for(k in seq_len(K_retain)) {
    pb$tick()
    formula_k <- formula(
      paste0("score_", k, " ~ ", rhs_formula)
    )
    lme_fit_list[[k]] <- lmer(formula = formula_k, 
                              data = df_scores,
                              control = control,
                              REML = REML)
  }
  
  # Step 4 - Return:
  list(poly_basis = poly_basis, lme_fit_list = lme_fit_list, K_retain = K_retain)
}


predict_scores_poly <- function(fit_poly_object, newdata) {
  # Step 1 -- Unpack objects
  poly_basis <- fit_poly_object$poly_basis
  lme_fit_list <- fit_poly_object$lme_fit_list
  K_retain <- fit_poly_object$K_retain
  # routine checks:
  stopifnot(sapply(lme_fit_list, function(x) {
    class(x) == "lmerMod"
  }))
  stopifnot(length(lme_fit_list) == K_retain)
  
  # Step 2 - Do predictions
  newdata_add_poly <- add_poly_to_df(df = newdata, poly_object = poly_basis)
  predictions_matrix <- sapply(lme_fit_list, predict, newdata = newdata_add_poly)
  
  # Step 3 - check and name output
  stopifnot(dim(predictions_matrix) == c(nrow(newdata), K_retain))
  colnames(predictions_matrix) <- paste0("pred_score_", seq_len(K_retain))
  predictions_matrix
}


predict_fd_poly <- function(fit_poly_object, newdata, pca_fd_obj) {
  
  # Extract K from poly fit object:
  K_retain <- fit_poly_object$K_retain
  
  # Predict scores:
  scores_predictions_matrix <- predict_scores_poly(fit_poly_object, newdata)
  
  # Use to predict functional data object:
  construct_fd_from_scores(
    pca_fd_obj = pca_fd_obj,
    scores_matrix = scores_predictions_matrix,
    K = K_retain)
}



extract_intercept_fd_poly <- function(fit_poly_object, pca_fd_obj) {
  # Extracts multivariate intercept function at a grid of 101 equally-spaced points
  # on the longitudinal domain [0,1].
  poly_basis <- fit_poly_object$poly_basis
  poly_design_matrix <- cbind(1, poly_basis)
  poly_basis_coefficients <- extract_fixef_coef(lme_list_object = fit_poly_object$lme_fit_list,
                                                K_retain = fit_poly_object$K_retain)
  poly_basis_coefficients <- poly_basis_coefficients[c("(Intercept)", paste0("poly_", seq_len(ncol(poly_basis)))), ]
  stopifnot(ncol(poly_design_matrix) == nrow(poly_basis_coefficients))
  poly_scores <- poly_design_matrix %*% poly_basis_coefficients
  
  intercept_on_longitudinal_grid <- construct_fd_from_scores(pca_fd_obj = pca_fd_obj,
                                                             scores_matrix = poly_scores,
                                                             K = fit_poly_object$K_retain, 
                                                             add_back_mean = TRUE)
  intercept_on_longitudinal_grid
}









