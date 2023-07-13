bootstrap_of_subjects <- function(df_for_bootstrap,
                                  k_retain,
                                  model = "spline_subject_ri_side",
                                  fixef_formula = "sex + speed_cent",
                                  df, # NB this is spline degrees of freedom
                                  control = lmerControl(optCtrl = list(xtol_rel=0,
                                                                             xtol_abs=1e-10,
                                                                             ftol_rel=0, 
                                                                             ftol_abs=1e-10)),
                                  REML = TRUE,
                                  diagonal_covariance = TRUE,
                                  model_function,
                                  B = 1000,
                                  par_mc = TRUE,
                                  n_cores = 4) {
  require(parallel)
  # checks input
  stopifnot("subject_id" %in% names(df_for_bootstrap))
  stopifnot(paste0("score_", seq_len(k_retain)) %in% names(df_for_bootstrap))
  # only support natural spline models for now
  if(model == "naive_spline_intercept") {
    model_fit_function <- fit_naive_spline_intercept
  } else if(model == "spline_subject_ri_side") {
    model_fit_function <- fit_spline_subject_ri_side
  } else if(model == "spline") {
    model_fit_function <- fit_spline
  } else {
    stop("model must be one of naive_spline_intercept, spline_subject_ri_side, spline.")
  }
  
  # Generate resampled indices before doing parralelisation to avoid RNG problems
  resampled_ids_list <- lapply(seq_len(B), function(b) {
    sample(unique(df_for_bootstrap$subject_id), replace = TRUE)
  })
  
  # do bootstrapping in parallel
  if(!par_mc) {
    bootstrap_list <- lapply(resampled_ids_list, function(x) {
      # create new data frame based on reasmpled ids:
      df_b <- purrr::map_dfr(.x = seq_along(x), .f = function(i) {
        df_for_bootstrap[df_for_bootstrap$subject_id == x[i], ]
      }, .id = "subject_id_b")
      # new pseudo-ID (this is important!)
      df_b$subject_id <- factor(df_b$subject_id_b)
      # fit models
      fit_b <- model_fit_function(df_scores = df_b, 
                                 K_retain = k_retain,
                                 df = df,
                                 fixef_formula = fixef_formula,
                                 control = control,
                                 REML = REML, 
                                 diagonal_covariance = diagonal_covariance,
                                 parallel = FALSE, 
                                 ncores = ncores)
      # extract fixed effects
      fixef_b <- extract_fixef_coef(fit_b$lme_fit_list, K_retain = k_retain)
      # return
      list(fixef = fixef_b, singular = sapply(fit_b$lme_fit_list, isSingular))
    })
  }
  
  if(par_mc) {
    bootstrap_list <- mclapply(resampled_ids_list, function(x) {
      # create new data frame based on resampled ids:
      df_b <- purrr::map_dfr(.x = seq_along(x), .f = function(i) {
        df_for_bootstrap[df_for_bootstrap$subject_id == x[i], ]
      }, .id = "subject_id_b")
      # new pseudo-ID (this is important!)
      df_b$subject_id <- factor(df_b$subject_id_b)
      # fit models
      fit_b <- model_fit_function(df_scores = df_b, 
                                  K_retain = k_retain,
                                  df = df,
                                  fixef_formula = fixef_formula,
                                  control = control,
                                  REML = REML, 
                                  diagonal_covariance = diagonal_covariance,
                                  parallel = FALSE, 
                                  ncores = ncores)
      # extract fixed effects
      fixef_b <- extract_fixef_coef(fit_b$lme_fit_list, K_retain = k_retain)
      # return
      list(fixef = fixef_b, singular = sapply(fit_b$lme_fit_list, isSingular))
    },
    mc.silent = FALSE,
    mc.cores = n_cores)
  }
  
  return(bootstrap_list)
}
