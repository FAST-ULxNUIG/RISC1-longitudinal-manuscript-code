extract_diag_var <- function(lme_fit_list) {
  var_mat <- sapply(lme_fit_list, function(x) {
    diag(vcov(x))
  })
  stopifnot(nrow(var_mat) == length(fixef(lme_fit_list[[1]])))
  stopifnot(ncol(var_mat) == length(lme_fit_list))
  rownames(var_mat) <- rownames(vcov(lme_fit_list[[1]]))
  var_mat
}

mvn_sim <- function(coef_point_est, # Vector of coefficients for the parameter function (length = k)
                    coef_covar_mat, # Estimated covariance matrix for coef_point_est (dimension = k x k)
                    Psi_basis, # Basis to which the coefficieients in coef_point_est correspond (dimension = npts x k)
                    N_simulation_mvn = 10000, # Number of samples to draw from MVN distribution.
                    coverage_level = 0.95) { # coverage level for the two-sided CI. 
  
  k <- length(coef_point_est)
  
  # Checks:
  stopifnot(is.vector(coef_point_est))
  stopifnot(k == ncol(Psi_basis))
  stopifnot(dim(coef_covar_mat) == c(k, k))
  # Check this is a proper covariance matrix being supplied:
  stopifnot(isSymmetric(coef_covar_mat))
  stopifnot(matrixcalc::is.positive.definite(coef_covar_mat))
  stopifnot((coverage_level > 0) & (coverage_level < 1))
  
  npts <- nrow(Psi_basis)
  
  # Get pointwise point estimate:
  fun_point_est <- (coef_point_est %*% t(Psi_basis))[1, ]
  
  # and pointwise standard error:
  fun_se <- sqrt(diag(
    Psi_basis %*% coef_covar_mat %*% t(Psi_basis)
  ))
  
  # Simulate from distribution of coefficients:
  coefs_samples <- mvtnorm::rmvnorm(n = N_simulation_mvn,
                                    mean = coef_point_est, 
                                    sigma = coef_covar_mat)
  # Convert simulated vectors to functions # by multiplying them by basis functions:
  fun_samples <- coefs_samples %*% t(Psi_basis)
  
  # Center samples on fixed effects point estimate
  fun_samples_cent <- sweep(fun_samples,
                            MARGIN = 2, 
                            STATS = fun_point_est,
                            FUN = "-",
                            check.margin = TRUE)
  
  # Now also scale these by pointwise standard error,
  # to give stanardised curves:
  fun_samples_stand <- sweep(fun_samples_cent,
                             MARGIN = 2, 
                             STATS = fun_se,
                             FUN = "/",
                             check.margin = TRUE)
  
  # calculate the maximum T statistic
  maxT_dist <- apply(X = fun_samples_stand,
                     MARGIN = 1,
                     FUN = function(x) max(abs(x)))
  
  # Obtain quantile of distribution:
  maxT_q <- quantile(maxT_dist, probs = coverage_level)
  
  if(maxT_q < 2) warning("Something is weird -- Simultaneous intervals narrower than pointwise!")
  
  lower_sim <- fun_point_est - maxT_q * fun_se
  upper_sim <- fun_point_est + maxT_q * fun_se
  
  lower_pw <- fun_point_est - 2 * fun_se
  upper_pw <- fun_point_est + 2 * fun_se
  
  list(point_est = fun_point_est,
    sim = list(lower = lower_sim, upper = upper_sim),
    pw = list(lower = lower_pw, upper = upper_pw),
    q = maxT_q,
    fun_se = fun_se)
}



