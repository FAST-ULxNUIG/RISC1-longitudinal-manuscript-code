generate_polynomial_model_basis_coefficient <- function(design_df,
                                                        N,
                                                        beta_poly_1_k,
                                                        beta_poly_2_k,
                                                        beta_poly_3_k,
                                                        beta_sex_k,
                                                        beta_speed_cent_k,
                                                        Q_star_k,
                                                        R_star_k,
                                                        s_k) {
  require(lme4)
  require(mvtnorm)
  # Check inputs of function:
  if(!setequal(x = names(design_df),
               y = c("subject_id", "side", "time", "stride_ind", "poly_1", "poly_2", "sex", "speed_cent"))) {
    stop("design_df must contain the following variables:
         subject_id, side, time, stride_ind, poly_1, poly_2, sex, speed_cent")
  }
  
  if(!(
    is.numeric(beta_poly_1_k) & is.numeric(beta_poly_2_k) &
    is.numeric(beta_poly_3_k) & is.numeric(beta_sex_k) &
    is.numeric(beta_speed_cent_k)
  )) stop("All fixed-effects basis coefficients must be numeric")
  # if they are, make them into vector:
  beta_vec <- c(beta_poly_1_k, beta_poly_2_k, beta_poly_3_k, beta_sex_k, beta_speed_cent_k)  
  
  if(!(is.numeric(Q_star_k) & all(dim(Q_star_k) == c(3, 3)))) {
    stop("Q_star_k must be a 3 x 3 matrix.")
  }
  
  if(!(matrixcalc::is.positive.definite(Q_star_k))) {
    stop("Q_star_k must be positive definite.")
  }
  
  if(!(is.numeric(R_star_k) & all(dim(R_star_k) == c(3, 3)))) {
    stop("R_star_k must be a 3 x 3 matrix.")
  }
  
  if(!(matrixcalc::is.positive.definite(R_star_k))) {
    stop("R_star_k must be positive definite.")
  }
  
  if(s_k <= 0) stop("s_k must be positive.")
  
  if(!(length(unique(design_df$subject_id)) == N)) stop("N must equal the number of unique subjects in df!")

  N_total <- nrow(design_df)
  # Start: Generate lme model object. ---------------------------------------
  design_df$score_k <- 1
  l_formula <- lFormula(
    formula = score_k ~ poly_1 + poly_2 + sex + speed_cent +
      (poly_1 + poly_2||subject_id) + 
      (poly_1 + poly_2||subject_id:side),
    data = design_df)
  
  # Extract fixed-effects design matrix: ------------------------------------
  # and create XB fixed effects contribution:
  X <- l_formula$X
  XB <- X %*% matrix(beta_vec, nrow = 5, ncol = 1)
  

  # Extract random effects design matrix: -----------------------------------
  Zt_list <- l_formula$reTrms$Ztlist

  ## Extract subject-level parameters first: --------------------------------
  Z_subject_xi_1 <- t(Zt_list[["1 | subject_id"]])
  Z_subject_xi_2 <- t(Zt_list[["0 + poly_1 | subject_id"]])
  Z_subject_xi_3 <- t(Zt_list[["0 + poly_2 | subject_id"]])
  # routine check:
  stopifnot(ncol(Z_subject_xi_1) == N)
  stopifnot(ncol(Z_subject_xi_2) == N)
  stopifnot(ncol(Z_subject_xi_3) == N)
  # Draw subject-level random effects:
  U <- rmvnorm(n = N, sigma = Q_star_k)
  U1 <- U[,1, drop = FALSE]
  U2 <- U[,2, drop = FALSE]
  U3 <- U[,3, drop = FALSE]
  ZU_subject <- Z_subject_xi_1 %*% U1 + Z_subject_xi_2 %*% U2 + Z_subject_xi_3 %*% U3
  

  ## Now subject and side-level params -------------------------------------
  Z_subject_side_xi_1 <- t(Zt_list[["1 | subject_id:side"]])
  Z_subject_side_xi_2 <- t(Zt_list[["0 + poly_1 | subject_id:side"]])
  Z_subject_side_xi_3 <- t(Zt_list[["0 + poly_2 | subject_id:side"]])
  # routine check:
  stopifnot(ncol(Z_subject_side_xi_1) == 2 * N)
  stopifnot(ncol(Z_subject_side_xi_2) == 2 * N)
  stopifnot(ncol(Z_subject_side_xi_3) == 2 * N)
  # Draw subject-level random effects:
  V <- rmvnorm(n = 2 * N, sigma = R_star_k)
  V1 <- V[,1, drop = FALSE]
  V2 <- V[,2, drop = FALSE]
  V3 <- V[,3, drop = FALSE]
  ZV_subject_side <- Z_subject_side_xi_1 %*% V1 + Z_subject_side_xi_2 %*% V2 + Z_subject_side_xi_3 %*% V3
  
  # Routine check:
  stopifnot(all(dim(ZU_subject) == dim(ZV_subject_side)))
  ZU_total <- ZU_subject + ZV_subject_side
  
  # Next, random error: -----------------------------------------------------
  var_epsilon_star_k <- rnorm(n = N_total, mean = 0, sd = sqrt(s_k))
  
  # Put everything together: ------------------------------------------------
  stopifnot(dim(XB) == c(N_total, 1))
  stopifnot(dim(ZU_total) == c(N_total, 1))
  XB <- XB[,1, drop = TRUE]
  ZU_total <- ZU_total[,1, drop = TRUE]
  
  score_k <- XB + ZU_total + var_epsilon_star_k

  stopifnot(is.numeric(score_k) & length(score_k) == N_total)
  
  score_k
}


