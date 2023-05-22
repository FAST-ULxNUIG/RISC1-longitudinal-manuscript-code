generate_basis_coefficient_matrix <- function(design_df,
                                              N,
                                              K = 10,
                                              beta_poly_1_vec,
                                              beta_poly_2_vec,
                                              beta_poly_3_vec,
                                              beta_sex_vec,
                                              beta_speed_cent_vec,
                                              Q_star_array,
                                              R_star_array,
                                              s_vec) {
  # Routine checks for inputs: ----------------------------------------------
  stopifnot(K >= 1)
  if(!setequal(x = names(design_df),
               y = c("subject_id", "side", "time", "stride_ind", "poly_1", "poly_2", "sex", "speed_cent"))) {
    stop("design_df must contain the following variables:
         subject_id, side, time, stride_ind, poly_1, poly_2, sex, speed_cent")
  }
  if(!(
    is.numeric(beta_poly_1_vec) & is.numeric(beta_poly_2_vec) &
    is.numeric(beta_poly_3_vec) & is.numeric(beta_sex_vec) &
    is.numeric(beta_speed_cent_vec)
  )) stop("All fixed-effects basis coefficients must be numeric")
  if(!is.numeric(s_vec)) stop("s_vec must be numeric")
  
  if(!(
    (length(beta_poly_1_vec) == K) & (length(beta_poly_2_vec) == K) &
    (length(beta_poly_3_vec) == K) & (length(beta_sex_vec) == K) &
    (length(beta_speed_cent_vec) == K)
  )) stop("All fixed-effects basis coefficients vectors must be of length K")
  if(!(length(s_vec) == k)) stop("s_vec must be of length K")
  
  if(!is.array(Q_star_array)) stop("Q_star_array must be an array")
  if(!all(dim(Q_star_array) == c(3, 3, K))) stop("Q_star_array must be of dimensions 3 x 3 x K")
  
  if(!is.array(R_star_array)) stop("R_star_array must be an array")
  if(!all(dim(R_star_array) == c(3, 3, K))) stop("R_star_array must be of dimensions 3 x 3 x K")
  

  # Now, start data generation  --------------------------------------------
  design_df_new <- design_df
  for(k in seq_len(K)) {
    score_name <- paste0("score_", k)
    design_df_new[, score_name] <- 
      generate_polynomial_model_basis_coefficient(
      design_df = design_df,
      N = N,
      beta_poly_1_k = beta_poly_1_vec[k],
      beta_poly_2_k = beta_poly_2_vec[k],
      beta_poly_3_k = beta_poly_3_vec[k],
      beta_sex_k = beta_sex_vec[k],
      beta_speed_cent_k = beta_speed_cent_vec[k],
      Q_star_k = Q_star_array[,,k],
      R_star_k = R_star_array[,,k],
      s_k = s_vec[k]
     )
  }
  # Check dimensions of outputted object: -----------------------------------
  stopifnot(dim(design_df_new) == c(N, ncol(design_df) + K))

  # Return: -----------------------------------------------------------------
  design_df_new
}



