# Generate simulated multilevel multivariate functional data.
generate_simulated_mvmllfd <- function(N,
                                       n_i,
                                       speed_sd,
                                       K = 10,
                                       beta_poly_1_vec,
                                       beta_poly_2_vec,
                                       beta_poly_3_vec,
                                       beta_sex_vec,
                                       beta_speed_cent_vec,
                                       Q_star_array,
                                       R_star_array,
                                       s_vec) {

  # Routine input checks: ---------------------------------------------------
  stopifnot(N >= 2)
  stopifnot(n_i >= 2)
  stopifnot(speed_sd > 0)
  
  # Generate design data.frame() --------------------------------------------
  design_dataset_df <- generate_design_multiple_subjects(N = N, n_i = n_i, speed_sd = speed_sd)
  

  # Use design to generate scores from scalar mixed model: ------------------
  scores_dataset_df <- generate_basis_coefficient_matrix(design_df = design_dataset_df, 
                                                       N = N,
                                                       K = K,
                                                       beta_poly_1_vec = beta_poly_1_vec, 
                                                       beta_poly_2_vec = beta_poly_2_vec,
                                                       beta_poly_3_vec = beta_poly_3_vec,
                                                       beta_sex_vec = beta_sex_vec,
                                                       beta_speed_cent_vec = beta_speed_cent_vec,
                                                       Q_star_array = Q_star_array,
                                                       R_star_array = R_star_array, 
                                                       s_vec = s_vec)
  
  
  
}