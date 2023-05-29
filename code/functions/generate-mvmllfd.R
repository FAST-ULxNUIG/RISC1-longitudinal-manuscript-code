# Generate simulated multilevel multivariate functional data.
generate_mvmllfd <- function(N,
                             pca_fd_obj,
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
  

  # Add polynomials to df: --------------------------------------------------
  longitudinal_grid <- seq(0, 1, by = 0.01)
  poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)
  design_and_poly_dataset_df <- add_poly_to_df(df = design_dataset_df, poly_object = poly_basis)
  
  # Use design to generate scores from scalar mixed model: ------------------
  scores_dataset_df <- generate_basis_coefficient_matrix(design_df = design_and_poly_dataset_df, 
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


  # Extract Generated Scores: -----------------------------------------------
  scores_df <- scores_dataset_df[, paste0("score_", seq_len(K))]
  stopifnot(nrow(scores_df) == 2 * N * n_i)
  scores_matrix <- as.matrix.data.frame(scores_df)
  

  # Generate functional data object from scores: ----------------------------
  mvmllfd_obj <- construct_fd_from_scores(pca_fd_obj = pca_fd_obj,
                                          scores_matrix = scores_matrix,
                                          K = K)
  # quick check on dimensions of objects being returned:
  stopifnot(ncol(mvmllfd_obj$coefs) == nrow(scores_dataset_df))
  list(df = design_dataset_df, fd_obj = mvmllfd_obj)
  
  }




