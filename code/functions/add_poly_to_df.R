add_poly_to_df <- function(df, poly_object) {
  # add orthonothormal polynomials to df, based on existing poly object.
  # check inputs:
  if(!(is.data.frame(df))) stop("df must be a data frame.")
  if(!("time" %in% names(df))) stop("df requires a 'time' variable.")
  if(!(is.numeric(df$time))) stop("time must be numeric.")
  if(!("poly" %in% class(poly_object))) stop("poly_object must be of class poly, created by the function poly()")
  # extract inputs
  poly_degree <- max(attr(poly_object, "degree"))
  poly_coefs_list <- attr(poly_object, "coefs")
  time <- df$time
  poly_matrix <- poly(x = time, degree = poly_degree, coefs = poly_coefs_list, raw = FALSE)
  colnames(poly_matrix) <- paste0("poly_", seq_len(poly_degree))
  # and return original df plus added poly columns:
  cbind(df, poly_matrix)
}