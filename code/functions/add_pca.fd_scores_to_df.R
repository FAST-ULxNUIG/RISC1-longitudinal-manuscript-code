add_scores_to_df <- function(pca.fd_obj, K_retain, df) {
  
  if (!(class(pca.fd_obj) == "pca.fd")) 
    stop("First argument is not a pca.fd object.")
  coefs <- pca.fd_obj$harmonics$coefs
  coefs_dim <- dim(coefs)
  ndim <- length(coefs_dim)
  
  if(!(is.integer(K_retain) & K_retain > 0)) stop("K_retain must be a positive integer")
  if(!(is.data.frame(df))) stop("df must be a data.frame()")
  if(!(nrow(pca.fd_obj$scores) == nrow(df))) stop("Number of observations in pca.fd_obj and df must match.")
  N <- nrow(df)
  if(!(ncol(pca.fd_obj$scores) == K_retain)) stop("Number of FPCs in pca.fd_obj and value of K_retain supplied must match.")
  
  # Extract pca scores from object:
  # if multivariate object, need to sum over dimensions
  if(ndim == 2) {
    scores <- pca.fd_obj$scores
  } else if(ndim == 3) {
    scores <- apply(pca.fd_obj$scores, c(1, 2), sum)
  }
  stopifnot(dim(scores) == c(N, K_retain))
  colnames(scores) <- paste0("score_", seq_len(K_retain))
  cbind(df, scores)
}








