construct_fd_from_scores <- function(pca_fd_obj, scores_matrix, K, add_back_mean = TRUE) {
  # needs fda package
  require(fda)
  
  # Routine checks for input:
  if (!(class(pca_fd_obj) == "pca.fd")) 
    stop("First argument is not a pca.fd object.")
  dimensionality_pca <- dim(pca_fd_obj$harmonics$coefs)
  mean_fd <- pca_fd_obj$meanfd
  if(dimensionality_pca[2] != K) {
    stop("Number of FPCs in pca_fd_obj is not equal to K.")
  }
  if(!(is.matrix(scores_matrix) & is.numeric(scores_matrix))) {
    stop("scores_matrix must be a numeric matrix.")
  }
  if(ncol(scores_matrix) != K) {
    stop("Number of columns in scores_matrix and K do not match.")
  }
  nrep <- nrow(scores_matrix)
  
  # -------------------------------------------------------------------------
  
  if (length(dimensionality_pca) == 2) {
    reconstructions_centered_coef <- t(scores_matrix %*% t(pca_fd_obj$harmonics$coefs))
  } else {
    nvar <- dimensionality_pca[3]
    nbasis <- dim(pca_fd_obj$harmonics$coefs)[1]
    reconstructions_centered_coef <- array(data = NA, dim = c(nbasis, nrep, nvar))
    for (j in 1:nvar) {
      reconstructions_centered_coef[,,j] <- t(scores_matrix %*% t(pca_fd_obj$harmonics$coefs[,,j]))
    }
  }
  reconstructed_fd <- fd(coef = reconstructions_centered_coef,
                                  basisobj = pca_fd_obj$harmonics$basis)
  
  if(add_back_mean) {
    # Step 4 - 'de-center' the observations -----------------------------------
    mean_fd_obj <- pca_fd_obj$meanfd
    reconstructed_fd <- decenter_fd_around_new_mean(fdobj = reconstructed_fd,
                                                    mean.fd.obj = mean_fd_obj)
    
  }
  
  # Some routine checks on the output:
  stopifnot(ncol(reconstructed_fd$coefs) == nrow(scores_matrix))
  stopifnot(nrow(reconstructed_fd$coefs) == dimensionality_pca[1])
  stopifnot(length(dim(reconstructed_fd$coefs)) == length(dimensionality_pca))
  
  
  # return:
  reconstructed_fd
}



