variance_explained_reconstruction <- function(true_fd_obj, reconstructed_fd_obj) {

  # Routine checks of inputs: -----------------------------------------------
  if (!(is.fd(true_fd_obj) || is.fdPar(true_fd_obj))) 
    stop("First argument is neither an fd or an fdPar object.")
  if (is.fdPar(true_fd_obj)) fdobvj = true_fd_obj$fd
  
  if (!(is.fd(reconstructed_fd_obj) || is.fdPar(reconstructed_fd_obj))) 
    stop("Second argument is neither an fd or an fdPar object.")
  if (is.fdPar(reconstructed_fd_obj)) fdobvj = reconstructed_fd_obj$fd
  
  if(!all(dim(true_fd_obj$coefs) == dim(reconstructed_fd_obj$coefs))) {
    stop("Dimensions of first and second arguments do not match.")
  }

  dimensionality_fd_obj <- dim(true_fd_obj$coefs)
  
  # Calculate PC of variance explained: -------------------------------------
  reconstruction_error <- reconstructed_fd_obj - true_fd_obj
  
  # have to do this separately based on whether fdobj is univariate or
  # multivariate
  
  # Numerator:
  if(length(dimensionality_fd_obj) == 2) {
    reconstruction_error_inprod <- sum(diag(inprod(reconstruction_error, reconstruction_error)))
  } else if(length(dimensionality_fd_obj) == 3) {
    reconstruction_error_inprod <- c(inprod(reconstruction_error, reconstruction_error))
  }

  # Total Var
  if(length(dimensionality_fd_obj) == 2) {
    total_var_inprod <- sum(diag(inprod(center.fd(true_fd_obj), center.fd(true_fd_obj))))
  } else if(length(dimensionality_fd_obj) == 3) {
    total_var_inprod <- c(inprod(center.fd(true_fd_obj), center.fd(true_fd_obj)))
  }
  
  vne <- reconstruction_error_inprod / total_var_inprod # variance not explained
  stopifnot(vne >= 0 & vne <= 1) # routine check
  ve <- 100 * (1 - vne) # convert to variance explained and to a percentage
  ve
  }


