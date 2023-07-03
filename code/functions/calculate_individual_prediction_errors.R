calculate_individual_prediction_errors <- function(pred_fd_obj, true_fd_obj, scaling_factor = 100) {
  
  # Routine checks of inputs: -----------------------------------------------
  if (!(is.fd(pred_fd_obj) || is.fdPar(pred_fd_obj))) 
    stop("First argument is neither an fd or an fdPar object.")
  if (is.fdPar(pred_fd_obj)) fdobvj = pred_fd_obj$fd
  
  if (!(is.fd(true_fd_obj) || is.fdPar(true_fd_obj))) 
    stop("Second argument is neither an fd or an fdPar object.")
  if (is.fdPar(true_fd_obj)) fdobvj = true_fd_obj$fd
  
  if(!all(dim(true_fd_obj$coefs) == dim(pred_fd_obj$coefs))) {
    stop("Dimensions of first and second arguments do not match.")
  }
  
  dimensionality_fd_obj <- dim(true_fd_obj$coefs)
  nrep <- dimensionality_fd_obj[2]
  ndim <- length(dimensionality_fd_obj)
  # Calculate PC of variance explained: -------------------------------------
  prediction_error <- pred_fd_obj - true_fd_obj
  
  # have to do this separately based on whether fdobj is univariate or
  # multivariate
  
  prediction_error_vec <- vector(mode = "numeric", length = nrep)
  
  for(i in seq_len(nrep)) {
    if(ndim == 2) {
      prediction_error_vec[i] <-c(inprod(prediction_error[i], prediction_error[i]))
    } else if(ndim == 3) {
      prediction_error_vec[i] <- sum(diag(inprod(prediction_error[i], prediction_error[i])))
    }
  }
  
  (1/scaling_factor) * prediction_error_vec
  
}