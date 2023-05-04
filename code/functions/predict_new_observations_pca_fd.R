predict_new_observations_pca.fd <- function(new_fd_obj, pca_fd_obj) {

  # ------------------------------------------------------------------------#
  # This function predicts the value of new functional observations from an
  # existing pca.fd object. In summary, this function first centers the new 
  # data around the mean used to center the data for the FPCA. Then it projects
  # the centered new data onto the fpcs to obtain the fpc scores. Next, the FPC
  # scores rotated by the basis coefficients of the FPCs, resulting in a 
  # basis representation (in terms of the original basis) of the (centered)
  # fpc-reconstructed functions. Finally, the functions are "de-centered", i.e.,
  # the mean used to center the data for the FPCA is added back on.
  # ------------------------------------------------------------------------#
  
  # Checks for inputs: ------------------------------------------------------
  if (!(is.fd(new_fd_obj) || is.fdPar(new_fd_obj))) 
    stop("First argument is neither an fd or an fdPar object.")
  if (is.fdPar(new_fd_obj)) fdobvj = new_fd_obj$fd
  
  if (!(class(pca_fd_obj) == "pca.fd")) 
    stop("Second argument is not a pca.fd object.")
  
  dimensionality_fd <- dim(new_fd_obj$coefs)
  dimensionality_pca <- dim(pca_fd_obj$harmonics$coefs)
  if(!(length(dimensionality_pca) == length(dimensionality_fd))) {
    stop("Dimensions of first and second arguments don't match.")
  }
  
  if(length(dimensionality_fd) > 2) {
    if(!(all(dimensionality_pca[c(1, 3)] == dimensionality_fd[c(1, 3)]))) {
      stop("Dimensions of first and second arguments don't match.")
    }
  }
  
  if(!(new_fd_obj$basis == pca_fd_obj$harmonics$basis)) {
    stop("For now, functional PCs and data must be represented on same basis.")
  }
  
  

  # Step 1 - Center the functions around the mean  --------------------------
  mean_fd_obj <- pca_fd_obj$meanfd
  new_fd_obj_cent <- center_fd_around_new_mean(fdobj = new_fd_obj,
                                               mean.fd.obj = mean_fd_obj)
  

  # Step 2 - Project the centered data onto the fpcs ------------------------
  scores_new_fd_obj_cent <- project_data_onto_fpcs(fdobj = new_fd_obj_cent, 
                         pca.fd_obj = pca_fd_obj)
  
  scores_new_fd_obj_cent <- apply(scores_new_fd_obj_cent, c(1, 2), sum)
  

  # Step 3 -- Reconstruct centered observations ---------------------------
  if (length(dimensionality_fd) == 2) {
    reconstructions_centered_coef <- t(scores_new_fd_obj_cent %*% t(pca_fd_obj$harmonics$coefs))
  } else {
    nvar <- dimensionality_fd[3]
    nrep <- dimensionality_fd[2]
    nbasis <- dim(pca_fd_obj$harmonics$coefs)[1]
    reconstructions_centered_coef <- array(data = NA, dim = c(nbasis, nrep, nvar))
    for (j in 1:nvar) {
      reconstructions_centered_coef[,,j] <- t(scores_new_fd_obj_cent %*% t(pca_fd_obj$harmonics$coefs[,,j]))
    }
  }
  reconstructed_centered_fd <- fd(coef = reconstructions_centered_coef,
                                  basisobj = pca_fd_obj$harmonics$basis)
  

  # Step 4 - 'de-center' the observations -----------------------------------
  reconstructed_fd <- decenter_fd_around_new_mean(fdobj = reconstructed_centered_fd,
                                                mean.fd.obj = mean_fd_obj)
  
  stopifnot(dim(reconstructed_fd$coefs) == dim(new_fd_obj$coefs))
  

  # Step 6 - Return reconstructed functions and pc scores -------------------
  list(
    reconstructed_fd = reconstructed_fd,
    scores = scores_new_fd_obj_cent
    )
}

