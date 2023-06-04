pca.fd_pve_cutoff <- function(fd_obj, pve, K_init) {

  # function to do (multivariate) fpca using the pca.fd() function
  # with k chosen based on proportion of variance explained (pve):

  # Checks for inputs: ------------------------------------------------------
  if (!(is.fd(fd_obj) || is.fdPar(fd_obj))) 
    stop("First argument is neither an fd or an fdPar object.")
  if (is.fdPar(fd_obj)) fdobvj = fd_obj$fd
  
  coefs <-fd_obj$coefs
  coefs_dim <- dim(coefs)
  ndim <- length(coefs_dim)
  
  if(!(pve > 0 & pve < 1)) stop("proportion of variance explained (pve) must be between 0 and 1.")
  
  if(!(is.integer(K_init) & K_init > 0)) stop("K_init must be a positive integer")
  

  # Initial PCA: ------------------------------------------------------------
  pca_init <- pca.fd(fdobj = fd_obj, nharm = K_init)
  # find k such that variance explained ~ 99.5%
  K_retain <- min(which(cumsum(pca_init$varprop) > pve))
  
  # Check we actually have used a large enough K_init and if not, then 
  # choose larger value
  if(identical(integer(0), K_retain)) {
    warning("Specififed pve not explained by first K_init FPCs -- choosing K_init bigger")
    if(ndim == 2) {
      K_init <- coefs_dim[1]
    } else if(ndim == 3) {
      K_init <- coefs_dim[1] * coefs_dim[3]
    }
    pca_init <- pca.fd(fdobj = fd_obj, nharm = K_init)
    K_retain <- min(which(cumsum(pca_init$varprop) > pve))
    }
  
  # re-extract the FPCs with new k
  k_seq <- seq_len(K_retain)

  # No only return specified number of components.
  pca_init$harmonics <- pca_init$harmonics[k_seq, ]
  if(ndim == 2) {
    pca_init$scores <- pca_init$scores[, k_seq, drop = FALSE]
  } else if(ndim == 3) {
    pca_init$scores <- pca_init$scores[, k_seq, , drop = FALSE]
  }
  stopifnot(ncol(pca_init$scores) == K_retain)
  pca_init$varprop <-  pca_init$varprop[k_seq]
  list(pca.fd = pca_init, K_retain = K_retain)
}
