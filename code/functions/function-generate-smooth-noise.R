generate_smooth_noise <- function(n, sigma, l) {
  # function to generate smooth noise from a Gaussian kernel on [0, 100]
  C <- function(t, tprime) {
    sigma^2 * dnorm(x = l * (t-tprime))
  }
  C_grid <- outer(X = 0:100, Y = 0:100, FUN = C)
  rmvnorm(n = n, sigma = C_grid)
}

add_smooth_noise_to_fd_obj <- function(fd_obj, sigma, l) {
    if (!(is.fd(fd_obj) || is.fdPar(fd_obj))) 
      stop("First argument is neither an fd or an fdPar object.")
    if (is.fdPar(fd_obj)) fdobvj = fd_obj$fd
    
    coef     <- as.array(fd_obj$coefs)
    coefd    <- dim(coef)
    nrep     <- coefd[2]
    ndim     <- length(coefd)
    basisobj <- fd_obj$basis
    nbasis   <- basisobj$nbasis
    
    if(!all(basisobj$rangeval == c(0, 100))) stop("Currently only implemented for functional data on [0, 100]")
    
    if(ndim == 2) {
      epsilon <- t(generate_smooth_noise(n = nrep, sigma = sigma, l = l))
    } else if(ndim == 3) {
      nvar <- coefd[3]
      epsilon <- array(NA, dim = c(101, nrep, nvar))
      for(j in seq_len(nvar)) {
        epsilon[,, j] <- t(generate_smooth_noise(n = nrep, sigma = sigma, l = l))
      }
    }
    epsilon_fd_obj <- Data2fd(argvals = 0:100, y = epsilon, basisobj = basisobj)
    
    fd_obj + epsilon_fd_obj
    
}
