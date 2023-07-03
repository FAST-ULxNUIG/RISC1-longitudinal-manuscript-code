# Function to calculate the integrated squared first derivative.  ---------
# Used to dtermine changes over the longitudinal domain.
calculate_rate_of_change <- function(arg_vals, y_vals, lambda_stable = 10^-10) {
  if(!(is.numeric(arg_vals) & is.numeric(y_vals))) stop("arg_vals and y_vals must be numeric")
  if(!(length(arg_vals) == length(y_vals))) stop("arg_vals and y_vals must be the same length")
  n_i <- length(arg_vals)
  bspline_basis <- create.bspline.basis(rangeval = range(arg_vals), nbasis = n_i, norder = 4)
  fd_par <- fdPar(fdobj = bspline_basis, Lfdobj = 2, lambda = lambda_stable)
  fd_obj <- smooth.basis(argvals = arg_vals, y = y_vals, fdParobj = fd_par)$fd
  deriv_fd_obj <- deriv.fd(expr = fd_obj, Lfdobj = 1)
  result <- inprod(deriv_fd_obj, deriv_fd_obj)
  stopifnot(dim(result) == c(1, 1))
  result[1,1]
}


