library(fda) # CRAN v5.5.1
source("code/functions/calculate_individual_prediction_errors.R")
# Simple test of function:
# on [0, 2 \pi] the true functions are: sin(t), cos(t)
# Reconstructions are: 0.75 sin(t), 0.75 cos(t) 
# It can be shown with simple algebra that the ve
# prediction eerors are
# both = 0.03125
# i.e., check
# https://www.symbolab.com/solver/integral-calculator/%5Cfrac%7B0.25%5E%7B2%7D%5Cint_%7B0%7D%5E%7B2%5Cpi%20%7D%20sin%5E%7B2%7D%5Cleft(x%5Cright)dx%7D%7B2%20%5Cpi%7D?or=input

grid_pts <- seq(0, 2 * pi, length = 100)
y <- cbind(sin(grid_pts), cos(grid_pts))
y_fd <- Data2fd(argvals = grid_pts, y = y)
y_recon <- 0.75 * cbind(sin(grid_pts), cos(grid_pts))
y_recon_fd <- Data2fd(argvals = grid_pts, y = y_recon)

testthat::expect_equal(
  calculate_individual_prediction_errors(pred_fd_obj = y_recon_fd,
                                         true_fd_obj = y_fd,
                                         scaling_factor = 2 * pi),
  c(0.03125, 0.03125),
  tolerance = 10^-7
)

# now let's simply reshape and put into a bivariate fd object with 2 observations
y_array <- y_recon_array <- array(NA, dim = c(length(grid_pts), 2, 2))
y_array[,1,] <- y
y_array[,2,] <- y[,c(2,1)]
y_bfd <- Data2fd(argvals = grid_pts, y = y_array)
y_recon_array[,1, ] <- y_recon 
y_recon_array[,2, ] <- y_recon[,c(2,1)]
y_recon_bfd <- Data2fd(argvals = grid_pts, y = y_recon_array)
testthat::expect_equal(
  calculate_individual_prediction_errors(pred_fd_obj = y_recon_bfd, true_fd_obj = y_bfd, scaling_factor = 2 * pi),
  rep(0.03125, 2) * 2,
  tolerance = 10^-7
)


# Now let's do another small change: --------------------------------------

y_array[,1,2] <- grid_pts
y_recon_array[,1,2] <- 0.5 * grid_pts

y_array[,2,2] <- 1
y_recon_array[,2,2] <- 0.5

y_bfd <- Data2fd(argvals = grid_pts, y = y_array)
y_recon_bfd <- Data2fd(argvals = grid_pts, y = y_recon_array)

calculate_individual_prediction_errors(pred_fd_obj = y_recon_bfd,
                                       true_fd_obj = y_bfd,
                                       scaling_factor = 2 * pi)

# for obs 1: dim 1 error = 0.03125 as before
# obs 1: dim 2 numerator = \int_0^{2\pi}(0.25t)^2 dt = 0.25 * (2 * pi)^3 / 3 - 20.67085
# obs 1: dim 2 denominator = 2\pi
# obs 1 dim 2 answer = 3.289868
# obs 1 answer = 3.289868 + 0.03125 = 3.321118

# for obs 2: dim 1 error = 0.03125 as before
# obs 2: dim 2 numerator = 0.5^2 * 2 \pi
# obs 2: dim 2 denominator = 2\pi
# obs 2 dim 2 answer = 0.25
# obs 2 answer = 0.25 + 0.03125 = 0.28125
testthat::expect_equal(
  calculate_individual_prediction_errors(pred_fd_obj = y_recon_bfd,
                                         true_fd_obj = y_bfd, 
                                         scaling_factor = 2 * pi),
  c(3.321118, 0.28125),
  tolerance = 10^-5 # close enough.
)



    