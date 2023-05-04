library(fda) # CRAN v5.5.1
source("code/functions/variance_explained_reconstruction.R")
# Simple test of function:
# on [0, 2 \pi] the true functions are: sin(t), cos(t)
# Reconstructions are: 0.75 sin(t), 0.75 cos(t) 
# It can be shown with simple algebra that the ve
# reconstruction is
# is = 100 * (1-(0.0625 * 2pi) / pi) = 87.5
grid_pts <- seq(0, 2 * pi, length = 100)
y <- cbind(sin(grid_pts), cos(grid_pts))
y_fd <- Data2fd(argvals = grid_pts, y = y)
y_recon <- 0.75 * cbind(sin(grid_pts), cos(grid_pts))
y_recon_fd <- Data2fd(argvals = grid_pts, y = y_recon)

testthat::expect_equal(
  variance_explained_reconstruction(true_fd_obj = y_fd, reconstructed_fd_obj = y_recon_fd),
  87.5,
  tolerance = 10^-6 # really a tolerance of 10^-8 on calculations because pve = (100 * prop ve)
)


