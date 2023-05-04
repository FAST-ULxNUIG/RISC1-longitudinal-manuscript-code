library(fda) # CRAN v5.5.1
source(file = "code/functions/center_fd_around_new_mean.R")
grid_pts <- seq(0, 2 * pi, length = 100)
y <- sin(grid_pts)
y_mean <- cos(grid_pts)
y_center <- sin(grid_pts) - cos(grid_pts)
y_fd <- Data2fd(argvals = grid_pts, y = y)
y_mean_fd <- Data2fd(argvals = grid_pts, y = y_mean)
y_center_fd <- Data2fd(argvals = grid_pts, y = y_center)
y_center_fd$fdnames[[3]] <- "Centered value"
# Centering using function:
y_center_using_func_fd <- center_fd_around_new_mean(fdobj = y_fd, mean.fd.obj = y_mean_fd)
# Test:
testthat::expect_equal(y_center_using_func_fd, y_center_fd)
