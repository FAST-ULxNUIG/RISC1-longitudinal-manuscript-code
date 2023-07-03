# ------------------------------------------------------------------------#
#  Some simple unit tests for calculating rate of change function.
# -----------------------------------------------------------------------#
source(here::here("code", "functions", "calculate_rate_of_change.R"))
library(testthat) # CRAN v3.1.4
library(fda)      # CRAN v5.5.1

t_grid <- seq(0, 1, length.out = 101)
y_test <- rep(1, length(t_grid))
testthat::expect_equal(calculate_rate_of_change(arg_vals = t_grid, y_vals = y_test),
                       0, tolerance = 10^-10)

y_test <- t_grid
testthat::expect_equal(calculate_rate_of_change(arg_vals = t_grid, y_vals = y_test),
                       1, tolerance = 10^-10)

y_test <- t_grid^2
testthat::expect_equal(calculate_rate_of_change(arg_vals = t_grid, y_vals = y_test, lambda_stable = 10^-12), 
                       4/3, 
                       tolerance = 10^-3) # this was best toilerance possible                   


y_test <- t_grid^3
testthat::expect_equal(calculate_rate_of_change(arg_vals = t_grid, y_vals = y_test, lambda_stable = 10^-12), 
                       9/5, 
                       tolerance = 10^-2) # again, tolerance less

# check it gives error if arguments are of different length:
testthat::expect_error(calculate_rate_of_change(arg_vals = t_grid, y_vals = 0:99))
