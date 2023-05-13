source(here::here("code", "functions", "add_poly_to_df.R"))
source(here::here("code", "functions", "generate_design.R"))

# Generate data anad add polynomials:
set.seed(1)
test <- generate_design_multiple_subjects(N = 2, n_i = 10, speed_sd = 1)
longitudinal_grid <- seq(0, 1, by = 0.01)
poly_basis <- poly(longitudinal_grid, degree = 2, raw = FALSE)
test_add_poly <- add_poly_to_df(df = test, poly_object = poly_basis)

# Check it doesn't change anything in df:
testthat::expect_equal(object = test,
                       test_add_poly[, (names(test_add_poly) != c("poly_1", "poly_2"))])


# Check that polynomials are calculated correctly:
par(mfrow = c(1, 2))
plot(poly_1 ~ time, data = test_add_poly)
lines(longitudinal_grid, poly_basis[,1])
plot(poly_2 ~ time, data = test_add_poly)
lines(longitudinal_grid, poly_basis[,2])
dev.off()

# Check if it it gives error if time is not numeric:
test$time <- factor(test_add_poly$time)
testthat::expect_error(add_poly_to_df(df = test, poly_object = poly_basis))



