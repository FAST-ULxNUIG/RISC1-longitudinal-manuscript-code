source(here::here("code", "functions", "generate_design.R"))
library(data.table) # CRAN v1.14.2
set.seed(1996)
test <- data.table(generate_design_multiple_subjects(N = 100000, n_i = 10, speed_sd = 1))

# Check every subject has 10 strides bilaterally:
stopifnot(test[, uniqueN(stride_ind), by = .(subject_id, side)][, all(V1 == 10)])
# Check that sex and speed don't vary within subjects:
stopifnot(test[, uniqueN(sex) ==1, by = .(subject_id)][, all(V1)])
stopifnot(test[, uniqueN(speed_cent) ==1, by = .(subject_id)][, all(V1)])
# Check that sex has the correct levels:
testthat::expect_equal(levels(test$sex), c("male", "female"))


# And check that distributions of variables are approximately correct:
testthat::expect_equal(mean(test$sex == "male"), 0.5, tolerance = 0.01)
testthat::expect_equal(test[, .(speed_cent = speed_cent[1]), by = subject_id][, sd(speed_cent)],
                       1,
                       tolerance = 0.01)
testthat::expect_equal(mean(test$speed_cent), 0, tolerance = 0.01)

# Check if it gives error of replicate strides not being produced for each subject:
testthat::expect_error(object = generate_design_multiple_subjects(N = 100, n_i = 1, speed_sd = 2))
# and if a negative standard deviation is being produced
testthat::expect_error(object = generate_design_multiple_subjects(N = 100, n_i = 10, speed_sd = -2))

