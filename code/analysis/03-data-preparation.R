# ------------------------------------------------------------------------#
# Prepare data for a multivariate functional longitudinal model:
# ------------------------------------------------------------------------#

# Packages and functions needed -------------------------------------------
library(data.table) # CRAN v1.14.0
library(fda)        # CRAN v5.5.1

# Helper functions:
# Custom functions to turn columns/ rows of matrix into elements of a list
source(file = here::here(
  "code",
  "functions",
  "functions-helper-smoothing.R"
))

outputs_path <- here::here("outputs")


# Import Data -------------------------------------------------------------
# Path to data:
# contains coefficients for registered and unregistered functions
data_path <- here::here(
  "data",
  "risc1-registered-basis-coef.rds"
)

# Read it in:
risc1_basis_coef <- readRDS(file = data_path)



# Define basis to combine with stored coefficients
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)




# Create subset of the data set to work with ------------------------------
# We are working with sagittal plane (i.e., flexion) kinematics
risc1_basis_coef_subset <- risc1_basis_coef[
  subject_id != "P_4092" & # This subject only has strides for one side
    plane_of_motion == "fle" &
    location %in% c("Hip", "Knee", "Ankle") &
    !is.na(retrospective_injury_status)] # Want this covariate non-missing

# How many different subjects?
risc1_basis_coef_subset[, uniqueN(subject_id)]
risc1_basis_coef_subset[subject_id == "P_4019"]
# [1] 288

# Remove subjects who have very sparse measurements at longitudinal level
# for now:
risc1_basis_coef_subset[,.(unique_strides = uniqueN(stride_num)),
                        by = .(subject_id, side)][, hist(unique_strides)]

risc1_basis_coef_subset[subject_id == "P_4019", .(unique_strides = uniqueN(stride_num)),
                        by = .(subject_id, side)]

(subject_ids_remove <- risc1_basis_coef_subset[,.(unique_strides = uniqueN(stride_num)),
                         by = .(subject_id, side)][unique_strides < 20, unique(subject_id)])
# [1] "P_4189" "P_4251" "P_4295" "P_4312"

risc1_basis_coef_subset[,.(unique_strides = uniqueN(stride_num)),
                        by = .(subject_id, side)][rev(order(unique_strides))]

risc1_basis_coef_subset <- risc1_basis_coef_subset[!(subject_id %in% subject_ids_remove)]
# 284

run_duration_dt <- fread("/Users/edwardgunning/Dropbox/Eddie Gunning/thesis-chapt-3/outputs/data/run_duration.csv")
fwrite(run_duration_dt, file = file.path(outputs_path, "run_duration_dt.csv"))

unique_strides_dt <- risc1_basis_coef_subset[, .(unique_strides = uniqueN(stride_num)), by = .(subject_id, side)]
fwrite(unique_strides_dt, file.path(outputs_path, "unique_strides_dt.csv"))




# Calculate the test proportion
risc1_basis_coef_subset[,.(test = 10 / uniqueN(stride_num)),
                        by = .(subject_id, side)][, mean(test)]
risc1_basis_coef_subset[,.(test = 10 / uniqueN(stride_num)),
                        by = .(subject_id, side)][, range(test)]



subset_coef_hip <- risc1_basis_coef_subset[location == "Hip"]
subset_coef_knee <- risc1_basis_coef_subset[location == "Knee"]
subset_coef_ankle <- risc1_basis_coef_subset[location == "Ankle"]

stopifnot(subset_coef_knee$subject_id == subset_coef_hip$subject_id)
stopifnot(subset_coef_knee$side == subset_coef_hip$side)
stopifnot(subset_coef_knee$trial_id== subset_coef_hip$trial_id)
stopifnot(subset_coef_ankle$subject_id == subset_coef_hip$subject_id)
stopifnot(subset_coef_ankle$side == subset_coef_hip$side)
stopifnot(subset_coef_ankle$trial_id== subset_coef_hip$trial_id)

subset_fd_knee <- fd(
  coef = coef_to_mat(subset_coef_knee[, paste0("lm_coef_", 1:80)]),
  basisobj = bspl80)

subset_fd_hip <- fd(
  coef = coef_to_mat(subset_coef_hip[, paste0("lm_coef_", 1:80)]),
  basisobj = bspl80)

subset_fd_ankle <- fd(
  coef = coef_to_mat(subset_coef_ankle[, paste0("lm_coef_", 1:80)]),
  basisobj = bspl80)

sample_for_plot <- sample(seq_len(ncol(subset_fd_hip$coefs)), size = 500, replace = FALSE)
par(mfrow = c(1, 3))
plot(subset_fd_hip[sample_for_plot, ])
title("Hip")
plot(subset_fd_knee[sample_for_plot, ])
title("Knee")
plot(subset_fd_ankle[sample_for_plot, ])
title("Ankle")

# -------------------------------------------------------------------------
# Make multivariate FDA object:
# -------------------------------------------------------------------------
coef_array <- array(data = NA, dim = c(dim(subset_fd_hip$coefs), 3))
coef_array[,, 1] <- subset_fd_hip$coefs
coef_array[,, 2] <- subset_fd_knee$coefs
coef_array[,, 3] <- subset_fd_ankle$coefs
mfd_obj <- fd(coef = coef_array, basisobj = bspl80)

objects_list <- list(
  mfd_obj = mfd_obj,
  subset_coef_hip = subset_coef_hip
)
saveRDS(object = objects_list, file = file.path(outputs_path, "data-objects.rds"))
