# ------------------------------------------------------------------------#
# Split the data into test and training sets for analysis.
# We leave out ten strides for each subject at random points along
# the treadmill run!
# ------------------------------------------------------------------------#

# Packages and functions needed -------------------------------------------
library(data.table) # CRAN v1.14.0
library(fda)        # CRAN v5.5.1


# -------------------------------------------------------------------------
outputs_path <- here::here("outputs")
data_objects_list <- readRDS(file = file.path(outputs_path, "data-objects.rds"))
source(file = here::here(
  "code",
  "functions",
  "functions-helper-smoothing.R"
))

# unpack:
mfd_obj <- data_objects_list$mfd_obj
subset_coef_hip <- data_objects_list$subset_coef_hip


# -------------------------------------------------------------------------

bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)

# -------------------------------------------------------------------------

subset_info_dt <- subset_coef_hip
subset_info_dt[, `:=`(
  height_cm_cent = scale(height_cm, center = TRUE, scale = FALSE),
  weight_kg_cent = scale(weight_kg, center = TRUE, scale = FALSE),
  age_cent = scale(age, center = TRUE, scale = FALSE),
  speed_cent = scale(self_selected_speed_kmph, center = TRUE, scale = FALSE))]

covariates_dt <- subset_info_dt[, .(trial_id,
                                    ris = factor(retrospective_injury_status, ordered = F),
                                    sex,
                                    side,
                                    stride_num,
                                    speed_cent,
                                    retrospective_injury_status,
                                    age_cent,
                                    height_cm_cent,
                                    weight_kg_cent,
                                    long_time,
                                    subject_id = factor(subject_id))]



# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Split into test and training sets


# -------------------------------------------------------------------------
test_strides <- covariates_dt[, .(test_strides = sample(unique(trial_id), size = 10)),
                                          by = .(subject_id, side)][, test_strides]

train_inds <- which(!(covariates_dt$trial_id %in% test_strides))
test_inds <- which(covariates_dt$trial_id %in% test_strides)

mfd_obj_train <- mfd_obj[train_inds,, ]
covariates_dt_train <- covariates_dt[train_inds]
stopifnot(ncol(mfd_obj_train$coefs) == covariates_dt_train[, .N])

mfd_obj_test <- mfd_obj[test_inds,, ]
covariates_dt_test <- covariates_dt[test_inds]
stopifnot(ncol(mfd_obj_test$coefs) == covariates_dt_test[, .N])

par(mfrow = c(1, 3))
plot(mfd_obj_test[1:100,,])

test_train_objects_list <- list(
  mfd_obj_train = mfd_obj_train,
  covariates_dt_train = covariates_dt_train,
  mfd_obj_test = mfd_obj_test,
  covariates_dt_test = covariates_dt_test
)

saveRDS(object = test_train_objects_list,
        file = file.path(outputs_path, "test-train-objects.rds"))

