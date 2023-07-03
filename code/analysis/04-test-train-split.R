# ------------------------------------------------------------------------#
# Split the data into test and training sets for analysis.
# We leave out ten strides for each subject at random points along
# the treadmill run!
# ------------------------------------------------------------------------#

# Packages and functions needed -------------------------------------------
library(data.table)   # CRAN v1.14.0
library(fda)          # CRAN v5.5.1
library(stringr)      # CRAN v1.4.0
library(modelsummary) # CRAN v1.4.1

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



                      # Create table 1 of characteristics: -------------------------------------
subset_info_dt <- subset_coef_hip
table1_dt <- subset_info_dt[, .(
  ris = retrospective_injury_status[1],
  speed = self_selected_speed_kmph[1],
  sex = sex[1],
  age = age[1],
  weight = weight_kg[1],
  height = height_cm[1]),
  by = subject_id]
table1_dt[, ris := factor(ris, levels = c("never_injured",
                                          "injured_greater_than_2_yr",
                                          "injured_1_to_2_yr", 
                                          "injured_less_than_1_yr"),
                          labels = c("Never Injured",
                                     "Injured $>2$ yr. ago",
                                     "Injured $1-2$ yr. ago",
                                     "Injured $<1$ yr. ago"))]
table1_dt[, sex := factor(sex, levels = c("male", "female"), labels = c("Male", "Female"))]
setnames(table1_dt,
         old = c("ris", "speed", "sex", "age", "weight", "height"),
         new = c("Retrospective Injury Status", "Speed (kmph)", "Sex", "Age (years)", "Weight (kg)", "Height (cm)"))
table1 <- datasummary_balance(~ 1, data = table1_dt[, - c("subject_id")],
                              output = "latex", escape = FALSE)
table1 <- stringr::str_replace(table1,
                               pattern = "& N & Pct.\\\\\\\\\n",
                               replacement = "& N & Pct.\\\\\\\\\n\\\\midrule\n")
table1 <- stringr::str_replace(table1, pattern = "N & Pct.", 
                               replacement = "\\\\textbf{N} & $\\\\mathbf{\\\\mathbf{(\\\\%)}}$")
table1 <- stringr::str_replace(table1, pattern = "Mean & Std. Dev.", 
                               replacement = "\\\\textbf{Mean} & \\\\textbf{Std. Dev.}")
table1 <- stringr::str_replace(table1, pattern = "\\\\end\\{tabular\\}\n\\\\end\\{table\\}",
                               replacement = "\\\\end\\{tabular\\}\\\n\\\\caption\\{Summary characteristics of the dataset used in the analysis.\\}\n\\\\label\\{tab:tab1.\\}\n\\\\end\\{table\\}")
writeLines(text = table1, con = file.path(outputs_path, "tables", "table_1.tex"))



# Scale variables: --------------------------------------------------------
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

stopifnot(length(unique(covariates_dt$trial_id))  == length(test_inds) + length(train_inds))
length(unique(covariates_dt$trial_id)) 
# 47150
length(test_inds)
#  5680
length(train_inds)
#  41470


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

