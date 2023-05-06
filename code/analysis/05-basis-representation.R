# ------------------------------------------------------------------------#
# Perform the mv-FPCA decomposition of the multivariate functional data.
# We use the `pca.fd()` function from the 'fda' package to calculate the 
# mv-FPCs from the univariate basis representations.
# ------------------------------------------------------------------------#

# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1 
library(data.table) # CRAN v1.14.2 
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1 
library(xtable)     # CRAN v1.8-4 

# Load all helper functions: ----------------------------------------------
# all custom written functions for analysis:
source(here::here("code/functions/center_fd_around_new_mean.R"))
source(here::here("code/functions/decenter_fd_around_new_mean.R"))
source(here::here("code/functions/project_data_onto_fpcs.R"))
source(here::here("code/functions/predict_new_observations_pca_fd.R"))
source(here::here("code/functions/variance_explained_reconstruction.R"))
source(here::here("code/functions/theme_gunning.R"))


# Some settings for the plots: --------------------------------------------
theme_gunning() # use default theme
theme_update(
  axis.text = element_text(size = 10.5),
  axis.title = element_text(size = 10.5),
  strip.text = element_text(size = 10.5),
  plot.title = element_text(size = 11.5),
  legend.text = element_text(size = 10.5)
)
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Path to save the outputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
outputs_path <- here::here("outputs")

# Read in data - the test and training fd objects: ------------------------
test_train_objects <- readRDS(file = file.path(outputs_path,
                                               "test-train-objects.rds"))
# and unpack:
mfd_obj_train <- test_train_objects$mfd_obj_train
covariates_dt_train <- test_train_objects$covariates_dt_train
mfd_obj_test <- test_train_objects$mfd_obj_test
covariates_dt_test <- test_train_objects$covariates_dt_test

# Caclulate total amount of test and training data:
N_train <- ncol(mfd_obj_train$coefs)
N_test <- ncol(mfd_obj_test$coefs)
stopifnot(N_train == covariates_dt_train[, .N])
stopifnot(N_test == covariates_dt_test[, .N])

# plot 3 samples random samples of size 5 for a test:
par(mfrow = c(1, 3))
plot(mfd_obj_train[sample(seq_len(N_train), size = 5, replace = FALSE)])
plot(mfd_obj_train[sample(seq_len(N_train), size = 5, replace = FALSE)])
plot(mfd_obj_train[sample(seq_len(N_train), size = 5, replace = FALSE)])

# Perform Multivariate Functional PCA (mv-FPCA) ---------------------------
# do mv-PCA:
mfpca_time <- system.time(
 {
   # do initial pca, retaining a large k
   mfpca_init <- pca.fd(fdobj = mfd_obj_train,
                       nharm = 100,
                       harmfdPar = fdPar(fdobj = mfd_obj_train$basis, 
                                         Lfdobj = 2,
                                         lambda = 10^-12))
  # find k such that variance explained ~ 99.5%
  k_retain <- min(which(cumsum(mfpca_init$varprop) > 0.995))
  # re-extract the FPCs with new k
  mfpca <- pca.fd(fdobj = mfd_obj_train,
                  nharm = k_retain,
                  harmfdPar = fdPar(fdobj = mfd_obj_train$basis,
                                    Lfdobj = 2, lambda = 10^-12))
 #calculate the mv-fpc scores by summing over the dimensions
 scores <- apply(mfpca$scores, c(1,2), sum)
 mean_fd <- mfpca$meanfd # extract the mean function used to centre the data
 }
)

                    # Create Scree Plot and Variance Explained Plot: --------------------------
# Use ggplot() to make publishable plot, involves reshaping data:
mfpc_dt <- data.table(
  fpc_number = seq_len(k_retain),
  eigen_value = mfpca$values[seq_len(k_retain)],
  var_pc = 100 * mfpca$varprop[seq_len(k_retain)])
mfpc_dt[, var_pc_cum := cumsum(var_pc)]  

p1 <- ggplot(data = mfpc_dt) +
  aes(x = fpc_number, y = eigen_value) +
  geom_line() +
  theme(legend.title = element_blank()) +
  geom_point(size = 0.3) +
  labs(y = "Eigenvalue ($\\lambda_k$)",
       x= "FPC Number ($k$)",
       title = "\\textbf{(a)} Scree Plot")

p2 <- ggplot(data = mfpc_dt) +
  aes(x = fpc_number, y = var_pc_cum) +
  geom_line() +
  theme(legend.title = element_blank()) +
  geom_point(size = 0.3) +
  labs(x= "FPC Number ($k$)",
       y = "$(\\%)$ of Variance Explained",
       title = "\\textbf{(b)} Explained Variance")  


(variance_explained_plot <- ggpubr::ggarrange(plotlist = list(p1, p2),
                                              ncol = 2, nrow = 1,
                                              common.legend = TRUE, 
                                              legend = "bottom")) 



# Reconstruct the test set observations using the mv-FPCA: ---------------
# Use custom function to do this
# (see  predict_new_observations_pca_fd.R)
reconstructions_test <- predict_new_observations_pca.fd(
  new_fd_obj = mfd_obj_test, 
  pca_fd_obj = mfpca)
# Extract the reconstructed fd object:
mfd_obj_reconstructions_test <- reconstructions_test$reconstructed_fd


# Quick plot:
set.seed(1996)
test_sample <- sample(1:N_test, size = 5)
par(mfrow = c(1, 1))
plot(mfd_obj_test[test_sample,1], lty = 2)
lines(mfd_obj_reconstructions_test[test_sample, 1], lty = 3)
plot(mfd_obj_test[test_sample,2], lty = 2)
lines(mfd_obj_reconstructions_test[test_sample, 2], lty = 3)
plot(mfd_obj_test[test_sample,3], lty = 2)
lines(mfd_obj_reconstructions_test[test_sample, 3], lty = 3)




# Plot reconstruction: ----------------------------------------------------
reconstruction_eval <- eval.fd(evalarg = 0:100, fdobj = mfd_obj_reconstructions_test[test_sample])
reconstruction_sample <- rbind(
  reconstruction_eval[,,1],
  reconstruction_eval[,,2],
  reconstruction_eval[,,3]
)

true_eval <- eval.fd(evalarg = 0:100, fdobj = mfd_obj_test[test_sample])
true_sample <- rbind(
  true_eval[,,1],
  true_eval[,,2],
  true_eval[,,3]
)

colnames(true_sample) <- colnames(reconstruction_sample) <- paste0("test_rep_", test_sample)

reconstruction_dt <- data.table(
  t = rep(0:100, 3),
  dimension = rep(c("Hip", "Knee", "Ankle"), each = 101),
  reconstruction_sample
)

truth_dt <- data.table(
  t = rep(0:100, 3),
  dimension = rep(c("Hip", "Knee",  "Ankle"), each = 101),
  true_sample
)

reconstruction_dt$reconstruct <- "Reconstructed"
truth_dt$reconstruct <- "True"

plot_reconstruction_dt <- rbind(
  reconstruction_dt, truth_dt
)

plot_reconstruction_dt_lng <- melt.data.table(
  data = plot_reconstruction_dt, 
  id.vars = c("t", "dimension", "reconstruct"),
  value.name = "angle",
  variable.name = "rep",
  measure.vars =  paste0("test_rep_", test_sample),
  variable.factor = FALSE, 
  value.factor = FALSE, 
  verbose = TRUE
)
plot_reconstruction_dt_lng[, rep:= factor(stringr::str_remove(rep, "rep"))]

plot_reconstruction_dt_lng[, dimension := factor(dimension,
                                                levels = c("Hip", "Knee", "Ankle"))]

(reconstruction_plot <- ggplot(data = plot_reconstruction_dt_lng) +
    aes(x = t,
        y = angle,
        colour = rep, 
        lty = reconstruct,
        alpha = reconstruct,
        group = interaction(reconstruct, rep)) +
    facet_wrap( ~ dimension, scales = "free_y") +
    scale_alpha_manual(values = c(0.75,1)) +
    geom_line(linewidth = 0.5) +
    labs(x = "Normalised Time ($\\%$ of Stride)", y = "Angle ($^{\\circ}$)",
         title = "\\textbf{(c)} Reconstruction of Test-Set Observations") +
    guides(colour ="none") +
    theme(legend.title = element_blank(),
          legend.margin = margin(b = -1),
          legend.position = "bottom"))



# Combine the reconstruction and variance explained plots: ----------------

# increase widespace to the sides of VE plot:
variance_explained_plot_narrow <- variance_explained_plot +
  theme(plot.margin = margin(0, 90, 0, 90))

# And combine two plots, stacking them:
basis_transform_plot <- ggpubr::ggarrange(variance_explained_plot_narrow,
                                          reconstruction_plot, 
                                          nrow = 2, 
                                          ncol = 1, 
                                          heights = c(0.445, 0.555))

tikz(file.path(plots_path, "basis-transform-plot.tex"),
     width = (1.25) * doc_width_inches, 
     height = 0.9375 * doc_width_inches)
print(basis_transform_plot)
dev.off()


# Cross-validation approaches for reconstruction error: -------------------
# (loso-cv can be time consuming!)
subject_id_dt_train <- covariates_dt_train[, .(subject_id = unique(subject_id))]
subject_id_dt_train[, ind := .I]
subject_id_dt_train[, fold := cut(ind, breaks = 10, labels = FALSE)]

ve_cv <- vector("numeric", length = 10)

for(k in 1:10) {
  print(paste0("fold", k))
  subject_ids_out <- subject_id_dt_train[fold == k, subject_id]
  subject_ids_in <- subject_id_dt_train[fold != k, subject_id]
  
  inds_out <- which(covariates_dt_train$subject_id %in% subject_ids_out)
  inds_in <- which(covariates_dt_train$subject_id %in% subject_ids_in)
  
  mfpca_in <- pca.fd(fdobj = mfd_obj_train[inds_in],
                     nharm = k_retain,
                     harmfdPar = fdPar(fdobj = mfd_obj_train$basis,
                                       Lfdobj = 2, lambda = 10^-12))
  predictions_out <- predict_new_observations_pca.fd(
    new_fd_obj = mfd_obj_train[inds_out],
    pca_fd_obj = mfpca_in)
  
  reconstructions_out <- predictions_out$reconstructed_fd
  
  ve_cv[k] <- variance_explained_reconstruction(
    true_fd_obj = mfd_obj_train[inds_out],
    reconstructed_fd_obj = reconstructions_out)
}

ten_fold_cv_dt <- data.table(
  fold = 1:10,
  ve_ten_fold_cv = ve_cv
) 


ten_fold_cv_dt[, mean(ve_ten_fold_cv)]

# -------------------------------------------------------------------------
N_sub <- uniqueN(subject_id_dt_train$subject_id)
ve_loso_cv <- vector("numeric", length = N_sub)
for(k in 1:N_sub) {
  print(paste0("fold", k))
  subject_ids_out <- subject_id_dt_train[ind == k, subject_id]
  subject_ids_in <- subject_id_dt_train[ind != k, subject_id]
  
  inds_out <- which(covariates_dt_train$subject_id %in% subject_ids_out)
  inds_in <- which(covariates_dt_train$subject_id %in% subject_ids_in)
  
  mfpca_in <- pca.fd(fdobj = mfd_obj_train[inds_in],
                     nharm = k_retain,
                     harmfdPar = fdPar(fdobj = mfd_obj_train$basis,
                                       Lfdobj = 2, lambda = 10^-12))
  predictions_out <- predict_new_observations_pca.fd(
    new_fd_obj = mfd_obj_train[inds_out],
    pca_fd_obj = mfpca_in)
  
  reconstructions_out <- predictions_out$reconstructed_fd
  
  ve_loso_cv[k] <- variance_explained_reconstruction(
    true_fd_obj = mfd_obj_train[inds_out],
    reconstructed_fd_obj = reconstructions_out)
}


hist(ve_loso_cv) # skewed. look at mean as well as median
median(ve_loso_cv)
mean(ve_loso_cv)

loso_cv_dt <- data.table(ind = 1:N_sub,
                         ve_loso_cv = ve_loso_cv)

cv_dt <- merge.data.table(x = subject_id_dt_train,
                          y = loso_cv_dt,
                          by = "ind")

cv_dt <- merge.data.table(x = cv_dt,
                          y = ten_fold_cv_dt,
                          by = "fold", all.x = TRUE)

fwrite(cv_dt, file = file.path(outputs_path, "tables", "cv-table.csv"))
fwrite(mfpc_dt, file = file.path(outputs_path, "tables", "mfpc-table.csv"))
# -------------------------------------------------------------------------


save_list <- list( 
  mfpca = mfpca,
  k_retain = k_retain,
  covariates_dt_test = covariates_dt_test,
  covariates_dt_train = covariates_dt_train,
  mfd_obj_test = mfd_obj_test,
  mfd_obj_train = mfd_obj_train,
  mfd_obj_reconstructions_test = mfd_obj_reconstructions_test,
  cv_dt = cv_dt,
  test_sample = test_sample,
  N_test = N_test,
  N_train = N_train,
  mfpc_dt = mfpc_dt
  )

saveRDS(save_list, file = "outputs/basis-transformation-results.rds")




