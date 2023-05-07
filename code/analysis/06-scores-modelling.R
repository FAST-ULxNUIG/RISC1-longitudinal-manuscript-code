# Load packages: ----------------------------------------------------------
library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155

# Read in results and unpack: ---------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
mfpca <- basis_transformation_results$mfpca
k_retain <- basis_transformation_results$k_retain
covariates_dt_train <- basis_transformation_results$covariates_dt_train
N_train <- basis_transformation_results$N_train

                    # Create scores: ----------------------------------------------------------
# scores for mv-fpca are in a 
# N_{Total} \times K \times P matrix
# Sum over the hip knee and ankle to get overall scores:
scores_train <- apply(mfpca$scores, c(1, 2), sum)
colnames(scores_train) <- paste0("score_", seq_len(k_retain))

# Some basic checks:
stopifnot(dim(scores_train) == c(N_train, k_retain))
stopifnot(nrow(scores_train) == nrow(covariates_dt_train))

# And join back up into data.table for modelling:
covariates_and_scores_dt_train <- cbind(covariates_dt_train, scores_train)



# Modelling ---------------------------------------------------------------





