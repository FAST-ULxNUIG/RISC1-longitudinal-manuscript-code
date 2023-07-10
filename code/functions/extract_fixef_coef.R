require(lme4)

# Functions to extract fixed-effects coefficient functions from a  --------
extract_fixef_coef <- function(lme_list_object, K_retain) {
  # routine input checks:
  if(!is.list(lme_list_object)) stop("lme_list_object must be a list!")
  if(!(length(lme_list_object) == K_retain)) stop("lme_list_object must be of length K_retain")
  if(!all(sapply(lme_list_object, function(x){class(x) == "lmerMod"}))) stop("Each element of lme_list_object must be a lmerMod object")
  # extraction:
  fixef_mat <- sapply(lme_list_object, fixef)
  stopifnot(ncol(fixef_mat) == K_retain)
  colnames(fixef_mat) == paste0("score_", K_retain)
  fixef_mat
}

compute_fixef_functions <- function(fit_object, pca_fd_obj, vars = NULL) {
  
  # vars = character vector containing the names of the fixed effects to extract.
  
  # Extract K from naive fit object:
  K_retain <- fit_object$K_retain
  fixef_coef_matrix <- extract_fixef_coef(lme_list_object = fit_object$lme_fit_list, K_retain = K_retain)
  
  if(!is.null(vars)) {
    if(!all(vars %in% rownames(fixef_coef_matrix))) {
      stop("vars must be a list of strings that correspond to names of fixed effects")
    }
    fixef_coef_matrix <- fixef_coef_matrix[vars, ]
  }
  
  # Use to predict functional data object:
  fd_return <- construct_fd_from_scores(
    pca_fd_obj = pca_fd_obj,
    scores_matrix = fixef_coef_matrix,
    K = K_retain, 
    add_back_mean = FALSE)
  
  # add names of variables:
  fd_return$fdnames$reps <- rownames(fixef_coef_matrix)
  fd_return
}

