# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice) # CRAN v0.12.3.1
library(xtable)
source(here::here("code","functions", "functions-helper-smoothing.R"))

# Path: -------------------------------------------------------------------
simulation_results_path <- "/Users/edwardgunning/Dropbox/simulation"
plots_path <- here::here("outputs", "figures")
outputs_path <- here::here("outputs")


# Read in: ----------------------------------------------------------------
## Baseline ---------------------------------------------------------------
N_sub_280_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_1.rds"
  )
)
## Sample Size: -----------------------------------------------------------
N_sub_500_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_500_prop_missing_0.1_long_strength_1.rds"
  )
)

N_sub_1000_prop_missing_0.1_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_1000_prop_missing_0.1_long_strength_1.rds"
  )
)

N_sub_1000_prop_missing_0.1_long_strength_1_part_2 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_1000_prop_missing_0.1_long_strength_part_21.rds"
  )
)

## Long Strength: --------------------------------------------------------
N_sub_280_prop_missing_0.1_long_strength_2 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_2.rds"
  )
)

N_sub_280_prop_missing_0.1_long_strength_3 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.1_long_strength_3.rds"
  )
)

## Prop Missing ------------------------------------------------------------
N_sub_280_prop_missing_0.2_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.2_long_strength_1.rds"
  )
)

N_sub_280_prop_missing_0.5_long_strength_1 <- readRDS(
  file.path(
    simulation_results_path,
    "N_sub_280_prop_missing_0.5_long_strength_1.rds"
  )
)






# Extract and Name: -------------------------------------------------------
results_280_0.1_1 <- N_sub_280_prop_missing_0.1_long_strength_1[1:500]
results_500_0.1_1 <- N_sub_500_prop_missing_0.1_long_strength_1[1:500]
results_1000_0.1_1 <- N_sub_1000_prop_missing_0.1_long_strength_1[1:250]
results_1000_0.1_1[251:500] <- N_sub_1000_prop_missing_0.1_long_strength_1[1:250]
results_280_0.1_2 <- N_sub_280_prop_missing_0.1_long_strength_2[1:500]
results_280_0.1_3 <- N_sub_280_prop_missing_0.1_long_strength_3[1:500]
results_280_0.2_1 <- N_sub_280_prop_missing_0.2_long_strength_1[1:500]
results_280_0.5_1 <- N_sub_280_prop_missing_0.5_long_strength_1[1:500]





format_singularity <- function(x) {
  stopifnot(ncol(x) == 4)
  npc <- length(x[["naive"]])
  fpc_num <- seq_len(12)
  if(is.null(x[["fpca"]])) x[["fpca"]] <- rep(NA, npc)
  
  
  if(npc < 12) {
    df <- data.frame(fpc_num = fpc_num,
                     poly = c(x[["poly"]], rep(NA, 12 - npc)),
                     naive = c(x[["naive"]], rep(NA, 12 - npc)),
                     spline_subject_ri_side = c(x[["spline_subject_ri_side"]], rep(NA, 12 - npc)),
                     fpca = c(x[["fpca"]], rep(NA, 12 - npc)))
                     
  } else if(npc == 12) {
    df <- data.frame(fpc_num = fpc_num,
                     poly = x[["poly"]],
                     naive = x[["naive"]],
                     spline_subject_ri_side = x[["spline_subject_ri_side"]],
                     fpca = x[["fpca"]])
  }
  
  df
  
  }



singularity280_0.1_1 <- lapply(results_280_0.1_1, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity280_0.1_1_df <- purrr::map_dfr(.x = singularity280_0.1_1, format_singularity, .id = "sim_rep")
singularity280_0.1_1_df$N <- 280
singularity280_0.1_1_df$prop_missing <- 0.1
singularity280_0.1_1_df$long_strength <- 1


singularity500_0.1_1 <- lapply(results_500_0.1_1, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity500_0.1_1_df <- purrr::map_dfr(.x = singularity500_0.1_1, format_singularity, .id = "sim_rep")
singularity500_0.1_1_df$N <- 500
singularity500_0.1_1_df$prop_missing <- 0.1
singularity500_0.1_1_df$long_strength <- 1

singularity1000_0.1_1 <- lapply(results_1000_0.1_1, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity1000_0.1_1_df <- purrr::map_dfr(.x = singularity1000_0.1_1, format_singularity, .id = "sim_rep")
singularity1000_0.1_1_df$N <- 1000
singularity1000_0.1_1_df$prop_missing <- 0.1
singularity1000_0.1_1_df$long_strength <- 1


singularity280_0.2_1 <- lapply(results_280_0.2_1, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity280_0.2_1_df <- purrr::map_dfr(.x = singularity280_0.2_1, format_singularity, .id = "sim_rep")
singularity280_0.2_1_df$N <- 280
singularity280_0.2_1_df$prop_missing <- 0.2
singularity280_0.2_1_df$long_strength <- 1

singularity280_0.5_1 <- lapply(results_280_0.5_1, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity280_0.5_1_df <- purrr::map_dfr(.x = singularity280_0.5_1, format_singularity, .id = "sim_rep")
singularity280_0.5_1_df$N <- 280
singularity280_0.5_1_df$prop_missing <- 0.5
singularity280_0.5_1_df$long_strength <- 1

singularity280_0.1_2 <- lapply(results_280_0.1_2, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity280_0.1_2_df <- purrr::map_dfr(.x = singularity280_0.1_2, format_singularity, .id = "sim_rep")
singularity280_0.1_2_df$N <- 280
singularity280_0.1_2_df$prop_missing <- 0.1
singularity280_0.1_2_df$long_strength <- 2

singularity280_0.1_3 <- lapply(results_280_0.1_3, function(y) {
  lapply(y$models, function(x) x$model_checks$singular_fit)
})
singularity280_0.1_3_df <- purrr::map_dfr(.x = singularity280_0.1_3, format_singularity, .id = "sim_rep")
singularity280_0.1_3_df$N <- 280
singularity280_0.1_3_df$prop_missing <- 0.1
singularity280_0.1_3_df$long_strength <- 3


singularity_df <- rbind(
  singularity280_0.1_1_df,
  singularity500_0.1_1_df,
  singularity1000_0.1_1_df,
  singularity280_0.2_1_df,
  singularity280_0.5_1_df,
  singularity280_0.1_2_df,
  singularity280_0.1_3_df
)



singularity_dt <- as.data.table(singularity_df)
singularity_dt_summary <- singularity_dt[fpc_num %in% 1:10,
                                         as.list(round(apply(.SD, 2, mean, na.rm = TRUE), 3)),
                                         by = .(N, prop_missing, long_strength),
                                         .SDcols = c("poly", "naive", "spline_subject_ri_side", "fpca")]

singularity_dt_summary <- singularity_dt_summary[, -c("naive")]



binomial_se <- function(phat, n) {
  sqrt((phat * (1-phat)) / n)
}
      
singularity_dt_summary[, c("poly", "spline_subject_ri_side", "fpca") := {
  my_list <- get_list_of_columns(apply(.SD,
                                       c(1,2),
                                       function(x) {
                                         if(x != 0) {
                                           
                                           paste0(x, " (", ifelse(round(binomial_se(x, n = 500 * 10), 3) > 0, round(binomial_se(x, n = 500 * 10), 3), "$<$ 0.001"), ")")
                                         } else if(x == 0) {paste0(x)}
                                       }))
  names(my_list) <- c("poly", "spline_subject_ri_side", "fpca")
  my_list
}, .SDcols = c("poly", "spline_subject_ri_side", "fpca")]



setnames(singularity_dt_summary, 
         old = c("N", "prop_missing", "long_strength", "poly", "spline_subject_ri_side", 
                 "fpca"),
         new = c("$N$", "Pr. Missing", "Lon. Strength", "Polynomial", "Spline", 
                 "ml-FPCA"))

bold <- function(x) {
  paste0("{\\bfseries \\small ", x, "}") 
}
Singularity_table <- xtable(singularity_dt_summary, 
                     label = "tab:singularity-table",
                     caption = "Proportion of singular fit warnings from the model fits. In cases where the proportion is non-zero, a Monte Carlo standard error estimate for the true proportion is reported in brackets to convey uncertainty due to the finite number of simulations.")
align(Singularity_table)[1] <- "l"
print(Singularity_table, 
      include.rownames = FALSE,
      file = file.path(outputs_path, "tables", "Singularity-simulation.tex"),
      sanitize.text.function = function(x){x},
      sanitize.colnames.function = bold,
      booktabs = TRUE)
fwrite(singularity_dt_summary, file = file.path(outputs_path, "tables", "Singularity-simulation.csv"))

