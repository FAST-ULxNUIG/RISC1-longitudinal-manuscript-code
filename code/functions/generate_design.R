# Function to generate mutliple-subject design in which an equal number of left and 
# right side observations are observed for each subject at alternating points on
# an equally-spaced grid of points.
generate_design_single_subject <- function(n_i, speed_sd) {
  stopifnot(is.numeric(speed_sd) & speed_sd > 0)
  stopifnot(n_i >= 2)
  t_grid <- seq(0, 1, length.out = 2 * n_i)
  df_subject <- data.frame(side = rep(c("left", "right"), each =  2 * n_i),
                           time = rep(t_grid, times = 2),
                           stride_ind = rep(seq_len(2 * n_i), times = 2))
  df_subject <- df_subject[(df_subject$side == "left" & (df_subject$stride_ind %% 2 == 0)) |
                             (df_subject$side == "right" & (df_subject$stride_ind %% 2 != 0)), ]
  df_subject <- df_subject[order(df_subject$stride_ind), ]
  df_subject$side <- factor(df_subject$side, levels = c("left", "right"))
  df_subject$speed_cent <- rnorm(n = 1, mean = 0, sd = speed_sd)
  df_subject$sex <- runif(n = 1, min = 0, max = 1)
  df_subject$sex <- ifelse(df_subject$sex <= 0.5, "male", "female")
  df_subject
}

generate_design_multiple_subjects <- function(N, n_i, speed_sd) {
  subject_inds <- seq_len(N)
  names(subject_inds) <- subject_inds
  df <- purrr::map_dfr(.x = subject_inds, 
                 .f = function(x) {generate_design_single_subject(n_i = n_i, speed_sd = speed_sd)},
                 .id = "subject_id")
  df$subject_id <- factor(df$subject_id, levels = paste(seq_len(N)))
  df$sex <- factor(df$sex, levels = c("male", "female"))
  df
}




