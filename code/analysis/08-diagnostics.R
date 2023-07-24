# ------------------------------------------------------------------------#
# Diagnostics from Mixed Model Fit:
# ------------------------------------------------------------------------#

# Load packages: ----------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(lme4)       # CRAN v1.1-30
library(nlme)       # CRAN v3.1-155
library(splines)
library(ggplot2)    # CRAN v3.4.2
library(tikzDevice)
library(ggpubr)
# source(here::here("code", "functions",
#                   "source_all_analysis_functions.R"))

# Graphics Settings: ------------------------------------------------------
# rough guide for sizing of plot outputs:
diagnostics_path <- here::here("outputs", "figures", "diagnostics")
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning()
theme_update(
  axis.text = element_text(size = 10.5),
  axis.title = element_text(size = 10.5),
  strip.text = element_text(size = 10.5),
  plot.title = element_text(size = 11.5),
  legend.text = element_text(size = 10.5),
  legend.title = element_text(size = 11)
)

# Load Results: -----------------------------------------------------------
basis_transformation_results <- readRDS(
  here::here("outputs","basis-transformation-results.rds"))
model_fit_results <- readRDS(
  here::here("outputs", "model-fit-results.rds"))
covariates_and_scores_dt_train <- model_fit_results$covariates_and_scores_dt_train

# -------------------------------------------------------------------------
spline_ri_model <- model_fit_results$spline_ri_model
naive_model <- model_fit_results$naive_model


# Loacations Scale Plot: --------------------------------------------------
scale_loc_plot <- function(m, line.col = "red", line.lty = 1,
                           line.lwd = 2) {
  # From here:
  # https://www.rdocumentation.org/packages/lme4/versions/1.1-33/topics/plot.merMod
  plot(m, sqrt(abs(resid(.))) ~ fitted(.),
       type = c("p", "smooth"),
       main = "Scale-Location Plot",
       par.settings = list(plot.line = list(alpha=1, col = line.col,
                                            lty = line.lty, lwd = line.lwd)))
}


# Now do our own with ggplot2
# extract data:
df_scale_loc_plot <- purrr::map_dfr(
  .x = spline_ri_model$lme_fit_list[1:12],
  .f = function(x) {
    data.frame(y = sqrt(abs(resid(x))), x = fitted(x))
    }, .id = "score")

(p <- ggplot(df_scale_loc_plot) +
  aes(x = x, y = y) +
  facet_wrap(~ factor(paste("mv-FPC", score), levels = paste("mv-FPC", 1:12)), scales = "free") +
  geom_point(alpha = 0.25) +
  geom_smooth(colour = "red") +
  labs(x = "Fitted", y = "$\\sqrt{|Resid|}$",
       title = "Location-Scale Plot for the First 12 mv-FPC Scores"))

tikz(file.path(diagnostics_path, "location-scale.tex"),
     width = doc_width_inches, 
     height = 0.75 * doc_width_inches, 
     standAlone = TRUE)
print(p)
dev.off()

tinytex::lualatex(file.path(diagnostics_path, "location-scale.tex"))

# ------------------------------------------------------------------------

plot_within_group_acf <- function(list.resids, list.inds) {
  long.resids <- c()
  for (i in 2:length(list.resids)) {
    
    # assumes data are ordered
    ind_vec <- list.inds[[i]]
    full_ind_vec <- seq(min(ind_vec), max(ind_vec), by = 1)
    resid_vec <- vector(mode = "numeric", length = length(full_ind_vec))
    resid_vec[full_ind_vec %in% ind_vec] <- list.resids[[i]]
    resid_vec[!(full_ind_vec %in% ind_vec)] <- NA
    
    #add 30 missing values to the vector of residuals
    if(i >= 2) {
      long.resids <- c(long.resids, rep(NA,30))
    }
    
    #append the residuals from the next core
    long.resids <- c(long.resids, resid_vec)
  }
  acf(long.resids, na.action=na.pass, lag.max=30, main=list('Within-Group Residual Autocorrelation Plot', cex=1))
}

get_within_group_acf <- function(list.resids, list.inds) {
  long.resids <- c()
  for (i in 2:length(list.resids)) {
    
    # assumes data are ordered
    ind_vec <- list.inds[[i]]
    full_ind_vec <- seq(min(ind_vec), max(ind_vec), by = 1)
    resid_vec <- vector(mode = "numeric", length = length(full_ind_vec))
    resid_vec[full_ind_vec %in% ind_vec] <- list.resids[[i]]
    resid_vec[!(full_ind_vec %in% ind_vec)] <- NA
    
    #add 30 missing values to the vector of residuals
    if(i >= 2) {
      long.resids <- c(long.resids, rep(NA,30))
    }
    
    #append the residuals from the next core
    long.resids <- c(long.resids, resid_vec)
  }
  acf(long.resids, na.action=na.pass, lag.max=30, plot = FALSE)
}




mind <- 1



for(mind in seq_len(27)) {
  print(paste0("mvfpc", mind))
  resid_df <- data.frame(resid = resid(spline_ri_model$lme_fit_list[[mind]]))
  resid_plot <- ggplot(data = resid_df) +
    aes(sample = resid) +
    stat_qq(size = 0.5) +
    stat_qq_line() +
    labs(x = "Theoretical", y = "Sample", title = "(c) Residuals") +
    scale_x_continuous(breaks = -3:3)
  
  subject_ri_plot <- ggplot(data = data.frame(ri = ranef(spline_ri_model$lme_fit_list[[mind]])[["subject_id"]][,1])) +
    aes(sample = ri) +
    stat_qq(size = 0.5) +
    stat_qq_line() +
    labs(x = "Theoretical", y = "Sample", title = "(a) Subject-Level RI") +
    scale_x_continuous(breaks = -3:3)
  
  subject_and_side_ri_plot <- ggplot(data = data.frame(ri = ranef(spline_ri_model$lme_fit_list[[mind]])[["subject_id:side" ]][,1])) +
    aes(sample = ri) +
    stat_qq(size = 0.5) +
    stat_qq_line() +
    labs(x = "Theoretical", y = "Sample", 
         title = "(b) Subject-and-Side-Level RI") +
    scale_x_continuous(breaks = -3:3)
  
  
  
  acf_naive <- get_within_group_acf(list.resids = split(residuals(naive_model$lme_fit_list[[mind]]),
                                                        covariates_and_scores_dt_train[, interaction(subject_id, side)]),
                                    list.inds = split(covariates_and_scores_dt_train$stride_num,
                                                      covariates_and_scores_dt_train[, interaction(subject_id, side)]))
  
  acf_model <- get_within_group_acf(list.resids = split(residuals(spline_ri_model$lme_fit_list[[mind]]),
                                                        covariates_and_scores_dt_train[, interaction(subject_id, side)]),
                                    list.inds = split(covariates_and_scores_dt_train$stride_num,
                                                      covariates_and_scores_dt_train[, interaction(subject_id, side)]))
  
  acf_df <- data.frame(lag = 0:30, naive = c(acf_naive$acf), spline = c(acf_model$acf))
  (acf_plot <- ggplot(data = acf_df) +
      aes(x = lag) +
      geom_point(aes(y = naive, colour = "Naive")) +
      geom_segment(aes(y = naive, colour = "Naive", xend = lag, yend = 0), alpha =0.5) +
      geom_line(aes(y = naive, colour = "Naive")) +
      geom_point(aes(y = spline, colour = "Spline")) +
      geom_segment(aes(y = spline, colour = "Spline", xend = lag, yend = 0), alpha =0.5) +
      geom_line(aes(y = spline, colour = "Spline")) +
      geom_hline(yintercept = 0) +
      labs(y = "ACF", x = "Lag", title = "(d) Residual ACF") +
      theme(legend.position = c(0.75, 0.75),
            legend.title = element_blank()))
  
  combined_plot <- ggarrange(subject_ri_plot, subject_and_side_ri_plot, resid_plot, acf_plot, 
                             ncol = 2, nrow = 2)
  
  tikz(file.path(diagnostics_path, paste0("diagnostics-mvfpc-", mind, ".tex")),
       width = doc_width_inches, 
       height = 1 * doc_width_inches, 
       standAlone = TRUE)
  print(combined_plot)
  dev.off()
  tinytex::lualatex(file.path(diagnostics_path, paste0("diagnostics-mvfpc-", mind, ".tex")))
}
