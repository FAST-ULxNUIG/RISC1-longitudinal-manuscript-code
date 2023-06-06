plot_list <- readRDS(here::here("outputs", "simulation", "simulated-data-plot.rds"))

library(fda)        # CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1

source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning()
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
plots_path <- here::here("outputs", "figures")


true_fd <- plot_list$mfd_obj_test
true_array <- eval.fd(evalarg = 0:100, fdobj = true_fd)
true_mat <- rbind(true_array[,,1],
                  true_array[,,2],
                  true_array[,,3])
colnames(true_mat) <- paste0("obs_", seq_len(200))
true_dt <- data.table(t_grid = rep(0:100, times = 3), dim = rep(c("Hip", "Knee", "Ankle"), each = 101), true_mat)
true_dt_lng <- melt.data.table(data = true_dt,
                               id.vars = c("t_grid", "dim"),
                               measure.vars =  paste0("obs_", seq_len(200)),
                               variable.name = "obs", 
                               value.name = "angle", 
                               variable.factor = FALSE, 
                               value.factor = FALSE)

simulated_fd <- plot_list$demo_fd_obj_test
simulated_array <- eval.fd(evalarg = 0:100, fdobj = simulated_fd)
simulated_mat <- rbind(simulated_array[,,1],
                  simulated_array[,,2],
                  simulated_array[,,3])
colnames(simulated_mat) <- paste0("obs_", seq_len(200))
simulated_dt <- data.table(t_grid = rep(0:100, times = 3), dim = rep(c("Hip", "Knee", "Ankle"), each = 101), simulated_mat)
simulated_dt_lng <- melt.data.table(data = simulated_dt,
                               id.vars = c("t_grid", "dim"),
                               measure.vars =  paste0("obs_", seq_len(200)),
                               variable.name = "obs", 
                               value.name = "angle", 
                               variable.factor = FALSE, 
                               value.factor = FALSE)


simulated_dt_lng$type <- "Simulated"
true_dt_lng$type <- "True"

plot_dt <- rbind(simulated_dt_lng, true_dt_lng)
plot_dt[, dim := factor(dim, levels = c("Hip", "Knee", "Ankle"))]

# -------------------------------------------------------------------------


p <- ggplot(data = plot_dt) +
  facet_grid(dim ~ type, scales = "free_y") +
  aes(x = t_grid, y = angle, group = obs) +
  geom_line(alpha = 0.25) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Angle ($^{\\circ}$)", x = "Normalised Time ($\\%$ of Stride)")

tikz(file.path(plots_path, "simulated-vs-true.tex"),
     width = 0.85 * doc_width_inches, 
     height = 1 * doc_width_inches, 
     standAlone = TRUE)
print(p)
dev.off()

tinytex::lualatex(file = file.path(plots_path, "simulated-vs-true.tex"))



