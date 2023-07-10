library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice)
outputs_path <- here::here("outputs")
functions_path <- here::here("code", "functions")
source(file.path(functions_path, "theme_gunning.R"))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
options(scipen = 999)
theme_gunning()

unique_strides_dt <- fread(file = file.path(outputs_path, "unique_strides_dt.csv"))
run_duration_dt <- fread(file = file.path(outputs_path, "run_duration_dt.csv"))

run_duration_dt[subject_id %in% unique_strides_dt$subject_id][, mean(run_duration>50 & run_duration<70)]

p1 <- ggplot(data = run_duration_dt[subject_id %in% unique_strides_dt$subject_id]) +
  aes(x = run_duration) +
  geom_histogram(colour = "darkgreen", fill = "grey", binwidth = 2) +
  scale_x_continuous(breaks = c(20, 40, 60, 80)) +
  geom_vline(xintercept = 60, linetype = 2) +
  labs(y = "Count",
       x = "Duration of Treadmill Run (seconds)", 
       title = "(a) Treadmill Run Duration") +
  geom_vline(xintercept = c(50, 70), linetype = 3)

p2 <- ggplot(unique_strides_dt) +
  aes(x = unique_strides) +
  geom_histogram(colour = "darkgreen", fill = "grey") +
  geom_vline(xintercept = 80, linetype = 2) +
  labs(y = "Count", 
       title = "(b) Number of Strides in the Dataset",
       x = "Number of Strides per Subject and Side")

combined_plot <- ggarrange(p1, p2, ncol = 2, nrow = 1)

tikz(here::here("outputs", "figures", "data-prep-plot.tex"),
     width  = 1 * doc_width_inches, 
     height = 0.5 * doc_width_inches,
     standAlone = TRUE)
combined_plot
dev.off()

tinytex::lualatex(here::here("outputs", "figures", "data-prep-plot.tex"))




