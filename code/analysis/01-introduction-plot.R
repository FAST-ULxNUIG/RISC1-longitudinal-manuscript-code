library(data.table) # Extension of `data.frame`, CRAN v1.14.2
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics, CRAN v3.4.0
library(tikzDevice) # R Graphics Output in LaTeX Format, CRAN v0.12.3.1

# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
options(scipen = 999)

# -------------------------------------------------------------------------
subset_dt_lng <- readRDS(file = here::here("data/data-for-intro-plot.rds"))
head(subset_dt_lng)
sample_hz <- 200
second_per_frame <- 1/sample_hz
subset_dt_lng[, seconds := second_per_frame * frame_new]
segment_times <- c(unique(subset_dt_lng[frame == 0, seconds]), max(subset_dt_lng$seconds))

theme_set(theme_bw())
theme_update(strip.text = element_text(size = 9),
             text = element_text(size = 9),
             panel.grid.minor = element_blank(),
             axis.title = element_text(size = 9),
             legend.text = element_text(size = 9, hjust = 0.5),
             axis.text = element_text(size = 9),
             legend.position = "bottom",
             legend.title = element_text(size = 9, hjust = 0.5),
             plot.title = element_text(hjust = 0.5, size = 11))
p <- ggplot(data = subset_dt_lng) +
  aes(x = seconds, y = angle, colour = factor(stride_num)) +
  facet_wrap(~ location, nrow = 3, ncol = 1, "free_y") +
  geom_line() +
  labs(colour = "Stride Number:") +
  geom_vline(xintercept = segment_times, col = "darkgrey", lty = 2) +
  labs(y = "Angle ($^{\\circ}$)", x = "Time (seconds)") +
  scale_x_continuous(expand = c(0.01, 0.01))
p

tikz(here::here("outputs", "figures", "subject-01-adjacent-strides.tex"),
     width  = 1 * doc_width_inches, 
     height = 0.75 * doc_width_inches)
p
dev.off()
