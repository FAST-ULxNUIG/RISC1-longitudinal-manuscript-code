# Load Packages Used ------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice) # R Graphics Output in LaTeX Format, CRAN v0.12.3.1


# Some plot Settings: -----------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
options(scipen = 999)
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

# Read in data: -----------------------------------------------------------
subset_dt <- readRDS(file = "data/data-for-strides-plot.rds")

# Basic settings from time variable: --------------------------------------
sample_hz <- 200
second_per_frame <- 1/sample_hz

# Subset the data: --------------------------------------------------------
# (demo will be on the knee angle)
# and re-shape from wide to long for plotting.
subset_knee_dt <- subset_dt[location=="Knee"]
subset_knee_dt_lng <- melt.data.table(data = subset_knee_dt,
                                      id.vars = c("location",
                                                  "stride_num",
                                                  "side",
                                                  "frame_num"),
                                      measure.vars = paste0("data_", 0:197),
                                      variable.name = "frame_new",
                                      value.name = "angle",
                                      variable.factor = FALSE,
                                      value.factor = FALSE,
                                      na.rm = TRUE,
                                      verbose = TRUE)
subset_knee_dt_lng[, frame_new := as.numeric(stringr::str_remove(frame_new, "data_"))]
subset_knee_dt_lng[, frame := frame_num + frame_new]
subset_knee_dt_lng[, frame := frame - min(frame)]
subset_knee_dt_lng[side == "right", stride_num := stride_num - 1]
subset_knee_dt_lng[, seconds := second_per_frame * frame]
subset_knee_dt_lng[,side := factor(side, 
                                   levels = c("left", "right"),
                                   labels = c("Left", "Right"))]


# Get stride start and end times: -----------------------------------------
segment_times <- unique(subset_knee_dt_lng[, .(side, frame_num)])
segment_times[, frame_num := frame_num - min(frame_num)]
segment_times[, seconds_time := second_per_frame * frame_num]



setorderv(subset_knee_dt_lng, c("location", "stride_num", "frame"))
head(subset_knee_dt_lng)


# Calculate midpoints of each stride: -------------------------------------
# (x-axis lables for each stride will be centred here) #
stride_midpoints <- subset_knee_dt_lng[, .(midpoint = mean(seconds)),
                                       by = .(side, stride_num)]

stride_midpoints[, label := paste("Stride", stride_num, side)]

# Create plot: ------------------------------------------------------------
# Create left (p1) and right-side (p2) plots separately and then combine


## Left: ------------------------------------------------------------------
p1 <- ggplot(data = subset_knee_dt_lng[side == "Left" &
                                         seconds <= 3.58 &
                                         seconds >= 0.365]) +
  geom_line(aes(x = seconds, y = angle), colour = "red3") +
  labs(y = "Angle ($^{\\circ}$)", x = "Time (seconds)") +
  scale_x_continuous(expand = c(0.01, 0.01),
                     breaks = stride_midpoints[side=="Left" & stride_num < 5]$midpoint,
                     labels = stride_midpoints[side=="Left" & stride_num < 5]$label) +
  geom_vline(data = segment_times[seconds_time <= 3.58 & seconds_time >= 0.365 & side == "Left"],
             mapping = aes(xintercept = seconds_time),
             colour = "red3",
             lty = 2, linewidth = 0.5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = "red3"),
        legend.position = "none",
        axis.text.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin = margin(0.5,1,0.5,1.2, "cm"),
        axis.title = element_blank()) +
  annotate(geom = "label", label = "$T=0$", x = 0.8, y = 25, size = 3) +
  annotate("segment", x = 0.675, y = 20, xend = 0.4, yend = 0,
           arrow = arrow(type = "closed", length = unit(0.05, "npc")))

p1


## Right ------------------------------------------------------------------
p2 <- ggplot(data = subset_knee_dt_lng[side == "Right" &
                                         seconds <= 3.58 & 
                                         seconds >= 0.365]) +
  geom_line(aes(x = seconds, y = angle), colour = "#619CFF") +
  labs(y = "Angle ($^{\\circ}$)", x = "Time (seconds)") +
  scale_x_continuous(expand = c(0.01, 0.01),
                     position = "top",
                     breaks = stride_midpoints[side=="Right" & stride_num > 0]$midpoint,
                     labels = stride_midpoints[side=="Right" & stride_num > 0]$label) +
  geom_vline(data = segment_times[seconds_time <= 3.58 &
                                    seconds_time >= 0.365 &
                                    side == "Right"],
             mapping = aes(xintercept = seconds_time),
             colour = "#619CFF",
             lty = 2, linewidth = 0.5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = "#619CFF"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        plot.margin = margin(0.5,1,0.5,1.2, "cm"),
        axis.title = element_blank()) +
  annotate(geom = "label", label = "$T_{i, \\mathrm{right}, 1}$", x = 1.11, y = 30, size = 3) +
  annotate("segment", x = 0.96, y = 25, xend = 0.72, yend = 5,
           arrow = arrow(type = "closed", length = unit(0.05, "npc")))

p <- ggarrange(p1, p2, ncol = 1, nrow = 2)


# Write as TikZ file: -----------------------------------------------------
tikz(here::here("outputs", 
                "plots",
                "subject-01-adjacent-strides-timings.tex"),
     width  = 1 * doc_width_inches, 
     height = 0.5 * doc_width_inches)
p
dev.off()

