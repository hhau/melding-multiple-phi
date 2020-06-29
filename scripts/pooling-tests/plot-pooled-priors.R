library(ggplot2)

source("scripts/common/plot-settings.R")

plot_data <- readRDS(
  file = "rds/pooling-tests/pooling-types-weights-results.rds"
)

p_1 <- ggplot(plot_data, aes(x = x, y = y, z = f_val)) +
  geom_contour(aes(colour = base_dist), bins = 7, alpha = 0.8) +
  facet_grid(
    rows = vars(weight_case), cols = vars(pooling_type),
    labeller = label_parsed
  ) +
  xlab(expression(phi["1:2"])) + 
  ylab(expression(phi["2:3"])) +
  scale_x_continuous(expand = c(0, 0), limits = c(-5, 5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 5)) +
  theme(
    legend.text.align = 0,
    panel.spacing.x = unit(1.5, "lines"),
    panel.spacing.y = unit(1, "lines"),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  ) +
  scale_colour_manual(
    labels = c(
      "normal" = "Gaussian",
      "student_t" = expression("Student" ~ "-" ~ italic("t")[4])
    ),
    values = c(
      "normal" = blues[2],
      "student_t" = highlight_col
    )
  ) +
  labs(
   colour = "Density" 
  ) +
  NULL

ggsave_fullpage(
  filename = "plots/pooling-tests/pooled-densities-2d.pdf",
  plot = p_1
)
