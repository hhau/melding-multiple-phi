library(dplyr)

source("scripts/common/plot-settings.R")

submodel_one_plot_data <- readRDS("rds/surv-example/plot-submodel-1-data.rds")
submodel_three_plot_data <- readRDS("rds/surv-example/plot-submodel-3-data.rds")

interesting_patients <- c(17, 22, 32)

plot_tbl <- bind_rows(submodel_one_plot_data, submodel_three_plot_data) %>%
  filter(patient_id %in% interesting_patients) %>%
  mutate(
    model_name = factor(
      x = model_name,
      levels = c(
        "event_time",
        "longitudinal"
      ),
      labels = c(
        "italic(z) ~ '--' ~ 'p'[1]",
        "italic(x) ~ '--' ~ 'p'[3]"
      )
    )
  )

outer_threshold_value <- plot_tbl %>%
  select(threshold_value) %>%
  filter(!is.na(threshold_value)) %>%
  distinct() %>%
  pull()

p1 <- ggplot(plot_tbl) +
  geom_point(aes(x = time, y = y_point)) +
  geom_line(
    data = plot_tbl %>% filter(!is.na(y_mean)),
    mapping = aes(x = time, y = y_mean, linetype = event_indicator, colour = model_name)
  ) +
  geom_ribbon(
    data = plot_tbl %>% filter(!is.na(y_lower)),
    aes(x = time, ymin = y_lower, ymax = y_upper), 
    alpha = 0.2
  ) +
  geom_hline(
    mapping = aes(yintercept = threshold_value),
    alpha = 0.5,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      y = outer_threshold_value,
      yend = outer_threshold_value,
      x = interval_lower,
      xend = interval_upper
    ),
    col = greens[3],
    size = rel(4),
    alpha = 0.5
  ) +
  facet_grid(
    cols = vars(patient_id),
    rows = vars(model_name),
    scales = "free",
    labeller = label_parsed
  ) +
  theme(
    strip.text.x = element_text(margin = margin(t = 2, b = 2)),
    strip.text.y = element_text(
      angle = 0,
      margin = margin(l = 2, r = 2)
    ),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  ) +
  scale_linetype_manual(
    values = c(
      "1" = "solid",
      "0" = "4212"
    )
  ) +
  scale_colour_manual(
    values = c(
      "italic(z) ~ '--' ~ 'p'[1]" = highlight_col,
      "italic(x) ~ '--' ~ 'p'[3]" = blues[2]
    )
  ) +
  xlab(expression(italic(t))) 

p1

ggsave_halfheight(
  filename = "plots/surv-example/both-longitudinal-submodels.pdf",
  plot = p1
)
