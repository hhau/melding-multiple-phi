library(dplyr)

source("scripts/common/plot-settings.R")

submodel_one_plot_tbl <- readRDS("rds/surv-example/plot-submodel-1-data.rds")
sim_settings <- readRDS("rds/surv-example/simulation-settings-and-joint-data.rds")

p1 <- ggplot(submodel_one_plot_tbl) +
  geom_point(aes(x = time, y = y_point)) +
  geom_line(
    data = submodel_one_plot_tbl %>% filter(!is.na(y_mean)),
    mapping = aes(x = time, y = y_mean, linetype = event_indicator)
  ) +
  geom_ribbon(
    data = submodel_one_plot_tbl %>% filter(!is.na(y_lower)),
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
      y = sim_settings$y_threshold,
      yend = sim_settings$y_threshold,
      x = interval_lower,
      xend = interval_upper
    ),
    col = blues[2],
    size = rel(2),
    alpha = 0.5
  ) +
  facet_wrap(vars(patient_id))

ggsave_fullpage(
  filename = "plots/surv-example/submodel-one-posterior.pdf",
  plot = p1
)
