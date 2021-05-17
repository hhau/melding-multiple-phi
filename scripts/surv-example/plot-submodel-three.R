library(dplyr)

source("scripts/common/plot-settings.R")

submodel_three_plot_tbl <- readRDS("rds/surv-example/plot-submodel-3-data.rds")

p1 <- ggplot(submodel_three_plot_tbl) +
  geom_point(aes(x = time, y = y_point)) +
  geom_line(
    data = submodel_three_plot_tbl %>% filter(!is.na(y_mean)),
    mapping = aes(x = time, y = y_mean, linetype = event_indicator)
  ) +
  geom_ribbon(
    data = submodel_three_plot_tbl %>% filter(!is.na(y_lower)),
    aes(x = time, ymin = y_lower, ymax = y_upper), 
    alpha = 0.2
  ) +
  facet_wrap(vars(patient_id))

ggsave_fullpage(
  filename = "plots/surv-example/submodel-three-posterior.pdf",
  plot = p1
)
