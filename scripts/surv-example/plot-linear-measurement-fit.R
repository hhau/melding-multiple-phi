library(ggplot2)
library(tidybayes)
library(dplyr)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

# read in data
simulated_data <- readRDS(
  file = "rds/surv-example/submodel-one-simulated-data.rds"
)
submodel_one_settings <- readRDS(
  file = "rds/surv-example/submodel-one-simulation-settings.rds"
)

# read in model output
model_output <- readRDS(
  file = "rds/surv-example/submodel-one-output.rds"
)

# build the first data plot
base_plot <- ggplot(simulated_data, aes(x = time, y = measurement)) +
  geom_point() +
  facet_wrap(~ patient_id) +
  bayesplot:::force_x_axis_in_facets() +
  theme(panel.spacing.x = unit(0.8, "lines"))

# add the posterior mean + 95% quantile interval
# this is a little slow, could consider caching?
# really should store all intermediary data
res <- model_output$samples %>%
  array_to_mcmc_list() %>%  
  spread_draws(plot_mu[patient_id, t_plot])

res$t_plot <- submodel_one_settings$t_plot[res$t_plot]

posterior_plot_data <- res %>% 
  point_interval(
    .width = 0.8,
    .point = mean,
    .interval = qi
  )

with_post_mean <- base_plot + 
  geom_line(
    data = posterior_plot_data, 
    aes(x = t_plot, y = plot_mu),
    inherit.aes = FALSE,
    col = highlight_col
  ) +
  geom_ribbon(
    data = posterior_plot_data, 
    aes(x = t_plot, ymin = .lower, ymax = .upper),
    inherit.aes = FALSE,
    alpha = 0.2
  )

# add the threshold
with_threshold <- with_post_mean +
  geom_abline(
    intercept = submodel_one_settings$y_threshold,
    slope = 0,
    alpha = 0.5,
    lty = "dashed"
  )

# add the event times density (use ggdist?)
event_res <- model_output$samples %>%
  array_to_mcmc_list() %>%  
  spread_draws(event_time[patient_id])
  
event_data <- event_res %>% 
  point_interval(.width = 0.8)

trimmed_event_data <- event_data %>% 
  rowwise() %>% 
  mutate(
    across(.lower, function(x) ifelse(between(x, 0, 1), x, 0)),
    across(.upper, function(y) ifelse(between(y, 0, 1), y, 1))
  ) %>% 
  filter(
    !(.lower == 0 & .upper == 1)
  )

final_version <- with_threshold + 
  geom_segment(
    data = trimmed_event_data,
    aes(
      y = submodel_one_settings$y_threshold,
      yend = submodel_one_settings$y_threshold,
      x = .lower,
      xend = .upper
    ),
    col = blues[2],
    size = rel(2),
    alpha = 0.6
  ) + 
  scale_x_continuous(
    limits = c(
      min(submodel_one_settings$t_plot),
      max(submodel_one_settings$t_plot)
    )
  )

ggsave_fullpage(
  filename = "plots/surv-example/submodel-one-posterior.pdf",
  plot = final_version
)
