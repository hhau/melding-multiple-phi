library(ggplot2)
library(tidybayes)
library(dplyr)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

# read in data
simulated_data <- readRDS(
  file = "rds/surv-example/submodel-one-simulated-data.rds"
)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

# read in model output
model_output <- readRDS(
  file = "rds/surv-example/submodel-one-output.rds"
)

# build the first data plot
base_plot <- ggplot(simulated_data, aes(x = time, y = measurement)) +
  geom_point() +
  facet_wrap(~ patient_id) +
  ylab(expression(italic(z))) +
  xlab(expression(italic(t))) +
  scale_x_continuous(
    breaks = seq(from = 0, to = 1, length.out = 3),
    labels = function(x) round(x, digits = 2),
    limits = c(0, 1)
  ) +
  # bayesplot:::force_x_axis_in_facets() +
  theme(
    panel.spacing.x = unit(0.3, "lines"),
    axis.text = element_text(size = rel(0.8)),
    strip.text = element_text(size = rel(0.8))
  )

# add the posterior mean + 95% quantile interval
# this is a little slow, could consider caching?
# really should store all intermediary data
res <- model_output$samples %>%
  array_to_mcmc_list() %>%  
  spread_draws(plot_mu[patient_id, plot_x])

res$plot_x <- sim_settings$x_plot[res$plot_x]

posterior_plot_data <- res %>% 
  point_interval(
    .width = 0.8,
    .point = mean,
    .interval = qi
  )

event_df <- tibble(
  patient_id = 1 : sim_settings$n_patients,
  event_indicator = sim_settings$event_indicator
)

posterior_plot_data <- posterior_plot_data %>%
  left_join(event_df, by = c("patient_id" = "patient_id")) %>% 
  mutate(event_indicator = as.factor(event_indicator))

with_post_mean <- base_plot + 
  geom_line(
    data = posterior_plot_data, 
    aes(x = plot_x, y = plot_mu, lty = event_indicator),
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = posterior_plot_data, 
    aes(x = plot_x, ymin = .lower, ymax = .upper),
    inherit.aes = FALSE,
    alpha = 0.2
  ) + 
  theme(
    legend.position = "none"
  ) +
  scale_linetype_manual(
    values = c(
      "1" = "solid",
      "0" = "dotted"
    )
  )

# add the threshold
with_threshold <- with_post_mean +
  geom_abline(
    intercept = sim_settings$y_threshold,
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
      y = sim_settings$y_threshold,
      yend = sim_settings$y_threshold,
      x = .lower,
      xend = .upper
    ),
    col = blues[2],
    size = rel(2),
    alpha = 0.6
  )

ggsave_fullpage(
  filename = "plots/surv-example/submodel-one-posterior.pdf",
  plot = final_version
)
