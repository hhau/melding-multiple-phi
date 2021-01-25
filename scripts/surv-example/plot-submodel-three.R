library(ggplot2)
library(tidybayes)
library(dplyr)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

# read in data
simulated_data <- readRDS(
  file = "rds/surv-example/submodel-three-simulated-data.rds"
)
submodel_three_settings <- readRDS(
  file = "rds/surv-example/submodel-three-simulation-settings.rds"
)
submodel_one_settings <- readRDS(
  file = "rds/surv-example/submodel-one-simulation-settings.rds"
)

# read in model output
model_output <- readRDS(
  file = "rds/surv-example/submodel-three-output.rds"
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
  spread_draws(plot_mu[patient_id, plot_x])

res$plot_x <- submodel_three_settings$x_plot[res$plot_x]

posterior_plot_data <- res %>% 
  point_interval(
    .width = 0.8,
    .point = mean,
    .interval = qi,
    .exclude = c(".chain", ".iteration", ".draw", ".row", "event_indicator")
  )

event_df <- tibble(
  patient_id = 1 : submodel_one_settings$n_patients,
  event_indicator = submodel_one_settings$event_indicator
)  

posterior_plot_data <- posterior_plot_data %>%
  left_join(event_df, by = c("patient_id" = "patient_id")) %>% 
  mutate(event_indicator = as.factor(event_indicator))

with_post_mean <- base_plot + 
  geom_line(
    data = posterior_plot_data, 
    aes(x = plot_x, y = plot_mu, col = event_indicator),
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
  )

ggsave_fullpage(
  filename = "plots/surv-example/submodel-three-posterior.pdf",
  plot = with_post_mean
)
