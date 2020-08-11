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
    .interval = qi
  )

with_post_mean <- base_plot + 
  geom_line(
    data = posterior_plot_data, 
    aes(x = plot_x, y = plot_mu),
    inherit.aes = FALSE,
    col = highlight_col
  ) +
  geom_ribbon(
    data = posterior_plot_data, 
    aes(x = plot_x, ymin = .lower, ymax = .upper),
    inherit.aes = FALSE,
    alpha = 0.2
  )

with_post_mean

ggsave_fullpage(
  filename = "plots/surv-example/submodel-three-posterior.pdf",
  plot = with_post_mean

)
