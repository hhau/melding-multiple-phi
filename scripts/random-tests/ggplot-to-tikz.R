library(ggplot2)
library(tikzDevice)
library(dplyr)

source('scripts/common/plot-settings.R')
source('scripts/mimic-example/GLOBALS.R')

scaling_list <- readRDS('rds/mimic-example/submodel-1-pf-data-list-format.rds')
scaling_data <- readRDS('rds/mimic-example/submodel-1-pf-data-stan-format.rds')
param_samples <- readRDS('rds/mimic-example/submodel-1-pf-samples-long.rds')
event_time_samples <- readRDS('rds/mimic-example/submodel-1-event-times-samples-long.rds')
fluid_param_samples <- readRDS('rds/mimic-example/submodel-3-fluid-samples-long.rds')
fluid_plot_mu_samples <- readRDS('rds/mimic-example/fluid-data-piecewise-plot-mu-tbl.rds')

indiv <- 6
draw <- 1034
chain <- 3

orig_scale_plot_dfs <- function(indiv_id, draw, chain) {
  local_list <- scaling_list[[indiv_id]]
  local_draw <- param_samples %>% 
    filter(
      i == indiv_id,
       .iteration == draw,
      .chain == chain
    )
  
  scale_mu <- local_list$y_obs_mean
  scale_sd <- local_list$y_obs_sd
  
  spline_coef <- local_draw %>% 
    filter(.variable == 'spline_coef') %>% 
    pull(.value)
  
  beta_zero <- local_draw %>% 
    filter(.variable == 'beta_zero') %>% 
    pull(.value)
  
  beta_zero_orig <- (beta_zero * scale_sd) + scale_mu
  
  event_time <- event_time_samples %>% 
    filter(
      i == indiv_id,
      .iteration == draw,
      .chain == chain,
      .variable == 'event_time'
    ) %>% 
    pull(.value)
  
  
  point_df <- tibble(
    x = c(0, event_time),
    y = c(beta_zero_orig, GLOBAL_ARDS_THRESHOLD),
    type = c('beta_zero', 'event_time')
  )
  
  spline_val <- as.numeric(local_list$x_plot_mat %*% spline_coef)
  spline_val_orig <- ((beta_zero + spline_val) * scale_sd) + scale_mu
  spline_df <- tibble(
    x = local_list$x_plot_seq,
    y = spline_val_orig,
    type = 'spline_orig_scale'
  )
  
  res <- list(point_df = point_df, spline_df = spline_df)
  return(res)
}

list_of_dfs <- orig_scale_plot_dfs(indiv, draw, chain)

tikz(
  file = 'tex-input/mimic-example/0090-pf-schematic.tex',
  width = display_settings$half_page_plot_width / 2.54,
  height = (display_settings$half_page_plot_height / (2.54 * 2))
)

p1 <- ggplot(mapping = aes(x = x, y = y)) +
  geom_point(
    data = list_of_dfs$point_df
  ) +
  geom_line(
    data = list_of_dfs$spline_df
  ) +
  geom_abline(
    intercept = GLOBAL_ARDS_THRESHOLD,
    slope = 0,
    linetype = 'dashed'
  ) +
  geom_segment(
    data = tibble(
      x = list_of_dfs$point_df$x[2], y = 210,
      xend = list_of_dfs$point_df$x[2], yend = GLOBAL_ARDS_THRESHOLD
    ),
    mapping = aes(x = x, xend = xend, y = y, yend = yend),
    linetype = 'dotted'
  ) +
  xlab(r"{$t$}") +
  ylab("P/F") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = c(210, NA))

print(p1)
dev.off()

fluid_plot_line_tbl <- tibble(
  x = seq(from = 0, to = max(list_of_dfs$spline_df$x), length.out = 300)
)

kappa <- 15
eta_before <- 3.5
eta_after <- 1.5
eta_zero_raw <- 0.2
fluid_vals <- array(dim = nrow(fluid_plot_line_tbl))

for (ii in seq_len(nrow(fluid_plot_line_tbl))) {
  x <- fluid_plot_line_tbl$x[ii]
  if (x < kappa) {
    res <- eta_zero_raw + eta_before * kappa + eta_before * (x - kappa)
  } else {
    res <- eta_zero_raw + eta_before * kappa + eta_after * (x - kappa)
  }
  fluid_vals[ii] <- res
}

fluid_plot_line_tbl$y <- fluid_vals

ggplot(fluid_plot_line_tbl, aes(x = x , y = y)) +
  geom_line() +
  xlab(r"{$t$}") +
  ylab("Fluid") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
