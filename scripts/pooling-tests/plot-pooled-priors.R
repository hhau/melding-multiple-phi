library(ggplot2)
library(purrr)
library(gridExtra)

source("scripts/common/plot-settings.R")
source("scripts/pooling-tests/sub-plot-maker.R")

n_contour_breaks <- 10
# magical numbers ahead:
contour_breaks <- exp(
  seq(from = -10, to = -2.25, length.out = n_contour_breaks)
)

n_plot_points <- 100
phi_12_plot_limit <- 8
phi_23_plot_limit <- 8

mean_p1 <- -2.5
mean_p3 <- 2.5
sigma_mat <- matrix(c(1, 0.8, 0.8, 1), ncol = 2, nrow = 2)

n_lambda <- 5
pooling_methods <- c("logarithmic", "linear")
lambda_one_values <- seq(from = 0, to = 0.5, length.out = n_lambda )
control_df <- expand.grid(
  pooling_method = pooling_methods, 
  lambda_one_value = lambda_one_values
)

res <- mapply(
  make_patchwork_output_plot,
  control_df$pooling_method,
  control_df$lambda_one_value,
  SIMPLIFY = FALSE
)

plot_output <- wrap_plots(
  res,
  ncol = 2,
  nrow = n_lambda,
  byrow = TRUE
)

ggsave_fullpage(
  filename = "plots/pooling-tests/version-two.pdf",
  plot = plot_output
)
