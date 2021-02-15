library(splines2)
library(gtools)
library(ggplot2)
library(parallel)
library(dplyr)

source("scripts/common/plot-settings.R")

n_point <- 100
x_seq <- seq(from = 0, to = 1, length.out = n_point)
n_internal_knot <- 1

X_mat <- iSpline(
  x = x_seq,
  knots = seq(
    from = 0, 
    to = 1, 
    length.out = n_internal_knot + 2
  )[-c(1, n_internal_knot + 2)],
  Boundary.knots = c(0, 1)
)

X_prime_mat <- iSpline(
  x = x_seq,
  knots = seq(
    from = 0, 
    to = 1, 
    length.out = n_internal_knot + 2
  )[-c(1, n_internal_knot + 2)],
  Boundary.knots = c(0, 1),
  derivs = 1
)

matplot(X_mat, type = "l")
matplot(X_prime_mat, type = "l")

flat_const <- 1
down_const <- 1 / 3

n_mc <- 1000
res <- mclapply(1 : n_mc, mc.cores = 5, function(iter_id) {
  flat_vec <- -abs(rt(n = ncol(X_mat), df = 30) + 1) / flat_const
  down_vec <- -abs(rt(n = ncol(X_mat), df = 10) + 1) / down_const
  
  flat_rand_intercept <- rnorm(n = 1, mean = 0, sd = 1)
  down_rand_intercept <- rnorm(n = 1, mean = 0, sd = 1)
  
  plot_df <- data.frame(
    x = rep(x_seq, 2),
    y = c(
      flat_rand_intercept + X_mat %*% flat_vec,
      down_rand_intercept + X_mat %*% down_vec
    ),
    col = rep(c("flat", "down"), each = n_point),
    iter = iter_id
  )
}) %>%
  bind_rows() %>%
  as_tibble()

n_traj <- 150
rand_iter <- sample(1 : n_mc, size = n_traj)

lines_df <- res %>%
  filter(iter %in% rand_iter)

plot_df <- res %>% 
  group_by(x, col) %>%
  mutate(
    mean = mean(y),
    lower = quantile(y, 0.025),
    upper = quantile(y, 0.975)
  ) %>%
  select(x, col, mean, lower, upper) %>%
  distinct()

ggplot(plot_df, aes(x = x, col = col)) +
  geom_line(aes(y = mean), size = rel(1.5)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(
    data = lines_df, 
    aes(y = y, group = interaction(iter, col)),
    # lty = "dashed",
    alpha = 0.2
  ) + 
  geom_line(
    data = tibble(x = c(0, 1), y = c(-10, -10)),
    aes(x = x , y = y),
    inherit.aes = FALSE
  )

# numerically find event time 
down_vec <- -abs(rt(n = ncol(X_mat), df = 10) + 1) / down_const
down_rand_intercept <- rnorm(n = 1, mean = 0, sd = 1)
flat_vec <- -abs(rt(n = ncol(X_mat), df = 30) + 1) / flat_const
flat_rand_intercept <- rnorm(n = 1, mean = 0, sd = 1)

threshold <- -10

fixed_basis <- function(x) {
  iSpline(
    x = x,
    knots = seq(
      from = 0, 
      to = 1, 
      length.out = n_internal_knot + 2
    )[-c(1, n_internal_knot + 2)],
    Boundary.knots = c(0, 1)
  )
}

optm_function_down <- function(t) {
  spline_val_at_t <- down_rand_intercept + fixed_basis(t) %*% down_vec
  sq_dist <- (threshold - spline_val_at_t)^2  
}

optm_function_flat <- function(t) {
  spline_val_at_t <- flat_rand_intercept + fixed_basis(t) %*% flat_vec
  sq_dist <- (threshold - spline_val_at_t)^2  
}

plot_spline <- function(intercept, coefs) {
  x_plot <- seq(from = 0, to = 1, length.out = 250)
  vals <- intercept + fixed_basis(x_plot) %*% coefs
  plot(x = x_plot, y = vals, type = "l")
  abline(h = threshold, lty = "33")
}

res_down <- optim(
  fn = optm_function_down,
  # gr = function(x) deriv_basis(x) %*% down_vec,
  par = 0.5,
  method = "L-BFGS-B",
  lower = 0, 
  upper = 1
)

plot_spline(down_rand_intercept, down_vec)
abline(v = res_down$par, col = "red")

res_flat <- optim(
  fn = optm_function_flat,
  # gr = function(x) deriv_basis(x) %*% flat_vec,
  par = 0.5,
  method = "L-BFGS-B",
  lower = 0, 
  upper = 1
)

plot_spline(flat_rand_intercept, flat_vec)
abline(v = res_flat$par, col = "red")
