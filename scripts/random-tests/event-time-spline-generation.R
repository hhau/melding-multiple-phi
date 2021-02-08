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

matplot(X_mat, type = "l")

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

# generalise

