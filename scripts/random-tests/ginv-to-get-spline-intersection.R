library(splines2)
library(gtools)
library(ggplot2)
library(parallel)
library(dplyr)
library(matlib)

# generalised inverse testing

flat_vec <- -abs(rnorm(n = ncol(X_mat), mean = 1 / flat_const, sd = 2)) 
down_vec <- -abs(rnorm(n = ncol(X_mat), mean = 1 / down_const, sd = 2))

flat_rand_intercept <- rnorm(n = 1, mean = 0, sd = 1)
down_rand_intercept <- rnorm(n = 1, mean = 0, sd = 1)

plot_df <- data.frame(
  x = rep(x_seq, 2),
  y = c(
    flat_rand_intercept + X_mat %*% flat_vec,
    down_rand_intercept + X_mat %*% down_vec
  ),
  col = rep(c("flat", "down"), each = n_point)
)

ggplot(plot_df, aes(x = x, y = y, col = col)) +
  geom_line()

y_star <- -10
target <- y_star - down_rand_intercept
g_beta_inv <- MASS::ginv(
  matrix(down_vec, ncol = 1)
)

lines_df <- data.frame(
  x = c(
    mean(target %*% g_beta_inv),
    gm_mean(target %*% g_beta_inv),
    hm_mean(target %*% g_beta_inv)
  ),
  point = c(
    "mean",
    "geo mean",
    "hm mean"
  )
)

ggplot(plot_df, aes(x = x, y = y, col = col)) +
  geom_line() +
  geom_vline(
    data = lines_df, 
    aes(xintercept = x, lty = point)
  )

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

hm_mean <- function(x) {
  1 / (mean(1 / x))
}

mean(target %*% g_beta_inv)
gm_mean(target %*% g_beta_inv)
hm_mean(target %*% g_beta_inv)

# nothing here is sensible nor should be used -- the only way to do this is 
# numerically
