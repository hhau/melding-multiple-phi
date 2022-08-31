library(tidyverse)

log_pooling_lambda <- rep(1 / 2, times = 3)

# phi_vec should contain (phi_{1 \cap 2}, phi_{2 \cap 3}) =
# (alpha_0, alpha_2, rho)

q_pool_log <- function(phi_vec, log = TRUE) {
  log_res <- sum(dnorm(phi_vec[1 : 2], sd = 100, log = TRUE)) * 
    sum(log_pooling_lambda[1 : 2])
  
  log_res <- log_res + dunif(phi_vec[3], min = 0, max = 10, log = TRUE) *
    sum(log_pooling_lambda[2 : 3])
  
  if (log) log_res else exp(log_res)
}

target_log_pooling <- function(phi_vec, log = TRUE) {
  log_res <- q_pool_log(phi_vec)
  log_res <- log_res - 2 * sum(dnorm(phi_vec[1 : 2], sd = 100, log = TRUE))
  log_res <- log_res - 2 * dunif(phi_vec[3], min = 0, max = 10, log = TRUE)
  
  if (log) log_res else exp(log_res)
}

n_plot <- 20

plot_df <- expand.grid(
  alpha_0 = seq(from = -4, to = -2, length.out = n_plot),
  alpha_2 = seq(from = 1.5, to = 3.5, length.out = n_plot),
  rho = seq(from = 2, to = 3, length.out = n_plot)
) %>% 
  mutate(
    val = apply(cbind(alpha_0, alpha_2, rho), 1, target_log_pooling)
  )

exp(-diff(range(plot_df$val)))


# the uniform prior is irrelevant in the linear pooling case, actually I guess
# it was also irrelevant in the log pooling case as well.
# in this example, linear pooling is completely independent of the pooling
# weights, which is kinda interesting.
# with log pooling weights lambda_1 = lambda_2 = 1 / 2, this is exactly the 
# same as log pooling
q_pool_lin <- function(phi_vec, log = TRUE) {

}

target_lin_pooling <- function(phi_vec, log = TRUE) {
  log_res <- -sum(dnorm(phi_vec[1 : 2], mean = 0, sd = 100, log = TRUE))
  if (log) log_res else exp(log_res)
}

plot_df_lin <- expand.grid(
  alpha_0 = seq(from = -4, to = -2, length.out = n_plot),
  alpha_2 = seq(from = 1.5, to = 3.5, length.out = n_plot),
  rho = seq(from = 2, to = 3, length.out = n_plot)
) %>% 
  mutate(
    val = apply(cbind(alpha_0, alpha_2, rho), 1, target_lin_pooling)
  )

exp(-diff(range(plot_df_lin$val)))
hist(plot_df_lin$val)

ggplot(plot_df_lin, aes(x = alpha_0, y = alpha_2, z = val)) +
  geom_contour(aes(colour = after_stat(level)))
