library(tibble)
library(parallel)
library(dplyr)
library(stringr)

source("scripts/pooling-tests/density-functions.R")

# parameters for plotting grid
n_grid_in_dir <- 100

phi_12_vals <- seq(from = -5, to = 5, length.out = n_grid_in_dir)
phi_23_vals <- seq(from = -5, to = 5, length.out = n_grid_in_dir)
plot_tbl <- as_tibble(expand.grid(phi_12 = phi_12_vals, phi_23 = phi_23_vals))

# density function parameters
mean_p1 <- c(-2)
mean_p2 <- matrix(c(0, 0), nrow = 2)
mean_p3 <- c(2)

sigma_mat <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
t_df <- 4

weight_cases <- list(
  "lambda[1]~'='~0.00" = c(outer_lambda_1 = 0.00, outer_lambda_2 = 1 - 2 * 0.00),
  "lambda[1]~'='~0.10" = c(outer_lambda_1 = 0.10, outer_lambda_2 = 1 - 2 * 0.10),
  "lambda[1]~'='~0.20" = c(outer_lambda_1 = 0.20, outer_lambda_2 = 1 - 2 * 0.20),
  "lambda[1]~'='~0.33" = c(outer_lambda_1 = 0.33333, outer_lambda_2 = 1 - 2 * 0.33333),
  "lambda[1]~'='~0.45" = c(outer_lambda_1 = 0.45, outer_lambda_2 = 1 - 2 * 0.45),
  "lambda[1]~'='~0.50" = c(outer_lambda_1 = 0.50, outer_lambda_2 = 1 - 2 * 0.50)
)

base_dists <- list(
  "normal",
  "student_t"
)

pooling_types <- list(
  "log",
  "linear"
)

res <- mclapply(X = base_dists, mc.cores = 2, function(a_base_dist) {
  mclapply(X = names(weight_cases), mc.cores = 2, function(a_weight_case) {
    mclapply(X = pooling_types, mc.cores = 2, function(a_pooling_type) {
      generator_name <- paste(
        c("f", a_base_dist, a_pooling_type, "generator"),
        collapse = "_"
      )
      density_fuction <- do.call(generator_name, args = as.list(weight_cases[[a_weight_case]]))
      f_val <- do.call(density_fuction, list(phi = rbind(plot_tbl$phi_12, plot_tbl$phi_23)))
      sub_res <- tibble(
        x = plot_tbl$phi_12,
        y = plot_tbl$phi_23,
        weight_case = a_weight_case,
        pooling_type = str_to_title(a_pooling_type),
        base_dist = a_base_dist,
        f_val = f_val
      )
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

saveRDS(
  file = "rds/pooling-tests/pooling-types-weights-results.rds",
  object = res
)
