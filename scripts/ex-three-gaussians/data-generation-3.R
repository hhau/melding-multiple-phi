# Model 3

n_data <- 100

phi_23 <- rnorm(n = 1, mean = 2, sd = 1)
sigma_y_3 <- abs(rnorm(n = 1, mean = 0, sd = 1))
y_3 <- rnorm(n = n_data, mean = phi_23, sd = sigma_y_3)

data_list <- list(
  model_index = 3,
  phi_23 = phi_23,
  y_3 = y_3
)

saveRDS(
  file = "rds/ex-three-gaussians/model-3-data.rds",
  object = data_list
)