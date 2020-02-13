# Model 1

n_data <- 100

phi_12 <- rnorm(n = 1, mean = 1, sd = 1)
sigma_y_1 <- abs(rnorm(n = 1, mean = 0, sd = 1))
y_1 <- rnorm(n = n_data, mean = phi_12, sd = sigma_y_1)

data_list <- list(
  model_index = 1,
  phi_12 = phi_12,
  y_1 = y_1
)

saveRDS(
  file = "rds/ex-three-gaussians/model-1-data.rds",
  object = data_list
)