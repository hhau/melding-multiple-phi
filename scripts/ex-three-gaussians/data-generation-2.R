# Model 2

n_data <- 100

model_1_data <- readRDS(file = "rds/ex-three-gaussians/model-1-data.rds")
model_3_data <- readRDS(file = "rds/ex-three-gaussians/model-3-data.rds")

# need to read from 
phi_12 <- model_1_data$phi_12
phi_23 <- model_3_data$phi_23

sigma_y_2 <- abs(rnorm(n = 1, mean = 0, sd = 1))
y_2 <- rnorm(n = n_data, mean = phi_12 + phi_23, sd = sigma_y_2)

data_list <- list(
  model_index = 2,
  y_2 = y_2
)

saveRDS(
  file = "rds/ex-three-gaussians/model-2-data.rds",
  object = data_list
)