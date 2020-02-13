library(rstan)

# data 
data_model_1 <- readRDS(
  file = "rds/ex-three-gaussians/model-1-data.rds"
)

data_model_2 <- readRDS(
  file = "rds/ex-three-gaussians/model-2-data.rds"
)

data_model_3 <- readRDS(
  file = "rds/ex-three-gaussians/model-3-data.rds"
)

stan_prefit <- stan_model("scripts/stan-files/ex-three-gaussians-melded-joint.stan")

stan_data <- list(
  n_1 = length(data_model_1$y_1),
  n_2 = length(data_model_2$y_2),
  n_3 = length(data_model_3$y_3),
  y_1 = data_model_1$y_1,
  y_2 = data_model_2$y_2,
  y_3 = data_model_3$y_3
)

model_fit <- sampling(
  stan_prefit, 
  data = stan_data
)

model_fit
