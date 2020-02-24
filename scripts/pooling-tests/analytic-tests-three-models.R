# combine the two independent components into a single MVN
mean_component_one <- c(-2, 2)
cov_component_one <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

# leave the correlated, 2d component as it is 
mean_component_two <- c(0, 0)
cov_component_two <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)

cov_combined <- solve(solve(cov_component_one) + solve(cov_component_two))
mean_combined <- cov_combined %*% (
  solve(cov_component_one, mean_component_one) +
  solve(cov_component_two, mean_component_two)  
)

cor_combined <- cov_combined / cov_combined[1, 1]
