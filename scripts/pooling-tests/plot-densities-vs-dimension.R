library(mvtnorm)
library(parallel)
library(tidyverse)

.extend_list <- function(...) {
  lists <- list(...)
  output <- lists[[1]]
  for (value in lists[2 : length(lists)]) {
    for (name in names(value)) {
      output[[name]] <- value[[name]]
    }
  }
  return(output)
}

generate_plot_value <- function(point, dimension, dens_obj) {
  f_lower_dimension_args <- .extend_list(
    dens_obj$dens_f_args,
    list(x = rep(point, dimension))
  )
  lower_dimension_value <- do.call(dens_obj$dens_f_name, f_lower_dimension_args)
  f_higher_dimension_args <- .extend_list(
    dens_obj$dens_f_args,
    list(x = rep(point, dimension * 2))
  )
  higher_dimension_value <- do.call(dens_obj$dens_f_name, f_higher_dimension_args)
  
  res <- tibble(
    d = dimension,
    ratio_val = higher_dimension_value / lower_dimension_value,
    point = point,
    density_type = dens_obj$name
  )
  
  return(res)
  
}

dimensions <- 1 : 50
test_points <- c(0, 2)
density_functions <- list(
  list(
    name = "Gaussian",
    dens_f_name = "dmvnorm",
    dens_f_args = list(log = FALSE)
  ),
  list(
    name = "Student-t",
    dens_f_name = "dmvt",
    dens_f_args = list(df = 5, log = FALSE)
  )
)

res <- lapply(density_functions, function(dens_obj) {
  lapply(test_points, function(a_point) {
    mclapply(dimensions, mc.cores = 5, function(a_dim) {
      generate_plot_value(
        point = a_point,
        dimension = a_dim,
        dens_obj = dens_obj
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()

ggplot(data = res, aes(x = d, y = ratio_val, lty = as_factor(point), colour = density_type)) + 
  # geom_point() +
  geom_line() +
  scale_y_log10()
