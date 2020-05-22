library(mvtnorm)
library(parallel)
library(tidyverse)

source("scripts/common/plot-settings.R")

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
  f_args <- .extend_list(
    dens_obj$density_f_args,
    list(x = rep(point, dimension))
  )
  value <- do.call(dens_obj$density_f_name, f_args)
  res <- tibble(
    d = dimension,
    density_f_val = value,
    point = point,
    density_type = dens_obj$name
  )

  return(res)
  
}

# plot arguments
dimensions <- 1 : 50
test_points <- c(0, 2)
density_functions <- list(
  list(
    name = "Gaussian",
    density_f_name = "dmvnorm",
    density_f_args = list(log = FALSE)
  ),
  list(
    name = "Student-t",
    density_f_name = "dmvt",
    density_f_args = list(df = 4, log = FALSE)
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

p1 <- ggplot(
  data = res,
  aes(x = d, y = density_f_val, lty = as_factor(point), colour = density_type)
) + 
  geom_line() +
  # geom_point() +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "Gaussian" = blues[2],
      "Student-t" = highlight_col
    ),
    labels = c( 
      "Gaussian" = "Gaussian",
      "Student-t" = expression("Student" ~ "-" ~ italic("t")[4])
    )
  ) +
  ylab(expression("log"(italic("f")('x'["d"])))) +
  labs(
    lty = expression("x"["d"]),
    colour = "Density"
  ) +
  theme(
    legend.text.align = 0,
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

ggsave_halfheight(
  file = "plots/pooling-tests/densities-vs-dimension.pdf",
  plot = p1
)