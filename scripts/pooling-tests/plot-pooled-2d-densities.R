library(tibble)
library(mvtnorm)
library(cubature)
library(pbapply)
library(parallel)
library(dplyr)

source("scripts/common/plot-settings.R")
source("scripts/pooling-tests/density-functions.R")

# going to do 2D
phi_12_vals <- seq(from = -5, to = 5, length.out = 250)
phi_23_vals <- seq(from = -5, to = 5, length.out = 250)
plot_tbl <- as_tibble(expand.grid(phi_12 = phi_12_vals, phi_23 = phi_23_vals))

weight_cases <- list(
  equal = c(outer_w_1 = 1/3, outer_w_2 = 1/3, outer_w_3 = 1/3),
  unequal = c(outer_w_1 = 1/9, outer_w_2 = 7/9, outer_w_3 = 1/9)
)

base_dists <- list(
  "normal",
  "student_t"
)

res <- mclapply(X = base_dists, mc.cores = 2, function(a_base_dist) {
  mclapply(X = names(weight_cases), mc.cores = 2, function(a_weight_case) {
    generator_name <- paste0("f_", a_base_dist, "_log_generator")
    density_fuction <- do.call(generator_name, args = as.list(weight_cases[[a_weight_case]]))
    f_val <- apply(cbind(plot_tbl$phi_12, plot_tbl$phi_23), 1, density_fuction)
    sub_res <- tibble(
      x = plot_tbl$phi_12,
      y = plot_tbl$phi_23,
      weight_case = a_weight_case,
      pooling_type = "logarithmic",
      base_dist = a_base_dist,
      f_val = f_val
    )
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

# facisnating! 
p_1 <- ggplot(res, aes(x = x, y = y, z = f_val)) +
  geom_raster(aes(fill = f_val)) +
  scale_fill_gradientn(colours = c(blues[1], blues[2], blues[3])) +
  geom_contour(colour = "white", alpha = 0.2, bins = 8) +
  facet_grid(rows = vars(weight_case), cols = vars(base_dist)) +
  xlab(expression(phi["1:2"])) + 
  ylab(expression(phi["2:3"])) +
  geom_hline(yintercept = 0, colour = highlight_col, alpha = 0.5) +  
  geom_vline(xintercept = 0, colour = highlight_col, alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(fill = expression(italic("f"))) +
  NULL

ggsave_halfheight(
  filename = "plots/pooling-tests/pooled-densities-2d.pdf",
  plot = p_1
)
