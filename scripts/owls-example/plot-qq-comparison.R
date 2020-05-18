source("scripts/common/plot-settings.R")

library(dplyr)
library(ggplot2)
library(tibble)

original_model_samples <- readRDS("rds/owls-example/original-ipm-samples.rds")
stage_two_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")

vars <- c("fec", "v[1]", "v[2]")
plot_quantiles <- seq(from = 0.001, to = 1 - 0.001, by = 0.001)
a_recode_vector <- c(
  "fec" = 'rho',
  "v[1]" = 'alpha[0]',
  "v[2]" = 'alpha[2]'
)

orig_samples_slim <- array(
  original_model_samples[, , vars],
  dim = c(
    dim(original_model_samples)[1] * dim(original_model_samples)[2],
    3
  ), 
  dimnames = list(
    NULL,
    vars
  )
)

meld_samples_slim <- array(
  stage_two_samples[, , vars],
  dim = c(
    dim(stage_two_samples)[1] * dim(stage_two_samples)[2],
    3
  ), 
  dimnames = list(
    NULL,
    vars
  )
)

res <- lapply(vars, function(a_var) {
  tibble(
    quantile = plot_quantiles,
    orig_quantile = quantile(orig_samples_slim[, a_var], probs = plot_quantiles),
    meld_quantile = quantile(meld_samples_slim[, a_var], probs = plot_quantiles),
    param = a_var
  )
}) %>% 
  bind_rows() %>% 
  group_by(param) %>% 
  mutate(
    standardised_orig_quantile = scale(orig_quantile),
    standardised_meld_quantile = scale(meld_quantile)
  )

res$param <- res$param %>% 
  recode(!!!a_recode_vector) 

plot_list <- lapply(unique(res$param), function(a_param) {
  plot_data <- filter(res, param == a_param)
  plot_limit_min <- min(
    min(plot_data$meld_quantile),
    min(plot_data$orig_quantile)
  )
  plot_limit_max <- max(
    max(plot_data$meld_quantile),
    max(plot_data$orig_quantile)
  )
  
  if (a_param == 'rho') {
    labs_and_theme <<- theme(
      aspect.ratio = 1,
      axis.title = element_blank()
    ) 
  } else if (a_param == 'alpha[0]') {
    labs_and_theme <<- theme(
      axis.title.x = element_blank(),
      aspect.ratio = 1
    ) 
  } else if (a_param == 'alpha[2]') {
    labs_and_theme <<- theme(
      aspect.ratio = 1,
      axis.title.y = element_blank()
    )
  }
  
  p1 <- ggplot(data = plot_data, aes(x = meld_quantile, y = orig_quantile)) +
    geom_line(colour = highlight_col, alpha = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(
      vars(param),
      ncol = 1,
      nrow = 1,
      labeller = label_parsed
    ) +
    scale_x_continuous(limits = c(plot_limit_min, plot_limit_max)) +
    scale_y_continuous(limits = c(plot_limit_min, plot_limit_max)) +
    coord_fixed() +
    labs_and_theme
  
  return(p1)
})

plot_1 <- plot_list[[1]] 
plot_2 <- plot_list[[2]] + ylab('IPM quantile')
plot_3 <- plot_list[[3]] + xlab("Melded quantile")

# load this afterwards because it redefines the + operator
library(patchwork)

combined_plot <- plot_2 + plot_3 + plot_1

ggsave_fullpage(
  filename = "plots/owls-example/orig-meld-qq-compare.pdf",
  plot = combined_plot, 
  adjust_height = -14
)
