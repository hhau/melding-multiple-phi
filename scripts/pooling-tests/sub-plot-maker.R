library(mvtnorm)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(latex2exp)

source("scripts/pooling-tests/density-functions.R")

make_patchwork_output_plot <- function(pooling_type, lambda_one_value) {
  pooled_joint_function <- do.call(
    paste0("f_normal_", pooling_type, "_generator"),
    args = list(
      outer_lambda_1 = lambda_one_value,
      outer_lambda_2 = 1 - 2 * (lambda_one_value)
    )
  )

  plot_df <- expand.grid(
    phi_12 = seq(
      from = -phi_12_plot_limit,
      to = phi_12_plot_limit,
      length.out = n_plot_points
    ),
    phi_23 = seq(
      from = -phi_23_plot_limit,
      to = phi_23_plot_limit,
      length.out = n_plot_points
    )
  ) %>% 
    mutate(
      twod_marginal = f_normal_twod_marginal(cbind(phi_12, phi_23)),
      d_type = "twod_marginal",
      phi_12_marginal = f_normal_phi_12_marginal(phi_12),
      phi_23_marginal = f_normal_phi_23_marginal(phi_23)
    )

  pooled_df <- plot_df %>% 
    select(phi_12, phi_23) %>% 
    mutate(
      twod_marginal = as.numeric(pooled_joint_function(
        rbind(phi_12, phi_23)
      )),
      d_type = "pooled_joint",
      phi_12_marginal = NA,
      phi_23_marginal = NA
    )

  full_df <- bind_rows(plot_df, pooled_df) %>%
    as_tibble()

  # contour_breaks <- quantile(
  #   full_df$twod_marginal,
  #   probs = seq(
  #     from = 0.8,
  #     to = 1,
  #     length.out = n_contour_breaks + 2
  #   )[2 : (n_contour_breaks + 1)]
  # )

  base_plot <- full_df %>% 
    select(-c(phi_12_marginal, phi_23_marginal)) %>%
    ggplot(aes(x = phi_12, y = phi_23, z = twod_marginal)) +
    geom_contour(
      aes(group = d_type, colour = d_type, alpha =..level..),
      breaks = contour_breaks,
    ) +
    scale_x_continuous(limits = c(-phi_12_plot_limit, phi_12_plot_limit)) +
    scale_y_continuous(limits = c(-phi_23_plot_limit, phi_23_plot_limit)) + 
    scale_alpha(range = c(0.5, 1)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = rel(0.75))
    ) +
    scale_discrete_manual(
      "color",
      values = c(
        "pooled_joint" = highlight_col,
        "twod_marginal" = blues[2]
      )
    ) +
    ggtitle(bquote(
      .(str_to_title(pooling_type)) ~ "pooling," ~ lambda[1] == ~ .(round(lambda_one_value, 3))
    ))

  phi_12_marginal_plot <- plot_df %>% 
    distinct(phi_12, phi_12_marginal) %>% 
    ggplot(aes(x = phi_12, y = phi_12_marginal)) +
    geom_line(colour = blues[2]) +
    scale_x_continuous(limits = c(-phi_12_plot_limit, phi_12_plot_limit)) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) 
  
  if (lambda_one_value == 0.5) {
    phi_12_marginal_plot <- phi_12_marginal_plot + xlab(TeX("$\\phi_{1 \\bigcap 2}$"))
  } else {
    phi_12_marginal_plot <- phi_12_marginal_plot + theme(axis.title.x = element_blank())
  }  

  phi_23_marginal_plot <- plot_df %>% 
    distinct(phi_23, phi_23_marginal) %>% 
    ggplot(aes(x = phi_23, y = phi_23_marginal)) +
    geom_line(colour = blues[2]) +
    scale_x_continuous(limits = c(-phi_23_plot_limit, phi_23_plot_limit)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    ) +
    coord_flip()
  
  if (pooling_type == "logarithmic") {
    phi_23_marginal_plot <- phi_23_marginal_plot + xlab(TeX("$\\phi_{2 \\bigcap 3}$"))
  } else {
    phi_23_marginal_plot <- phi_23_marginal_plot + theme(axis.title.y = element_blank())
  }

  p1 <- (phi_23_marginal_plot + base_plot + plot_spacer() + phi_12_marginal_plot +
    plot_layout(
      ncol = 2,
      nrow = 2,
      widths = c(1, 9),
      heights = c(9, 1)
    )) & 
    theme(plot.margin = unit(c(0,0,0,0), units = "cm"))

  return(p1)

}
