library(ggplot2)
library(bayesplot)

# Theme settings
theme_set(theme_classic())
theme_replace(
  panel.grid.major = element_line(),
  panel.grid.minor = element_line(linetype = "dashed", size = rel(2/3)),
  legend.text = element_text(size = rel(1.1)),
  legend.title = element_text(size = rel(1.1))
)

bayesplot_theme_set(theme_classic())
bayesplot_theme_replace(
  panel.grid.major = element_line(),
  panel.grid.minor = element_line(linetype = "dashed", size = rel(2/3)),
  legend.text = element_text(size = rel(1.1)),
  legend.title = element_text(size = rel(1.1))
)

# Colours - should all be darkest[1] to lightest[n] 
blues <- c(
  "#00214F",
  "#2C7FB8",
  "#A6E9FF"
) 

greens <- c(
  "#364723",
  "#557438",
  "#749C4D",
  "#98CD65"
)

## Burgundy(ish) red highlight
highlight_col <- wesanderson::wes_palette("FantasticFox1")[5]

# ggplot saving settings - @mbertolacci
display_settings <- list(
  full_page_plot_width = 15,
  full_page_plot_height = 21,
  half_page_plot_width = 7,
  half_page_plot_height = 10,
  png_plot_dpi = 300,
  highlight_colour = highlight_col
)

ggsave_base <- function(filename, plot, bg = 'transparent', ...) {
  ggsave(
    filename,
    plot,
    units = 'cm',
    dpi = display_settings$png_plot_dpi,
    bg = bg,
    ...
  )
}

ggsave_fullwidth <- function(filename, plot, ...) {
  ggsave_base(
    filename,
    plot,
    width = display_settings$full_page_plot_width,
    ...
  )
}

ggsave_fullpage <- function(filename, plot, adjust_height = 0, ...) {
  ggsave_base(
    filename,
    plot,
    width = display_settings$full_page_plot_width,
    height = display_settings$full_page_plot_height + adjust_height,
    ...
  )
}

ggsave_halfheight <- function(filename, plot, ...) {
  ggsave_base(
    filename,
    plot,
    width = display_settings$full_page_plot_width,
    height = display_settings$half_page_plot_height,
    ...
  )
}

parsed_map <- function(map) {
  function(x) parse(text = map[x])
}
