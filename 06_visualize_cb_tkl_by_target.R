#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Visualize CB and TKL by Target Value  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Factor effects for the cb and tkl error methods (both target CFI) at each
# target level (0.99-0.95), using error_results.RData (from 01_generate-error.R).
# Same layout as figure 3 of 04_visualize_wb_by_target.R, except that the CFI
# row is replaced by RMSEA: cb and tkl hit their CFI target (essentially)
# exactly, so a CFI row would be flat by construction.
# Writes one ggplot2 figure per method to ./figures/cfi_by_target/

# Load packages
library(ggplot2); library(dplyr); library(tidyr)

# Create directory
dir.create("./figures/cfi_by_target", recursive = TRUE, showWarnings = FALSE)

# Load results
load("./error_results.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared setup ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Surface and ink
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink_2 <- "#52514e"; grid_col <- "#e6e5e1"

# Method colors (as in 02_visualize_error.R)
method_colors <- c(cb = "#2a78d6", tkl = "#1baf7a")

# Light -> dark steps of a method's color for the ordered target
# (CFI 0.99 -> 0.95 = least -> most misfit), built like the amber steps in
# 04_visualize_wb_by_target.R
color_steps <- function(base, n = 5){
  mix <- function(color, with, weight){
    rgb(t((1 - weight) * col2rgb(color) + weight * col2rgb(with)), maxColorValue = 255)
  }
  colorRampPalette(c(mix(base, "white", 0.45), base, mix(base, "black", 0.55)))(n)
}

# Theme
theme_cfi <- function(){
  theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = surface, color = NA),
      panel.background = element_rect(fill = surface, color = NA),
      panel.grid.major = element_line(color = grid_col, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      text = element_text(color = ink),
      axis.text = element_text(color = ink_2),
      axis.title = element_text(color = ink_2),
      strip.text = element_text(color = ink, face = "bold", hjust = 0),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = ink_2),
      plot.title.position = "plot",
      legend.position = "top",
      legend.justification = "left",
      plot.margin = margin(12, 16, 10, 12)
    )
}

# Save helper
save_plot <- function(plot, name, width, height){
  ggsave(
    file.path("figures", "cfi_by_target", name), plot, width = width, height = height,
    dpi = 200, bg = surface
  )
}

# Design factors and the outcomes shown
factor_names <- c(
  BLOCKS = "Blocks", NODES = "Nodes", WITHIN = "Within",
  BETWEEN = "Between", DIFFUSION = "Diffusion", N = "N"
)
outcomes <- c("RMSEA", "SRMR", "Max residual")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Factor effects at each target ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Median outcome at each factor level, one line per target
factor_effects <- function(method, title){

  # Method rows only, with labeled target (0.99 first = least misfit)
  data <- error_results %>%
    filter(method == !!method) %>%
    mutate(target = as.numeric(value))
  target_levels <- sort(unique(data$target), decreasing = TRUE)
  data <- data %>%
    mutate(target_f = factor(
      target, levels = target_levels, labels = paste("CFI =", target_levels)
    ))

  # Level order across factors
  level_order <- as.character(sort(unique(unlist(data[names(factor_names)]))))

  # Summarize
  factor_df <- bind_rows(lapply(names(factor_names), function(f){
    data %>%
      group_by(target_f, level = as.character(.data[[f]])) %>%
      summarize(
        RMSEA = median(RMSEA), SRMR = median(SRMR),
        `Max residual` = median(max_residual), .groups = "drop"
      ) %>%
      mutate(factor = factor_names[[f]])
  })) %>%
    pivot_longer(all_of(outcomes), names_to = "index", values_to = "median") %>%
    mutate(
      level = factor(level, levels = level_order),
      factor = factor(factor, levels = factor_names),
      index = factor(index, levels = outcomes)
    )

  # Plot
  ggplot(factor_df, aes(x = level, y = median, color = target_f, group = target_f)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2, shape = 21, fill = surface, stroke = 1.1) +
    facet_grid(index ~ factor, scales = "free", switch = "y") +
    scale_color_manual(values = color_steps(method_colors[[method]]), name = NULL) +
    labs(
      title = title,
      subtitle = "Median outcome at each factor level, one line per CFI target (rows share a y-scale)",
      x = NULL, y = NULL
    ) +
    theme_cfi() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

}

save_plot(
  factor_effects(
    "cb", "cb: at a fixed CFI, more blocks means less misfit; within and diffusion stay flat"
  ),
  "01_factor_effects_cb.png", 12, 6.5
)

save_plot(
  factor_effects(
    "tkl", "tkl: at a fixed CFI, more blocks means less misfit; within and diffusion stay flat"
  ),
  "02_factor_effects_tkl.png", 12, 6.5
)

cat("Saved figures to ./figures/cfi_by_target/\n")
