#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Visualize WB by Target Value  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Breakdown of the wb error method (targets RMSEA) at each target level
# (0.01-0.05), using error_results.RData (from 01_generate_error.R).
# Complements 03_visualize_wb.R, which mostly pools over targets.
# Writes five ggplot2 figures to ./figures/wb_by_target/

# Load packages
library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

# Create directory
dir.create("./figures/wb_by_target", recursive = TRUE, showWarnings = FALSE)

# Load results
load("./error_results.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared setup ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Surface and ink
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink_2 <- "#52514e"; grid_col <- "#e6e5e1"

# wb's color (amber) in light -> dark steps: used for the ordered target
# (0.01 -> 0.05) and for ordered design factors
wb_color <- "#eda100"
target_colors <- c("#f3c969", "#eeb22e", "#d99200", "#a86f00", "#6b4600")
shades <- c("#f3c969", "#eda100", "#8a5d00")

# Theme
theme_wb <- function(){
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
    file.path("figures", "wb_by_target", name), plot, width = width, height = height,
    dpi = 200, bg = surface
  )
}

# wb rows only, with labeled target
target_levels <- sort(unique(as.numeric(error_results$value[error_results$method == "wb"])))
wb <- error_results %>%
  filter(method == "wb") %>%
  mutate(
    target = as.numeric(value),
    target_f = factor(target, levels = target_levels, labels = paste("RMSEA =", target_levels))
  )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 1. Full distributions by target ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Empirical CDF of each outcome, one line per target
dist_long <- wb %>%
  select(target_f, CFI, SRMR, `Max residual` = max_residual) %>%
  pivot_longer(c(CFI, SRMR, `Max residual`), names_to = "index", values_to = "estimate") %>%
  mutate(index = factor(index, levels = c("CFI", "SRMR", "Max residual")))

# CFI cutoffs (CFI panel only)
cutoffs <- data.frame(
  index = factor("CFI", levels = levels(dist_long$index)), x = c(0.90, 0.95)
)

p1 <- ggplot(dist_long, aes(x = estimate, color = target_f)) +
  geom_vline(data = cutoffs, aes(xintercept = x), linetype = "dashed", color = ink_2) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_wrap(~ index, scales = "free_x") +
  scale_color_manual(values = target_colors, name = NULL) +
  scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Higher RMSEA targets shift and widen every outcome distribution",
    subtitle = "Cumulative share of replications at or below each value; dashed lines = CFI of 0.90 and 0.95",
    x = NULL, y = "Cumulative share"
  ) +
  theme_wb() +
  theme(panel.spacing.x = unit(1.4, "lines"))

save_plot(p1, "01_distributions_by_target.png", 11, 4.6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 2. Where CFI breaks down ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Median CFI in each blocks x nodes cell at each target
cell_df <- wb %>%
  group_by(target_f, BLOCKS, NODES) %>%
  summarize(CFI = median(CFI), .groups = "drop") %>%
  mutate(
    BLOCKS = factor(BLOCKS), NODES = factor(NODES),
    label = sprintf("%.2f", CFI)
  )

p2 <- ggplot(cell_df, aes(x = NODES, y = BLOCKS, fill = CFI)) +
  geom_tile(color = surface, linewidth = 1.2) +
  geom_text(aes(label = label), color = ink, size = 3.3) +
  facet_wrap(~ target_f, nrow = 1) +
  scale_fill_gradient(
    low = "#c98500", high = "#fdf3d6", name = "Median CFI",
    guide = guide_colorbar(barwidth = unit(10, "lines"))
  ) +
  labs(
    title = "Low CFI concentrates in many-block, many-node cells as the target rises",
    subtitle = "Median CFI by blocks and nodes per block, for each RMSEA target",
    x = "Nodes per block", y = "Blocks"
  ) +
  theme_wb() +
  theme(panel.grid = element_blank())

save_plot(p2, "02_cfi_cells_by_target.png", 13, 4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 3. Factor effects at each target ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Median outcome at each factor level, one line per target
factor_names <- c(
  BLOCKS = "Blocks", NODES = "Nodes", WITHIN = "Within",
  BETWEEN = "Between", DIFFUSION = "Diffusion", N = "N"
)
level_order <- as.character(sort(unique(unlist(wb[names(factor_names)]))))

factor_df <- bind_rows(lapply(names(factor_names), function(f){
  wb %>%
    group_by(target_f, level = as.character(.data[[f]])) %>%
    summarize(
      CFI = median(CFI), SRMR = median(SRMR),
      `Max residual` = median(max_residual), .groups = "drop"
    ) %>%
    mutate(factor = factor_names[[f]])
})) %>%
  pivot_longer(c(CFI, SRMR, `Max residual`), names_to = "index", values_to = "median") %>%
  mutate(
    level = factor(level, levels = level_order),
    factor = factor(factor, levels = factor_names),
    index = factor(index, levels = c("CFI", "SRMR", "Max residual"))
  )

p3 <- ggplot(factor_df, aes(x = level, y = median, color = target_f, group = target_f)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2, shape = 21, fill = surface, stroke = 1.1) +
  facet_grid(index ~ factor, scales = "free", switch = "y") +
  scale_color_manual(values = target_colors, name = NULL) +
  labs(
    title = "Blocks matter most, and more at higher targets; within and diffusion stay flat",
    subtitle = "Median outcome at each factor level, one line per RMSEA target (rows share a y-scale)",
    x = NULL, y = NULL
  ) +
  theme_wb() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

save_plot(p3, "03_factor_effects_by_target.png", 12, 6.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 4. Recovery misses by target ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Walktrap misses by nodes (its driver); PRISM misses by blocks (its driver)
walktrap_df <- wb %>%
  group_by(target, NODES = factor(NODES, labels = paste(sort(unique(NODES)), "nodes"))) %>%
  summarize(miss = mean(!walktrap), .groups = "drop")

prism_df <- wb %>%
  group_by(target, BLOCKS = factor(BLOCKS, labels = paste(sort(unique(BLOCKS)), "blocks"))) %>%
  summarize(miss = mean(!prism), .groups = "drop")

recovery_panel <- function(data, group, title){
  ggplot(data, aes(x = target, y = miss, color = .data[[group]])) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.4, shape = 21, fill = surface, stroke = 1.2) +
    scale_color_manual(values = shades, name = NULL) +
    scale_x_continuous(breaks = target_levels) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      limits = c(0, NA)
    ) +
    labs(title = title, x = "Target RMSEA", y = "Miss rate") +
    theme_wb() +
    theme(plot.title = element_text(face = "bold", size = 12))
}

p4 <- (recovery_panel(walktrap_df, "NODES", "Walktrap: flat in the target, driven by nodes") |
  recovery_panel(prism_df, "BLOCKS", "PRISM: rises with the target, driven by blocks")) +
  plot_annotation(
    title = "Community recovery misses by RMSEA target",
    subtitle = "Share of replications failing to recover the true communities (Louvain missed none)",
    theme = theme_wb()
  )

save_plot(p4, "04_recovery_misses_by_target.png", 11, 4.6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 5. Global vs. local misfit by target ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Sample of replications per target (full data is too dense to read)
set.seed(1)
sample_df <- wb %>%
  group_by(target_f) %>%
  slice_sample(n = 2000) %>%
  ungroup()

p5 <- ggplot(sample_df, aes(x = CFI, y = max_residual)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = ink_2) +
  geom_point(color = wb_color, alpha = 0.3, size = 1) +
  facet_wrap(~ target_f, nrow = 1) +
  labs(
    title = "At a fixed RMSEA, CFI and max residual still vary widely, more so at higher targets",
    subtitle = "2,000 sampled replications per target; dashed line = CFI of 0.95",
    x = "Achieved CFI", y = "Max residual"
  ) +
  theme_wb()

save_plot(p5, "05_cfi_vs_residual_by_target.png", 13, 4.2)

cat("Saved figures to ./figures/wb_by_target/\n")
