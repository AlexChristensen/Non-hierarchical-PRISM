#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Visualize the WB Method       ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Detailed breakdown of the wb error method (targets RMSEA at 0.01-0.05)
# using error_results.RData (from 01_generate_error.R).
# Writes five ggplot2 figures to ./figures/wb/
#
# wb hits its RMSEA target essentially exactly (|error| < 1e-11), so the
# interesting questions are what it does to the *other* outcomes and which
# design factors drive that.

# Load packages
library(ggplot2); library(dplyr); library(tidyr)

# Create directory
dir.create("./figures/wb", recursive = TRUE, showWarnings = FALSE)

# Load results
load("./error_results.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared setup ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Surface and ink
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink_2 <- "#52514e"; grid_col <- "#e6e5e1"

# wb's categorical color (amber) with three lightness steps for ordered
# factors (blocks / nodes): light -> dark = low -> high
wb_color <- "#eda100"
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
    file.path("figures", "wb", name), plot, width = width, height = height,
    dpi = 200, bg = surface
  )
}

# wb rows only
wb <- error_results %>%
  filter(method == "wb") %>%
  mutate(
    target = as.numeric(value),
    BLOCKS_f = factor(BLOCKS, labels = paste(sort(unique(BLOCKS)), "blocks")),
    NODES_f = factor(NODES, labels = paste(sort(unique(NODES)), "nodes"))
  )

# Sanity check: wb hits the RMSEA target
stopifnot(max(abs(wb$RMSEA - wb$target)) < 1e-9)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 1. Implied CFI by blocks and nodes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Precomputed box statistics (so hidden outliers do not stretch the scale)
cfi_box <- wb %>%
  group_by(BLOCKS_f, NODES_f, target) %>%
  summarize(stats = list(boxplot.stats(CFI)$stats), .groups = "drop") %>%
  mutate(
    ymin = sapply(stats, `[`, 1), lower = sapply(stats, `[`, 2),
    middle = sapply(stats, `[`, 3), upper = sapply(stats, `[`, 4),
    ymax = sapply(stats, `[`, 5),
    target = factor(target)
  )

p1 <- ggplot(cfi_box, aes(x = target, fill = NODES_f, color = NODES_f)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = ink_2) +
  geom_boxplot(
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.35, linewidth = 0.5, width = 0.8,
    position = position_dodge(width = 0.85)
  ) +
  facet_wrap(~ BLOCKS_f) +
  scale_fill_manual(values = shades) +
  scale_color_manual(values = shades) +
  labs(
    title = "The CFI that wb implies falls sharply with more blocks",
    subtitle = "Achieved CFI by RMSEA target (outliers hidden); dashed line = CFI of 0.95",
    x = "Target RMSEA", y = "Achieved CFI", fill = NULL, color = NULL
  ) +
  theme_wb()

save_plot(p1, "01_cfi_by_blocks_nodes.png", 11, 4.8)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 2. How often CFI is acceptable ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Share of replications at or above two common CFI cutoffs
cutoff_df <- wb %>%
  group_by(BLOCKS_f, NODES_f, target) %>%
  summarize(
    `CFI >= .95` = mean(CFI >= 0.95), `CFI >= .90` = mean(CFI >= 0.90),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(`CFI >= .95`, `CFI >= .90`), names_to = "cutoff", values_to = "share"
  ) %>%
  mutate(cutoff = factor(cutoff, levels = c("CFI >= .95", "CFI >= .90")))

p2 <- ggplot(cutoff_df, aes(x = target, y = share, color = NODES_f)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3, shape = 21, fill = surface, stroke = 1.2) +
  facet_grid(cutoff ~ BLOCKS_f) +
  scale_color_manual(values = shades) +
  scale_x_continuous(breaks = sort(unique(wb$target))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "With 8 blocks, an RMSEA of 0.03 or higher often fails CFI cutoffs",
    subtitle = "Share of replications meeting each CFI cutoff, by RMSEA target",
    x = "Target RMSEA", y = "Share of replications", color = NULL
  ) +
  theme_wb()

save_plot(p2, "02_cfi_cutoff_share.png", 10, 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 3. SRMR and max residual ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Median and middle 80% by target, blocks, and nodes
local_df <- wb %>%
  select(BLOCKS_f, NODES_f, target, SRMR, `Max residual` = max_residual) %>%
  pivot_longer(c(SRMR, `Max residual`), names_to = "index", values_to = "estimate") %>%
  group_by(index, BLOCKS_f, NODES_f, target) %>%
  summarize(
    median = median(estimate),
    lower = quantile(estimate, 0.10), upper = quantile(estimate, 0.90),
    .groups = "drop"
  ) %>%
  mutate(index = factor(index, levels = c("SRMR", "Max residual")))

p3 <- ggplot(local_df, aes(x = target, y = median, color = BLOCKS_f, fill = BLOCKS_f)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2, shape = 21, fill = surface, stroke = 1.2) +
  facet_grid(index ~ NODES_f, scales = "free_y") +
  scale_color_manual(values = shades) +
  scale_fill_manual(values = shades) +
  scale_x_continuous(breaks = sort(unique(wb$target))) +
  labs(
    title = "Local misfit grows with the RMSEA target and, more weakly, with blocks",
    subtitle = "Median and middle 80% of SRMR and max residual by RMSEA target",
    x = "Target RMSEA", y = NULL, color = NULL, fill = NULL
  ) +
  theme_wb()

save_plot(p3, "03_srmr_residual_by_design.png", 11, 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 4. Which design factors matter ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Median outcome at each level of each factor, pooled over all RMSEA targets
# (the design is balanced, so levels differ only in the factor itself)
factor_names <- c(
  BLOCKS = "Blocks", NODES = "Nodes", WITHIN = "Within",
  BETWEEN = "Between", DIFFUSION = "Diffusion", N = "N"
)
# Numeric order of every level so each panel sorts correctly
level_order <- as.character(sort(unique(unlist(wb[names(factor_names)]))))

factor_df <- bind_rows(lapply(names(factor_names), function(f){
  wb %>%
    group_by(level = as.character(.data[[f]])) %>%
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

p4 <- ggplot(factor_df, aes(x = level, y = median, group = 1)) +
  geom_line(color = wb_color, linewidth = 0.8) +
  geom_point(size = 2.8, shape = 21, color = wb_color, fill = surface, stroke = 1.3) +
  facet_grid(index ~ factor, scales = "free", switch = "y") +
  labs(
    title = "Blocks drive wb's side effects most; nodes and between add some",
    subtitle = paste0(
      "Median outcome at each factor level, pooled over RMSEA targets; ",
      "within, diffusion, and N barely matter (y-axes do not start at zero)"
    ),
    x = NULL, y = NULL
  ) +
  theme_wb() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

save_plot(p4, "04_factor_effects.png", 12, 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 5. Where community recovery fails ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Miss rate of each algorithm (Louvain never missed) by blocks, nodes, and between
miss_df <- wb %>%
  group_by(BLOCKS, NODES, BETWEEN) %>%
  summarize(Walktrap = mean(!walktrap), PRISM = mean(!prism), .groups = "drop") %>%
  pivot_longer(c(Walktrap, PRISM), names_to = "algorithm", values_to = "miss") %>%
  mutate(
    algorithm = factor(algorithm, levels = c("Walktrap", "PRISM")),
    BLOCKS = factor(BLOCKS), NODES = factor(NODES),
    BETWEEN = factor(BETWEEN, labels = paste("Between =", sort(unique(BETWEEN)))),
    label = sprintf("%.1f%%", 100 * miss)
  )

p5 <- ggplot(miss_df, aes(x = NODES, y = BLOCKS, fill = miss)) +
  geom_tile(color = surface, linewidth = 1.2) +
  geom_text(aes(label = label), color = ink, size = 3.6) +
  facet_grid(algorithm ~ BETWEEN) +
  scale_fill_gradient(
    low = "#fdf3d6", high = "#d99200",
    labels = scales::percent_format(accuracy = 1), name = "Miss rate"
  ) +
  labs(
    title = "Walktrap misses 4-node blocks; PRISM's rare misses cluster at 8 blocks",
    subtitle = paste0(
      "Share of replications failing to recover the true communities ",
      "(Louvain missed none)"
    ),
    x = "Nodes per block", y = "Blocks"
  ) +
  theme_wb() +
  theme(panel.grid = element_blank(), legend.position = "right")

save_plot(p5, "05_recovery_misses.png", 10, 5)

cat("Saved figures to ./figures/wb/\n")
