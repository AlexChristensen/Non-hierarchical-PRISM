#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Visualize Error Results       ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Reads error_results.RData (from 01_generate_error.R) and writes
# four ggplot2 figures to ./figures/
#
# Design of the study (from 01_generate_error.R):
#   - Four error-generation methods, each aimed at one fit index:
#       cb, tkl -> CFI | wb -> RMSEA | hb -> SRMR
#   - `value` is the *target* level of that index
#   - CFI, TLI, RMSEA, SRMR, max_residual are what was actually achieved
#   - walktrap / louvain / prism are TRUE when the algorithm recovered
#     the true community structure

# Load packages
library(ggplot2); library(dplyr); library(tidyr)

# Create directory
dir.create("./figures", showWarnings = FALSE)

# Load results
load("./error_results.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared setup ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Categorical colors, assigned in fixed order (never cycled)
# (blue, orange, aqua, yellow: first four slots of the validated reference palette)
method_colors <- c(cb = "#2a78d6", hb = "#eb6834", tkl = "#1baf7a", wb = "#eda100")
algo_colors <- c(Walktrap = "#2a78d6", Louvain = "#eb6834", PRISM = "#1baf7a")

# Surface and ink
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink_2 <- "#52514e"; grid_col <- "#e6e5e1"

# Readable labels
metric_labels <- c(cfi = "CFI", rmsea = "RMSEA", srmr = "SRMR")
method_labels <- c(cb = "cb (CFI)", hb = "hb (SRMR)", tkl = "tkl (CFI)", wb = "wb (RMSEA)")

# Theme
theme_error <- function(){
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
      plot.caption = element_text(color = ink_2, hjust = 0),
      plot.title.position = "plot",
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      plot.margin = margin(12, 16, 10, 12)
    )
}

# Prepare data
results <- error_results %>%
  mutate(
    target = as.numeric(value),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    method = factor(method, levels = names(method_colors)),
    achieved = case_when(
      metric == "CFI" ~ CFI,
      metric == "RMSEA" ~ RMSEA,
      metric == "SRMR" ~ SRMR
    )
  )

# Save helper
save_plot <- function(plot, name, width, height){
  ggsave(
    file.path("figures", name), plot, width = width, height = height,
    dpi = 200, bg = surface
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 1. Achieved vs. target fit ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Mean achieved index at each target (dashed line = perfect agreement)
achieved_df <- results %>%
  group_by(metric, method, target) %>%
  summarize(
    mean = mean(achieved), sd = sd(achieved), .groups = "drop"
  )

p1 <- ggplot(achieved_df, aes(x = target, y = mean, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = ink_2) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd), width = 0, linewidth = 0.6
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = surface, stroke = 1.2) +
  facet_wrap(~ metric, scales = "free") +
  scale_color_manual(values = method_colors, labels = method_labels) +
  labs(
    title = "Error generation hits its target index, except hb for SRMR",
    subtitle = paste0(
      "Mean achieved value (±1 SD) by target; dashed line = perfect agreement. ",
      "cb and tkl overlap exactly"
    ),
    x = "Target value", y = "Achieved value"
  ) +
  theme_error()

save_plot(p1, "01_achieved_vs_target.png", 11, 4.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 2. Local misfit at equal global fit ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Median and middle 80% of the largest residual correlation
residual_df <- results %>%
  group_by(metric, method, target) %>%
  summarize(
    median = median(max_residual),
    lower = quantile(max_residual, 0.10),
    upper = quantile(max_residual, 0.90),
    .groups = "drop"
  )

p2 <- ggplot(residual_df, aes(x = target, y = median, color = method, fill = method)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = surface, stroke = 1.2) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(
    title = "Methods that reach the same global fit leave very different local misfit",
    subtitle = "Largest absolute residual by target; line = median, band = middle 80%",
    x = "Target value", y = "Max residual"
  ) +
  theme_error()

save_plot(p2, "02_max_residual_by_target.png", 11, 4.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 3. Other indices under the CFI methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# cb and tkl both target CFI: what happens to the other indices?
cfi_long <- results %>%
  filter(metric == "CFI") %>%
  select(method, target, TLI, RMSEA, SRMR) %>%
  pivot_longer(c(TLI, RMSEA, SRMR), names_to = "index", values_to = "estimate") %>%
  mutate(
    index = factor(index, levels = c("TLI", "RMSEA", "SRMR")),
    target = factor(target)
  )

p3 <- ggplot(cfi_long, aes(x = target, y = estimate, color = method, fill = method)) +
  geom_boxplot(
    alpha = 0.2, outlier.shape = NA, linewidth = 0.5, width = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ index, scales = "free_y") +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(
    title = "At the same CFI, cb and tkl differ mainly in the SRMR they imply",
    subtitle = paste0(
      "Distribution of other fit indices across CFI targets (outliers hidden); ",
      "TLI and RMSEA are nearly identical"
    ),
    x = "Target CFI", y = "Achieved value"
  ) +
  theme_error()

save_plot(p3, "03_other_indices_cfi_methods.png", 11, 4.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 4. Community recovery ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Recovery rate of each algorithm by nodes per block
recovery_df <- results %>%
  group_by(NODES, WITHIN, BETWEEN) %>%
  summarize(
    Walktrap = mean(walktrap), Louvain = mean(louvain), PRISM = mean(prism),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(Walktrap, Louvain, PRISM), names_to = "algorithm", values_to = "recovery"
  ) %>%
  mutate(
    algorithm = factor(algorithm, levels = names(algo_colors)),
    WITHIN = factor(WITHIN, labels = paste("Within =", sort(unique(WITHIN)))),
    BETWEEN = factor(BETWEEN, labels = paste("Between =", sort(unique(BETWEEN))))
  )

p4 <- ggplot(recovery_df, aes(x = NODES, y = recovery, color = algorithm)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = surface, stroke = 1.2) +
  facet_grid(WITHIN ~ BETWEEN) +
  scale_x_continuous(breaks = sort(unique(results$NODES))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = algo_colors) +
  labs(
    title = "Walktrap struggles with 4-node blocks; Louvain and PRISM stay near-perfect",
    subtitle = "Share of replications recovering the true communities (y-axis does not start at zero)",
    x = "Nodes per block", y = "Recovery rate"
  ) +
  theme_error()

save_plot(p4, "04_community_recovery.png", 10, 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 5. All four methods on common axes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Methods target different indices at different levels, so compare them on
# what they all produce: the max residual at a given *achieved* index value.
# One point per method x target (median of both axes); band = middle 80%
# of max residual.
common_df <- results %>%
  group_by(method, target) %>%
  summarize(
    `1 - CFI` = 1 - median(CFI), RMSEA = median(RMSEA), SRMR = median(SRMR),
    median = median(max_residual),
    lower = quantile(max_residual, 0.10),
    upper = quantile(max_residual, 0.90),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(`1 - CFI`, RMSEA, SRMR), names_to = "index", values_to = "achieved"
  ) %>%
  mutate(index = factor(index, levels = c("1 - CFI", "RMSEA", "SRMR")))

p5 <- ggplot(
  common_df, aes(x = achieved, y = median, color = method, fill = method)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = surface, stroke = 1.2) +
  facet_wrap(~ index, scales = "free_x") +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(
    title = "At equal global fit, cb leaves the least local misfit and hb the most",
    subtitle = paste0(
      "Max residual against the achieved index (median per target; band = middle 80%). ",
      "Worse fit is to the right in every panel"
    ),
    x = "Achieved index value", y = "Max residual"
  ) +
  theme_error()

save_plot(p5, "05_methods_common_axes.png", 11, 4.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 6. Community recovery by method ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Recovery rate by method, algorithm, and nodes per block
method_recovery_df <- results %>%
  group_by(method, NODES) %>%
  summarize(
    Walktrap = mean(walktrap), Louvain = mean(louvain), PRISM = mean(prism),
    .groups = "drop"
  ) %>%
  pivot_longer(
    c(Walktrap, Louvain, PRISM), names_to = "algorithm", values_to = "recovery"
  ) %>%
  mutate(
    algorithm = factor(algorithm, levels = names(algo_colors)),
    NODES = factor(NODES, labels = paste(sort(unique(NODES)), "nodes per block"))
  )

p6 <- ggplot(method_recovery_df, aes(x = method, y = recovery, color = algorithm)) +
  geom_point(
    size = 3, shape = 21, fill = surface, stroke = 1.2,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ NODES) +
  scale_x_discrete(labels = method_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = algo_colors) +
  labs(
    title = "Community recovery does not depend on the error method",
    subtitle = "Share of replications recovering the true communities (y-axis does not start at zero)",
    x = NULL, y = "Recovery rate"
  ) +
  theme_error() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(p6, "06_recovery_by_method.png", 10, 4.8)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 7. Every index under every method ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Like figure 3, but for all four methods. Each column is a method on its own
# target axis; each row is one achieved index on a y-scale shared across methods
# (so panels in a row are directly comparable). The row matching a method's own
# target index is the target itself (near-constant within each level).
index_labels <- c(
  CFI = "CFI", TLI = "TLI", RMSEA = "RMSEA", SRMR = "SRMR",
  max_residual = "Max residual"
)

all_long <- results %>%
  select(method, target, CFI, TLI, RMSEA, SRMR, max_residual) %>%
  pivot_longer(
    c(CFI, TLI, RMSEA, SRMR, max_residual),
    names_to = "index", values_to = "estimate"
  ) %>%
  mutate(
    index = factor(index, levels = names(index_labels), labels = index_labels),
    target = factor(target)
  )

# Precompute box statistics so hidden outliers do not stretch the y-scales
box_df <- all_long %>%
  group_by(method, index, target) %>%
  summarize(
    stats = list(boxplot.stats(estimate)$stats), .groups = "drop"
  ) %>%
  mutate(
    ymin = sapply(stats, `[`, 1), lower = sapply(stats, `[`, 2),
    middle = sapply(stats, `[`, 3), upper = sapply(stats, `[`, 4),
    ymax = sapply(stats, `[`, 5)
  )

p7 <- ggplot(box_df, aes(x = target, color = method, fill = method)) +
  geom_boxplot(
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.2, linewidth = 0.5, width = 0.7
  ) +
  facet_grid(
    index ~ method, scales = "free", switch = "y",
    labeller = labeller(method = method_labels)
  ) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(
    title = "How each error method shapes every fit index",
    subtitle = paste0(
      "Achieved values by the method's own target (outliers hidden); ",
      "rows share a y-scale, x-axes differ by method"
    ),
    x = "Target value (of the method's own index)", y = NULL
  ) +
  theme_error() +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    panel.spacing = unit(0.8, "lines")
  )

save_plot(p7, "07_all_indices_by_method.png", 12, 10)

cat("Saved figures to ./figures/\n")
