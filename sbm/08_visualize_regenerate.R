#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Visualize Base-Network Redraws ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Reads error_results.RData (from 01_generate-error.R) and writes two ggplot2
# figures to ./figures/regenerate/
#
# `regenerate` counts how many base networks were discarded before one was
# recoverable (or simulation succeeded) for a replication. It is a property of
# the base network, so it is identical for all 20 rows (4 methods x 5 targets)
# of a replication: keep one row per replication (7,200 base networks).

# Load packages
library(ggplot2); library(dplyr); library(tidyr)

# Create directory
dir.create("./figures/regenerate", recursive = TRUE, showWarnings = FALSE)

# Load results
load("./error_results.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared setup ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Surface and ink
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink_2 <- "#52514e"; grid_col <- "#e6e5e1"

# Purple: not one of the method or algorithm colors in 02_visualize_error.R
regen_color <- "#7a5fc4"

# Theme
theme_regen <- function(){
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
    file.path("figures", "regenerate", name), plot, width = width, height = height,
    dpi = 200, bg = surface
  )
}

# One row per replication (regenerate is constant across its 20 rows)
base <- error_results %>%
  filter(method == "cb", value == "0.99") %>%
  mutate(regenerate = as.numeric(regenerate))

stopifnot(nrow(base) == nrow(error_results) / 20)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 1. Redraws by design factor ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Share needing at least one redraw and mean redraws, at each factor level
factor_names <- c(
  BLOCKS = "Blocks", NODES = "Nodes", WITHIN = "Within",
  BETWEEN = "Between", DIFFUSION = "Diffusion"
)
level_order <- as.character(sort(unique(unlist(base[names(factor_names)]))))

factor_df <- bind_rows(lapply(names(factor_names), function(f){
  base %>%
    group_by(level = as.character(.data[[f]])) %>%
    summarize(
      `Share with 1+ redraws` = mean(regenerate > 0),
      `Mean redraws` = mean(regenerate), .groups = "drop"
    ) %>%
    mutate(factor = factor_names[[f]])
})) %>%
  pivot_longer(
    c(`Share with 1+ redraws`, `Mean redraws`), names_to = "index", values_to = "value"
  ) %>%
  mutate(
    level = factor(level, levels = level_order),
    factor = factor(factor, levels = factor_names),
    index = factor(index, levels = c("Share with 1+ redraws", "Mean redraws"))
  )

p1 <- ggplot(factor_df, aes(x = level, y = value, group = 1)) +
  geom_line(color = regen_color, linewidth = 0.8) +
  geom_point(size = 2.8, shape = 21, color = regen_color, fill = surface, stroke = 1.3) +
  facet_grid(index ~ factor, scales = "free", switch = "y") +
  labs(
    title = "Redraws concentrate in 4-node blocks and in weaker block structure",
    subtitle = paste0(
      "Base networks discarded before a usable one was found, by design factor ",
      "(7,200 base networks; y-axes start at zero)"
    ),
    x = NULL, y = NULL
  ) +
  expand_limits(y = 0) +
  theme_regen() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

save_plot(p1, "01_redraws_by_factor.png", 12, 5.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### 2. Redraws in every design cell ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Mean redraws in each of the 72 conditions: blocks x nodes inside each
# within x diffusion (rows) by between (columns) panel
cell_df <- base %>%
  group_by(BLOCKS, NODES, WITHIN, BETWEEN, DIFFUSION) %>%
  summarize(mean = mean(regenerate), .groups = "drop") %>%
  mutate(
    BLOCKS = factor(BLOCKS), NODES = factor(NODES),
    BETWEEN = factor(BETWEEN, labels = paste("Between =", sort(unique(BETWEEN)))),
    row = factor(
      paste0("Within = ", WITHIN, ", Diffusion = ", DIFFUSION),
      levels = paste0(
        "Within = ", rep(sort(unique(WITHIN)), each = 2),
        ", Diffusion = ", rep(sort(unique(DIFFUSION)), 2)
      )
    ),
    label = ifelse(mean == 0, "0", sprintf("%.2f", mean))
  )

p2 <- ggplot(cell_df, aes(x = NODES, y = BLOCKS, fill = mean)) +
  geom_tile(color = surface, linewidth = 1.2) +
  geom_text(aes(label = label), color = ink, size = 3.3) +
  facet_grid(row ~ BETWEEN) +
  scale_fill_gradient(low = "#f1edfa", high = "#9a82d9", name = "Mean redraws") +
  labs(
    title = "Redraws peak at many blocks of 4 nodes, and weak structure makes it worse",
    subtitle = "Mean number of discarded base networks per replication in each of the 72 conditions",
    x = "Nodes per block", y = "Blocks"
  ) +
  theme_regen() +
  theme(
    panel.grid = element_blank(), legend.position = "right",
    strip.text.y = element_text(angle = 0)
  )

save_plot(p2, "02_redraws_by_condition.png", 10, 8)

cat("Saved figures to ./figures/regenerate/\n")
