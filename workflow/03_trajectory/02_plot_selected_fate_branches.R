#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(patchwork)
  library(tidydr)
})

base_dir <- Sys.getenv(
  "SCRNA_MONOCLE_DIR",
  unset = "/public3/DSC/single_cell/Result/figer_new/monocle3"
)
cds <- readRDS(Sys.getenv(
  "SCRNA_PSEUDOTIME_RDS",
  unset = file.path(base_dir, "data_pseudotime.rds")
))

# These columns preserve the published interactive selection. The frozen
# data_pseudotime object is required because a later cds_MM_foam checkpoint
# contains a different Fate 2 selection.
metadata <- as.data.frame(SummarizedExperiment::colData(cds))
required_fate_columns <- c("subset_1", "subset_2")
if (!all(required_fate_columns %in% colnames(metadata))) {
  stop("The saved cds_MM_foam object does not contain subset_1/subset_2 labels.")
}

fate1_n <- sum(as.character(metadata$subset_1) == "Fate 1", na.rm = TRUE)
fate2_n <- sum(as.character(metadata$subset_2) == "Fate 2", na.rm = TRUE)
if (fate1_n == 0L || fate2_n == 0L) {
  stop("The saved subset_1/subset_2 columns do not contain Fate 1/Fate 2 cells.")
}
message(sprintf("Using original saved selections: Fate 1 = %d; Fate 2 = %d", fate1_n, fate2_n))

fate_colors <- c(
  "Fate 1" = "blue",
  "Fate 2" = "red",
  "Unselected" = "gray80"
)

branch_theme <- theme_dr() +
  theme(
    aspect.ratio = 0.90,
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    text = element_text(family = "Arial"),
    legend.position = "none"
  )

p1 <- plot_cells(
  cds,
  color_cells_by = "subset_1",
  label_cell_groups = FALSE,
  cell_size = 1,
  trajectory_graph_segment_size = 0.75
) +
  coord_flip() +
  scale_color_manual(
    name = "Cell fate",
    values = fate_colors,
    limits = names(fate_colors),
    breaks = names(fate_colors),
    drop = FALSE
  ) +
  branch_theme +
  ggtitle("Cell fate 1") +
  annotate(
    "label",
    x = -8, y = -3,
    label = "Cell fate 1",
    hjust = 1.2, vjust = -1,
    size = 6, color = "black", fill = "grey90",
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.05, "lines")
  )

p2 <- plot_cells(
  cds,
  color_cells_by = "subset_2",
  label_cell_groups = FALSE,
  cell_size = 1,
  trajectory_graph_segment_size = 0.75
) +
  coord_flip() +
  scale_color_manual(
    name = "Cell fate",
    values = fate_colors,
    limits = names(fate_colors),
    breaks = names(fate_colors),
    drop = FALSE
  ) +
  branch_theme +
  ggtitle("Cell fate 2") +
  annotate(
    "label",
    x = -5.5, y = 2.5,
    label = "Cell fate 2",
    hjust = -0.5, vjust = 2,
    size = 6, color = "black", fill = "grey90",
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.05, "lines")
  )

legend_data <- data.frame(
  Cell_fate = factor(names(fate_colors), levels = names(fate_colors)),
  x = 1,
  y = seq_along(fate_colors)
)
legend_plot <- ggplot(
  legend_data,
  aes(x = x, y = y, color = Cell_fate)
) +
  geom_point(size = 4) +
  scale_color_manual(
    name = "Cell fate",
    values = fate_colors,
    limits = names(fate_colors),
    drop = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_void() +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA)
  )
shared_legend <- cowplot::get_legend(legend_plot)

p_selected_cells <- cowplot::plot_grid(
  p1,
  p2,
  shared_legend,
  nrow = 1,
  rel_widths = c(1, 1, 0.18),
  align = "h",
  axis = "tb"
)

ggsave(
  file.path(base_dir, "F2.7_MM_selected_cells.pdf"),
  p_selected_cells,
  width = 14,
  height = 7,
  device = cairo_pdf
)
ggsave(
  file.path(base_dir, "F2.7_MM_selected_cells.png"),
  p_selected_cells,
  width = 18,
  height = 8,
  dpi = 300,
  bg = "white",
  device = png,
  type = "cairo"
)

message("Saved fate-branch panels with trajectory_graph_segment_size = 0.75")
