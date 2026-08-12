#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
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

plot_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.75,
  cell_size = 1,
  alpha = 0.8
) +
  coord_flip() +
  facet_wrap(~Sample_Type, nrow = 1, scales = "free_y") +
  scale_color_gradientn(
    colours = c("blue", "cyan", "green", "yellow", "orange", "red"),
    name = "Pseudotime",
    guide = guide_colorbar(barwidth = 1.5, title.position = "top")
  ) +
  theme_dr() +
  theme(
    aspect.ratio = 0.90,
    strip.text = element_text(size = 20),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(
      size = 15,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    text = element_text(family = "Arial")
  ) +
  annotate(
    "label",
    x = -5.5, y = 2.5,
    label = "Cell fate 2",
    hjust = -0.5, vjust = 2,
    size = 6, color = "black", fill = "grey90",
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.05, "lines")
  ) +
  annotate(
    "label",
    x = -8, y = -3,
    label = "Cell fate 1",
    hjust = 1.2, vjust = -1,
    size = 6, color = "black", fill = "grey90",
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.05, "lines")
  )

ggsave(
  file.path(base_dir, "F2.5_MM_Pseudotim.pdf"),
  plot_pt,
  width = 14,
  height = 7,
  device = cairo_pdf
)
ggsave(
  file.path(base_dir, "F2.5_MM_Pseudotim.png"),
  plot_pt,
  width = 18,
  height = 8,
  dpi = 300,
  bg = "white",
  device = png,
  type = "cairo"
)

message("Saved global pseudotime panels with aspect.ratio = 0.90")
