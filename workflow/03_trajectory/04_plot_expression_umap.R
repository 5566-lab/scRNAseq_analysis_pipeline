#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
})

base_dir <- Sys.getenv(
  "SCRNA_MONOCLE_DIR",
  unset = "/public3/DSC/single_cell/Result/figer_new/monocle3"
)
out_dir <- file.path(base_dir, "Expression_UMAP_no_trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(Sys.getenv(
  "SCRNA_PSEUDOTIME_RDS",
  unset = file.path(base_dir, "data_pseudotime.rds")
))

top_20_all <- c(
  "CD36", "FABP5", "PFKP", "STAT1", "ISG15",
  "APOBEC3A", "CCL2", "TNF", "S100A6", "LGALS3",
  "CLEC10A", "CD302", "AXL", "AOAH", "CH25H",
  "PHACTR1", "FGL2", "THBS1", "F13A1", "AREG"
)

expression_scale <- scale_color_gradientn(
  colors = c("gray90", "red2", "red3"),
  values = c(0, 0.5, 1)
)

p_top20 <- plot_cells(
  cds,
  genes = top_20_all,
  label_groups_by_cluster = FALSE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  cell_size = 1,
  show_trajectory_graph = FALSE
) +
  coord_flip() +
  facet_wrap(~feature_label, ncol = 5, scales = "free") +
  expression_scale +
  labs(color = "Expression") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    legend.position = "right",
    text = element_text(family = "Arial")
  )

ggsave(file.path(out_dir, "S2.2_Top20_expression_UMAP_no_trajectory.pdf"),
       p_top20, width = 20, height = 16, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_dir, "S2.2_Top20_expression_UMAP_no_trajectory.png"),
       p_top20, width = 20, height = 16, dpi = 300, limitsize = FALSE)

p_apobec3a <- plot_cells(
  cds,
  genes = "APOBEC3A",
  label_groups_by_cluster = FALSE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  cell_size = 1,
  show_trajectory_graph = FALSE
) +
  coord_flip() +
  facet_wrap(~Sample_Type, ncol = 2, scales = "free") +
  scale_color_gradientn(
    colors = c("gray90", "red2", "red3"),
    values = c(0, 0.25, 1)
  ) +
  labs(title = "APOBEC3A Single-Cell Expression Map", color = "APOBEC3A") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 22, face = "bold"),
    strip.background = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    text = element_text(family = "Arial")
  )

ggsave(file.path(out_dir, "F2.5_APOBEC3A_AC_PA_no_trajectory.pdf"),
       p_apobec3a, width = 15, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir, "F2.5_APOBEC3A_AC_PA_no_trajectory.png"),
       p_apobec3a, width = 15, height = 6, dpi = 400)

message("Saved expression UMAPs without pseudotime trajectories to: ", out_dir)
