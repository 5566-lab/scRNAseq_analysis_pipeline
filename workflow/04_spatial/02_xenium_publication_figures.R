#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(FNN)
  library(yaml)
})

ROOT <- Sys.getenv(
  "AST_ROOT",
  unset = normalizePath("workflow/04_spatial", mustWork = TRUE)
)
source(file.path(ROOT, "R", "utils.R"))
cfg <- yaml::read_yaml(file.path(ROOT, "config", "config.yaml"))
fig_dir <- cfg$outputs$figures
tab_dir <- cfg$outputs$tables
ensure_dir(fig_dir)
ensure_dir(tab_dir)

message("Loading Xenium object")
obj <- readRDS(cfg$inputs$xenium_rds)
obj$severity <- factor(obj$disease, levels = severity_levels)
DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
obj <- add_signature_scores(obj, assay = DefaultAssay(obj), prefix = "sig")

score_cols <- grep("^sig_", colnames(obj@meta.data), value = TRUE)
write.csv(
  kruskal_table(obj@meta.data, score_cols, "severity"),
  file.path(tab_dir, "xenium_signature_severity_kruskal.csv"),
  row.names = FALSE
)

cell_tab <- obj@meta.data %>%
  count(severity, name = "n_cells")
write.csv(cell_tab, file.path(tab_dir, "xenium_cell_overview.csv"), row.names = FALSE)

p_overview <- ggplot(cell_tab, aes(x = severity, y = n_cells, fill = severity)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = severity_cols) +
  labs(x = NULL, y = "Cells", title = "Xenium cell-level dataset overview") +
  theme_pub() +
  theme(legend.position = "none")
save_figure(p_overview, file.path(fig_dir, "Fig4A_Xenium_dataset_overview.pdf"), width = 5.2, height = 3.7)

comp <- composition_table(obj@meta.data, "severity", "predicted.id")
write.csv(comp, file.path(tab_dir, "xenium_celltype_composition_by_severity.csv"), row.names = FALSE)
p_comp <- plot_composition(comp, "Xenium cell type composition by plaque severity")
save_figure(p_comp, file.path(fig_dir, "Fig4B_Xenium_celltype_composition_by_severity.pdf"), width = 7.2, height = 4.8)

p_scores <- plot_score_by_severity(obj@meta.data, score_cols, "severity", "Xenium cell-level programs by plaque severity")
save_figure(p_scores, file.path(fig_dir, "Fig5A_Xenium_signature_scores_by_severity.pdf"), width = 10.8, height = 6.4)

myeloid_scores <- obj@meta.data %>%
  filter(predicted.id == "Myeloid") %>%
  select(severity, predicted.id, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "signature", values_to = "score") %>%
  mutate(signature = sub("^sig_", "", signature))
write.csv(myeloid_scores, file.path(tab_dir, "xenium_myeloid_signature_scores_long.csv"), row.names = FALSE)

p_myeloid <- ggplot(myeloid_scores, aes(x = severity, y = score, fill = severity)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.12, alpha = 0.82) +
  geom_boxplot(width = 0.13, outlier.size = 0.1, linewidth = 0.18, alpha = 0.9) +
  facet_wrap(~ signature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = severity_cols) +
  labs(x = NULL, y = "Z-scored module expression", title = "Xenium myeloid cells carry severity-associated plaque programs") +
  theme_pub(9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_figure(p_myeloid, file.path(fig_dir, "Fig5B_Xenium_myeloid_signature_scores_by_severity.pdf"), width = 10.8, height = 6.4)

marker_genes <- c("APOBEC3A", "LGALS3", "PLA2G7", "SPP1", "TREM2", "LPL", "APOE", "GPNMB", "IL1B", "ISG15", "MX1", "LYVE1", "MRC1", "CD68", "CD14", "FCN1", "ACTA2", "TAGLN", "PECAM1", "VWF", "CD3D")
marker_genes <- present_genes(obj, marker_genes)
marker_expr <- fetch_data_layer(obj, marker_genes, layer = "data")
marker_summary <- bind_cols(obj@meta.data %>% select(severity, predicted.id), as_tibble(marker_expr)) %>%
  pivot_longer(all_of(marker_genes), names_to = "gene", values_to = "expr") %>%
  group_by(severity, predicted.id, gene) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), pct_expr = mean(expr > 0, na.rm = TRUE), .groups = "drop")
write.csv(marker_summary, file.path(tab_dir, "xenium_marker_expression_by_celltype_severity.csv"), row.names = FALSE)

p_marker <- marker_summary %>%
  filter(predicted.id %in% c("Myeloid", "SMCPericyte", "ModSMC", "Endothelium", "TCells", "Fibroblast1")) %>%
  mutate(severity = factor(severity, levels = severity_levels)) %>%
  ggplot(aes(x = severity, y = gene, size = pct_expr, color = mean_expr)) +
  geom_point(alpha = 0.9) +
  facet_wrap(~ predicted.id, ncol = 3) +
  scale_color_viridis_c(option = "magma") +
  scale_size_continuous(labels = scales::percent_format(accuracy = 1), range = c(0.7, 5.8)) +
  labs(x = NULL, y = NULL, color = "Mean expression", size = "Expressing cells", title = "Xenium marker programs by cell type and plaque severity") +
  theme_pub(8.5)
save_figure(p_marker, file.path(fig_dir, "Fig5C_Xenium_marker_dotplot_by_celltype_severity.pdf"), width = 10.8, height = 8.2)

coord_df <- collect_tissue_coordinates(obj) %>%
  mutate(cell = as.character(cell)) %>%
  bind_cols(obj@meta.data[match(.$cell, rownames(obj@meta.data)), c("severity", "predicted.id", score_cols), drop = FALSE]) %>%
  filter(!is.na(x), !is.na(y), !is.na(severity), !is.na(predicted.id))

set.seed(1)
rep_fovs <- coord_df %>%
  count(severity, image_name, name = "n") %>%
  group_by(severity) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()
write.csv(rep_fovs, file.path(tab_dir, "xenium_representative_fovs_for_spatial_maps.csv"), row.names = FALSE)

plot_cells <- coord_df %>%
  semi_join(rep_fovs, by = c("severity", "image_name")) %>%
  mutate(facet_label = paste0(image_name, " (", severity, ")")) %>%
  group_by(facet_label, predicted.id) %>%
  mutate(.rand = runif(n())) %>%
  arrange(.rand, .by_group = TRUE) %>%
  slice_head(n = 2500) %>%
  ungroup()
p_space_type <- ggplot(plot_cells, aes(x = x, y = -y, color = predicted.id)) +
  geom_point(size = 0.08, alpha = 0.55) +
  facet_wrap(~ facet_label, nrow = 1) +
  coord_equal() +
  scale_color_manual(values = celltype_cols, na.value = "grey75") +
  labs(x = NULL, y = NULL, color = "Cell type", title = "Xenium cell-level spatial organization across plaque severity") +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_figure(p_space_type, file.path(fig_dir, "Fig6A_Xenium_spatial_celltype_severity_sampled.pdf"), width = 12.2, height = 4.7)

score_space <- coord_df %>%
  semi_join(rep_fovs, by = c("severity", "image_name")) %>%
  filter(predicted.id == "Myeloid") %>%
  mutate(facet_label = paste0(image_name, " (", severity, ")")) %>%
  group_by(facet_label) %>%
  mutate(.rand = runif(n())) %>%
  arrange(.rand, .by_group = TRUE) %>%
  slice_head(n = 5000) %>%
  ungroup() %>%
  select(severity, facet_label, x, y, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "signature", values_to = "score") %>%
  mutate(signature = sub("^sig_", "", signature)) %>%
  filter(signature %in% c("APOBEC3A_axis", "Foam_LAM", "ISG_myeloid", "Inflammatory_myeloid", "Resident_LYVE1_TRM"))
p_space_scores <- ggplot(score_space, aes(x = x, y = -y, color = score)) +
  geom_point(size = 0.08, alpha = 0.65) +
  facet_grid(signature ~ facet_label) +
  coord_equal() +
  scale_color_viridis_c(option = "magma") +
  labs(x = NULL, y = NULL, color = "Score", title = "Xenium spatial localization of myeloid disease programs") +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(face = "bold"))
save_figure(p_space_scores, file.path(fig_dir, "Fig6B_Xenium_myeloid_spatial_signature_scores_sampled.pdf"), width = 11, height = 8.8)

message("Running nearest-neighbor enrichment")
nn_df <- coord_df %>% select(cell, image_name, x, y, severity, predicted.id)
neighbor_tab <- bind_rows(lapply(split(nn_df, nn_df$image_name), function(df) {
  if (nrow(df) < 25 || !any(df$predicted.id == "Myeloid")) return(tibble())
  k_use <- min(20, nrow(df) - 1)
  knn <- FNN::get.knn(as.matrix(df[, c("x", "y")]), k = k_use)
  myeloid_idx <- which(df$predicted.id == "Myeloid")
  neighbor_types <- matrix(df$predicted.id[knn$nn.index[myeloid_idx, ]], nrow = length(myeloid_idx))
  bind_rows(lapply(seq_along(myeloid_idx), function(i) {
    tab <- prop.table(table(factor(neighbor_types[i, ], levels = names(celltype_cols))))
    tibble(
      focal_cell = df$cell[myeloid_idx[i]],
      image_name = df$image_name[myeloid_idx[i]],
      severity = df$severity[myeloid_idx[i]],
      neighbor_type = names(tab),
      neighbor_fraction = as.numeric(tab)
    )
  }))
}))
write.csv(neighbor_tab, file.path(tab_dir, "xenium_myeloid_neighbor_fractions.csv"), row.names = FALSE)

neighbor_summary <- neighbor_tab %>%
  group_by(severity, neighbor_type) %>%
  summarise(mean_fraction = mean(neighbor_fraction), .groups = "drop") %>%
  filter(neighbor_type %in% c("Myeloid", "SMCPericyte", "ModSMC", "Endothelium", "TCells", "Fibroblast1", "Fibroblast2"))
write.csv(neighbor_summary, file.path(tab_dir, "xenium_myeloid_neighbor_summary.csv"), row.names = FALSE)

p_neighbor <- ggplot(neighbor_summary, aes(x = severity, y = mean_fraction, fill = severity)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.15) +
  facet_wrap(~ neighbor_type, ncol = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = severity_cols) +
  labs(x = NULL, y = "Mean neighbor fraction around myeloid cells", title = "Xenium myeloid-cell neighborhoods remodel with plaque severity") +
  theme_pub(9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_figure(p_neighbor, file.path(fig_dir, "Fig6C_Xenium_myeloid_neighborhood_by_severity.pdf"), width = 9.5, height = 5.6)

saveRDS(obj@meta.data[, c("severity", "predicted.id", score_cols), drop = FALSE], file.path(tab_dir, "xenium_metadata_with_signature_scores.rds"))
message("Xenium figures complete")
