#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
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

message("Loading Visium object")
obj <- readRDS(cfg$inputs$visium_rds)
obj$severity <- factor(obj$category, levels = severity_levels)
DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
obj <- add_signature_scores(obj, assay = DefaultAssay(obj), prefix = "sig")

score_cols <- grep("^sig_", colnames(obj@meta.data), value = TRUE)
write.csv(
  kruskal_table(obj@meta.data, score_cols, "severity"),
  file.path(tab_dir, "visium_signature_severity_kruskal.csv"),
  row.names = FALSE
)

sample_tab <- obj@meta.data %>%
  distinct(sample, severity) %>%
  count(severity, name = "n_samples")
spot_tab <- obj@meta.data %>%
  count(severity, name = "n_spots")
overview_tab <- full_join(sample_tab, spot_tab, by = "severity")
write.csv(overview_tab, file.path(tab_dir, "visium_sample_spot_overview.csv"), row.names = FALSE)

p_overview <- overview_tab %>%
  pivot_longer(c(n_samples, n_spots), names_to = "metric", values_to = "n") %>%
  ggplot(aes(x = severity, y = n, fill = severity)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = severity_cols) +
  labs(x = NULL, y = "Count", title = "Visium FFPE dataset overview") +
  theme_pub() +
  theme(legend.position = "none")
save_figure(p_overview, file.path(fig_dir, "Fig1A_Visium_dataset_overview.pdf"), width = 6.4, height = 3.4)

comp <- composition_table(obj@meta.data, "severity", "cell.type.max")
write.csv(comp, file.path(tab_dir, "visium_celltype_composition_by_severity.csv"), row.names = FALSE)
p_comp <- plot_composition(comp, "Visium dominant spot cell type by plaque severity")
save_figure(p_comp, file.path(fig_dir, "Fig1B_Visium_celltype_composition_by_severity.pdf"), width = 7.2, height = 4.8)

p_scores <- plot_score_by_severity(obj@meta.data, score_cols, "severity", "Visium spatial transcriptomic programs by plaque severity")
save_figure(p_scores, file.path(fig_dir, "Fig2A_Visium_signature_scores_by_severity.pdf"), width = 10.8, height = 6.4)

marker_genes <- c("APOBEC3A", "LGALS3", "PLA2G7", "SPP1", "TREM2", "LPL", "APOE", "FABP5", "GPNMB", "IL1B", "S100A8", "S100A9", "IFIT1", "ISG15", "MX1", "LYVE1", "SELENOP")
marker_genes <- present_genes(obj, marker_genes)
marker_expr <- fetch_data_layer(obj, marker_genes, layer = "data")
marker_summary <- bind_cols(obj@meta.data %>% select(severity, sample, cell.type.max), as_tibble(marker_expr)) %>%
  pivot_longer(all_of(marker_genes), names_to = "gene", values_to = "expr") %>%
  group_by(severity, gene) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), pct_expr = mean(expr > 0, na.rm = TRUE), .groups = "drop")
write.csv(marker_summary, file.path(tab_dir, "visium_marker_expression_by_severity.csv"), row.names = FALSE)

top_marker_plot <- marker_summary %>%
  mutate(severity = factor(severity, levels = severity_levels)) %>%
  ggplot(aes(x = severity, y = gene, size = pct_expr, color = mean_expr)) +
  geom_point(alpha = 0.9) +
  scale_color_viridis_c(option = "magma") +
  scale_size_continuous(labels = scales::percent_format(accuracy = 1), range = c(1.3, 6.5)) +
  labs(x = NULL, y = NULL, color = "Mean expression", size = "Expressing spots", title = "Visium marker expression across plaque severity") +
  theme_pub()
save_figure(top_marker_plot, file.path(fig_dir, "Fig2B_Visium_marker_dotplot_by_severity.pdf"), width = 7.2, height = 5.6)

coord_df <- collect_tissue_coordinates(obj) %>%
  rename(barcode = cell) %>%
  bind_cols(obj@meta.data[match(.$barcode, rownames(obj@meta.data)), c("sample", "severity", "cell.type.max", score_cols), drop = FALSE])

rep_samples <- obj@meta.data %>%
  count(severity, sample, sort = TRUE) %>%
  group_by(severity) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pull(sample)
write.csv(tibble(representative_sample = rep_samples), file.path(tab_dir, "visium_representative_samples_for_spatial_maps.csv"), row.names = FALSE)

coord_rep <- coord_df %>%
  filter(sample %in% rep_samples) %>%
  mutate(facet_label = paste0(sample, " (", severity, ")"))
p_space_type <- ggplot(coord_rep, aes(x = imagecol, y = -imagerow, color = cell.type.max)) +
  geom_point(size = 0.42, alpha = 0.9) +
  facet_wrap(~ facet_label, nrow = 1) +
  coord_equal() +
  scale_color_manual(values = celltype_cols, na.value = "grey75") +
  labs(x = NULL, y = NULL, color = "Dominant cell type", title = "Representative Visium spatial domains") +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_figure(p_space_type, file.path(fig_dir, "Fig3A_Visium_representative_spatial_celltypes.pdf"), width = 12, height = 4.3)

space_scores <- coord_rep %>%
  select(sample, severity, facet_label, imagecol, imagerow, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "signature", values_to = "score") %>%
  mutate(signature = sub("^sig_", "", signature)) %>%
  filter(signature %in% c("APOBEC3A_axis", "Foam_LAM", "ISG_myeloid", "Inflammatory_myeloid", "Resident_LYVE1_TRM"))
p_space_scores <- ggplot(space_scores, aes(x = imagecol, y = -imagerow, color = score)) +
  geom_point(size = 0.34, alpha = 0.9) +
  facet_grid(signature ~ facet_label) +
  coord_equal() +
  scale_color_viridis_c(option = "magma") +
  labs(x = NULL, y = NULL, color = "Score", title = "Visium spatial localization of myeloid and APOBEC3A-related programs") +
  theme_void(base_size = 8) +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))
save_figure(p_space_scores, file.path(fig_dir, "Fig3B_Visium_representative_spatial_signature_scores.pdf"), width = 12.5, height = 8.8)

saveRDS(obj@meta.data[, c("sample", "severity", "cell.type.max", score_cols), drop = FALSE], file.path(tab_dir, "visium_metadata_with_signature_scores.rds"))
message("Visium figures complete")
