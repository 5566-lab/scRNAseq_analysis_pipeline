#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
input_rds <- project_path(cfg, cfg$outputs$macspectrum_rds)
if (!file.exists(input_rds)) {
  input_rds <- project_path(cfg, cfg$outputs$hdwgcnna_rds)
}
message_step("Loading object for Monocle3: ", input_rds)
obj <- readRDS(input_rds)

expr <- GetAssayData(obj, assay = "RNA", slot = "counts")
cell_metadata <- obj@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expr),
  row.names = rownames(expr)
)

cds <- new_cell_data_set(expr, cell_metadata = cell_metadata, gene_metadata = gene_metadata)

modules <- GetModules(obj)
selected_modules <- unlist(cfg$monocle3$selected_modules)
selected_genes <- modules |>
  filter(color %in% selected_modules) |>
  pull(gene_name) |>
  unique()
selected_genes <- selected_genes[selected_genes %in% rowData(cds)$gene_short_name]

if (length(selected_genes) > 0) {
  cds <- preprocess_cds(cds, num_dim = cfg$monocle3$num_dim, use_genes = selected_genes)
} else {
  cds <- preprocess_cds(cds, num_dim = cfg$monocle3$num_dim)
}

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.n_neighbors = cfg$monocle3$umap_neighbors,
  umap.min_dist = cfg$monocle3$umap_min_dist
)

if ("umap" %in% names(obj@reductions)) {
  seurat_umap <- Embeddings(obj, reduction = "umap")
  common_cells <- intersect(colnames(cds), rownames(seurat_umap))
  reducedDims(cds)$UMAP[common_cells, ] <- seurat_umap[common_cells, ]
}

cds <- cluster_cells(cds, resolution = cfg$monocle3$cluster_resolution)
cds <- learn_graph(
  cds,
  close_loop = FALSE,
  learn_graph_control = list(
    minimal_branch_len = cfg$monocle3$minimal_branch_len,
    prune_graph = TRUE
  )
)
cds <- order_cells(cds)

obj$pseudotime <- pseudotime(cds)[colnames(obj)]
out_rds <- project_path(cfg, cfg$outputs$pseudotime_rds)
safe_save_rds(obj, out_rds)

out_dir <- project_path(cfg, "results/monocle3")
ensure_dir(out_dir)
saveRDS(cds, file.path(out_dir, "cds_monocyte_macrophage.rds"))

cell_col <- cfg$cell_types$cell_type_column
p_celltype <- plot_cells(cds, color_cells_by = cell_col, show_trajectory_graph = TRUE)
ggsave(file.path(out_dir, "F2.5_monocle_celltype_cluster.pdf"), p_celltype, width = 18, height = 8)

p_pt <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
ggsave(file.path(out_dir, "F2.5_MM_Pseudotime.pdf"), p_pt, width = 18, height = 8)

if ("APOBEC3A" %in% rowData(cds)$gene_short_name) {
  p_a3a <- plot_cells(cds, genes = "APOBEC3A", show_trajectory_graph = TRUE)
  ggsave(file.path(out_dir, "F2.5_monocle_APOBEC3A.pdf"), p_a3a, width = 18, height = 8)
}

trace_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = cfg$pathway$ncores)
safe_write_csv(trace_genes, file.path(out_dir, "trace_genes.csv"), row.names = TRUE)

message_step("Saved Monocle3 outputs: ", out_dir)

