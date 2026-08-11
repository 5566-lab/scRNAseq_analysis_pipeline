#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
set.seed(cfg$project$seed)

existing <- cfg$inputs$existing_integrated_rds
out_integrated <- project_path(cfg, cfg$outputs$integrated_rds)
out_annotated <- project_path(cfg, cfg$outputs$annotated_rds)
ensure_parent_dir(out_integrated)
ensure_parent_dir(out_annotated)

if (file.exists(existing)) {
  message_step("Using existing integrated Seurat object: ", existing)
  obj <- readRDS(existing)
} else {
  prepared <- project_path(cfg, cfg$outputs$prepared_objects_rds)
  message_step("Building integrated object from prepared list: ", prepared)
  object_list <- readRDS(prepared)
  object_list <- lapply(object_list, standard_qc_filter, cfg = cfg)
  obj <- merge(object_list[[1]], y = object_list[-1], add.cell.ids = names(object_list))
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = cfg$integration$dims)

  if (tolower(cfg$integration$method) == "harmony") {
    obj <- RunHarmony(obj, group.by.vars = cfg$integration$batch_column)
    reduction <- "harmony"
  } else {
    reduction <- "pca"
  }

  obj <- FindNeighbors(obj, reduction = reduction, dims = seq_len(cfg$integration$dims))
  obj <- FindClusters(obj, resolution = cfg$integration$resolution)
  obj <- RunUMAP(obj, reduction = reduction, dims = seq_len(cfg$integration$dims))
}

safe_save_rds(obj, out_integrated)

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
safe_write_csv(markers, project_path(cfg, "results/tables/all_cluster_markers.csv"))

plot_dir <- project_path(cfg, "results/figures/integration")
ensure_dir(plot_dir)
if ("umap" %in% names(obj@reductions)) {
  p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  ggsave(file.path(plot_dir, "umap_clusters.pdf"), p1, width = 8, height = 6)
  if ("Source_GSE" %in% colnames(obj@meta.data)) {
    p2 <- DimPlot(obj, reduction = "umap", group.by = "Source_GSE")
    ggsave(file.path(plot_dir, "umap_source_gse.pdf"), p2, width = 8, height = 6)
  }
}

safe_save_rds(obj, out_annotated)
message_step("Saved integrated object: ", out_integrated)

