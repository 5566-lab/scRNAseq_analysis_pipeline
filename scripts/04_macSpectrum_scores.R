#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(macSpectrum)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
input_rds <- project_path(cfg, cfg$outputs$hdwgcnna_rds)
message_step("Loading hdWGCNA object: ", input_rds)
obj <- readRDS(input_rds)

expr <- GetAssayData(obj, assay = cfg$macspectrum$assay, slot = cfg$macspectrum$slot)
mac_mtx <- as.data.frame(as.matrix(expr)) |>
  tibble::rownames_to_column("geneid")

ensembl_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = mac_mtx$geneid,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)
mac_mtx$geneid <- ensembl_ids
mac_mtx <- na.omit(mac_mtx)

cell_col <- cfg$cell_types$cell_type_column
obj@meta.data[[cell_col]] <- normalize_foam_labels(obj@meta.data[[cell_col]])
feature <- factor(obj@meta.data[[cell_col]])
names(feature) <- rownames(obj@meta.data)

score <- macspec(mac_mtx, feature, select_hu_mo = cfg$macspectrum$species)
obj@meta.data <- merge(
  obj@meta.data,
  score[, c("MPI", "AMDI")],
  by.x = "row.names",
  by.y = "row.names",
  all.x = TRUE
) |>
  tibble::column_to_rownames("Row.names")

out_rds <- project_path(cfg, cfg$outputs$macspectrum_rds)
safe_save_rds(obj, out_rds)

plot_dir <- project_path(cfg, "results/figures/macSpectrum")
ensure_dir(plot_dir)
ggsave(file.path(plot_dir, "FeaturePlot_MPI.pdf"), FeaturePlot(obj, features = "MPI"), width = 8, height = 6)
ggsave(file.path(plot_dir, "FeaturePlot_AMDI.pdf"), FeaturePlot(obj, features = "AMDI"), width = 8, height = 6)
ggsave(file.path(plot_dir, "VlnPlot_MPI_AMDI.pdf"), VlnPlot(obj, features = c("MPI", "AMDI"), group.by = cell_col, pt.size = 0), width = 10, height = 6)

message_step("Saved macSpectrum-scored object: ", out_rds)

