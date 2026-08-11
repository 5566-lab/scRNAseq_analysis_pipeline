#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
set.seed(cfg$project$seed)

input_rds <- project_path(cfg, cfg$outputs$annotated_rds)
if (!file.exists(input_rds)) {
  input_rds <- project_path(cfg, cfg$outputs$integrated_rds)
}
message_step("Loading Seurat object: ", input_rds)
obj <- readRDS(input_rds)

cell_col <- cfg$cell_types$cell_type_column
labels <- unlist(cfg$cell_types$monocyte_macrophage_labels)
if (cell_col %in% colnames(obj@meta.data)) {
  obj@meta.data[[cell_col]] <- normalize_foam_labels(obj@meta.data[[cell_col]])
  obj <- subset(obj, cells = rownames(obj@meta.data)[obj@meta.data[[cell_col]] %in% normalize_foam_labels(labels)])
}

DefaultAssay(obj) <- "RNA"
obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = cfg$hdwgcnna$wgcna_name
)
obj <- MetacellsByGroups(
  obj,
  group.by = cfg$hdwgcnna$group_by,
  k = cfg$hdwgcnna$metacell_k,
  max_shared = 10,
  ident.group = cfg$hdwgcnna$group_by
)
obj <- NormalizeMetacells(obj)
obj <- SetDatExpr(obj, group_name = NULL, group.by = cfg$hdwgcnna$group_by, assay = "RNA", slot = "data")
obj <- ConstructNetwork(
  obj,
  soft_power = cfg$hdwgcnna$soft_power,
  minModuleSize = cfg$hdwgcnna$min_module_size,
  tom_name = cfg$hdwgcnna$wgcna_name
)
obj <- ModuleEigengenes(obj)
obj <- ModuleConnectivity(obj)
obj <- ResetModuleNames(obj, new_name = "MM")

out_dir <- project_path(cfg, "results/hdWGCNA/Mo_Ma")
ensure_dir(out_dir)
modules <- GetModules(obj)
safe_write_csv(modules, file.path(out_dir, "data_modules.csv"), row.names = FALSE)

hub_df <- GetHubGenes(obj, n_hubs = 25)
safe_write_csv(hub_df, file.path(out_dir, "data_genes.csv"), row.names = FALSE)
safe_save_rds(obj, project_path(cfg, cfg$outputs$hdwgcnna_rds))

pdf(file.path(out_dir, "data_Dendrogram.pdf"), width = 10, height = 8)
PlotDendrogram(obj, main = "Mo/Ma hdWGCNA Dendrogram")
dev.off()

pdf(file.path(out_dir, "data_ModuleFeaturePlot.pdf"), width = 15, height = 12)
print(ModuleFeaturePlot(obj, features = "hMEs", order = TRUE))
dev.off()

message_step("Saved hdWGCNA outputs: ", out_dir)

