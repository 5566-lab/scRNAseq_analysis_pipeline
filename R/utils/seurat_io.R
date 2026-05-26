suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

create_seurat_from_counts <- function(counts, project, cfg) {
  CreateSeuratObject(
    counts = counts,
    project = project,
    min.cells = cfg$qc$min_cells,
    min.features = cfg$qc$min_features
  )
}

add_percent_mt <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

standard_qc_filter <- function(obj, cfg) {
  add_percent_mt(obj)
  subset(
    obj,
    subset = nFeature_RNA >= cfg$qc$min_features &
      nFeature_RNA <= cfg$qc$max_features &
      percent.mt <= cfg$qc$max_percent_mt
  )
}

read_10x_if_present <- function(path) {
  if (!dir.exists(path)) {
    stop("10X directory does not exist: ", path)
  }
  Read10X(data.dir = path)
}

safe_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  invisible(path)
}

safe_write_csv <- function(object, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(object, path, row.names = row.names)
  invisible(path)
}

normalize_foam_labels <- function(x) {
  x <- as.character(x)
  x[x == "FOAM_cells1"] <- "Foam cells1"
  x[x == "FOAM_cells2"] <- "Foam cells2"
  x
}

