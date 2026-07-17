#!/usr/bin/env Rscript
repo <- normalizePath(getwd(), mustWork = TRUE)
files <- list.files(repo, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("archive/original_scripts", files)]
for (file in files) {
  parse(file)
  message("OK ", file)
}

required <- c(
  "README.md",
  "configs/config.yaml",
  "metadata/sample_manifest.csv",
  "docs/code_function_archive.md",
  "docs/external_consensus_mmi_method.md",
  "data/gene_sets/external_consensus_mmi_no_c1q.csv",
  "data/gene_sets/external_consensus_mmi_no_c1q_sources.csv"
)
missing <- required[!file.exists(file.path(repo, required))]
if (length(missing) > 0) {
  stop("Missing required files: ", paste(missing, collapse = ", "), call. = FALSE)
}

mmi <- read.csv(
  file.path(repo, "data/gene_sets/external_consensus_mmi_no_c1q.csv"),
  stringsAsFactors = FALSE
)
required_mmi_columns <- c("direction", "gene_symbol")
if (!all(required_mmi_columns %in% colnames(mmi))) {
  stop("Frozen MMI table is missing required columns", call. = FALSE)
}
mature <- unique(mmi$gene_symbol[mmi$direction == "mature_positive"])
immature <- unique(mmi$gene_symbol[mmi$direction == "immature_negative"])
if (length(mature) != 476L || length(immature) != 476L) {
  stop("Frozen MMI must contain 476 genes in each direction", call. = FALSE)
}
if (length(intersect(mature, immature)) > 0L) {
  stop("Frozen MMI directions must not overlap", call. = FALSE)
}

mmi_sources <- read.csv(
  file.path(repo, "data/gene_sets/external_consensus_mmi_no_c1q_sources.csv"),
  stringsAsFactors = FALSE
)
expected_sources <- c(
  "GSE5099_Macrophage", "GSE5099_Monocyte",
  "GSE11864_Macrophage", "GSE11864_Monocyte",
  "HPCA_Macrophage", "HPCA_Monocyte"
)
if (!setequal(unique(mmi_sources$source_set), expected_sources)) {
  stop("External MMI source membership is incomplete or contains an unexpected module", call. = FALSE)
}

message("Repository validation passed")
