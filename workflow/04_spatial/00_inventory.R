#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
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
xenium <- readRDS(cfg$inputs$xenium_rds)
message("Loading Visium object")
visium <- readRDS(cfg$inputs$visium_rds)

dataset_inventory <- tibble(
  dataset = c("GSE315246_Xenium", "GSE314851_Visium_FFPE"),
  platform = c("Xenium", "10x Visium FFPE"),
  n_features = c(nrow(xenium), nrow(visium)),
  n_observations = c(ncol(xenium), ncol(visium)),
  assays = c(paste(Assays(xenium), collapse = ";"), paste(Assays(visium), collapse = ";")),
  reductions = c(paste(Reductions(xenium), collapse = ";"), paste(Reductions(visium), collapse = ";")),
  images = c(paste(Images(xenium), collapse = ";"), paste(Images(visium), collapse = ";"))
)
write.csv(dataset_inventory, file.path(tab_dir, "dataset_inventory.csv"), row.names = FALSE)

xenium_meta_summary <- xenium@meta.data %>%
  mutate(disease = factor(disease, levels = severity_levels)) %>%
  count(disease, predicted.id, name = "n") %>%
  group_by(disease) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
write.csv(xenium_meta_summary, file.path(tab_dir, "xenium_celltype_by_disease.csv"), row.names = FALSE)

visium_meta_summary <- visium@meta.data %>%
  mutate(category = factor(category, levels = severity_levels)) %>%
  count(category, sample, cell.type.max, name = "n") %>%
  group_by(category, sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
write.csv(visium_meta_summary, file.path(tab_dir, "visium_celltype_by_sample_category.csv"), row.names = FALSE)

signature_availability <- bind_rows(
  available_signature_table(xenium, "GSE315246_Xenium"),
  available_signature_table(visium, "GSE314851_Visium_FFPE")
)
write.csv(signature_availability, file.path(tab_dir, "signature_gene_availability.csv"), row.names = FALSE)

raw_archives <- tibble(
  archive = c("GSE277441_RAW.tar", "GSE277170_RAW.tar", "GSE283269_RAW.tar"),
  platform_hint = c("Xenium cell-level raw files", "NanoString GeoMx DSP DCC files", "10x Visium/CytAssist raw files"),
  path = file.path(cfg$inputs$spatial_root, archive)
) %>%
  rowwise() %>%
  mutate(exists = file.exists(path), size_gb = ifelse(exists, round(file.info(path)$size / 1024^3, 3), NA_real_)) %>%
  ungroup()
write.csv(raw_archives, file.path(tab_dir, "raw_archive_inventory.csv"), row.names = FALSE)

print(dataset_inventory)
print(signature_availability)
