#!/usr/bin/env Rscript
source("R/utils/config.R")

cfg_path <- parse_config_arg()
steps <- c(
  "scripts/01_prepare_seurat_objects.R",
  "scripts/02_integrate_cluster.R",
  "scripts/03_hdWGCNA_mo_ma.R",
  "scripts/04_macSpectrum_scores.R",
  "scripts/04b_auc_macrophage_signatures.R",
  "scripts/05_monocle_pseudotime.R",
  "scripts/06_gsea_gsva.R"
)

for (step in steps) {
  message_step("Running ", step)
  status <- system2("Rscript", c(step, "--config", cfg_path))
  if (!identical(status, 0L)) {
    stop("Pipeline failed at ", step, call. = FALSE)
  }
}

message_step("Pipeline completed")
