#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

spatial_root <- Sys.getenv(
  "AST_SPATIAL_ROOT",
  unset = "/public3/DSC/single_cell/spatial"
)
repo_root <- Sys.getenv(
  "AST_ROOT",
  unset = normalizePath("workflow/04_spatial", mustWork = TRUE)
)
output_root <- Sys.getenv("AST_OUTPUT_ROOT", unset = file.path(repo_root, "results"))
tab_dir <- file.path(output_root, "tables")

co_detection <- function(mat, metadata, dataset, group_name = NA_character_, min_gene_pct = 0.002) {
  target <- "APOBEC3A"
  if (!target %in% rownames(mat)) {
    stop(target, " not present in matrix")
  }
  target_expr <- as.numeric(mat[target, ])
  target_pos <- target_expr > 0
  n_target_pos <- sum(target_pos)
  n_target_neg <- sum(!target_pos)
  binary <- mat > 0
  gene_detected <- Matrix::rowSums(binary)
  gene_pct <- gene_detected / ncol(mat)
  keep <- which(gene_pct >= min_gene_pct & rownames(mat) != target)

  both <- Matrix::rowSums(binary[keep, target_pos, drop = FALSE])
  gene_pos_target_neg <- gene_detected[keep] - both
  target_pos_gene_neg <- n_target_pos - both
  neither <- n_target_neg - gene_pos_target_neg
  mean_pos <- Matrix::rowMeans(mat[keep, target_pos, drop = FALSE])
  mean_neg <- Matrix::rowMeans(mat[keep, !target_pos, drop = FALSE])
  jaccard <- both / (n_target_pos + gene_detected[keep] - both)

  fisher_p <- vapply(seq_along(keep), function(i) {
    suppressWarnings(fisher.test(
      matrix(
        c(both[i], target_pos_gene_neg[i], gene_pos_target_neg[i], neither[i]),
        nrow = 2
      ),
      alternative = "greater"
    )$p.value)
  }, numeric(1))

  out <- tibble(
    dataset = dataset,
    group_name = group_name,
    gene = rownames(mat)[keep],
    n_units = ncol(mat),
    apobec3a_positive_units = n_target_pos,
    apobec3a_detection_rate = n_target_pos / ncol(mat),
    gene_detection_rate = as.numeric(gene_pct[keep]),
    both_detected = as.numeric(both),
    jaccard = as.numeric(jaccard),
    odds_ratio = as.numeric((both + 0.5) * (neither + 0.5) /
      ((target_pos_gene_neg + 0.5) * (gene_pos_target_neg + 0.5))),
    mean_expr_in_APOBEC3A_pos = as.numeric(mean_pos),
    mean_expr_in_APOBEC3A_neg = as.numeric(mean_neg),
    log2_fc_APOBEC3A_pos_vs_neg = log2((as.numeric(mean_pos) + 1e-6) / (as.numeric(mean_neg) + 1e-6)),
    fisher_p = fisher_p,
    p_adj = p.adjust(fisher_p, method = "BH")
  ) %>%
    arrange(p_adj, desc(odds_ratio), desc(log2_fc_APOBEC3A_pos_vs_neg))

  if ("disease" %in% colnames(metadata)) {
    det_by_severity <- metadata %>%
      mutate(APOBEC3A_detected = target_pos) %>%
      group_by(disease) %>%
      summarise(
        n_units = n(),
        apobec3a_positive_units = sum(APOBEC3A_detected),
        apobec3a_detection_rate = mean(APOBEC3A_detected),
        .groups = "drop"
      )
    write.csv(
      det_by_severity,
      file.path(tab_dir, paste0("apobec3a_detection_by_severity_", dataset, "_", group_name, ".csv")),
      row.names = FALSE
    )
  }
  out
}

message("Loading Visium")
vis <- readRDS(file.path(spatial_root, "GSE314851_Visium_FFPE_integrated.rds"))
vis_mat <- GetAssayData(vis, assay = "Spatial", layer = "data")
vis_meta <- vis@meta.data %>% mutate(disease = category)
vis_res <- co_detection(vis_mat, vis_meta, "Visium_spot", "all_spots", min_gene_pct = 0.002)
write.csv(vis_res, file.path(tab_dir, "apobec3a_visium_spot_coexpression.csv"), row.names = FALSE)

message("Loading Xenium")
xen <- readRDS(file.path(spatial_root, "GSE315246_xenium.obj.integrated.rds"))
xen_mat <- GetAssayData(xen, assay = "Xenium", layer = "counts")
xen_meta <- xen@meta.data
xen_res <- co_detection(xen_mat, xen_meta, "Xenium_cell", "all_cells", min_gene_pct = 0.002)
write.csv(xen_res, file.path(tab_dir, "apobec3a_xenium_cell_coexpression.csv"), row.names = FALSE)

myeloid_idx <- which(xen_meta$predicted.id == "Myeloid")
myeloid_res <- co_detection(
  xen_mat[, myeloid_idx, drop = FALSE],
  xen_meta[myeloid_idx, , drop = FALSE],
  "Xenium_cell",
  "myeloid_cells",
  min_gene_pct = 0.002
)
write.csv(myeloid_res, file.path(tab_dir, "apobec3a_xenium_myeloid_coexpression.csv"), row.names = FALSE)

message("Top Visium co-detected genes")
print(head(vis_res, 20))
message("Top Xenium all-cell co-detected genes")
print(head(xen_res, 20))
message("Top Xenium myeloid co-detected genes")
print(head(myeloid_res, 20))
