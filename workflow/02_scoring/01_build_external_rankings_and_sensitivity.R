#!/usr/bin/env Rscript

# Compare external MMI signatures built from the top N genes in the complete
# GSE5099 and GSE11864 macrophage-versus-monocyte differential rankings.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(Seurat)
  library(AUCell)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

repo_root <- Sys.getenv(
  "PIPELINE_REPO_ROOT",
  unset = normalizePath(".", mustWork = TRUE)
)
workspace_dir <- Sys.getenv(
  "SCRNA_WORKSPACE_ROOT",
  unset = "/public3/DSC/single_cell"
)
input_rds <- Sys.getenv(
  "SCRNA_SCORED_RDS",
  unset = file.path(workspace_dir, "Result", "figer_new", "hdWGCNA", "Mo_Ma", "data_scored.rds")
)
out_dir <- Sys.getenv(
  "SCORING_OUTPUT_DIR",
  unset = file.path(repo_root, "results", "scoring")
)
geo_cache <- Sys.getenv(
  "SCORING_GEO_CACHE",
  unset = file.path(workspace_dir, "GEOquery_cache")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_cache, recursive = TRUE, showWarnings = FALSE)

top_n_values <- c(100L, 200L, 250L, 300L)
sources <- c("GSE5099", "GSE11864")

data <- readRDS(input_rds)
rna_features <- rownames(data[["RNA"]])
auc_max_rank <- ceiling(0.05 * length(rna_features))

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[grepl("///|//|;", x)] <- NA_character_
  x[x == ""] <- NA_character_
  toupper(x)
}

collapse_expression <- function(expression_matrix, probe_ids, gene_symbols) {
  gene_symbols <- clean_symbols(gene_symbols)
  matched_rows <- match(probe_ids, rownames(expression_matrix))
  keep <- !is.na(matched_rows) & !is.na(gene_symbols)
  expression_matrix <- expression_matrix[matched_rows[keep], , drop = FALSE]
  gene_symbols <- gene_symbols[keep]
  summed <- rowsum(expression_matrix, group = gene_symbols, reorder = FALSE)
  counts <- as.numeric(table(factor(gene_symbols, levels = rownames(summed))))
  summed / counts
}

fit_paired_effect <- function(expression_matrix, condition, donor) {
  condition <- factor(condition, levels = c("Monocyte", "Macrophage"))
  donor <- factor(donor)
  design <- model.matrix(~ donor + condition)
  fit <- eBayes(lmFit(expression_matrix, design), trend = TRUE)
  coef_name <- "conditionMacrophage"
  tibble(
    gene_symbol = rownames(expression_matrix),
    logFC = fit$coefficients[, coef_name],
    moderated_t = fit$t[, coef_name],
    p_value = fit$p.value[, coef_name],
    p_adjusted_BH = p.adjust(fit$p.value[, coef_name], method = "BH")
  )
}

message("Reconstructing the complete GSE5099 differential ranking...")
gse5099 <- getGEO(
  "GSE5099", GSEMatrix = TRUE, getGPL = FALSE, destdir = geo_cache
)
gse5099_matrices <- lapply(gse5099, function(eset) {
  pd <- pData(eset)
  keep <- grepl(
    "Monocyte at T0|Macrophage at 7 days", pd$title, ignore.case = TRUE
  )
  eset <- eset[, keep]
  pd <- pData(eset)
  condition <- ifelse(
    grepl("Macrophage at 7 days", pd$title, ignore.case = TRUE),
    "Macrophage", "Monocyte"
  )
  donor <- sub(".*rep([0-9]+).*", "\\1", pd$title, ignore.case = TRUE)
  logical_name <- paste(condition, donor, sep = "_rep")

  gpl_id <- annotation(eset)
  gpl <- getGEO(gpl_id, destdir = geo_cache, AnnotGPL = TRUE)
  annot <- Table(gpl)
  symbol_column <- intersect(
    c("Gene symbol", "Gene Symbol", "GENE_SYMBOL"), colnames(annot)
  )[1]
  if (is.na(symbol_column)) {
    stop("No gene-symbol annotation column found for ", gpl_id, ".")
  }
  collapsed <- collapse_expression(exprs(eset), annot$ID, annot[[symbol_column]])
  colnames(collapsed) <- logical_name
  collapsed[, order(colnames(collapsed)), drop = FALSE]
})

gse5099_expression <- do.call(rbind, gse5099_matrices)
gse5099_expression <- rowsum(
  gse5099_expression,
  group = rownames(gse5099_expression),
  reorder = FALSE
) / as.numeric(table(factor(
  rownames(gse5099_expression),
  levels = unique(rownames(gse5099_expression))
)))
gse5099_columns <- colnames(gse5099_expression)
gse5099_effects <- fit_paired_effect(
  gse5099_expression,
  condition = sub("_rep.*", "", gse5099_columns),
  donor = sub(".*_rep", "", gse5099_columns)
) %>%
  filter(
    gene_symbol %in% rna_features,
    is.finite(moderated_t),
    is.finite(logFC)
  ) %>%
  arrange(desc(moderated_t)) %>%
  mutate(
    rank_macrophage = rank(-moderated_t, ties.method = "first"),
    rank_monocyte = rank(moderated_t, ties.method = "first")
  )

message("Reconstructing the complete paired GSE11864 differential ranking...")
gse11864 <- getGEO(
  "GSE11864", GSEMatrix = TRUE, getGPL = FALSE, destdir = geo_cache
)[[1]]
gse11864_pd <- pData(gse11864)
gse11864_keep <- grepl(
  "fresh monocytes$|, M-CSF$", gse11864_pd$title, ignore.case = TRUE
)
gse11864 <- gse11864[, gse11864_keep]
gse11864_pd <- pData(gse11864)
gse11864_expression <- exprs(gse11864)
if (max(gse11864_expression, na.rm = TRUE) > 100) {
  gse11864_expression <- log2(gse11864_expression)
}
gse11864_expression <- normalizeBetweenArrays(
  gse11864_expression, method = "quantile"
)
probe_symbols <- AnnotationDbi::mapIds(
  hgu133plus2.db,
  keys = rownames(gse11864_expression),
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)
gse11864_expression <- collapse_expression(
  gse11864_expression,
  probe_ids = names(probe_symbols),
  gene_symbols = unname(probe_symbols)
)
gse11864_condition <- ifelse(
  grepl("fresh monocytes$", gse11864_pd$title, ignore.case = TRUE),
  "Monocyte", "Macrophage"
)
gse11864_donor <- sub("Donor ([0-9]+).*", "\\1", gse11864_pd$title)
gse11864_effects <- fit_paired_effect(
  gse11864_expression,
  condition = gse11864_condition,
  donor = gse11864_donor
) %>%
  filter(
    gene_symbol %in% rna_features,
    is.finite(moderated_t),
    is.finite(logFC)
  ) %>%
  arrange(desc(moderated_t)) %>%
  mutate(
    rank_macrophage = rank(-moderated_t, ties.method = "first"),
    rank_monocyte = rank(moderated_t, ties.method = "first")
  )

write.csv(
  gse5099_effects,
  file.path(out_dir, "GSE5099_complete_external_effect_ranking.csv"),
  row.names = FALSE
)
write.csv(
  gse11864_effects,
  file.path(out_dir, "GSE11864_complete_external_effect_ranking.csv"),
  row.names = FALSE
)

effect_tables <- list(
  GSE5099 = gse5099_effects,
  GSE11864 = gse11864_effects
)

gene_sets <- list()
gene_set_membership <- list()
conflict_records <- list()
set_size_records <- list()
final_sets <- list()

for (top_n in top_n_values) {
  source_top <- lapply(sources, function(source) {
    effects <- effect_tables[[source]]
    list(
      macrophage = effects %>%
        arrange(desc(moderated_t)) %>%
        slice_head(n = top_n) %>%
        pull(gene_symbol),
      monocyte = effects %>%
        arrange(moderated_t) %>%
        slice_head(n = top_n) %>%
        pull(gene_symbol)
    )
  })
  names(source_top) <- sources

  macrophage_union <- unique(c(
    source_top$GSE5099$macrophage,
    source_top$GSE11864$macrophage
  ))
  monocyte_union <- unique(c(
    source_top$GSE5099$monocyte,
    source_top$GSE11864$monocyte
  ))
  direction_conflicts <- sort(intersect(macrophage_union, monocyte_union))
  macrophage_final <- sort(setdiff(macrophage_union, direction_conflicts))
  monocyte_final <- sort(setdiff(monocyte_union, direction_conflicts))

  set_label <- paste0("Top", top_n)
  final_sets[[set_label]] <- list(
    mature = macrophage_final,
    immature = monocyte_final
  )
  gene_sets[[paste0(set_label, "_Mature")]] <- macrophage_final
  gene_sets[[paste0(set_label, "_Immature")]] <- monocyte_final

  set_size_records[[set_label]] <- tibble(
    top_n_per_source_per_direction = top_n,
    GSE5099_macrophage = length(source_top$GSE5099$macrophage),
    GSE11864_macrophage = length(source_top$GSE11864$macrophage),
    macrophage_source_overlap = length(intersect(
      source_top$GSE5099$macrophage,
      source_top$GSE11864$macrophage
    )),
    macrophage_union_before_conflict_removal = length(macrophage_union),
    final_macrophage_genes = length(macrophage_final),
    GSE5099_monocyte = length(source_top$GSE5099$monocyte),
    GSE11864_monocyte = length(source_top$GSE11864$monocyte),
    monocyte_source_overlap = length(intersect(
      source_top$GSE5099$monocyte,
      source_top$GSE11864$monocyte
    )),
    monocyte_union_before_conflict_removal = length(monocyte_union),
    final_monocyte_genes = length(monocyte_final),
    direction_conflicts_removed = length(direction_conflicts)
  )

  if (length(direction_conflicts)) {
    conflict_records[[set_label]] <- tibble(
      top_n_per_source_per_direction = top_n,
      gene_symbol = direction_conflicts,
      GSE5099_macrophage = direction_conflicts %in%
        source_top$GSE5099$macrophage,
      GSE5099_monocyte = direction_conflicts %in%
        source_top$GSE5099$monocyte,
      GSE11864_macrophage = direction_conflicts %in%
        source_top$GSE11864$macrophage,
      GSE11864_monocyte = direction_conflicts %in%
        source_top$GSE11864$monocyte,
      action = "removed_from_both_directions"
    )
  }

  for (direction_label in c("mature_positive", "immature_negative")) {
    genes <- if (direction_label == "mature_positive") {
      macrophage_final
    } else {
      monocyte_final
    }
    source_direction <- if (direction_label == "mature_positive") {
      "macrophage"
    } else {
      "monocyte"
    }
    gene_set_membership[[paste(set_label, direction_label, sep = "_")]] <- tibble(
      top_n_per_source_per_direction = top_n,
      direction = direction_label,
      gene_symbol = genes,
      in_GSE5099_topN = genes %in%
        source_top$GSE5099[[source_direction]],
      in_GSE11864_topN = genes %in%
        source_top$GSE11864[[source_direction]],
      source_support_count =
        as.integer(in_GSE5099_topN) + as.integer(in_GSE11864_topN),
      GSE5099_rank_in_direction = if (direction_label == "mature_positive") {
        gse5099_effects$rank_macrophage[
          match(genes, gse5099_effects$gene_symbol)
        ]
      } else {
        gse5099_effects$rank_monocyte[
          match(genes, gse5099_effects$gene_symbol)
        ]
      },
      GSE11864_rank_in_direction = if (direction_label == "mature_positive") {
        gse11864_effects$rank_macrophage[
          match(genes, gse11864_effects$gene_symbol)
        ]
      } else {
        gse11864_effects$rank_monocyte[
          match(genes, gse11864_effects$gene_symbol)
        ]
      }
    )
  }
}

gene_set_membership <- bind_rows(gene_set_membership)
conflict_records <- bind_rows(conflict_records)
set_size_summary <- bind_rows(set_size_records)

write.csv(
  gene_set_membership,
  file.path(out_dir, "topN_deduplicated_gene_set_membership.csv"),
  row.names = FALSE
)
write.csv(
  conflict_records,
  file.path(out_dir, "topN_direction_conflicts_removed.csv"),
  row.names = FALSE
)
write.csv(
  set_size_summary,
  file.path(out_dir, "topN_gene_set_size_and_overlap_summary.csv"),
  row.names = FALSE
)

incremental_records <- list()
previous_sets <- NULL
for (top_n in top_n_values) {
  current_sets <- final_sets[[paste0("Top", top_n)]]
  for (direction_label in c("mature", "immature")) {
    current <- current_sets[[direction_label]]
    if (is.null(previous_sets)) {
      incremental_records[[paste(top_n, direction_label, "baseline", sep = "_")]] <-
        tibble(
          top_n_per_source_per_direction = top_n,
          direction = direction_label,
          comparison_to_previous_top_n = NA_integer_,
          change = "baseline_top100",
          gene_symbol = current
        )
    } else {
      previous_top_n <- top_n_values[match(top_n, top_n_values) - 1L]
      added <- setdiff(current, previous_sets[[direction_label]])
      removed <- setdiff(previous_sets[[direction_label]], current)
      retained <- intersect(current, previous_sets[[direction_label]])
      incremental_records[[paste(top_n, direction_label, "added", sep = "_")]] <-
        tibble(
          top_n_per_source_per_direction = top_n,
          direction = direction_label,
          comparison_to_previous_top_n = previous_top_n,
          change = "added",
          gene_symbol = added
        )
      incremental_records[[paste(top_n, direction_label, "removed", sep = "_")]] <-
        tibble(
          top_n_per_source_per_direction = top_n,
          direction = direction_label,
          comparison_to_previous_top_n = previous_top_n,
          change = "removed",
          gene_symbol = removed
        )
      incremental_records[[paste(top_n, direction_label, "retained", sep = "_")]] <-
        tibble(
          top_n_per_source_per_direction = top_n,
          direction = direction_label,
          comparison_to_previous_top_n = previous_top_n,
          change = "retained",
          gene_symbol = retained
        )
    }
  }
  previous_sets <- current_sets
}
incremental_records <- bind_rows(incremental_records)
write.csv(
  incremental_records,
  file.path(out_dir, "topN_incremental_gene_changes.csv"),
  row.names = FALSE
)

message("Building one shared RNA/counts ranking and calculating eight AUC sets...")
expression_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
rankings <- AUCell_buildRankings(
  expression_matrix,
  plotStats = FALSE,
  splitByBlocks = TRUE,
  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
  verbose = TRUE
)
auc <- AUCell_calcAUC(
  gene_sets,
  rankings,
  normAUC = TRUE,
  aucMaxRank = auc_max_rank,
  nCores = 1,
  verbose = TRUE
)
auc_scores <- as.data.frame(t(getAUC(auc)), check.names = FALSE)
auc_scores <- auc_scores[colnames(data), , drop = FALSE]

score_meta <- tibble(cell_barcode = colnames(data))
for (top_n in top_n_values) {
  set_label <- paste0("Top", top_n)
  mature_auc <- auc_scores[[paste0(set_label, "_Mature")]]
  immature_auc <- auc_scores[[paste0(set_label, "_Immature")]]
  score_meta[[paste0(set_label, "_Mature_AUCell")]] <- mature_auc
  score_meta[[paste0(set_label, "_Immature_AUCell")]] <- immature_auc
  score_meta[[paste0(set_label, "_MMI")]] <- mature_auc - immature_auc
}
saveRDS(score_meta, file.path(out_dir, "topN_MMI_AUCell_cell_scores.rds"))
write.csv(
  score_meta,
  gzfile(file.path(out_dir, "topN_MMI_AUCell_cell_scores.csv.gz")),
  row.names = FALSE
)

score_values <- as.data.frame(
  score_meta[, setdiff(colnames(score_meta), "cell_barcode"), drop = FALSE],
  check.names = FALSE
)
rownames(score_values) <- score_meta$cell_barcode
data <- AddMetaData(data, score_values)
rm(expression_matrix, rankings, auc, auc_scores)
invisible(gc())

subtype_summary <- bind_rows(lapply(top_n_values, function(top_n) {
  set_label <- paste0("Top", top_n)
  mature_column <- paste0(set_label, "_Mature_AUCell")
  immature_column <- paste0(set_label, "_Immature_AUCell")
  mmi_column <- paste0(set_label, "_MMI")
  data[[]] %>%
    group_by(Celltype_raw1) %>%
    summarise(
      top_n_per_source_per_direction = top_n,
      final_mature_genes = length(final_sets[[set_label]]$mature),
      final_immature_genes = length(final_sets[[set_label]]$immature),
      n_cells = n(),
      Mature_AUC_mean = mean(.data[[mature_column]], na.rm = TRUE),
      Immature_AUC_mean = mean(.data[[immature_column]], na.rm = TRUE),
      MMI_mean = mean(.data[[mmi_column]], na.rm = TRUE),
      MMI_median = median(.data[[mmi_column]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      MMI_mean_rank_high_to_low = rank(-MMI_mean, ties.method = "min"),
      MMI_median_rank_high_to_low = rank(-MMI_median, ties.method = "min")
    )
}))
write.csv(
  subtype_summary,
  file.path(out_dir, "topN_MMI_subtype_summary_and_ranks.csv"),
  row.names = FALSE
)

fc_order <- subtype_summary %>%
  filter(Celltype_raw1 %in% c("Foam cells1", "Foam cells2", "LAM")) %>%
  select(
    top_n_per_source_per_direction, final_mature_genes,
    final_immature_genes, Celltype_raw1, MMI_mean, MMI_median,
    MMI_mean_rank_high_to_low, MMI_median_rank_high_to_low
  )
write.csv(
  fc_order,
  file.path(out_dir, "topN_FC1_FC2_LAM_order.csv"),
  row.names = FALSE
)

rank_changes <- subtype_summary %>%
  select(
    top_n_per_source_per_direction, Celltype_raw1,
    MMI_mean, MMI_median,
    MMI_mean_rank_high_to_low, MMI_median_rank_high_to_low
  ) %>%
  arrange(Celltype_raw1, top_n_per_source_per_direction)
write.csv(
  rank_changes,
  file.path(out_dir, "topN_subtype_rank_changes.csv"),
  row.names = FALSE
)

sample_scores <- bind_rows(lapply(top_n_values, function(top_n) {
  mmi_column <- paste0("Top", top_n, "_MMI")
  data[[]] %>%
    group_by(Source_GSE, Patient_ID, Sample_Type, Celltype_raw1) %>%
    summarise(
      top_n_per_source_per_direction = top_n,
      n_cells = n(),
      MMI_median = median(.data[[mmi_column]], na.rm = TRUE),
      .groups = "drop"
    )
}))
write.csv(
  sample_scores,
  file.path(out_dir, "topN_sample_subtype_scores.csv"),
  row.names = FALSE
)

paired_fc2_lam <- sample_scores %>%
  filter(Celltype_raw1 %in% c("Foam cells2", "LAM")) %>%
  select(
    top_n_per_source_per_direction, Source_GSE, Patient_ID,
    Sample_Type, Celltype_raw1, MMI_median
  ) %>%
  pivot_wider(names_from = Celltype_raw1, values_from = MMI_median) %>%
  filter(!is.na(`Foam cells2`), !is.na(LAM)) %>%
  mutate(LAM_minus_FC2 = LAM - `Foam cells2`) %>%
  group_by(top_n_per_source_per_direction) %>%
  summarise(
    n_pairs = n(),
    median_LAM_minus_FC2 = median(LAM_minus_FC2),
    mean_LAM_minus_FC2 = mean(LAM_minus_FC2),
    n_LAM_higher = sum(LAM_minus_FC2 > 0),
    n_FC2_higher = sum(LAM_minus_FC2 < 0),
    p_paired = wilcox.test(
      LAM, `Foam cells2`, paired = TRUE, exact = FALSE
    )$p.value,
    .groups = "drop"
  )
write.csv(
  paired_fc2_lam,
  file.path(out_dir, "topN_FC2_LAM_paired_sample_tests.csv"),
  row.names = FALSE
)

save_plot <- function(filename, plot, width, height) {
  ggsave(
    file.path(out_dir, paste0(filename, ".pdf")),
    plot = plot, device = cairo_pdf, width = width, height = height
  )
  ggsave(
    file.path(out_dir, paste0(filename, ".png")),
    plot = plot, width = width, height = height, dpi = 300, bg = "white"
  )
}

plot_top_n_mmi <- function(top_n, shared_limits = NULL) {
  plot_data <- subtype_summary %>%
    filter(top_n_per_source_per_direction == .env$top_n) %>%
    arrange(MMI_mean) %>%
    mutate(Celltype_raw1 = factor(Celltype_raw1, levels = Celltype_raw1))
  mature_n <- unique(plot_data$final_mature_genes)
  immature_n <- unique(plot_data$final_immature_genes)

  plot <- ggplot(
    plot_data,
    aes(x = Celltype_raw1, y = MMI_mean, fill = MMI_mean)
  ) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
    theme_classic() +
    labs(
      x = "Cell Subtypes",
      y = "Macrophage Maturation Index (MMI)",
      title = paste0(
        "Top ", top_n, "/source/direction",
        " (deduplicated: ", mature_n, " mature, ", immature_n, " monocyte)"
      )
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, size = 12, face = "bold"
      ),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  if (!is.null(shared_limits)) {
    plot <- plot + coord_cartesian(ylim = shared_limits)
  }
  plot
}

all_limits <- range(subtype_summary$MMI_mean, na.rm = TRUE)
padding <- 0.06 * diff(all_limits)
shared_limits <- all_limits + c(-padding, padding)

plots <- lapply(top_n_values, function(top_n) {
  plot <- plot_top_n_mmi(top_n)
  save_plot(
    paste0("02_Top", top_n, "_Macrophage_Maturation_Index_MMI"),
    plot, 12, 6
  )
  plot_top_n_mmi(top_n, shared_limits = shared_limits)
})

combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
  plot_annotation(
    title = "GSE5099 + GSE11864 top-N external MMI comparison"
  ) &
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
save_plot(
  "02_TopN_Macrophage_Maturation_Index_MMI_combined",
  combined, 22, 14
)

subtype_palette <- c(
  "Classical Mono" = "#D55E00",
  "Non-classical Mono" = "#E69F00",
  "Inflammatory Mono" = "#CC79A7",
  "ISG+ Mono" = "#F0E442",
  "Transitional Mac" = "#999999",
  "CX3CR1+ TRM" = "#56B4E9",
  "LYVE1+ TRM" = "#009E73",
  "Foam cells1" = "#0072B2",
  "Foam cells2" = "#6A51A3",
  "LAM" = "#000000"
)

p_rank <- ggplot(
  rank_changes,
  aes(
    x = top_n_per_source_per_direction,
    y = MMI_mean_rank_high_to_low,
    color = Celltype_raw1,
    group = Celltype_raw1
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = subtype_palette) +
  scale_x_continuous(breaks = top_n_values) +
  scale_y_reverse(breaks = seq_len(10)) +
  theme_classic(base_size = 12) +
  labs(
    x = "Top genes per source and direction",
    y = "MMI mean rank (1 = highest)",
    color = "Cell subtype"
  ) +
  theme(legend.position = "right")
save_plot("TopN_subtype_MMI_rank_trajectory", p_rank, 10, 7)

p_mean <- ggplot(
  rank_changes,
  aes(
    x = top_n_per_source_per_direction,
    y = MMI_mean,
    color = Celltype_raw1,
    group = Celltype_raw1
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = subtype_palette) +
  scale_x_continuous(breaks = top_n_values) +
  theme_classic(base_size = 12) +
  labs(
    x = "Top genes per source and direction",
    y = "Mean MMI",
    color = "Cell subtype"
  ) +
  theme(legend.position = "right")
save_plot("TopN_subtype_MMI_mean_trajectory", p_mean, 10, 7)

write.csv(
  tibble(
    parameter = c(
      "external_sources", "top_n_per_source_per_direction",
      "deduplication", "direction_conflict_handling",
      "input_assay", "input_layer", "aucMaxRank", "normAUC",
      "selection_uses_carotid_study_data"
    ),
    value = c(
      "GSE5099;GSE11864", paste(top_n_values, collapse = ";"),
      "union_within_direction",
      "remove_conflicting_gene_from_both_directions",
      "RNA", "counts", auc_max_rank, "TRUE", "false"
    )
  ),
  file.path(out_dir, "topN_MMI_parameter_manifest.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Finished. Outputs written to: ", out_dir)
