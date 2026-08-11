#!/usr/bin/env Rscript

# Final publication-oriented MPI/MMI analysis.
# MMI uses exactly 200 source-balanced consensus genes per direction from
# GSE5099 and GSE11864. Plot geometry follows single_cell.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scales)
})

set.seed(20260730)

top_n <- 200L
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
comparison_dir <- Sys.getenv(
  "SCORING_OUTPUT_DIR",
  unset = file.path(repo_root, "results", "scoring")
)
gene_membership_file <- file.path(
  repo_root, "data", "gene_sets",
  "MMI_GSE5099_GSE11864_Top200_membership.csv"
)
conflict_file <- file.path(
  repo_root, "data", "gene_sets",
  "MMI_GSE5099_GSE11864_Top200_direction_conflicts.csv"
)
comparison_score_rds <- file.path(
  comparison_dir, "topN_MMI_AUCell_cell_scores.rds"
)
out_dir <- file.path(comparison_dir, "final_mpi_mmi")
score_rds <- file.path(out_dir, "ConsensusTop200_MMI_AUCell_cell_scores.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mature_column <- "Top200_Mature_AUCell"
immature_column <- "Top200_Monocyte_AUCell"
mmi_column <- "Top200_MMI"

# Display controls only; they do not alter scores or statistical tests.
mmi_color_quantiles <- c(0.02, 0.98)
density_trim_quantiles <- c(0.02, 0.98)
density_contour_breaks <- seq(0.2, 0.9, by = 0.1)

save_plot <- function(filename, plot, width, height) {
  ggsave(
    file.path(out_dir, paste0(filename, ".pdf")), plot = plot,
    device = cairo_pdf, width = width, height = height
  )
  ggsave(
    file.path(out_dir, paste0(filename, ".png")), plot = plot,
    width = width, height = height, dpi = 300, bg = "white"
  )
}

message("Preparing the fixed GSE5099 + GSE11864 Top200 gene sets...")
if (!file.exists(gene_membership_file)) {
  stop("Top-N gene membership file not found: ", gene_membership_file)
}
membership_all <- read.csv(gene_membership_file, stringsAsFactors = FALSE)
required_membership_columns <- c(
  "top_n_final_per_direction", "direction", "consensus_rank", "gene_symbol",
  "consensus_mean_percentile", "consensus_worst_source_percentile",
  "GSE5099_rank_in_direction", "GSE11864_rank_in_direction",
  "selected_in_GSE5099_individual_top200",
  "selected_in_GSE11864_individual_top200"
)
missing_membership_columns <- setdiff(
  required_membership_columns, colnames(membership_all)
)
if (length(missing_membership_columns)) {
  stop(
    "Gene membership file is missing columns: ",
    paste(missing_membership_columns, collapse = ", ")
  )
}

membership <- membership_all %>%
  filter(top_n_final_per_direction == top_n)
if (!nrow(membership)) {
  stop("No Top", top_n, " records found in the gene membership file.")
}
mature_genes <- sort(unique(
  membership$gene_symbol[membership$direction == "mature_positive"]
))
immature_genes <- sort(unique(
  membership$gene_symbol[membership$direction == "immature_negative"]
))
conflicting_genes <- if (file.exists(conflict_file)) {
  read.csv(conflict_file, stringsAsFactors = FALSE) %>%
    pull(gene_symbol) %>%
    unique() %>%
    sort()
} else {
  character()
}

if (length(intersect(mature_genes, immature_genes))) {
  stop("Positive and negative gene sets still overlap after conflict removal.")
}
if (length(mature_genes) != top_n || length(immature_genes) != top_n) {
  stop("The frozen MMI signatures must contain exactly Top", top_n, " genes per direction.")
}

write.csv(
  membership,
  file.path(out_dir, "Top200_MMI_gene_set_membership.csv"),
  row.names = FALSE
)
write.csv(
  tibble(gene_symbol = conflicting_genes, action = "removed_from_both_directions"),
  file.path(out_dir, "Top200_MMI_removed_direction_conflicts.csv"),
  row.names = FALSE
)

message(
  "Fixed Top200 gene sets: ", length(mature_genes), " mature-positive and ",
  length(immature_genes), " monocyte-high genes; ",
  length(conflicting_genes),
  " direction-discordant genes were excluded before consensus ranking."
)

message("Loading the original scored monocyte/macrophage object...")
data <- readRDS(input_rds)
if (!"RNA" %in% Assays(data) || !"counts" %in% Layers(data[["RNA"]])) {
  stop("The RNA assay must contain a counts layer for the AUCell workflow.")
}

mpi_columns <- c("Polarization_Index", "M1_Score", "M2_Score")
if (!all(mpi_columns %in% colnames(data[[]]))) {
  message("MPI columns are absent; calculating them from the frozen 145/165 signatures...")
  mpi_file <- file.path(repo_root, "data", "gene_sets", "mpi_signatures.csv")
  mpi_table <- read.csv(mpi_file, stringsAsFactors = FALSE)
  mpi_sets <- split(mpi_table$gene, mpi_table$gene_set)
  m1_genes <- intersect(mpi_sets[["MPI M1-like signature"]], rownames(data[["RNA"]]))
  m2_genes <- intersect(mpi_sets[["MPI M2-like signature"]], rownames(data[["RNA"]]))
  expression_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
  mpi_rankings <- AUCell_buildRankings(
    expression_matrix,
    plotStats = FALSE,
    splitByBlocks = TRUE,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    verbose = TRUE
  )
  mpi_auc <- AUCell_calcAUC(
    list(M1 = m1_genes, M2 = m2_genes),
    mpi_rankings,
    normAUC = TRUE,
    aucMaxRank = ceiling(0.05 * nrow(expression_matrix)),
    nCores = 1,
    verbose = TRUE
  )
  mpi_scores <- as.data.frame(t(getAUC(mpi_auc)), check.names = FALSE)
  mpi_scores <- mpi_scores[colnames(data), , drop = FALSE]
  data$M1_Score <- mpi_scores$M1
  data$M2_Score <- mpi_scores$M2
  data$Polarization_Index <- data$M1_Score - data$M2_Score
  rm(expression_matrix, mpi_rankings, mpi_auc, mpi_scores)
  invisible(gc())
}
required_meta <- c(
  "Celltype_raw1", "Sample_Type", "Patient_ID", "Source_GSE",
  "Polarization_Index", "M1_Score", "M2_Score"
)
missing_meta <- setdiff(required_meta, colnames(data[[]]))
if (length(missing_meta)) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

rna_features <- rownames(data[["RNA"]])
mature_genes <- intersect(mature_genes, rna_features)
immature_genes <- intersect(immature_genes, rna_features)
auc_max_rank <- ceiling(0.05 * length(rna_features))
if (length(intersect(mature_genes, immature_genes))) {
  stop("Detected positive and negative gene sets overlap unexpectedly.")
}

score_columns <- c(mature_column, immature_column, mmi_column)
if (file.exists(score_rds)) {
  message("Using the self-contained final Top200 AUCell score cache...")
  score_meta <- readRDS(score_rds)
} else if (file.exists(comparison_score_rds)) {
  message("Importing the verified Top200 scores from the completed Top-N comparison...")
  comparison_scores <- readRDS(comparison_score_rds)
  required_comparison_scores <- c(
    "cell_barcode", "Top200_Mature_AUCell",
    "Top200_Immature_AUCell", "Top200_MMI"
  )
  missing_comparison_scores <- setdiff(
    required_comparison_scores, colnames(comparison_scores)
  )
  if (length(missing_comparison_scores)) {
    stop(
      "Top-N score cache is missing columns: ",
      paste(missing_comparison_scores, collapse = ", ")
    )
  }
  score_meta <- comparison_scores %>%
    transmute(
      cell_barcode,
      !!mature_column := Top200_Mature_AUCell,
      !!immature_column := Top200_Immature_AUCell,
      !!mmi_column := Top200_MMI
    )
  rm(comparison_scores)
} else {
  message("Building per-cell gene rankings from RNA/counts...")
  expression_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
  cells_rankings <- AUCell_buildRankings(
    expression_matrix,
    plotStats = FALSE,
    splitByBlocks = TRUE,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    verbose = TRUE
  )

  gene_sets <- list(
    Top200_Mature = mature_genes,
    Top200_Monocyte = immature_genes
  )
  message("Calculating AUCell AUC values with aucMaxRank = ", auc_max_rank, "...")
  cells_auc <- AUCell_calcAUC(
    gene_sets,
    cells_rankings,
    normAUC = TRUE,
    aucMaxRank = auc_max_rank,
    nCores = 1,
    verbose = TRUE
  )
  auc_scores <- as.data.frame(t(getAUC(cells_auc)), check.names = FALSE)
  auc_scores <- auc_scores[colnames(data), , drop = FALSE]
  colnames(auc_scores) <- c(mature_column, immature_column)
  auc_scores[[mmi_column]] <-
    auc_scores[[mature_column]] - auc_scores[[immature_column]]
  score_meta <- auc_scores %>%
    rownames_to_column("cell_barcode")
  rm(expression_matrix, cells_rankings, cells_auc, auc_scores)
  invisible(gc())
}

score_index <- match(colnames(data), score_meta$cell_barcode)
if (
  anyNA(score_index) ||
    !all(score_columns %in% colnames(score_meta)) ||
    anyNA(score_meta[, score_columns, drop = FALSE])
) {
  stop("Top200 AUCell scores do not match the current Seurat object.")
}
new_scores <- as.data.frame(
  score_meta[score_index, score_columns, drop = FALSE],
  check.names = FALSE
)
rownames(new_scores) <- colnames(data)
data <- AddMetaData(data, new_scores)
saveRDS(score_meta, score_rds)
write.csv(
  score_meta,
  gzfile(file.path(out_dir, "ConsensusTop200_MMI_AUCell_cell_scores.csv.gz")),
  row.names = FALSE
)

data$Sample_Type <- factor(
  as.character(data$Sample_Type),
  levels = c("Atherosclerotic Core", "Proximal Adjacent")
)

format_p <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 0.001) return(sprintf("p = %.2e", p))
  sprintf("p = %.3f", p)
}

paired_or_unpaired_wilcox <- function(df, score_column) {
  paired <- df %>%
    select(Source_GSE, Patient_ID, Sample_Type, all_of(score_column)) %>%
    pivot_wider(names_from = Sample_Type, values_from = all_of(score_column)) %>%
    filter(!is.na(`Atherosclerotic Core`), !is.na(`Proximal Adjacent`))

  if (nrow(paired) >= 3L) {
    test_result <- wilcox.test(
      paired$`Atherosclerotic Core`, paired$`Proximal Adjacent`,
      paired = TRUE, exact = FALSE
    )
    return(tibble(
      score = score_column,
      test = "paired Wilcoxon signed-rank test",
      analysis_unit = "paired patient",
      n_core = nrow(paired),
      n_adjacent = nrow(paired),
      n_pairs = nrow(paired),
      p_value = test_result$p.value
    ))
  }

  complete <- df %>%
    filter(!is.na(Sample_Type), !is.na(.data[[score_column]]))
  test_result <- wilcox.test(
    complete[[score_column]] ~ complete$Sample_Type,
    paired = FALSE, exact = FALSE
  )
  tibble(
    score = score_column,
    test = "Wilcoxon rank-sum test",
    analysis_unit = "biological sample",
    n_core = sum(complete$Sample_Type == "Atherosclerotic Core"),
    n_adjacent = sum(complete$Sample_Type == "Proximal Adjacent"),
    n_pairs = nrow(paired),
    p_value = test_result$p.value
  )
}

add_single_bracket <- function(plot, values, label) {
  value_range <- diff(range(values, na.rm = TRUE))
  if (!is.finite(value_range) || value_range == 0) value_range <- 0.01
  y_max <- max(values, na.rm = TRUE)
  y_tick <- y_max + 0.08 * value_range
  y_bar <- y_max + 0.12 * value_range
  y_label <- y_max + 0.18 * value_range
  plot +
    annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, linewidth = 0.6) +
    annotate("segment", x = 1, xend = 1, y = y_tick, yend = y_bar, linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2, y = y_tick, yend = y_bar, linewidth = 0.6) +
    annotate("text", x = 1.5, y = y_label, label = label, size = 4.5)
}

add_multiple_brackets <- function(plot, values, tests, group_levels) {
  value_range <- diff(range(values, na.rm = TRUE))
  if (!is.finite(value_range) || value_range == 0) value_range <- 0.01
  y_max <- max(values, na.rm = TRUE)
  for (i in seq_len(nrow(tests))) {
    x1 <- match(tests$group1[i], group_levels)
    x2 <- match(tests$group2[i], group_levels)
    y_tick <- y_max + (0.08 + 0.11 * (i - 1)) * value_range
    y_bar <- y_max + (0.11 + 0.11 * (i - 1)) * value_range
    y_label <- y_max + (0.145 + 0.11 * (i - 1)) * value_range
    plot <- plot +
      annotate("segment", x = x1, xend = x2, y = y_bar, yend = y_bar, linewidth = 0.55) +
      annotate("segment", x = x1, xend = x1, y = y_tick, yend = y_bar, linewidth = 0.55) +
      annotate("segment", x = x2, xend = x2, y = y_tick, yend = y_bar, linewidth = 0.55) +
      annotate(
        "text", x = mean(c(x1, x2)), y = y_label,
        label = format_p(tests$p_adjusted_BH[i]), size = 4
      )
  }
  plot
}

sample_scores <- data[[]] %>%
  filter(!is.na(Sample_Type), !is.na(Patient_ID), !is.na(Source_GSE)) %>%
  group_by(Source_GSE, Patient_ID, Sample_Type) %>%
  summarise(
    n_cells = n(),
    MPI_median = median(Polarization_Index, na.rm = TRUE),
    MMI_median = median(.data[[mmi_column]], na.rm = TRUE),
    .groups = "drop"
  )

group_tests <- bind_rows(
  paired_or_unpaired_wilcox(sample_scores, "MPI_median"),
  paired_or_unpaired_wilcox(sample_scores, "MMI_median")
) %>%
  mutate(p_adjusted_BH = p.adjust(p_value, method = "BH"))
write.csv(
  sample_scores,
  file.path(out_dir, "Top200_MMI_MPI_biological_sample_scores.csv"),
  row.names = FALSE
)
write.csv(
  group_tests,
  file.path(out_dir, "Top200_MMI_MPI_sample_level_tests.csv"),
  row.names = FALSE
)

# Theme and plot geometry follow the corresponding blocks in single_cell.R.
theme_custom <- theme_bw(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

message("Generating original-code-style MPI and MMI figures...")

plot_data_polar <- data[[]] %>%
  group_by(Celltype_raw1) %>%
  summarise(Mean_Polar = mean(Polarization_Index, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Mean_Polar)) %>%
  mutate(Celltype_raw1 = factor(Celltype_raw1, levels = Celltype_raw1))

p_mpi_bar <- ggplot(
  plot_data_polar,
  aes(x = Celltype_raw1, y = Mean_Polar, fill = Mean_Polar > 0)
) +
  geom_bar(stat = "identity", color = "white", width = 0.8) +
  scale_fill_manual(
    values = c("TRUE" = "#BC3C29FF", "FALSE" = "#0072B5FF"),
    labels = c("TRUE" = "M1 Dominant", "FALSE" = "M2 Dominant")
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  theme_classic() +
  labs(x = "Cell Subtypes", y = "Macrophage Polarization Index") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )
save_plot("01_Macrophage_Polarization_Direction_MPI", p_mpi_bar, 12, 6)

plot_data_mmi <- data[[]] %>%
  group_by(Celltype_raw1) %>%
  summarise(Mean_MMI = mean(.data[[mmi_column]], na.rm = TRUE), .groups = "drop") %>%
  arrange(Mean_MMI) %>%
  mutate(Celltype_raw1 = factor(Celltype_raw1, levels = Celltype_raw1))

p_mmi_bar <- ggplot(
  plot_data_mmi,
  aes(x = Celltype_raw1, y = Mean_MMI, fill = Mean_MMI)
) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
  theme_classic() +
  labs(x = "Cell Subtypes", y = "Macrophage Maturation Index (MMI)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "none"
  )
save_plot("02_Macrophage_Maturation_Index_MMI", p_mmi_bar, 12, 6)

p_mpi_umap <- FeaturePlot(
  data,
  features = "Polarization_Index",
  pt.size = 0.8,
  reduction = "umap",
  max.cutoff = "q98",
  order = TRUE
) +
  scale_colour_gradient2(
    low = "#0072B5FF", mid = "lightgrey", high = "#BC3C29FF",
    midpoint = 0, name = "MPI"
  ) +
  ggtitle("Macrophage Polarization Index (MPI)") +
  coord_flip()
save_plot("03_FeaturePlot_MPI_coord_flipped", p_mpi_umap, 10, 8)

mmi_limits <- as.numeric(quantile(
  data[[mmi_column, drop = TRUE]], mmi_color_quantiles, na.rm = TRUE
))
p_mmi_umap <- FeaturePlot(
  data,
  features = mmi_column,
  pt.size = 0.8,
  reduction = "umap",
  order = TRUE
) +
  scale_colour_gradient2(
    low = "#0072B5FF", mid = "lightgrey", high = "#FF0000",
    midpoint = 0, limits = mmi_limits, oob = scales::squish, name = "MMI"
  ) +
  ggtitle("Macrophage Maturation Index (MMI)") +
  coord_flip()
save_plot("04_FeaturePlot_MMI_coord_flipped", p_mmi_umap, 10, 8)
save_plot("05_FeaturePlot_MPI_MMI_combined", p_mpi_umap | p_mmi_umap, 18, 8)

plot_data_m1_m2 <- data[[]] %>%
  group_by(Celltype_raw1) %>%
  summarise(
    Mean_M1 = mean(M1_Score, na.rm = TRUE),
    Mean_M2 = mean(M2_Score, na.rm = TRUE),
    .groups = "drop"
  )
plot_long_mpi <- plot_data_m1_m2 %>%
  pivot_longer(c(Mean_M1, Mean_M2), names_to = "Metric", values_to = "Score") %>%
  mutate(Score = ifelse(Metric == "Mean_M2", -Score, Score))

p_mpi_components <- ggplot(
  plot_long_mpi,
  aes(x = Celltype_raw1, y = Score, fill = Metric)
) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = c("Mean_M1" = "#E64B35", "Mean_M2" = "#4DBBD5"),
    labels = c("Mean_M1" = "M1 AUC Score", "Mean_M2" = "M2 AUC Score")
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_y_continuous(labels = abs) +
  theme_classic() +
  labs(x = "Cell Subtypes", y = "AUC Score") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )
save_plot("06_Macrophage_AUCell_M1_vs_M2_Scores", p_mpi_components, 12, 6)

plot_data_mmi_components <- data[[]] %>%
  group_by(Celltype_raw1) %>%
  summarise(
    Mature = mean(.data[[mature_column]], na.rm = TRUE),
    Immature = mean(.data[[immature_column]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Mature, Immature), names_to = "Metric", values_to = "Score") %>%
  mutate(Score = ifelse(Metric == "Immature", -Score, Score))

p_mmi_components <- ggplot(
  plot_data_mmi_components,
  aes(x = Celltype_raw1, y = Score, fill = Metric)
) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = c("Mature" = "#31A354", "Immature" = "#756BB1"),
    labels = c("Mature" = "Mature AUC Score", "Immature" = "Immature AUC Score")
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_y_continuous(labels = abs) +
  theme_classic() +
  labs(x = "Cell Subtypes", y = "AUC Score") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )
save_plot("07_Macrophage_AUCell_Mature_vs_Immature_Scores", p_mmi_components, 12, 6)

sample_type_comparisons <- list(c("Atherosclerotic Core", "Proximal Adjacent"))
p_mpi_violin <- VlnPlot(
  data,
  features = "Polarization_Index",
  group.by = "Sample_Type",
  pt.size = 0,
  cols = c("#D95F02", "#1B9E77")
) +
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  labs(x = "Macrophage", y = "MPI") +
  stat_compare_means(
    comparisons = sample_type_comparisons, method = "wilcox.test",
    label = "p.signif", size = 5, bracket.size = 0.6
  ) +
  ylim(NA, max(data$Polarization_Index, na.rm = TRUE) * 1.5) +
  ggtitle("MPI: After Downsampling") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

p_mmi_violin <- VlnPlot(
  data,
  features = mmi_column,
  group.by = "Sample_Type",
  pt.size = 0,
  cols = c("#D95F02", "#1B9E77")
) +
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  labs(x = "Macrophage", y = "MMI") +
  stat_compare_means(
    comparisons = sample_type_comparisons, method = "wilcox.test",
    label = "p.signif", size = 5, bracket.size = 0.6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  coord_cartesian(clip = "off") +
  ggtitle("MMI: After Downsampling") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.margin = margin(t = 14, r = 8, b = 8, l = 8)
  )

save_plot("08_MPI_MMI_VlnPlot_original_code", p_mpi_violin | p_mmi_violin, 12, 6)

target_subtypes <- c("Foam cells1", "Foam cells2", "LAM")
sub_data_3cell <- subset(data, subset = Celltype_raw1 %in% target_subtypes)
sub_data_3cell$Celltype_raw1 <- factor(
  sub_data_3cell$Celltype_raw1,
  levels = target_subtypes
)
comparisons_3cell <- list(
  c("Foam cells1", "Foam cells2"),
  c("Foam cells1", "LAM"),
  c("Foam cells2", "LAM")
)
colors_3cell <- c("#1B9E77", "#D95F02", "#7570B3")

sample_subtype_scores <- data[[]] %>%
  filter(
    Celltype_raw1 %in% target_subtypes,
    !is.na(Patient_ID), !is.na(Source_GSE), !is.na(Sample_Type)
  ) %>%
  group_by(Source_GSE, Patient_ID, Sample_Type, Celltype_raw1) %>%
  summarise(
    n_cells = n(),
    MPI_median = median(Polarization_Index, na.rm = TRUE),
    MMI_median = median(.data[[mmi_column]], na.rm = TRUE),
    .groups = "drop"
  )

paired_subtype_test <- function(df, score_column, comparison) {
  paired <- df %>%
    filter(Celltype_raw1 %in% comparison) %>%
    select(Source_GSE, Patient_ID, Sample_Type, Celltype_raw1, all_of(score_column)) %>%
    pivot_wider(names_from = Celltype_raw1, values_from = all_of(score_column)) %>%
    filter(!is.na(.data[[comparison[1]]]), !is.na(.data[[comparison[2]]]))
  p_value <- if (nrow(paired) >= 3L) {
    wilcox.test(
      paired[[comparison[1]]], paired[[comparison[2]]],
      paired = TRUE, exact = FALSE
    )$p.value
  } else {
    NA_real_
  }
  tibble(
    score = score_column,
    group1 = comparison[1],
    group2 = comparison[2],
    analysis_unit = "paired biological sample",
    n_pairs = nrow(paired),
    p_value = p_value
  )
}

subtype_tests <- bind_rows(lapply(
  c("MPI_median", "MMI_median"),
  function(score_name) bind_rows(lapply(
    comparisons_3cell,
    function(comparison) paired_subtype_test(
      sample_subtype_scores, score_name, comparison
    )
  ))
)) %>%
  group_by(score) %>%
  mutate(p_adjusted_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write.csv(
  sample_subtype_scores,
  file.path(out_dir, "Top200_MMI_MPI_sample_subtype_scores.csv"),
  row.names = FALSE
)
write.csv(
  subtype_tests,
  file.path(out_dir, "Top200_MMI_MPI_sample_subtype_tests.csv"),
  row.names = FALSE
)

p_mpi_3cell <- VlnPlot(
  sub_data_3cell,
  features = "Polarization_Index",
  group.by = "Celltype_raw1",
  pt.size = 0,
  cols = colors_3cell
) +
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(
    comparisons = comparisons_3cell, method = "wilcox.test",
    label = "p.signif", size = 5, bracket.size = 0.6
  ) +
  labs(x = NULL, y = "Macrophage Polarization Index (MPI)") +
  ggtitle(NULL) +
  NoLegend() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

p_mmi_3cell <- VlnPlot(
  sub_data_3cell,
  features = mmi_column,
  group.by = "Celltype_raw1",
  pt.size = 0,
  cols = colors_3cell
) +
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(
    comparisons = comparisons_3cell, method = "wilcox.test",
    label = "p.signif", size = 5, bracket.size = 0.6
  ) +
  labs(x = NULL, y = "Macrophage Maturation Index (MMI)") +
  ggtitle(NULL) +
  NoLegend() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
save_plot("09_S2.6_MPI_MMI_3cell_VlnPlot", p_mpi_3cell | p_mmi_3cell, 12, 6)

density_data <- data[[]] %>%
  filter(Celltype_raw1 %in% target_subtypes) %>%
  group_by(Celltype_raw1) %>%
  filter(
    between(
      Polarization_Index,
      quantile(Polarization_Index, density_trim_quantiles[1], na.rm = TRUE),
      quantile(Polarization_Index, density_trim_quantiles[2], na.rm = TRUE)
    ),
    between(
      .data[[mmi_column]],
      quantile(.data[[mmi_column]], density_trim_quantiles[1], na.rm = TRUE),
      quantile(.data[[mmi_column]], density_trim_quantiles[2], na.rm = TRUE)
    )
  ) %>%
  ungroup()

p_density <- ggplot(
  density_data,
  aes(x = Polarization_Index, y = .data[[mmi_column]], color = Celltype_raw1)
) +
  geom_density_2d(
    linewidth = 0.8, alpha = 0.7,
    contour_var = "ndensity", breaks = density_contour_breaks
  ) +
  scale_color_manual(
    values = c(
      "LAM" = "#1B9E77",
      "Foam cells1" = "#D95F02",
      "Foam cells2" = "#7570B3"
    )
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  theme_custom +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_blank()
  ) +
  labs(x = "Polarization Index (MPI)", y = "Maturation Index (MMI)")
save_plot("10_Trajectory_Density_MPI_MMI", p_density, 8, 6)

subtype_summary <- data[[]] %>%
  group_by(Celltype_raw1) %>%
  summarise(
    n_cells = n(),
    MPI_mean = mean(Polarization_Index, na.rm = TRUE),
    MPI_median = median(Polarization_Index, na.rm = TRUE),
    Mature_AUC_mean = mean(.data[[mature_column]], na.rm = TRUE),
    Immature_AUC_mean = mean(.data[[immature_column]], na.rm = TRUE),
    MMI_mean = mean(.data[[mmi_column]], na.rm = TRUE),
    MMI_median = median(.data[[mmi_column]], na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  subtype_summary,
  file.path(out_dir, "Top200_MMI_MPI_subtype_summary.csv"),
  row.names = FALSE
)

sample_summary <- data[[]] %>%
  group_by(Sample_Type) %>%
  summarise(
    n_cells = n(),
    MPI_mean = mean(Polarization_Index, na.rm = TRUE),
    MPI_median = median(Polarization_Index, na.rm = TRUE),
    MMI_mean = mean(.data[[mmi_column]], na.rm = TRUE),
    MMI_median = median(.data[[mmi_column]], na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  sample_summary,
  file.path(out_dir, "Top200_MMI_MPI_sample_type_summary.csv"),
  row.names = FALSE
)

parameter_manifest <- tibble(
  parameter = c(
    "input_assay", "input_layer", "aucMaxRank", "normAUC",
    "external_sources", "final_top_n_per_direction",
    "cross_source_consensus_ranking", "direction_conflict_handling",
    "mature_positive_genes", "immature_negative_genes",
    "mature_genes_in_both_individual_source_top200",
    "immature_genes_in_both_individual_source_top200",
    "direction_discordant_genes_excluded",
    "selection_uses_carotid_study_data",
    "MPI_definition",
    "group_test_unit", "multiple_testing"
  ),
  value = c(
    "RNA", "counts", as.character(auc_max_rank), "TRUE",
    "GSE5099;GSE11864", as.character(top_n),
    "mean_of_source_specific_directional_percentile_ranks",
    "exclude_genes_with_opposite_effect_signs_between_sources",
    as.character(length(mature_genes)), as.character(length(immature_genes)),
    as.character(sum(
      membership$direction == "mature_positive" &
        membership$selected_in_GSE5099_individual_top200 &
        membership$selected_in_GSE11864_individual_top200
    )),
    as.character(sum(
      membership$direction == "immature_negative" &
        membership$selected_in_GSE5099_individual_top200 &
        membership$selected_in_GSE11864_individual_top200
    )),
    as.character(length(conflicting_genes)),
    "FALSE",
    "AUCell(frozen 145-gene M1-like) minus AUCell(frozen 165-gene M2-like)",
    "biological sample or paired patient", "Benjamini-Hochberg"
  )
)
write.csv(
  parameter_manifest,
  file.path(out_dir, "Top200_MMI_MPI_parameter_manifest.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

writeLines(
  c(
    "FINAL ANALYSIS: GSE5099 + GSE11864 Top200 MPI/MMI.",
    "MMI definition: Top200_Mature_AUCell - Top200_Monocyte_AUCell.",
    "External sources: GSE5099 and GSE11864 only.",
    paste0(
      "The final gene sets contain exactly ", top_n,
      " mature-positive and ", top_n, " monocyte-high genes."
    ),
    paste0(
      "Only genes measured in both external sources and showing the same differential direction ",
      "were eligible; ", length(conflicting_genes),
      " genes with opposite effect signs were excluded."
    ),
    paste0(
      "Eligible genes were ranked in each source and direction, converted to percentile ranks, ",
      "and ordered by the mean of the two source-specific percentiles."
    ),
    "No carotid expression, carotid subtype label, or carotid MMI result was used for gene selection.",
    paste0("Final mature-positive genes: ", length(mature_genes), "."),
    paste0("Final monocyte-high genes: ", length(immature_genes), "."),
    "The full direction-discordant audit list is provided as a separate CSV file.",
    "AUCell workflow follows the main code: RNA/counts rankings, normalized AUC, top 5% aucMaxRank.",
    "MPI is Polarization_Index = M1_Score - M2_Score using the frozen 145/165 literature-derived signatures.",
    "Plot geometry follows GSE159677_Carotid_MainProject/single_cell.R.",
    "Violin annotations retain the legacy cell-level stat_compare_means code.",
    "Biological-sample and paired-patient tests are exported separately as CSV files.",
    "Benjamini-Hochberg adjustment is applied within each exported sample-level test family.",
    paste0("MMI FeaturePlot limits: quantiles ", paste(mmi_color_quantiles, collapse = "-"), "."),
    paste0("Density display trim by subtype: quantiles ", paste(density_trim_quantiles, collapse = "-"), "."),
    paste0("Density contour normalized-density breaks: ", paste(density_contour_breaks, collapse = ", "), ".")
  ),
  file.path(out_dir, "README.txt")
)

message("Finished. Outputs written to: ", out_dir)
