#!/usr/bin/env Rscript

# Score an external consensus MMI with the AUCell workflow used by single_cell.R.
# The exploratory C1Q source module is excluded; directional conflicts are removed.

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

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
set.seed(cfg$project$seed)

input_rds <- project_path(cfg, cfg$outputs$auc_scored_rds)
source_gene_set_file <- project_path(
  cfg, cfg$auc_signatures$external_mmi_source_gene_sets
)
out_dir <- project_path(cfg, cfg$auc_signatures$publication_output_dir)
score_rds <- file.path(out_dir, "MMI_AUCell_noC1Q_cell_scores.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mature_column <- "External_Mature_AUCell"
immature_column <- "External_Monocyte_Immaturity_AUCell"
mmi_column <- "MMI_AUCell_noC1Q"

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

message("Preparing the external consensus gene sets...")
source_sets <- read.csv(source_gene_set_file, stringsAsFactors = FALSE) %>%
  filter(source_set != "C1Q_TissueMac_Core")

expected_sources <- c(
  "GSE5099_Macrophage", "GSE5099_Monocyte",
  "GSE11864_Macrophage", "GSE11864_Monocyte",
  "HPCA_Macrophage", "HPCA_Monocyte"
)
missing_sources <- setdiff(expected_sources, unique(source_sets$source_set))
if (length(missing_sources)) {
  stop("Missing required source gene sets: ", paste(missing_sources, collapse = ", "))
}

positive_union <- unique(source_sets$gene_symbol[
  source_sets$direction == "mature_positive"
])
negative_union <- unique(source_sets$gene_symbol[
  source_sets$direction == "immature_negative"
])
conflicting_genes <- sort(intersect(positive_union, negative_union))
mature_genes <- sort(setdiff(positive_union, conflicting_genes))
immature_genes <- sort(setdiff(negative_union, conflicting_genes))

if (length(intersect(mature_genes, immature_genes))) {
  stop("Positive and negative gene sets still overlap after conflict removal.")
}
if (length(mature_genes) < 50L || length(immature_genes) < 50L) {
  stop("Too few genes remain after conflict removal.")
}

write.csv(
  bind_rows(
    tibble(direction = "mature_positive", gene_symbol = mature_genes),
    tibble(direction = "immature_negative", gene_symbol = immature_genes)
  ),
  file.path(out_dir, "MMI_AUCell_noC1Q_clean_gene_sets.csv"),
  row.names = FALSE
)
write.csv(
  tibble(gene_symbol = conflicting_genes, action = "removed_from_both_directions"),
  file.path(out_dir, "MMI_AUCell_noC1Q_removed_conflicting_genes.csv"),
  row.names = FALSE
)
write.csv(
  source_sets,
  file.path(out_dir, "MMI_AUCell_noC1Q_source_gene_sets.csv"),
  row.names = FALSE
)

message(
  "Clean gene sets: ", length(mature_genes), " mature-positive and ",
  length(immature_genes), " immature-negative genes; removed ",
  length(conflicting_genes), " conflicts."
)

message("Loading the original scored monocyte/macrophage object...")
data <- readRDS(input_rds)
required_meta <- c(
  "Celltype_raw1", "Sample_Type", "Patient_ID", "Source_GSE",
  "Polarization_Index", "M1_Score", "M2_Score"
)
missing_meta <- setdiff(required_meta, colnames(data[[]]))
if (length(missing_meta)) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}
if (!"RNA" %in% Assays(data) || !"counts" %in% Layers(data[["RNA"]])) {
  stop("The RNA assay must contain a counts layer for the original AUCell workflow.")
}

rna_features <- rownames(data[["RNA"]])
mature_genes <- intersect(mature_genes, rna_features)
immature_genes <- intersect(immature_genes, rna_features)
auc_max_rank <- ceiling(
  cfg$auc_signatures$auc_max_rank_fraction * length(rna_features)
)
if (length(intersect(mature_genes, immature_genes))) {
  stop("Detected positive and negative gene sets overlap unexpectedly.")
}

score_columns <- c(mature_column, immature_column, mmi_column)
if (all(score_columns %in% colnames(data[[]]))) {
  message("Using AUCell scores stored in the pipeline object...")
  score_meta <- data[[]] %>%
    rownames_to_column("cell_barcode") %>%
    select(cell_barcode, all_of(score_columns))
  saveRDS(score_meta, score_rds)
  write.csv(
    score_meta,
    gzfile(file.path(out_dir, "MMI_AUCell_noC1Q_cell_scores.csv.gz")),
    row.names = FALSE
  )
} else if (file.exists(score_rds)) {
  message("Using cached AUCell scores...")
  score_meta <- readRDS(score_rds)
  score_index <- match(colnames(data), score_meta$cell_barcode)
  if (anyNA(score_index) || !all(score_columns %in% colnames(score_meta))) {
    stop("Cached AUCell scores do not match the current Seurat object.")
  }
  new_scores <- score_meta[score_index, score_columns, drop = FALSE]
  rownames(new_scores) <- colnames(data)
  data <- AddMetaData(data, new_scores)
} else {
  message("Building per-cell gene rankings from RNA/counts...")
  expression_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
  cells_rankings <- AUCell_buildRankings(
    expression_matrix,
    plotStats = FALSE,
    splitByBlocks = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(
      workers = cfg$auc_signatures$ncores
    ),
    verbose = TRUE
  )

  gene_sets <- list(
    External_Mature = mature_genes,
    External_Monocyte_Immaturity = immature_genes
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
  data <- AddMetaData(data, auc_scores)

  score_meta <- auc_scores %>%
    rownames_to_column("cell_barcode")
  saveRDS(score_meta, score_rds)
  write.csv(
    score_meta,
    gzfile(file.path(out_dir, "MMI_AUCell_noC1Q_cell_scores.csv.gz")),
    row.names = FALSE
  )
  rm(expression_matrix, cells_rankings, cells_auc, auc_scores)
  invisible(gc())
}

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
  file.path(out_dir, "MMI_AUCell_noC1Q_biological_sample_scores.csv"),
  row.names = FALSE
)
write.csv(
  group_tests,
  file.path(out_dir, "MMI_AUCell_noC1Q_sample_level_tests.csv"),
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
  geom_point(
    data = sample_scores,
    aes(x = Sample_Type, y = MPI_median),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.08, height = 0),
    shape = 21, size = 2.2, stroke = 0.45, fill = "white", color = "black"
  ) +
  labs(x = "Macrophage", y = "MPI") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  ggtitle("MPI by lesion region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
p_mpi_violin <- add_single_bracket(
  p_mpi_violin,
  data$Polarization_Index,
  paste0(
    group_tests$analysis_unit[group_tests$score == "MPI_median"],
    " n=", group_tests$n_pairs[group_tests$score == "MPI_median"], ", ",
    format_p(group_tests$p_adjusted_BH[group_tests$score == "MPI_median"])
  )
)

p_mmi_violin <- VlnPlot(
  data,
  features = mmi_column,
  group.by = "Sample_Type",
  pt.size = 0,
  cols = c("#D95F02", "#1B9E77")
) +
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  geom_point(
    data = sample_scores,
    aes(x = Sample_Type, y = MMI_median),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.08, height = 0),
    shape = 21, size = 2.2, stroke = 0.45, fill = "white", color = "black"
  ) +
  labs(x = "Macrophage", y = "MMI") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  ggtitle("MMI by lesion region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
p_mmi_violin <- add_single_bracket(
  p_mmi_violin,
  data[[mmi_column, drop = TRUE]],
  paste0(
    group_tests$analysis_unit[group_tests$score == "MMI_median"],
    " n=", group_tests$n_pairs[group_tests$score == "MMI_median"], ", ",
    format_p(group_tests$p_adjusted_BH[group_tests$score == "MMI_median"])
  )
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
  file.path(out_dir, "MMI_AUCell_noC1Q_sample_subtype_scores.csv"),
  row.names = FALSE
)
write.csv(
  subtype_tests,
  file.path(out_dir, "MMI_AUCell_noC1Q_sample_subtype_tests.csv"),
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
  geom_point(
    data = sample_subtype_scores,
    aes(x = Celltype_raw1, y = MPI_median),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.08, height = 0),
    shape = 21, size = 1.8, stroke = 0.4, fill = "white", color = "black"
  ) +
  labs(x = NULL, y = "Macrophage Polarization Index (MPI)") +
  ggtitle(NULL) +
  NoLegend() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
p_mpi_3cell <- add_multiple_brackets(
  p_mpi_3cell,
  sub_data_3cell$Polarization_Index,
  subtype_tests %>% filter(score == "MPI_median"),
  target_subtypes
)

p_mmi_3cell <- VlnPlot(
  sub_data_3cell,
  features = mmi_column,
  group.by = "Celltype_raw1",
  pt.size = 0,
  cols = colors_3cell
) +
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  geom_point(
    data = sample_subtype_scores,
    aes(x = Celltype_raw1, y = MMI_median),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.08, height = 0),
    shape = 21, size = 1.8, stroke = 0.4, fill = "white", color = "black"
  ) +
  labs(x = NULL, y = "Macrophage Maturation Index (MMI)") +
  ggtitle(NULL) +
  NoLegend() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
p_mmi_3cell <- add_multiple_brackets(
  p_mmi_3cell,
  sub_data_3cell[[mmi_column, drop = TRUE]],
  subtype_tests %>% filter(score == "MMI_median"),
  target_subtypes
)
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
  file.path(out_dir, "MMI_AUCell_noC1Q_subtype_summary.csv"),
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
  file.path(out_dir, "MMI_AUCell_noC1Q_sample_type_summary.csv"),
  row.names = FALSE
)

parameter_manifest <- tibble(
  parameter = c(
    "input_assay", "input_layer", "aucMaxRank", "normAUC",
    "mature_positive_genes", "immature_negative_genes",
    "conflicting_genes_removed", "C1Q_source_module_included",
    "group_test_unit", "multiple_testing"
  ),
  value = c(
    "RNA", "counts", as.character(auc_max_rank), "TRUE",
    as.character(length(mature_genes)), as.character(length(immature_genes)),
    as.character(length(conflicting_genes)), "FALSE",
    "biological sample or paired patient", "Benjamini-Hochberg"
  )
)
write.csv(
  parameter_manifest,
  file.path(out_dir, "MMI_AUCell_noC1Q_parameter_manifest.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

writeLines(
  c(
    "MMI definition: External_Mature_AUCell - External_Monocyte_Immaturity_AUCell.",
    "External sources: GSE5099, GSE11864, and HPCA-derived fixed markers.",
    "The separately curated C1Q tissue-macrophage source module was excluded before union construction.",
    paste0("Mature-positive genes after cleaning: ", length(mature_genes), "."),
    paste0("Immature-negative genes after cleaning: ", length(immature_genes), "."),
    paste0("Removed from both directions: ", paste(conflicting_genes, collapse = ", "), "."),
    "AUCell workflow follows the main code: RNA/counts rankings, normalized AUC, top 5% aucMaxRank.",
    "MPI is the existing Polarization_Index and was not recalculated or modified.",
    "Plot geometry follows GSE159677_Carotid_MainProject/single_cell.R.",
    "Violin geometry shows cells; overlaid white points and inferential tests use biological-sample medians.",
    "Paired patient-level Wilcoxon tests are used when at least three complete pairs are available.",
    "Benjamini-Hochberg adjustment is applied within each displayed family of tests.",
    paste0("MMI FeaturePlot limits: quantiles ", paste(mmi_color_quantiles, collapse = "-"), "."),
    paste0("Density display trim by subtype: quantiles ", paste(density_trim_quantiles, collapse = "-"), "."),
    paste0("Density contour normalized-density breaks: ", paste(density_contour_breaks, collapse = ", "), ".")
  ),
  file.path(out_dir, "README.txt")
)

message("Finished. Outputs written to: ", out_dir)
