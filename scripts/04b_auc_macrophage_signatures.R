#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(AUCell)
  library(showtext)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()

m1_features <- c(
  "ACOD1", "AFDN", "AIM1", "AKAP13", "ALAS1", "ANKRD22", "APOBEC3A", "APOL2",
  "APOL3", "BST2", "C15orf48", "CALHM6", "CCL18", "CCL19", "CCL5", "CCR10",
  "CCR7", "CD1B", "CD274", "CD38", "CD40", "CD64", "CD74", "CD80", "CD86",
  "CDC42SE2", "COG6", "CRYBG1", "CSTF3", "CXCL10", "CXCL11", "CXCL8", "CXCL9",
  "CXCR10", "DEFA3", "EBI3", "ELOVL5", "ENSA", "EPSTI1", "FAM177A1", "FAM26F",
  "FAS", "FBP1", "FBXO6", "FCGR1A", "FCGR1B", "FCGR1C", "FDX1", "FUS", "GBP1",
  "GBP4", "GBP5", "GLS", "GPAT3", "HLA-A", "HLA-DMA", "HLA-DMB", "HLA-DRA",
  "HLA-DRB1", "HLA-DRB3", "HSD11B1", "HSP90AB4P", "IDO1", "IFI44", "IFIH1",
  "IFIT2", "IFIT3", "IFIT5", "IL12", "IL12A", "IL12B", "IL1A", "IL1B", "IL23",
  "IL23A", "IL6", "IL8", "INOS", "IRF1", "IRF5", "IRG1", "ISG15", "ISG20",
  "KMO", "KYNU", "LAMP3", "LGALS3BP", "LHFPL2", "LTF", "LY75", "MARCKSL1",
  "MFF", "MGEA5", "MGST1", "MHCII", "MLLT4", "MRAS", "MT2A", "MX1", "MX2",
  "NDRG2", "NFKB2", "NMES1", "NOS2", "NT5C3A", "NUB1", "OAS2", "OAS3", "OASL",
  "OGA", "P2RX7", "PLAUR", "PLD1", "PML", "PNPLA6", "PNPT1", "PTGES", "PTGS2",
  "PTX3", "RBM17", "RCN1", "RHOF", "RIPK2", "RSAD2", "SLAMF1", "SLAMF7",
  "SLC15A3", "SLC27A3", "SLC29A3", "SLC2A3", "SLC2A6", "SOAT1", "SPN", "STX11",
  "TAP1", "TAP2", "TAPBP", "TNF", "TNFA", "TNFAIP3", "TRAF1", "VAMP5", "WARS",
  "WARS1", "XIRP1"
)

m2_features <- c(
  "ABI3", "ACOT11", "ADAP1", "ADAP2", "ADORA3", "ALDH1A1", "ALOX15", "APPL2",
  "ARG1", "ARG2", "ARHGAP26", "ARHGAP4", "ARSA", "ARSB", "BABAM2", "BIN1",
  "BLVRB", "BRE", "CCL13", "CCL17", "CCL18", "CCL20", "CCL22", "CCL24", "CCL4",
  "CD14", "CD163", "CD163L1", "CD200R", "CD200R1", "CD206", "CD209", "CD23",
  "CD274", "CD276", "CD32", "CD36", "CHMP2A", "CLEC7A", "CNRIP1", "COMMD1",
  "CRYL1", "CSF1R", "CST3", "CTSA", "CTSB", "CTSC", "CTSD", "CUL4B", "CYB5R4",
  "DAB2", "DCD", "EGF", "EMB", "F13A1", "FAH", "FASL", "FASLG", "FCER2",
  "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGRT", "FIGF", "FN1", "FOLR2",
  "FUCA1", "GALE", "GAS7", "GATA3", "GATM", "GLMP", "GLUL", "GNPDA1", "GPR183",
  "GRAMD4", "HAVCR2", "HECTD3", "HEXA", "HEXB", "HMOX1", "IL10", "IL17RB",
  "IL1R2", "IL1RA", "IL1RN", "IL4R", "IL4RA", "IRF4", "ITSN1", "LACC1", "LGMN",
  "LRP1", "LYVE1", "MANBA", "MARCO", "ME1", "MGLL", "MMP1", "MMP12", "MMP14",
  "MMP19", "MMP9", "MPEG1", "MPI", "MRC1", "MSR1", "NAIP", "NAPRT", "NDUFA4",
  "NDUFB3", "NEU1", "NIF3L1", "NLN", "NPL", "NUBP1", "NUDT2", "P2RY11", "PARP1",
  "PDCD1LG2", "PDPK1", "PITHD1", "PLA2G15", "PLXDC2", "PMVK", "PREX1", "PRKCE",
  "QPRT", "RASA1", "RENBP", "RNASE6", "RNASET2", "SAR1B", "SDF4", "SDSL",
  "SERINC1", "SERPINB2", "SGPL1", "SLC9A9", "SLCO2B1", "SOCS1", "SOCS3", "STAB1",
  "SYPL1", "TANGO2", "TG", "TGFB1", "TGFB2", "TGFB3", "TGFBR2", "TGM2",
  "TMEM176B", "TNFSF12", "TNFSF8", "TRIM47", "VEGFA", "VEGFB", "VEGFC", "VEGFD",
  "VPS50", "VTCN1", "WDR64", "WDR81", "WNT7B"
)

load_external_mmi_gene_sets <- function(cfg, obj) {
  gene_set_file <- project_path(cfg, cfg$auc_signatures$external_mmi_gene_sets)
  if (!file.exists(gene_set_file)) {
    stop("External MMI gene-set file does not exist: ", gene_set_file, call. = FALSE)
  }
  gene_set_table <- read.csv(gene_set_file, stringsAsFactors = FALSE)
  required_columns <- c("direction", "gene_symbol")
  if (!all(required_columns %in% colnames(gene_set_table))) {
    stop("External MMI gene-set file must contain: ", paste(required_columns, collapse = ", "), call. = FALSE)
  }

  mature <- unique(gene_set_table$gene_symbol[
    gene_set_table$direction == "mature_positive"
  ])
  immature <- unique(gene_set_table$gene_symbol[
    gene_set_table$direction == "immature_negative"
  ])
  conflicts <- intersect(mature, immature)
  if (length(conflicts) > 0) {
    stop("External MMI directions overlap: ", paste(conflicts, collapse = ", "), call. = FALSE)
  }

  detected <- rownames(obj[["RNA"]])
  sets <- list(
    External_Mature_AUCell = intersect(mature, detected),
    External_Monocyte_Immaturity_AUCell = intersect(immature, detected)
  )
  if (any(lengths(sets) < 50L)) {
    stop("Fewer than 50 external MMI genes were detected in one direction.", call. = FALSE)
  }
  sets
}

input_rds <- project_path(cfg, cfg$outputs$hdwgcnna_rds)
message_step("Loading object for original MPI and external-consensus MMI scoring: ", input_rds)
obj <- readRDS(input_rds)

external_mmi_sets <- load_external_mmi_gene_sets(cfg, obj)
gene_sets <- c(
  list(M1_Score = m1_features, M2_Score = m2_features),
  external_mmi_sets
)

expr <- GetAssayData(obj, assay = "RNA", layer = "counts")
rankings <- AUCell_buildRankings(
  expr,
  plotStats = FALSE,
  splitByBlocks = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(
    workers = cfg$auc_signatures$ncores
  )
)
auc_max_rank <- ceiling(cfg$auc_signatures$auc_max_rank_fraction * nrow(rankings))
auc <- AUCell_calcAUC(
  gene_sets,
  rankings,
  normAUC = TRUE,
  aucMaxRank = auc_max_rank,
  nCores = 1
)
auc_scores <- as.data.frame(t(getAUC(auc)))

auc_scores <- auc_scores |>
  mutate(
    MMI_AUCell_noC1Q = External_Mature_AUCell - External_Monocyte_Immaturity_AUCell,
    AMDI_Index = MMI_AUCell_noC1Q,
    Polarization_Index = M1_Score - M2_Score
  )

obj <- AddMetaData(obj, auc_scores)
out_rds <- project_path(cfg, cfg$outputs$auc_scored_rds)
safe_save_rds(obj, out_rds)

plot_dir <- project_path(
  cfg,
  "results/figures/legacy_cell_level_auc_plots_NOT_FOR_INFERENCE"
)
ensure_dir(plot_dir)
warning(
  "Plots from 04b preserve the legacy cell-level display only. ",
  "Use outputs from 04c_publication_mpi_mmi_plots.R for inference."
)
cell_col <- cfg$cell_types$cell_type_column

setup_plot_font <- function(cfg) {
  font_family <- cfg$auc_signatures$font_family
  font_files <- unlist(cfg$auc_signatures$fonts)
  if (!is.null(font_family) && length(font_files) == 4 && all(file.exists(font_files))) {
    showtext::font_add(
      family = font_family,
      regular = font_files[["regular"]],
      bold = font_files[["bold"]],
      italic = font_files[["italic"]],
      bolditalic = font_files[["bolditalic"]]
    )
    showtext::showtext_auto()
    return(font_family)
  }
  "sans"
}

normalize_sample_type <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("AC", "Atherosclerotic Core", "Carotid Atherosclerotic Core") ~ "Atherosclerotic Core",
    x %in% c("PA", "Proximal Adjacent", "Carotid Proximal Adjacent") ~ "Proximal Adjacent",
    TRUE ~ x
  )
}

require_metadata <- function(obj, cols) {
  missing <- setdiff(cols, colnames(obj@meta.data))
  if (length(missing) > 0) {
    stop("AUCell scored object is missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

get_sig_label <- function(p) {
  if (is.na(p)) {
    return("ns")
  }
  if (p <= 0.0001) {
    return("****")
  }
  if (p <= 0.001) {
    return("***")
  }
  if (p <= 0.01) {
    return("**")
  }
  if (p <= 0.05) {
    return("*")
  }
  "ns"
}

safe_wilcox_label <- function(df, y_col, group_col) {
  complete <- df[!is.na(df[[y_col]]) & !is.na(df[[group_col]]), , drop = FALSE]
  if (length(unique(complete[[group_col]])) != 2) {
    return("ns")
  }
  get_sig_label(wilcox.test(complete[[y_col]] ~ complete[[group_col]])$p.value)
}

make_index_vln <- function(df, y_col, y_lab, font_family_use, theme_vln) {
  y_values <- df[[y_col]]
  y_min <- min(y_values, na.rm = TRUE)
  y_max <- max(y_values, na.rm = TRUE)
  y_range <- y_max - y_min
  if (y_range == 0) {
    y_range <- max(abs(y_max) * 0.1, 0.01)
  }

  y_bracket <- y_max + 0.12 * y_range
  y_tick <- y_max + 0.10 * y_range
  y_label <- y_max + 0.155 * y_range
  y_bottom <- y_min - 0.08 * y_range
  y_top <- y_max + 0.25 * y_range
  sig_label <- safe_wilcox_label(df, y_col, "Sample_Type")

  ggplot(df, aes(x = Sample_Type, y = .data[[y_col]], fill = Sample_Type)) +
    geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.45, alpha = 1) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.7, linewidth = 0.45) +
    scale_fill_manual(values = c("Atherosclerotic Core" = "#D95F02", "Proximal Adjacent" = "#1B9E77")) +
    geom_segment(aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    geom_segment(aes(x = 1, xend = 1, y = y_tick, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    geom_segment(aes(x = 2, xend = 2, y = y_tick, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    annotate("text", x = 1.5, y = y_label, label = sig_label, size = 6, family = font_family_use, fontface = "bold") +
    labs(x = NULL, y = y_lab) +
    coord_cartesian(ylim = c(y_bottom, y_top), clip = "off") +
    theme_vln
}

make_3cell_index_vln <- function(df, y_col, y_lab, comparisons, theme_vln) {
  y_min <- min(df[[y_col]], na.rm = TRUE)
  y_max <- max(df[[y_col]], na.rm = TRUE)
  y_range <- y_max - y_min
  if (y_range == 0) {
    y_range <- max(abs(y_max) * 0.1, 0.01)
  }

  present_groups <- as.character(unique(df$Celltype_raw1))
  valid_comparisons <- Filter(function(x) all(x %in% present_groups), comparisons)

  p <- ggplot(df, aes(x = Celltype_raw1, y = .data[[y_col]], fill = Celltype_raw1)) +
    geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.45, alpha = 1) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.7, linewidth = 0.45) +
    scale_fill_manual(values = c("Foam cells1" = "#1B9E77", "Foam cells2" = "#D95F02", "LAM" = "#7570B3"))

  if (length(valid_comparisons) > 0) {
    p <- p + stat_compare_means(
      comparisons = valid_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      size = 7.5,
      bracket.size = 0.6,
      tip.length = 0.02,
      y.position = y_max + seq(0.12, by = 0.10, length.out = length(valid_comparisons)) * y_range
    )
  }

  p +
    labs(x = NULL, y = y_lab) +
    coord_cartesian(ylim = c(y_min - 0.08 * y_range, y_max + 0.42 * y_range), clip = "off") +
    theme_vln
}

save_plot_pair <- function(plot, basename, width, height) {
  ggsave(file.path(plot_dir, paste0(basename, ".pdf")), plot = plot, device = cairo_pdf, width = width, height = height)
  ggsave(file.path(plot_dir, paste0(basename, ".png")), plot = plot, width = width, height = height, dpi = 300)
}

require_metadata(obj, c("Sample_Type", "Celltype_raw1", "Polarization_Index", "AMDI_Index"))
obj$Sample_Type <- factor(normalize_sample_type(obj$Sample_Type), levels = c("Atherosclerotic Core", "Proximal Adjacent"))
obj$Celltype_raw1 <- normalize_foam_labels(obj$Celltype_raw1)

font_family_use <- setup_plot_font(cfg)
theme_custom_plus2 <- theme_bw(base_family = font_family_use, base_size = 13) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
theme_vln_like_example <- theme_classic(base_family = font_family_use) +
  theme(
    text = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )
theme_vln_3cell <- theme_vln_like_example +
  theme(plot.margin = margin(t = 18, r = 10, b = 10, l = 10))

plot_df <- obj@meta.data |>
  filter(!is.na(Sample_Type)) |>
  select(Sample_Type, Polarization_Index, AMDI_Index)
if (nrow(plot_df) > 0) {
  p_mpi_amdi <- make_index_vln(plot_df, "Polarization_Index", "Macrophage Polarization Index (MPI)", font_family_use, theme_vln_like_example) |
    make_index_vln(plot_df, "AMDI_Index", "Macrophage Maturation Index (MMI)", font_family_use, theme_vln_like_example)
  save_plot_pair(p_mpi_amdi, "S2.6_MPI_AMDI_VlnPlot_exampleStyle_fontPlus2", 10, 6)
}

celltype_order_3 <- c("Foam cells1", "Foam cells2", "LAM")
celltype_comparisons_3 <- list(c("Foam cells1", "Foam cells2"), c("Foam cells2", "LAM"), c("Foam cells1", "LAM"))
plot_df_3cell <- obj@meta.data |>
  filter(Celltype_raw1 %in% celltype_order_3) |>
  select(Celltype_raw1, Polarization_Index, AMDI_Index) |>
  mutate(Celltype_raw1 = factor(Celltype_raw1, levels = celltype_order_3))
if (nrow(plot_df_3cell) > 0 && length(unique(plot_df_3cell$Celltype_raw1)) >= 2) {
  p_mpi_mdi_3cell <- make_3cell_index_vln(plot_df_3cell, "Polarization_Index", "Macrophage Polarization Index (MPI)", celltype_comparisons_3, theme_vln_3cell) |
    make_3cell_index_vln(plot_df_3cell, "AMDI_Index", "Macrophage Differentiation Index (MDI)", celltype_comparisons_3, theme_vln_3cell)
  save_plot_pair(p_mpi_mdi_3cell, "S2.6_MPI_MDI_3cell_VlnPlot_fontPlus2", 12, 6)
}

density_df <- obj@meta.data |>
  filter(Celltype_raw1 %in% c("LAM", "Foam cells1", "Foam cells2")) |>
  mutate(Celltype_raw1 = factor(Celltype_raw1, levels = c("LAM", "Foam cells1", "Foam cells2")))
if (nrow(density_df) > 0) {
  p_density <- ggplot(density_df, aes(x = Polarization_Index, y = AMDI_Index)) +
    geom_density_2d(aes(color = Celltype_raw1), linewidth = 0.8, alpha = 0.7) +
    scale_color_manual(values = c("LAM" = "#1B9E77", "Foam cells1" = "#D95F02", "Foam cells2" = "#7570B3")) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    theme_custom_plus2 +
    theme(legend.position = "right", legend.title = element_text(face = "bold"), legend.key = element_blank()) +
    labs(x = "Polarization Index (MPI)", y = "Maturation Index (AMDI)", color = "Cell type")
  save_plot_pair(p_density, "S2.6_Trajectory_Density_fontPlus2", 8, 6)
}

target_genes <- c("CXCL8", "IL1B", "TIMP1", "FOLR2", "TREM2", "C1QB", "APOE", "PLIN2")
target_genes <- target_genes[target_genes %in% rownames(obj)]
sub_data_foam <- subset(obj, subset = Celltype_raw1 %in% c("LAM", "Foam cells1", "Foam cells2"))
if (length(target_genes) > 0 && ncol(sub_data_foam) > 0) {
  sub_data_foam$Celltype_raw1 <- factor(sub_data_foam$Celltype_raw1, levels = c("LAM", "Foam cells1", "Foam cells2"))
  present_groups <- as.character(unique(sub_data_foam$Celltype_raw1))
  foam_comparisons <- list(c("LAM", "Foam cells1"), c("LAM", "Foam cells2"), c("Foam cells1", "Foam cells2"))
  foam_comparisons <- Filter(function(x) all(x %in% present_groups), foam_comparisons)
  plots_foam <- lapply(target_genes, function(gene) {
    p <- VlnPlot(
      sub_data_foam,
      features = gene,
      group.by = "Celltype_raw1",
      pt.size = 0,
      cols = c("#1B9E77", "#D95F02", "#7570B3")
    ) +
      theme_custom_plus2 +
      labs(y = "Expression Level", x = "")

    if (length(foam_comparisons) > 0) {
      p <- p + stat_compare_means(
        comparisons = foam_comparisons,
        method = "wilcox.test",
        label = "p.signif",
        bracket.size = 0.6,
        tip.length = 0.02,
        size = 7,
        vjust = 0.5
      )
    }

    p +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
      theme(
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5)
      )
  })
  combined_violin_plots <- wrap_plots(plots_foam, ncol = 3) & NoLegend()
  save_plot_pair(combined_violin_plots, "S2.6_Top10_Foam_Markers_VlnPlot_fontPlus2", 15, 15)
}

message_step("Saved original MPI and external-consensus MMI scores: ", out_rds)
