#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(AUCell)
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

load_mono_immaturity_genes <- function(cfg, obj) {
  marker_table <- project_path(cfg, cfg$auc_signatures$marker_table)
  if (file.exists(marker_table)) {
    markers <- read.csv(marker_table, check.names = FALSE)
    return(markers |>
      filter(cluster == cfg$auc_signatures$mono_cluster_label) |>
      arrange(p_val_adj, desc(avg_log2FC)) |>
      slice_head(n = cfg$auc_signatures$mono_top_n) |>
      pull(gene) |>
      unique())
  }

  cell_col <- cfg$cell_types$cell_type_column
  if (!cell_col %in% colnames(obj@meta.data)) {
    warning("No marker table or cell-type column found; Mono_Immaturity_Score will be skipped")
    return(character())
  }
  Idents(obj) <- cell_col
  markers <- FindMarkers(obj, ident.1 = cfg$auc_signatures$mono_cluster_label, only.pos = TRUE)
  markers |>
    tibble::rownames_to_column("gene") |>
    arrange(p_val_adj, desc(avg_log2FC)) |>
    slice_head(n = cfg$auc_signatures$mono_top_n) |>
    pull(gene)
}

input_rds <- project_path(cfg, cfg$outputs$macspectrum_rds)
if (!file.exists(input_rds)) {
  input_rds <- project_path(cfg, cfg$outputs$hdwgcnna_rds)
}
message_step("Loading object for custom AUCell scoring: ", input_rds)
obj <- readRDS(input_rds)

mono_features <- load_mono_immaturity_genes(cfg, obj)
gene_sets <- list(M1_Score = m1_features, M2_Score = m2_features)
if (length(mono_features) > 0) {
  gene_sets$Mono_Immaturity_Score <- mono_features
}

expr <- GetAssayData(obj, assay = "RNA", slot = "counts")
rankings <- AUCell_buildRankings(expr, nCores = cfg$auc_signatures$ncores, plotStats = FALSE)
auc <- AUCell_calcAUC(gene_sets, rankings)
auc_scores <- as.data.frame(t(getAUC(auc)))

auc_scores <- auc_scores |>
  mutate(
    AMDI_Index = if ("Mono_Immaturity_Score" %in% colnames(auc_scores)) -Mono_Immaturity_Score else NA_real_,
    Polarization_Index = M1_Score - M2_Score
  )

obj <- AddMetaData(obj, auc_scores)
out_rds <- project_path(cfg, cfg$outputs$auc_scored_rds)
safe_save_rds(obj, out_rds)

plot_dir <- project_path(cfg, "results/figures/auc_macrophage_signatures")
ensure_dir(plot_dir)
cell_col <- cfg$cell_types$cell_type_column

polar_summary <- obj@meta.data |>
  group_by(.data[[cell_col]]) |>
  summarise(Mean_Polar = mean(Polarization_Index, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(Mean_Polar))
p_polar <- ggplot(polar_summary, aes(x = reorder(.data[[cell_col]], Mean_Polar), y = Mean_Polar, fill = Mean_Polar > 0)) +
  geom_col(color = "white", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#BC3C29FF", "FALSE" = "#0072B5FF")) +
  labs(x = "Cell subtype", y = "Macrophage Polarization Index") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(plot_dir, "Macrophage_Polarization_Direction.pdf"), p_polar, width = 10, height = 6)

if (!all(is.na(obj$AMDI_Index))) {
  amdi_summary <- obj@meta.data |>
    group_by(.data[[cell_col]]) |>
    summarise(Mean_AMDI = mean(AMDI_Index, na.rm = TRUE), .groups = "drop") |>
    arrange(Mean_AMDI)
  p_amdi <- ggplot(amdi_summary, aes(x = reorder(.data[[cell_col]], Mean_AMDI), y = Mean_AMDI, fill = Mean_AMDI)) +
    geom_col(color = "black", width = 0.7) +
    coord_flip() +
    scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
    labs(x = "Cell subtype", y = "Macrophage Maturation Index") +
    theme_classic() +
    theme(legend.position = "none")
  ggsave(file.path(plot_dir, "Macrophage_Maturation_Index.pdf"), p_amdi, width = 10, height = 6)
}

ggsave(file.path(plot_dir, "FeaturePlot_Polarization_Index.pdf"), FeaturePlot(obj, features = "Polarization_Index"), width = 8, height = 6)
if (!all(is.na(obj$AMDI_Index))) {
  ggsave(file.path(plot_dir, "FeaturePlot_AMDI_Index.pdf"), FeaturePlot(obj, features = "AMDI_Index"), width = 8, height = 6)
}

message_step("Saved custom AUCell-scored object: ", out_rds)
