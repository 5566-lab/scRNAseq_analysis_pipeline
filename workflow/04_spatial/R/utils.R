suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
})

severity_levels <- c("Mild", "Moderate", "Severe")
severity_cols <- c(Mild = "#4C78A8", Moderate = "#F58518", Severe = "#B6424B")

celltype_cols <- c(
  Myeloid = "#7A3E9D",
  SMCPericyte = "#4C78A8",
  ModSMC = "#72B7B2",
  Endothelium = "#54A24B",
  TCells = "#E45756",
  BCells = "#F2CF5B",
  Fibroblast1 = "#9C755F",
  Fibroblast2 = "#BAB0AC",
  Lymphatic = "#88CCEE",
  Mast = "#CC6677",
  Glia = "#117733",
  PlasmaCells = "#AA4499",
  Proliferating = "#332288"
)

signature_list <- list(
  APOBEC3A_axis = c("APOBEC3A", "ISG15", "MX1", "IFIT1", "IFIT3", "OAS1", "OAS2"),
  Foam_LAM = c("APOE", "APOC1", "TREM2", "LPL", "SPP1", "GPNMB", "PLA2G7", "LGALS3", "FABP5"),
  ISG_myeloid = c("ISG15", "MX1", "IFIT1", "IFI6", "IFI44L", "OAS1", "OAS2", "STAT1"),
  Inflammatory_myeloid = c("IL1B", "S100A8", "S100A9", "CXCL8", "CCL2", "CCL3", "CCL4", "NFKBIA"),
  Resident_LYVE1_TRM = c("LYVE1", "SELENOP", "MRC1", "FOLR2", "F13A1", "C1QA", "C1QB", "C1QC"),
  SMC_ModSMC = c("ACTA2", "TAGLN", "MYH11", "CNN1", "LGALS3", "SPP1", "FN1", "COL1A1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ACKR1", "CLDN5", "RAMP2"),
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R", "CCR7", "NKG7")
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30")
    )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_figure <- function(plot, filename, width = 7, height = 5, dpi = 320) {
  ensure_dir(dirname(filename))
  ggsave(filename, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  png_file <- sub("\\.pdf$", ".png", filename)
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(c(filename, png_file))
}

present_genes <- function(object, genes) {
  intersect(unique(genes), rownames(object))
}

available_signature_table <- function(object, dataset) {
  bind_rows(lapply(names(signature_list), function(sig) {
    genes <- signature_list[[sig]]
    tibble(
      dataset = dataset,
      signature = sig,
      n_requested = length(unique(genes)),
      n_available = length(present_genes(object, genes)),
      available_genes = paste(present_genes(object, genes), collapse = ";")
    )
  }))
}

collect_tissue_coordinates <- function(object) {
  bind_rows(lapply(Images(object), function(img) {
    cc <- GetTissueCoordinates(object, image = img)
    cc$cell <- if ("cell" %in% colnames(cc)) as.character(cc$cell) else rownames(cc)
    cc$image_name <- img
    as_tibble(cc)
  }))
}

fetch_data_layer <- function(object, genes, layer = "data") {
  tryCatch(
    FetchData(object, vars = genes, layer = layer),
    error = function(e) {
      FetchData(object, vars = genes)
    }
  )
}

zscore_signature <- function(object, genes, assay = NULL, layer = "data") {
  genes <- present_genes(object, genes)
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(object)))
  }
  if (!is.null(assay)) DefaultAssay(object) <- assay
  expr <- fetch_data_layer(object, genes, layer = layer)
  expr <- as.matrix(expr)
  if (ncol(expr) == 1) {
    score <- as.numeric(scale(expr[, 1]))
  } else {
    z <- scale(expr)
    z[is.nan(z)] <- 0
    score <- rowMeans(z, na.rm = TRUE)
  }
  names(score) <- colnames(object)
  score
}

add_signature_scores <- function(object, assay = NULL, prefix = "sig") {
  for (sig in names(signature_list)) {
    object[[paste0(prefix, "_", sig)]] <- zscore_signature(object, signature_list[[sig]], assay = assay)
  }
  object
}

composition_table <- function(meta, group_col, type_col) {
  meta %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[type_col]])) %>%
    count(.data[[group_col]], .data[[type_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup() %>%
    rename(group = all_of(group_col), cell_type = all_of(type_col))
}

plot_composition <- function(tab, title) {
  ggplot(tab, aes(x = group, y = frac, fill = cell_type)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.15) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = celltype_cols, na.value = "grey70") +
    labs(x = NULL, y = "Composition", fill = "Cell type", title = title) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

plot_score_by_severity <- function(meta, score_cols, severity_col, title_prefix) {
  plot_df <- meta %>%
    select(all_of(c(severity_col, score_cols))) %>%
    pivot_longer(all_of(score_cols), names_to = "signature", values_to = "score") %>%
    mutate(
      severity = factor(.data[[severity_col]], levels = severity_levels),
      signature = sub("^sig_", "", signature)
    ) %>%
    filter(!is.na(severity), !is.na(score))

  ggplot(plot_df, aes(x = severity, y = score, fill = severity)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.15, alpha = 0.82) +
    geom_boxplot(width = 0.13, outlier.size = 0.12, linewidth = 0.2, alpha = 0.85) +
    facet_wrap(~ signature, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = severity_cols) +
    labs(x = NULL, y = "Z-scored module expression", title = title_prefix) +
    theme_pub(9) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
}

kruskal_table <- function(meta, score_cols, severity_col) {
  bind_rows(lapply(score_cols, function(score_col) {
    dat <- meta %>%
      select(all_of(c(severity_col, score_col))) %>%
      filter(!is.na(.data[[severity_col]]), !is.na(.data[[score_col]]))
    p <- tryCatch(kruskal.test(dat[[score_col]] ~ dat[[severity_col]])$p.value, error = function(e) NA_real_)
    means <- dat %>%
      group_by(.data[[severity_col]]) %>%
      summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(level = as.character(.data[[severity_col]]))
    tibble(
      signature = sub("^sig_", "", score_col),
      p_kruskal = p,
      mild_mean = means$mean_score[match("Mild", means$level)],
      moderate_mean = means$mean_score[match("Moderate", means$level)],
      severe_mean = means$mean_score[match("Severe", means$level)]
    )
  })) %>%
    mutate(p_adj = p.adjust(p_kruskal, method = "BH"))
}
