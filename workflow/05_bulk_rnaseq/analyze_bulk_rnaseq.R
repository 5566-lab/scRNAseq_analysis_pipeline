#!/usr/bin/env Rscript

# APOBEC3A-knockout bulk RNA-seq analysis. The clone-specific and combined
# contrasts share one manifest and one implementation. Combined inference uses
# raw integer counts with clone included in the DESeq2 design.

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(optparse)
  library(pheatmap)
  library(tibble)
  library(yaml)
})

option_list <- list(
  make_option(c("-c", "--config"), default = "configs/config.yaml"),
  make_option(c("-x", "--contrast"), default = "all",
              help = "clone13, clone37, combined, or all"),
  make_option(c("--validate-only"), action = "store_true", default = FALSE,
              dest = "validate_only")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg_path <- normalizePath(opt$config, mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(cfg_path), ".."), mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)$bulk_rnaseq
set.seed(123)

resolve_path <- function(path) {
  if (grepl("^/", path)) path else file.path(repo_root, path)
}

read_featurecounts <- function(path) {
  tab <- read.delim(
    path, header = TRUE, comment.char = "#", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(tab) < 7L) stop("featureCounts table has no sample columns: ", path)
  genes <- sub("\\.[0-9]+$", "", tab[[1]])
  counts <- as.matrix(tab[, -(1:6), drop = FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- genes
  if (anyDuplicated(genes)) {
    counts <- rowsum(counts, group = genes, reorder = FALSE)
    storage.mode(counts) <- "integer"
  }
  colnames(counts) <- sub(
    "\\.sorted\\.bam$", "", basename(colnames(counts))
  )
  counts
}

read_manifest <- function(path) {
  manifest <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Sample_Name", "Condition", "Batch")
  missing <- setdiff(required, colnames(manifest))
  if (length(missing)) {
    stop("Sample manifest is missing: ", paste(missing, collapse = ", "))
  }
  manifest %>%
    transmute(
      sample = Sample_Name,
      condition = factor(Condition, levels = c("WT", "KO")),
      clone = factor(Batch, levels = c("clone13", "clone37"))
    ) %>%
    distinct(sample, .keep_all = TRUE)
}

select_contrast <- function(counts, manifest, contrast_name) {
  selected <- switch(
    contrast_name,
    clone13 = filter(manifest, clone == "clone13"),
    clone37 = filter(manifest, clone == "clone37"),
    combined = manifest,
    stop("Unknown contrast: ", contrast_name)
  )
  missing <- setdiff(selected$sample, colnames(counts))
  if (length(missing)) {
    stop("Count table is missing manifest samples: ", paste(missing, collapse = ", "))
  }
  list(
    counts = counts[, selected$sample, drop = FALSE],
    col_data = column_to_rownames(as.data.frame(selected), "sample")
  )
}

map_symbols <- function(ensembl_ids) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(ensembl_ids)
  }
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  ifelse(is.na(symbols) | symbols == "", ensembl_ids, unname(symbols))
}

save_pca <- function(dds, out_dir, contrast_name) {
  vsd <- vst(dds, blind = FALSE)
  pca <- plotPCA(vsd, intgroup = c("condition", "clone"), returnData = TRUE)
  percent_var <- round(100 * attr(pca, "percentVar"))
  p <- ggplot(pca, aes(PC1, PC2, color = condition, shape = clone, label = name)) +
    geom_point(size = 3.2) +
    geom_text_repel(size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(WT = "#0072B2", KO = "#D55E00")) +
    labs(
      title = paste("APOBEC3A knockout RNA-seq:", contrast_name),
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    theme_classic(base_size = 11)
  ggsave(file.path(out_dir, "PCA.pdf"), p, width = 7, height = 6)
  ggsave(file.path(out_dir, "PCA.png"), p, width = 7, height = 6, dpi = 300)
  invisible(vsd)
}

save_volcano <- function(result_table, out_dir, padj_cutoff, lfc_cutoff) {
  dat <- result_table %>%
    mutate(
      status = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "Up",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      ),
      neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin))
    )
  labels <- dat %>%
    filter(status != "Not significant") %>%
    arrange(padj) %>%
    slice_head(n = 12)
  p <- ggplot(dat, aes(log2FoldChange, neg_log10_padj, color = status)) +
    geom_point(alpha = 0.65, size = 1.2) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = 2) +
    geom_text_repel(data = labels, aes(label = SYMBOL), size = 3, max.overlaps = Inf) +
    scale_color_manual(values = c(
      Down = "#0072B2", `Not significant` = "grey75", Up = "#D55E00"
    )) +
    labs(x = "log2 fold change (KO / WT)", y = "-log10 adjusted P", color = NULL) +
    theme_classic(base_size = 11)
  ggsave(file.path(out_dir, "Volcano.pdf"), p, width = 8, height = 7)
  ggsave(file.path(out_dir, "Volcano.png"), p, width = 8, height = 7, dpi = 300)
}

save_heatmap <- function(dds, result_table, out_dir) {
  selected <- result_table %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 50)
  if (!nrow(selected)) return(invisible(NULL))
  mat <- assay(vst(dds, blind = FALSE))[selected$ENSEMBL, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  rownames(mat) <- make.unique(selected$SYMBOL)
  annotation <- as.data.frame(colData(dds)[, c("condition", "clone"), drop = FALSE])
  pheatmap(
    mat,
    annotation_col = annotation,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    fontsize_row = 6,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    filename = file.path(out_dir, "Top50_DEG_heatmap.pdf"),
    width = 8,
    height = 10
  )
}

run_enrichment <- function(result_table, out_dir, padj_cutoff, lfc_cutoff) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("clusterProfiler/org.Hs.eg.db unavailable; enrichment skipped")
    return(invisible(NULL))
  }
  for (direction in c("Up", "Down")) {
    genes <- result_table %>%
      filter(
        !is.na(padj), padj < padj_cutoff,
        if (direction == "Up") log2FoldChange >= lfc_cutoff else log2FoldChange <= -lfc_cutoff
      ) %>%
      pull(SYMBOL) %>%
      unique()
    ids <- suppressMessages(clusterProfiler::bitr(
      genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
    ))
    if (!nrow(ids)) next
    go <- suppressMessages(clusterProfiler::enrichGO(
      ids$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP",
      pAdjustMethod = "BH", readable = TRUE
    ))
    write.csv(as.data.frame(go), file.path(out_dir, paste0("GO_BP_", direction, ".csv")), row.names = FALSE)
  }

  if (!requireNamespace("msigdbr", quietly = TRUE)) return(invisible(NULL))
  pathways <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2") %>%
    select(gs_name, entrez_gene) %>%
    filter(!is.na(entrez_gene))
  ranked <- result_table %>%
    filter(!is.na(stat), !is.na(SYMBOL)) %>%
    arrange(desc(abs(stat))) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  mapped <- suppressMessages(clusterProfiler::bitr(
    ranked$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )) %>%
    inner_join(select(ranked, SYMBOL, stat), by = "SYMBOL") %>%
    arrange(desc(abs(stat))) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  gene_list <- mapped$stat
  names(gene_list) <- mapped$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  gsea <- suppressMessages(clusterProfiler::GSEA(
    gene_list, TERM2GENE = pathways, pvalueCutoff = 1,
    pAdjustMethod = "BH", seed = TRUE, verbose = FALSE
  ))
  write.csv(as.data.frame(gsea), file.path(out_dir, "GSEA_C2_results.csv"), row.names = FALSE)
}

run_gsva <- function(dds, out_dir, contrast_name) {
  if (!requireNamespace("GSVA", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE) ||
      !requireNamespace("limma", quietly = TRUE)) {
    message("GSVA/msigdbr/limma unavailable; GSVA skipped")
    return(invisible(NULL))
  }
  sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(gene_symbol)), .groups = "drop")
  gene_sets <- setNames(sets$genes, sets$gs_name)
  expr <- assay(vst(dds, blind = FALSE))
  rownames(expr) <- map_symbols(rownames(expr))
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  scores <- tryCatch(
    GSVA::gsva(expr, gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE),
    error = function(e) {
      param <- GSVA::gsvaParam(expr, gene_sets, kcdf = "Gaussian")
      GSVA::gsva(param, verbose = FALSE)
    }
  )
  col_data <- as.data.frame(colData(dds))
  design <- if (contrast_name == "combined") {
    model.matrix(~ clone + condition, data = col_data)
  } else {
    model.matrix(~ condition, data = col_data)
  }
  fit <- limma::eBayes(limma::lmFit(scores, design))
  coef_name <- "conditionKO"
  tab <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("pathway")
  write.csv(tab, file.path(out_dir, "GSVA_differential_pathways.csv"), row.names = FALSE)
}

run_contrast <- function(counts, manifest, contrast_name) {
  selected <- select_contrast(counts, manifest, contrast_name)
  design_formula <- if (contrast_name == "combined") ~ clone + condition else ~ condition
  dds <- DESeqDataSetFromMatrix(
    countData = selected$counts,
    colData = selected$col_data,
    design = design_formula
  )
  keep <- rowSums(counts(dds) >= cfg$min_count) >= cfg$min_samples
  dds <- DESeq(dds[keep, ], quiet = TRUE)
  res <- results(dds, contrast = c("condition", "KO", "WT"))
  tab <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    mutate(SYMBOL = map_symbols(ENSEMBL)) %>%
    select(ENSEMBL, SYMBOL, everything()) %>%
    arrange(padj)

  out_dir <- file.path(resolve_path(cfg$output_dir), contrast_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(tab, file.path(out_dir, "DESeq2_full_results.csv"), row.names = FALSE)
  write.csv(
    as.data.frame(counts(dds, normalized = TRUE)) %>% rownames_to_column("ENSEMBL"),
    file.path(out_dir, "normalized_counts.csv"), row.names = FALSE
  )
  write.csv(as.data.frame(colData(dds)) %>% rownames_to_column("sample"),
            file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)
  save_pca(dds, out_dir, contrast_name)
  save_volcano(tab, out_dir, cfg$padj_cutoff, cfg$abs_log2fc_cutoff)
  save_heatmap(dds, tab, out_dir)
  run_enrichment(tab, out_dir, cfg$padj_cutoff, cfg$abs_log2fc_cutoff)
  run_gsva(dds, out_dir, contrast_name)
  saveRDS(dds, file.path(out_dir, "dds.rds"))
  message("Completed bulk RNA-seq contrast: ", contrast_name)
}

counts_matrix <- read_featurecounts(resolve_path(cfg$count_table))
sample_manifest <- read_manifest(resolve_path(cfg$sample_manifest))
requested <- if (opt$contrast == "all") unlist(cfg$contrasts) else opt$contrast
if (opt$validate_only) {
  for (contrast_name in requested) {
    selected <- select_contrast(counts_matrix, sample_manifest, contrast_name)
    design <- if (contrast_name == "combined") {
      model.matrix(~ clone + condition, data = selected$col_data)
    } else {
      model.matrix(~ condition, data = selected$col_data)
    }
    if (qr(design)$rank != ncol(design)) stop("Design is not full rank: ", contrast_name)
    cat(sprintf(
      "PASS  %s: %d genes, %d samples, design=%s\n",
      contrast_name, nrow(selected$counts), ncol(selected$counts),
      paste(colnames(design), collapse = "+")
    ))
  }
  quit(status = 0, save = "no")
}
for (contrast_name in requested) {
  run_contrast(counts_matrix, sample_manifest, contrast_name)
}
