#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(GSVA)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
set.seed(cfg$project$seed)

calculate_cpm <- function(count_matrix) {
  lib_sizes <- colSums(count_matrix)
  log2(t(t(count_matrix) / lib_sizes) * 1e6 + 1)
}

prepare_count_matrix <- function(count_table) {
  count_data <- read.table(count_table, header = TRUE, row.names = 1, check.names = FALSE)
  count_data$geneID <- sub("\\..*", "", rownames(count_data))
  count_data <- aggregate(. ~ geneID, data = count_data, FUN = mean)
  rownames(count_data) <- count_data$geneID
  count_data$geneID <- NULL
  id_map <- bitr(rownames(count_data), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  id_map <- id_map[!duplicated(id_map$ENSEMBL), ]
  count_data$ENTREZID <- id_map$ENTREZID[match(rownames(count_data), id_map$ENSEMBL)]
  count_data <- na.omit(count_data)
  count_data <- aggregate(. ~ ENTREZID, data = count_data, FUN = mean)
  rownames(count_data) <- count_data$ENTREZID
  count_data$ENTREZID <- NULL
  calculate_cpm(as.matrix(count_data))
}

run_gsea_from_deg <- function(deg_table, geneset, output_dir) {
  diff_data <- read.csv(deg_table, header = TRUE, row.names = 1, check.names = FALSE)
  genelist <- diff_data$log2FoldChange
  names(genelist) <- sub("\\..*", "", rownames(diff_data))
  genelist <- sort(na.omit(genelist), decreasing = TRUE)
  entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(genelist), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  names(genelist) <- entrez
  genelist <- sort(na.omit(genelist), decreasing = TRUE)
  res <- GSEA(genelist, TERM2GENE = geneset, pvalueCutoff = cfg$pathway$gsea_pvalue_cutoff, pAdjustMethod = "BH", eps = 0, seed = cfg$project$seed)
  safe_write_csv(as.data.frame(res), file.path(output_dir, "gsea_results.csv"))
  res
}

run_gsva_from_counts <- function(count_table, gene_sets, output_dir) {
  expr <- prepare_count_matrix(count_table)
  params <- gsvaParam(exprData = expr, geneSets = gene_sets, kcdf = "Poisson", absRanking = FALSE)
  score <- gsva(params, verbose = TRUE)
  group <- ifelse(grepl("WT", colnames(score), ignore.case = TRUE), "WT", "Mutant")
  design <- model.matrix(~0 + group)
  colnames(design) <- make.names(colnames(design))
  contrast <- makeContrasts(Mutant_vs_WT = groupMutant - groupWT, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score, design), contrast))
  diff <- topTable(fit, number = Inf, adjust.method = "BH")
  safe_write_csv(diff, file.path(output_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)
  diff
}

input_rds <- project_path(cfg, cfg$outputs$pseudotime_rds)
message_step("Loading pseudotime object: ", input_rds)
obj <- readRDS(input_rds)
out_dir <- project_path(cfg, cfg$pathway$output_dir)
ensure_dir(out_dir)

if ("cdsMM_sub1" %in% colnames(obj@meta.data)) {
  obj$foam_state <- ifelse(obj$cdsMM_sub1 == "Yes", "Foam_Cell", "Macrophage")
} else {
  cell_col <- cfg$cell_types$cell_type_column
  obj$foam_state <- ifelse(grepl("Foam", obj@meta.data[[cell_col]], ignore.case = TRUE), "Foam_Cell", "Macrophage")
}

Idents(obj) <- "foam_state"
markers <- FindMarkers(obj, ident.1 = "Foam_Cell", ident.2 = "Macrophage", logfc.threshold = 0.5, min.pct = 0.1)
safe_write_csv(markers, file.path(out_dir, "foam_vs_macrophage_markers.csv"), row.names = TRUE)

genelist <- markers$avg_log2FC
names(genelist) <- rownames(markers)
genelist <- sort(na.omit(genelist), decreasing = TRUE)
gene_map <- bitr(names(genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genelist <- genelist[gene_map$SYMBOL]
names(genelist) <- gene_map$ENTREZID

msig <- msigdbr(species = cfg$project$species, category = cfg$pathway$msig_category)
geneset <- msig |> select(gs_name, entrez_gene)
gene_sets <- split(msig$entrez_gene, msig$gs_name)

single_cell_gsea <- GSEA(genelist, TERM2GENE = geneset, pvalueCutoff = cfg$pathway$gsea_pvalue_cutoff, pAdjustMethod = "BH", eps = 0, seed = cfg$project$seed)
safe_write_csv(as.data.frame(single_cell_gsea), file.path(out_dir, "gsea_results.csv"))

expr <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "data"))
expr_df <- data.frame(SYMBOL = rownames(expr), expr, check.names = FALSE) |>
  inner_join(gene_map, by = "SYMBOL") |>
  select(-SYMBOL) |>
  aggregate(. ~ ENTREZID, data = _, FUN = mean)
rownames(expr_df) <- expr_df$ENTREZID
expr_matrix <- as.matrix(expr_df[, setdiff(colnames(expr_df), "ENTREZID")])
params <- gsvaParam(exprData = expr_matrix, geneSets = gene_sets, kcdf = "Poisson", absRanking = FALSE)
gsva_scores <- gsva(params, verbose = TRUE)
design <- model.matrix(~ obj$foam_state)
fit <- eBayes(lmFit(gsva_scores, design))
diff_pathways <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
safe_write_csv(diff_pathways, file.path(out_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)

volcano <- ggplot(diff_pathways, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < cfg$pathway$gsva_fdr_cutoff & abs(logFC) > cfg$pathway$gsva_logfc_cutoff), alpha = 0.7) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_hline(yintercept = -log10(cfg$pathway$gsva_fdr_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-cfg$pathway$gsva_logfc_cutoff, cfg$pathway$gsva_logfc_cutoff), linetype = "dashed") +
  labs(x = "Log2 fold change", y = "-log10 adjusted P")
ggsave(file.path(out_dir, "GSVA_volcano.pdf"), volcano, width = 8, height = 6)

for (clone in names(cfg$apobec3a_ko)) {
  clone_dir <- file.path(out_dir, "APOBEC3A_KO", clone)
  ensure_dir(clone_dir)
  run_gsea_from_deg(cfg$apobec3a_ko[[clone]]$deg_table, geneset, clone_dir)
  run_gsva_from_counts(cfg$apobec3a_ko[[clone]]$count_table, gene_sets, clone_dir)
}

message_step("Saved GSEA/GSVA outputs: ", out_dir)

