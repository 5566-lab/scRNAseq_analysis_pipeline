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

input_rds <- cfg$scrna$scored_rds
message_step("Loading scored monocyte/macrophage object: ", input_rds)
obj <- readRDS(input_rds)
out_dir <- project_path(cfg, cfg$scrna_pathway$output_dir)
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

msig <- msigdbr(species = cfg$project$species, category = cfg$scrna_pathway$msig_category)
geneset <- msig |> select(gs_name, entrez_gene)
gene_sets <- split(msig$entrez_gene, msig$gs_name)

single_cell_gsea <- GSEA(genelist, TERM2GENE = geneset, pvalueCutoff = cfg$scrna_pathway$gsea_pvalue_cutoff, pAdjustMethod = "BH", eps = 0, seed = cfg$project$seed)
safe_write_csv(as.data.frame(single_cell_gsea), file.path(out_dir, "gsea_results.csv"))

expr <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
expr_df <- data.frame(SYMBOL = rownames(expr), expr, check.names = FALSE) |>
  inner_join(gene_map, by = "SYMBOL") |>
  select(-SYMBOL) |>
  aggregate(. ~ ENTREZID, data = _, FUN = mean)
rownames(expr_df) <- expr_df$ENTREZID
expr_matrix <- as.matrix(expr_df[, setdiff(colnames(expr_df), "ENTREZID")])
params <- gsvaParam(exprData = expr_matrix, geneSets = gene_sets, kcdf = "Gaussian", absRanking = FALSE)
gsva_scores <- gsva(params, verbose = TRUE)
design <- model.matrix(~ obj$foam_state)
fit <- eBayes(lmFit(gsva_scores, design))
diff_pathways <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
safe_write_csv(diff_pathways, file.path(out_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)

volcano <- ggplot(diff_pathways, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < cfg$scrna_pathway$gsva_fdr_cutoff & abs(logFC) > cfg$scrna_pathway$gsva_logfc_cutoff), alpha = 0.7) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_hline(yintercept = -log10(cfg$scrna_pathway$gsva_fdr_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-cfg$scrna_pathway$gsva_logfc_cutoff, cfg$scrna_pathway$gsva_logfc_cutoff), linetype = "dashed") +
  labs(x = "Log2 fold change", y = "-log10 adjusted P")
ggsave(file.path(out_dir, "GSVA_volcano.pdf"), volcano, width = 8, height = 6)

message_step("Saved GSEA/GSVA outputs: ", out_dir)
