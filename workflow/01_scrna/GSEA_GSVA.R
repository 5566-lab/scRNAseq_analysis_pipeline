data<- readRDS("/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/data_pseudotime.rds")
subset_cells <- data
subset_cells@meta.data$cdsMM_sub1 <- ifelse(
  subset_cells@meta.data$cdsMM_sub1 == "Yes",
  "Foam_Cell",
  "Macrophage"
)
Idents(subset_cells) <- "cdsMM_sub1"



####代谢通路差异分析####
library(scMetabolism)
subset_cells[["RNA"]] <- as(subset_cells[["RNA"]], "Assay")
subset_cells <- sc.metabolism.Seurat(
  obj = subset_cells,
  method = "VISION",       # 可选 AUCell/ssGSEA/GSVA
  imputation = F,          # 关闭数据插补（加快速度）
  ncores = 4,              # 并行计算核数
  metabolism.type = "KEGG" # 使用 KEGG
)

# 提取代谢评分矩阵
metabolism_scores <- subset_cells@assays$METABOLISM$score

df <- subset_cells@meta.data

avg_df = aggregate(t(metabolism_scores),
                   list(df$Celltype_raw1),
                   mean)
rownames(avg_df) <- avg_df$Group.1
avg_df <- subset(avg_df, select = -Group.1)
rownames(avg_df)[rownames(avg_df) == "FOAM_cells1"] <- "Foam cells1"
rownames(avg_df)[rownames(avg_df) == "FOAM_cells2"] <- "Foam cells2"
custom_row_order <- c("Foam cells2", "Foam cells1", "CM", "cDC1", "Monocyte","Macrophage", "TrMs")
avg_df <- avg_df[custom_row_order, , drop=FALSE]

# 计算通路方差并筛选
pathway_variance <- apply(avg_df, 2, var)
top20_pathways <- names(sort(pathway_variance, decreasing = TRUE)[1:10])
avg_df_top20 <- avg_df[, top20_pathways, drop=FALSE]
avg_df_top20 <- avg_df_top20[custom_row_order, , drop=FALSE]

# 生成热图
pdf("/public3/DSC/single_cell/GSE159677/GSEA/metabolism_heatmap.pdf", width=12, height=8)
pheatmap(t(avg_df_top20),
         show_colnames = T,
         scale = 'row',
         cluster_rows = T,
         color = colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F,
         main = "",#Top 20 Differential Metabolic Pathways by Cell Type
         fontsize_row = 12,
         fontsize_col = 14,
         fontsize = 12,
         angle_col = 45) +
        theme(
          text = element_text(family = "Arial"))
dev.off()
saveRDS(subset_cells, file = "/public3/DSC/single_cell/GSE159677/GSEA/metabolism.rds")


top_n <- 5

# 初始化存储所有显著通路的列表
all_top_pathways <- list()

# 遍历每个细胞类型
for (cell_type in rownames(avg_df)) {
  # 获取当前细胞类型的通路活性值（转换为数值向量）
  cell_data <- as.numeric(avg_df[cell_type, ])
  names(cell_data) <- colnames(avg_df)  # 保持通路名称

  # 按活性值降序排列并选择前5（处理可能的NA值）
  sorted_pathways <- names(sort(cell_data, decreasing = TRUE, na.last = TRUE))[1:top_n]

  # 存储结果
  all_top_pathways[[cell_type]] <- sorted_pathways
}

# 合并所有通路并去重
selected_pathways <- unique(unlist(all_top_pathways))
# 提取目标通路数据
avg_df_selected <- avg_df[, selected_pathways, drop = FALSE]

# 按细胞类型顺序排列
avg_df_selected <- avg_df_selected[custom_row_order, , drop = FALSE]

pdf("/public3/DSC/single_cell/GSE159677/GSEA/Selected_Metabolic_heatmap.pdf", width=15, height=15)
pheatmap(t(avg_df),
         show_colnames = T,
         scale = 'row',
         cluster_rows = T,
         color = colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = T,
         main = "",#"Selected Metabolic Pathways Activity by Cell Type"
         fontsize_row = 12,
         fontsize_col = 14,
         fontsize = 12,
         angle_col = 45) +
  theme(
    text = element_text(family = "Arial"))
dev.off()





####GSEA和GSVA分析####

markers <- FindMarkers(
  object = subset_cells,
  ident.1 = "Foam_Cell",
  ident.2 = "Macrophage",
  group.by = "cdsMM_sub1",
  logfc.threshold = 0.5,     # 过滤低变化基因
  only.pos = FALSE,          # 包含上下调基因
  min.pct = 0.1,             # 基因至少表达于10%细胞
  return.thresh = 0.01       # 返回p_val_adj<0.01的基因
)
###GSEA###
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

output_dir <- "/public3/DSC/single_cell/GSE159677/GSEA"
# 转换基因ID为ENTREZID
genelist <- markers$avg_log2FC
names(genelist) <- rownames(markers)
genelist <- sort(genelist, decreasing = TRUE)

gene_map <- bitr(names(genelist),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
genelist <- genelist[gene_map$SYMBOL]
names(genelist) <- gene_map$ENTREZID

geneset <- msigdbr(species = "Homo sapiens",
                        category = "C2",
                        ) %>%
  dplyr::select(gs_name, entrez_gene)

gsea_res <- GSEA(genelist,
                 TERM2GENE = geneset,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 eps = 0,
                 seed = 123 )
write.csv(as.data.frame(gsea_res),
          file = file.path(output_dir, "gsea_results.csv"),
          row.names = FALSE)
library(enrichplot)
gseaplot2(gsea_res,
          geneSetID = 1:10,  # 选择前10条通路
          base_size = 12,
          title = "Top Enriched Pathways in Foam Cells")
#COATES_MACROPHAGE_M1_VS_M2_UP → 科茨巨噬细胞M1对比M2上调基因 泡沫细胞偏向M2
#WUNDER_INFLAMMATORY_RESPONSE_AND_CHOLESTEROL_UP → 温德炎症反应与胆固醇上调基因 泡沫细胞下调
#WINTER_HYPOXIA_METAGENE → 温特缺氧代谢基因特征 泡沫细胞上调
#HINATA_NFKB_TARGETS_FIBROBLAST_UP → 日向NF-κB靶标成纤维细胞上调基因 #泡沫细胞上调
#LENAOUR_DENDRITIC_CELL_MATURATION_UP → 勒努尔树突细胞成熟上调基因 泡沫细胞下调
#REACTOME_COMPLEMENT_CASCADE → Reactome 补体级联反应
#KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY_MEDIATED_TRANSCRIPTION → KEGG 病原体HTLV-1 Tax蛋白通过NFY介导转录
#WP_GLYCOLYSIS_AND_GLUCONEOGENESIS → 糖酵解与糖异生通路
#REACTOME_PD_1_SIGNALING → Reactome PD-1信号通路 #泡沫细胞下调
#MA_RAT_AGING_UP → 马氏大鼠衰老上调基因 泡沫细胞上调
#GROSS_HIF1A_TARGETS_DN → 格罗斯HIF1α靶标下调基因 泡沫细胞上调
#REACTOME_INNATE_IMMUNE_SYSTEM → Reactome 先天免疫系统
#KEGG_COMPLEMENT_AND_COAGULATION_CASCADES 补体与凝血级联通路 泡沫细胞下调
#WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1 → 王氏SIRT1抑制的成脂基因 泡沫细胞上调
#LEE_BMP2_TARGETS_UP → 李氏BMP2靶标上调基因
#MOOTHA_GLYCOLYSIS → 穆塔糖酵解基因集 泡沫细胞上调
#CHEBOTAEV_GR_TARGETS_DN → 切博塔耶夫GR靶标下调基因
#CADWELL_ATG16L1_TARGETS_UP → 卡德韦尔ATG16L1靶标上调基因
#REACTOME_INTERFERON_GAMMA_SIGNALING → Reactome γ干扰素信号通路 泡沫细胞下调
#LEE_TARGETS_OF_PTCH1_AND_SUFU_UP → 李氏PTCH1与SUFU靶标上调基因
#WP_ZINC_HOMEOSTASIS → 锌稳态通路 泡沫细胞上调
#BROWN_MYELOID_CELL_DEVELOPMENT_UP → 布朗髓系细胞发育上调基因 泡沫细胞上调
#ONDER_CDH1_TARGETS_2_DN → 昂德尔CDH1靶标2下调基因
#SENESE_HDAC1_AND_HDAC2_TARGETS_UP → 塞内塞HDAC1与HDAC2靶标上调基因
#WESTON_VEGFA_TARGETS → 韦斯顿VEGFA靶标基因集
#HALMOS_CEBPA_TARGETS_UP → 哈尔莫斯CEBPA靶标上调基因
#WONG_ADULT_TISSUE_STEM_MODULE → 黄氏成体组织干细胞模块 #泡沫细胞下调
#SERVITJA_ISLET_HNF1A_TARGETS_UP → 塞尔维察胰岛HNF1α靶标上调基因
#RUTELLA_RESPONSE_TO_CSF2RB_AND_IL4_UP → 鲁特拉响应CSF2RB与IL4上调基因
#BENPORATH_PRC2_TARGETS → 本波拉思PRC2靶标基因
#YOSHIMURA_MAPK8_TARGETS_UP → 吉村MAPK8靶标上调基因
#WP_PPAR_SIGNALING → PPAR信号通路 泡沫细胞上调
#WP_OXIDATIVE_DAMAGE_RESPONSE → 氧化损伤响应通路 #泡沫细胞下调
#IGLESIAS_E2F_TARGETS_UP → 伊格莱西亚斯E2F靶标上调基因 #泡沫细胞下调
#REACTOME_DAP12_INTERACTIONS → Reactome DAP12相互作用
#WP_COPPER_HOMEOSTASIS → 铜稳态通路
#REACTOME_G_ALPHA_I_SIGNALLING_EVENTS → Reactome Gαi信号事件 泡沫细胞下调
#MARZEC_IL2_SIGNALING_DN → 马尔泽克IL-2信号下调基因 泡沫细胞下调
#HECKER_IFNB1_TARGETS → 赫克IFNβ1靶标基因 泡沫细胞下调

target_pathways <- c("COATES_MACROPHAGE_M1_VS_M2_UP", "WP_PPAR_SIGNALING")

pdf(file = file.path(output_dir, "gsea_plot.pdf"), width = 10, height = 8)
gseaplot2(gsea_res,
          geneSetID = target_pathways,  # 指定目标通路
          base_size = 12,
          title = "GSEA of Macrophage Polarization and PPAR Signaling")
dev.off()

library(GseaVis)
gseaNb(gsea_res, geneSetID = "COATES_MACROPHAGE_M1_VS_M2_UP",subPlot = 3,
       addPval = T,
       pvalX = 0.85,pvalY = 0.75,
       nesDigit = 4,
       pDigit = 4)


dotplotGsea(data = gsea_res,
            topn= 10,
            str.width = 20 # 折叠通路名，我这里不能用...
)

dotplotGsea(data = gsea_res,
            topn= 10,
            order.by = "NES",
            add.seg = T,
            line.col = 'orange',
            line.type = 1
)


volcanoGsea(data = gsea_res,
            nudge.y = c(-0.8,0.8)
)

####GSVA###
library(GSVA)
library(GSEABase)
library(limma)
expr_matrix <- as.matrix(subset_cells@assays$RNA$data)

# 2.2 转换基因ID为ENTREZID（与GSEA保持一致）

expr_df <- data.frame(SYMBOL = rownames(expr_matrix), expr_matrix)
expr_df <- inner_join(expr_df, gene_map, by = "SYMBOL") %>%
  dplyr::select(-SYMBOL) %>%
  aggregate(. ~ ENTREZID, data = ., FUN = mean)

rownames(expr_df) <- expr_df$ENTREZID
expr_matrix <- as.matrix(expr_df[, -1])

kegg_df <- msigdbr(species = "Homo sapiens", category = "C2")
kegg_list <- split(kegg_df$entrez_gene, kegg_df$gs_name)

# 创建参数对象（新版GSVA语法）
params <- gsvaParam(
  exprData = expr_matrix,
                    geneSets = kegg_list,
                    kcdf = "Poisson",
                    absRanking = FALSE)
#执行GSVA计算
gsva_scores <- gsva(params, verbose = TRUE)

#差异通路分析

design <- model.matrix(~ subset_cells$cdsMM_sub1)
fit <- lmFit(gsva_scores, design)
fit <- eBayes(fit)
diff_pathways <- topTable(fit,
                          coef = 2,
                          number = Inf,
                          adjust.method = "BH")

sig_pathways <- diff_pathways %>%
  dplyr::filter(abs(logFC) > 0.5 & adj.P.Val < 0.01)
write.csv(diff_pathways, file = file.path(output_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)
diff_pathways <- read.csv(file.path(output_dir, "GSVA_diff_pathways_all.csv"))
#WP_HIF1A_AND_PPARG_REGULATION_OF_GLYCOLYSIS
#WikiPathways通路：HIF1A和PPARγ对糖酵解的调控 泡沫细胞上调
#WP_GLYCOLYSIS_IN_SENESCENCE
#WikiPathways通路：衰老中的糖酵解  泡沫细胞上调
#GOERING_BLOOD_HDL_CHOLESTEROL_QTL_CIS
#Goering研究：血液HDL胆固醇的顺式数量性状位点（QTL）泡沫细胞下调
#SANDERSON_PPARA_TARGETS
#Sanderson研究：PPARA靶基因 泡沫细胞上调
#WP_FATTY_ACID_TRANSPORTERS
#WikiPathways通路：脂肪酸转运蛋白 泡沫细胞下调
#WP_CCL18_SIGNALING
#WikiPathways通路：CCL18信号 上调
#SA_FAS_SIGNALING
#SA研究：FAS信号 下调
#REACTOME_SCAVENGING_BY_CLASS_B_RECEPTORS
#Reactome通路：B类受体的清除作用 上调
#KONDO_HYPOXIA
#Kondo研究：缺氧反应  上调
#REACTOME_HDL_REMODELING
#高密度脂蛋白（HDL）重塑 下调
#REACTOME_ENOS_ACTIVATION
#eNOS（内皮型一氧化氮合酶）激活 下调
#REACTOME_METABOLISM_OF_NITRIC_OXIDE_NOS3_ACTIVATION_AND_REGULATION
#一氧化氮代谢（NOS3激活与调控） 下调
#RAMJAUN_APOPTOSIS_BY_TGFB1_VIA_MAPK1_UP
#Ramjaun研究：TGFB1通过MAPK1诱导凋亡（上调） 下调

volcano_plot <- ggplot(diff_pathways,
                       aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.01 & abs(logFC) > 0.5,
                                ifelse(logFC > 0, "Up", "Down"), "NS")),
             alpha = 0.7) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = subset(diff_pathways, adj.P.Val < 0.01 & abs(logFC) > 0.5),
    aes(label = rownames(subset(diff_pathways, adj.P.Val < 0.01 & abs(logFC) > 0.5))),  # 关键修改
    size = 3,
    max.overlaps = 20) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adj.Pvalue)")
ggsave(file.path(output_dir,"GSVA_volcano.pdf"), plot = volcano_plot, width = 8, height = 6)

library(ggplot2)
sig_pathways$Direction <- ifelse(sig_pathways$logFC > 0, "Up", "Down")
top_pathways <- sig_pathways %>%
  tibble::rownames_to_column("Pathway") %>%
  group_by(Direction) %>%
  arrange(desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  tibble::column_to_rownames("Pathway")

ggplot(top_pathways, aes(x = reorder(rownames(top_pathways), logFC),
                         y = logFC, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#B83D3D", "#1A5592")) +
  coord_flip() +
  labs(x = "Pathway", y = "Log2 Fold Change") +
  theme_minimal()


####APOBEC3A_KO####
output_dir <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO"

diff_data <- read.csv(
  "/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/ko/3A_total.csv",
  header = TRUE,
  row.names = 1         # 第一列为基因ID（如ENSG00000183853.18_12）
)
rownames(diff_data) <- sub("\\..*", "", rownames(diff_data))
count_filtered <- count_data[rowSums(diff_data > 1) >= 3, ]



genelist <- diff_data$log2FoldChange
names(genelist) <- rownames(diff_data)
genelist <- na.omit(genelist)
genelist <- sort(genelist, decreasing = TRUE)

#去除版本号
names(genelist) <- sub("\\..*", "", names(genelist))

genelist_entrez <- mapIds(
  org.Hs.eg.db,
  keys = names(genelist),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
names(genelist) <- genelist_entrez
genelist <- na.omit(genelist)

gsea_res <- GSEA(
  geneList = genelist,
  TERM2GENE = geneset,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 0,
  seed = 123  # 保证可重复性
)
write.csv(as.data.frame(gsea_res),
          file = file.path(output_dir, "gsea_results.csv"),
          row.names = FALSE)

# 读取两个GSEA结果文件
file1 <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO/gsea_results.csv"
file2 <- "/public3/DSC/single_cell/GSE159677/GSEA/gsea_results.csv"

gsea_res1 <- read.csv(file1, check.names = FALSE)  # 文件1: A3A_KO结果
gsea_res2 <- read.csv(file2, check.names = FALSE)  #文件2: 单核细胞分化泡沫细胞

common_pathways <- intersect(gsea_res1$Description, gsea_res2$Description)

# 提取各文件的NES和p.adjust
result <- bind_rows(
  gsea_res1 %>%
    filter(Description %in% common_pathways) %>%
    select(Description, "NES_A3A_KO" = NES, "p.adjust_A3A_KO" = p.adjust),
  gsea_res2 %>%
    filter(Description %in% common_pathways) %>%
    select(Description, "NES_General" = NES, "p.adjust_General" = p.adjust)
) %>%
  group_by(Description) %>%
  summarise(
    NES_A3A_KO = first(na.omit(NES_A3A_KO)),  # 若存在多行取第一个
    p.adjust_A3A_KO = first(na.omit(p.adjust_A3A_KO)),
    NES_General = first(na.omit(NES_General)),
    p.adjust_General = first(na.omit(p.adjust_General))
  )
#APOBEC3A敲除后富集通路应与泡沫细胞分化通路相反 clone13 #第一次测序
#DAZARD_RESPONSE_TO_UV_NHEK_UP
#FULCHER_INFLAMMATORY_RESPONSE_LECTIN_VS_LPS_UP
#GHANDHI_BYSTANDER_IRRADIATION_UP
#HINATA_NFKB_TARGETS_FIBROBLAST_UP
#ONDER_CDH1_TARGETS_2_DN
#REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX
#REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION
#SENESE_HDAC1_AND_HDAC2_TARGETS_UP
#VERHAAK_AML_WITH_NPM1_MUTATED_UP




count_data <- read.table("/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/count.txt",
                         header = TRUE,
                         row.names = 1)

# 新建geneID列（值为当前行名）
count_data$geneID <- rownames(count_data)

# 对geneID列去版本号（保留点号前的部分）
count_data$geneID <- sub("\\..*", "", count_data$geneID)
count_data<- aggregate(. ~ geneID, data = count_data, FUN = mean)
# 将处理后的geneID赋给行名
rownames(count_data) <- count_data$geneID

# 删除临时列geneID
count_data$geneID <- NULL

ensembl_ids <- rownames(count_data)

id_map <- bitr(
  ensembl_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

id_map <- id_map[!duplicated(id_map$ENSEMBL), ]

# - 合并到表达矩阵
count_data$ENTREZID <- id_map$ENTREZID[match(rownames(count_data), id_map$ENSEMBL)]
count_data <- na.omit(count_data)  # 删除未匹配的基因
count_data<- aggregate(. ~ ENTREZID, data = count_data, FUN = mean)

# - 将ENTREZID设为行名
rownames(count_data) <- count_data$ENTREZID
count_data$ENTREZID <- NULL

calculate_cpm <- function(count_matrix) {
  # 1. 计算每列（样本）的总reads数
  lib_sizes <- colSums(count_matrix)

  # 2. 计算CPM = (count / lib_size) * 1e6
  cpm_matrix <- t(t(count_matrix) / lib_sizes) * 1e6

  # 3. 对结果取log2(CPM + 1)（避免log(0)）
  log2_cpm <- log2(cpm_matrix + 1)
  return(log2_cpm)
}

# 应用标准化（假设count_data已去除版本号并合并重复基因）
count_data_normalized <- calculate_cpm(count_data)

params <- gsvaParam(
  exprData = as.matrix(count_data_normalized),
  geneSets = kegg_list,
  kcdf = "Poisson",
  absRanking = FALSE)

gsva_result <- gsva(params, verbose = TRUE)

group <- ifelse(grepl("WT", colnames(gsva_result)), "WT", "Mutant")
design <- model.matrix(~0 + group)
colnames(design) <- c("Mutant", "WT")

# 构建对比矩阵（Mutant vs WT）
contrast_matrix <- makeContrasts(Mutant_vs_WT = Mutant - WT, levels = design)

# 线性模型拟合
fit <- lmFit(gsva_result, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取差异通路结果
diff_pathways <- topTable(fit2, number = Inf, adjust.method = "BH")
sig_pathways <- subset(diff_pathways, adj.P.Val < 0.01 & abs(logFC) > 0.5)

write.csv(sig_pathways, file = file.path(output_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)

# 读取两个GSVA结果文件
file1 <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO/GSVA_diff_pathways_all.csv"
file2 <- "/public3/DSC/single_cell/GSE159677/GSEA/GSVA_diff_pathways_all.csv"

gsva_res1 <- read.csv(file1, check.names = FALSE)  # 文件1: A3A_KO结果
gsva_res2 <- read.csv(file2, check.names = FALSE)  #文件2: 单核细胞分化泡沫细胞

common_pathways <- intersect(gsva_res1[,1], gsva_res2[,1])

colnames(gsva_res1)[1] <- "Pathway"
colnames(gsva_res2)[1] <- "Pathway"

# 提取各文件的logFC和p.adjust
result <- bind_rows(
  gsva_res1 %>%
    filter(.[[1]] %in% common_pathways) %>%
    select(1 , "logFC_A3A" = logFC, "p.adjust_A3A_KO" = adj.P.Val),
  gsva_res2 %>%
    filter(.[[1]]  %in% common_pathways) %>%
    select(1 , "logFC_General" = logFC, "p.adjust_General" = adj.P.Val)
) %>%
  group_by(Pathway) %>%
  summarise(
    logFC_A3A = first(na.omit(logFC_A3A)),  # 若存在多行取第一个
    p.adjust_A3A_KO = first(na.omit(p.adjust_A3A_KO)),
    logFC_General = first(na.omit(logFC_General)),
    p.adjust_General = first(na.omit(p.adjust_General))
  )


filtered_result <- result %>%
  filter(
    abs(logFC_A3A) > 0.5 &
      abs(logFC_General) > 0.5 &
      sign(logFC_A3A) != sign(logFC_General)  # 符号相反
  )


#BIOCARTA_EOSINOPHILS_PATHWAY
#GAURNIER_PSMD4_TARGETS
#GOUYER_TATI_TARGETS_DN
#VETTER_TARGETS_OF_PRKCA_AND_ETS1_DN
#WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE





####APOBEC3A#### clone37 第二次测序
library(enrichplot)

output_dir <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO/clone37"

diff_data <- read.csv(
  "/public3/DSC/single_cell/DEG_ALL.csv",
  header = TRUE,
  row.names = 1         # 第一列为基因ID（如ENSG00000183853.18_12）
)
rownames(diff_data) <- sub("\\..*", "", rownames(diff_data))
count_filtered <- count_data[rowSums(diff_data > 1) >= 3, ]



genelist <- diff_data$log2FoldChange
names(genelist) <- rownames(diff_data)
genelist <- na.omit(genelist)
genelist <- sort(genelist, decreasing = TRUE)

#去除版本号
names(genelist) <- sub("\\..*", "", names(genelist))

genelist_entrez <- mapIds(
  org.Hs.eg.db,
  keys = names(genelist),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
names(genelist) <- genelist_entrez
genelist <- na.omit(genelist)

gsea_res <- GSEA(
  geneList = genelist,
  TERM2GENE = geneset,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 0,
  seed = 123  # 保证可重复性
)
write.csv(as.data.frame(gsea_res),
          file = file.path(output_dir, "gsea_results.csv"),
          row.names = FALSE)

# 读取两个GSEA结果文件
file1 <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO/clone37/gsea_results.csv"
file2 <- "/public3/DSC/single_cell/GSE159677/GSEA/gsea_results.csv"

gsea_res1 <- read.csv(file1, check.names = FALSE)  # 文件1: A3A_KO结果
gsea_res2 <- read.csv(file2, check.names = FALSE)  #文件2: 单核细胞分化泡沫细胞

common_pathways <- intersect(gsea_res1$Description, gsea_res2$Description)

# 提取各文件的NES和p.adjust
result <- bind_rows(
  gsea_res1 %>%
    filter(Description %in% common_pathways) %>%
    select(Description, "NES_A3A_KO" = NES, "p.adjust_A3A_KO" = p.adjust),
  gsea_res2 %>%
    filter(Description %in% common_pathways) %>%
    select(Description, "NES_General" = NES, "p.adjust_General" = p.adjust)
) %>%
  group_by(Description) %>%
  summarise(
    NES_A3A_KO = first(na.omit(NES_A3A_KO)),  # 若存在多行取第一个
    p.adjust_A3A_KO = first(na.omit(p.adjust_A3A_KO)),
    NES_General = first(na.omit(NES_General)),
    p.adjust_General = first(na.omit(p.adjust_General))
  )
#APOBEC3A敲除后富集通路应与泡沫细胞分化通路相反 clone37
#KEGG系统性红斑狼疮
#RUTELLA_HGF响应(对比CSF2RB/IL4)上调
#WU细胞迁移
#YAMASHITA前列腺癌甲基化基因
#REACTOME细胞外基质组织
#ONDER_CDH1靶标(2型)下调
#WESTON_VEGFA靶标
#WINTER低氧代谢特征
#GHANDHI旁观者辐射效应上调
#VERHAAK_NPM1突变型AML上调
#KAN三氧化二砷响应
#SENGUPTA鼻咽癌(LMP1相关)下调
#DELYS甲状腺癌上调
#OISHI胆管瘤干细胞样下调
#NABA核心基质组
#SENGUPTA鼻咽癌下调
#REACTOME细胞外基质降解
#CROMER肿瘤发生上调
#LEE肝癌(DENA诱导)上调



count_data <- read.table("/dsk2/who/panxy/RNAediting/CAD/count.txt",
                         header = TRUE,
                         row.names = 1)

# 新建geneID列（值为当前行名）
count_data$geneID <- rownames(count_data)

# 对geneID列去版本号（保留点号前的部分）
count_data$geneID <- sub("\\..*", "", count_data$geneID)
count_data<- aggregate(. ~ geneID, data = count_data, FUN = mean)
# 将处理后的geneID赋给行名
rownames(count_data) <- count_data$geneID

# 删除临时列geneID
count_data$geneID <- NULL

ensembl_ids <- rownames(count_data)

id_map <- bitr(
  ensembl_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

id_map <- id_map[!duplicated(id_map$ENSEMBL), ]

# - 合并到表达矩阵
count_data$ENTREZID <- id_map$ENTREZID[match(rownames(count_data), id_map$ENSEMBL)]
count_data <- na.omit(count_data)  # 删除未匹配的基因
count_data<- aggregate(. ~ ENTREZID, data = count_data, FUN = mean)

# - 将ENTREZID设为行名
rownames(count_data) <- count_data$ENTREZID
count_data$ENTREZID <- NULL

calculate_cpm <- function(count_matrix) {
  # 1. 计算每列（样本）的总reads数
  lib_sizes <- colSums(count_matrix)

  # 2. 计算CPM = (count / lib_size) * 1e6
  cpm_matrix <- t(t(count_matrix) / lib_sizes) * 1e6

  # 3. 对结果取log2(CPM + 1)（避免log(0)）
  log2_cpm <- log2(cpm_matrix + 1)
  return(log2_cpm)
}

# 应用标准化（假设count_data已去除版本号并合并重复基因）
count_data_normalized <- calculate_cpm(count_data)

params <- gsvaParam(
  exprData = as.matrix(count_data_normalized),
  geneSets = kegg_list,
  kcdf = "Poisson",
  absRanking = FALSE)

gsva_result <- gsva(params, verbose = TRUE)

group <- ifelse(grepl("WT", colnames(gsva_result)), "WT", "Mutant")
design <- model.matrix(~0 + group)
colnames(design) <- c("Mutant", "WT")

# 构建对比矩阵（Mutant vs WT）
contrast_matrix <- makeContrasts(Mutant_vs_WT = Mutant - WT, levels = design)

# 线性模型拟合
fit <- lmFit(gsva_result, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取差异通路结果
diff_pathways <- topTable(fit2, number = Inf, adjust.method = "BH")
sig_pathways <- subset(diff_pathways, adj.P.Val < 0.01 & abs(logFC) > 0.5)

write.csv(sig_pathways, file = file.path(output_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)

# 读取两个GSVA结果文件
file1 <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO/clone13/GSVA_diff_pathways_all.csv"
file2 <- "/public3/DSC/single_cell/GSE159677/GSEA/GSVA_diff_pathways_all.csv"

gsva_res1 <- read.csv(file1, check.names = FALSE)  # 文件1: A3A_KO结果
gsva_res2 <- read.csv(file2, check.names = FALSE)  #文件2: 单核细胞分化泡沫细胞

common_pathways <- intersect(gsva_res1[,1], gsva_res2[,1])

colnames(gsva_res1)[1] <- "Pathway"
colnames(gsva_res2)[1] <- "Pathway"

# 提取各文件的logFC和p.adjust
result <- bind_rows(
  gsva_res1 %>%
    filter(.[[1]] %in% common_pathways) %>%
    select(1 , "logFC_A3A" = logFC, "p.adjust_A3A_KO" = adj.P.Val),
  gsva_res2 %>%
    filter(.[[1]]  %in% common_pathways) %>%
    select(1 , "logFC_General" = logFC, "p.adjust_General" = adj.P.Val)
) %>%
  group_by(Pathway) %>%
  summarise(
    logFC_A3A = first(na.omit(logFC_A3A)),  # 若存在多行取第一个
    p.adjust_A3A_KO = first(na.omit(p.adjust_A3A_KO)),
    logFC_General = first(na.omit(logFC_General)),
    p.adjust_General = first(na.omit(p.adjust_General))
  )


filtered_result <- result %>%
  filter(
    abs(logFC_A3A) > 0.5 &         # logFC_A3A绝对值>1
      abs(logFC_General) > 0.5 &        # logFC_General绝对值>1
      sign(logFC_A3A) != sign(logFC_General)  # 符号相反
  )


#BIOCARTA_EOSINOPHILS_PATHWAY
#FUJIWARA_PARK2_IN_LIVER_CANCER_UP
#HEBERT_MATRISOME_TNBC_BONE_METASTASIS_TUMOR_CELL_DERIVED
#KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY_MEDIATED_TRANSCRIPTION
#KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES
#KEGG_MEDICUS_REFERENCE_LYSOSOMAL_CA2_RELEASE
#REACTOME_COMMON_PATHWAY_OF_FIBRIN_CLOT_FORMATION
#REACTOME_SCAVENGING_BY_CLASS_B_RECEPTORS
#UNTERMAN_IPF_VS_CTRL_NK_CELL_UP
#UNTERMAN_PROGRESSIVE_VS_STABLE_IPF_NK_CELL_DN
#VETTER_TARGETS_OF_PRKCA_AND_ETS1_DN






#####APOBEC3A#####13+37
count_data1 <- read.table("/dsk2/who/panxy/RNAediting/CAD/count.txt",
                         header = TRUE,
                         row.names = 1)
count_data2 <- read.table("/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/count.txt",
                         header = TRUE,
                         row.names = 1)



# 检查行名是否一致
all(rownames(count_data1) == rownames(count_data2))  # 应该返回 TRUE

# 合并列（适用于不同样本）
count_data <- cbind(count_data1, count_data2)

# 新建geneID列（值为当前行名）
count_data$geneID <- rownames(count_data)

# 对geneID列去版本号（保留点号前的部分）
count_data$geneID <- sub("\\..*", "", count_data$geneID)
count_data<- aggregate(. ~ geneID, data = count_data, FUN = mean)
# 将处理后的geneID赋给行名
rownames(count_data) <- count_data$geneID

# 删除临时列geneID
count_data$geneID <- NULL

ensembl_ids <- rownames(count_data)

id_map <- bitr(
  ensembl_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

id_map <- id_map[!duplicated(id_map$ENSEMBL), ]

# - 合并到表达矩阵
count_data$ENTREZID <- id_map$ENTREZID[match(rownames(count_data), id_map$ENSEMBL)]
count_data <- na.omit(count_data)  # 删除未匹配的基因
count_data<- aggregate(. ~ ENTREZID, data = count_data, FUN = mean)

# - 将ENTREZID设为行名
rownames(count_data) <- count_data$ENTREZID
count_data$ENTREZID <- NULL

calculate_cpm <- function(count_matrix) {
  # 1. 计算每列（样本）的总reads数
  lib_sizes <- colSums(count_matrix)

  # 2. 计算CPM = (count / lib_size) * 1e6
  cpm_matrix <- t(t(count_matrix) / lib_sizes) * 1e6

  # 3. 对结果取log2(CPM + 1)（避免log(0)）
  log2_cpm <- log2(cpm_matrix + 1)
  return(log2_cpm)
}

# 应用标准化（假设count_data已去除版本号并合并重复基因）
count_data_normalized <- calculate_cpm(count_data)

params <- gsvaParam(
  exprData = as.matrix(count_data_normalized),
  geneSets = kegg_list,
  kcdf = "Poisson",
  absRanking = FALSE)

gsva_result <- gsva(params, verbose = TRUE)

group <- ifelse(grepl("WT", colnames(gsva_result)), "WT", "Mutant")
design <- model.matrix(~0 + group)
colnames(design) <- c("Mutant", "WT")

# 构建对比矩阵（Mutant vs WT）
contrast_matrix <- makeContrasts(Mutant_vs_WT = Mutant - WT, levels = design)

# 线性模型拟合
fit <- lmFit(gsva_result, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取差异通路结果
diff_pathways <- topTable(fit2, number = Inf, adjust.method = "BH")
sig_pathways <- subset(diff_pathways, adj.P.Val < 0.01 & abs(logFC) > 0.5)

write.csv(sig_pathways, file = file.path(output_dir, "GSVA_diff_pathways_all.csv"), row.names = TRUE)

# 读取两个GSVA结果文件
file1 <- "/public3/DSC/single_cell/GSE159677/GSEA/A3A_KO/GSVA_diff_pathways_all.csv"
file2 <- "/public3/DSC/single_cell/GSE159677/GSEA/GSVA_diff_pathways_all.csv"

gsva_res1 <- read.csv(file1, check.names = FALSE)  # 文件1: A3A_KO结果
gsva_res2 <- read.csv(file2, check.names = FALSE)  #文件2: 单核细胞分化泡沫细胞

common_pathways <- intersect(gsva_res1[,1], gsva_res2[,1])

colnames(gsva_res1)[1] <- "Pathway"
colnames(gsva_res2)[1] <- "Pathway"

# 提取各文件的logFC和p.adjust
result <- bind_rows(
  gsva_res1 %>%
    filter(.[[1]] %in% common_pathways) %>%
    select(1 , "logFC_A3A" = logFC, "p.adjust_A3A_KO" = adj.P.Val),
  gsva_res2 %>%
    filter(.[[1]]  %in% common_pathways) %>%
    select(1 , "logFC_General" = logFC, "p.adjust_General" = adj.P.Val)
) %>%
  group_by(Pathway) %>%
  summarise(
    logFC_A3A = first(na.omit(logFC_A3A)),  # 若存在多行取第一个
    p.adjust_A3A_KO = first(na.omit(p.adjust_A3A_KO)),
    logFC_General = first(na.omit(logFC_General)),
    p.adjust_General = first(na.omit(p.adjust_General))
  )


filtered_result <- result %>%
  filter(
    abs(logFC_A3A) > 0.5 &
      abs(logFC_General) > 0.5 &
      sign(logFC_A3A) != sign(logFC_General)  # 符号相反
  )
#Biocarta嗜酸性粒细胞通路
#Gaurnier的PSMD4靶基因
#Gouyer的TATI（肿瘤相关胰蛋白酶抑制剂）靶基因下调
#KORKOLA_CHORIOCARCINOMA_UP
#LINDGREN_BLADDER_CANCER_CLUSTER_2A_UP
#LOPEZ_MESOTHELIOMA_SURVIVAL_DN
#REACTOME_DEFECTIVE_CHST14_CAUSES_EDS_MUSCULOCONTRACTURAL_TYPE
#Vetter的PRKCA（蛋白激酶Cα）和ETS1靶基因下调
#WikiPathways的细胞因子与炎症反应通路

###gsea#####
diff_data <- read.csv(
  "/public3/DSC/single_cell/GSE159677/GSVA/Differential_Expression_Results_13+37_vs_WT.csv",
  header = TRUE,
  row.names = 1
)
rownames(diff_data) <- diff_data$symbol
genelist <- diff_data$log2FoldChange
names(genelist) <- rownames(diff_data)
genelist <- na.omit(genelist)
genelist <- sort(genelist, decreasing = TRUE)
genelist_entrez <- bitr(names(genelist),
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

# 更新 genelist 的命名
names(genelist) <- genelist_entrez$ENTREZID
genelist <- na.omit(genelist)
gsea_res <- GSEA(
  geneList = genelist,
  TERM2GENE = geneset,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 0,
  seed = 123  # 保证可重复性
)
write.csv(as.data.frame(gsea_res),
          file = file.path(output_dir, "13+37_gsea_results.csv"),
          row.names = FALSE)
file1 <- "/public3/DSC/single_cell/GSE159677/GSEA/13+37_gsea_results.csv"
file2 <- "/public3/DSC/single_cell/GSE159677/GSEA/gsea_results.csv"

gsea_res1 <- read.csv(file1, check.names = FALSE)  # 文件1: A3A_KO结果
gsea_res2 <- read.csv(file2, check.names = FALSE)  #文件2: 单核细胞分化泡沫细胞

common_pathways <- intersect(gsea_res1$Description, gsea_res2$Description)

# 提取各文件的NES和p.adjust
result <- full_join(
  gsea_res1 %>%
    filter(Description %in% common_pathways) %>%
    select(Description, "NES_A3A_KO" = NES, "p.adjust_A3A_KO" = p.adjust),
  gsea_res2 %>%
    filter(Description %in% common_pathways) %>%
    select(Description, "NES_General" = NES, "p.adjust_General" = p.adjust)
) %>%
  group_by(Description)
#	GHANDHI_旁观者辐射激活通路
#ONDER_CDH1靶标下调2
#HINATA_NFKB靶标（成纤维细胞激活）
#MARTENS_维甲酸响应激活通路
#VERHAAK_NPM1突变型AML激活通路
#BILD_HRAS致癌特征通路
#SENESE_HDAC1与HDAC2靶标激活通路
#RUTELLA_HGF响应vs CSF2RB与IL4联合响应激活通路
#DAZARD_紫外线响应（NHEK细胞）激活通路
#CROMER_肿瘤发生激活通路
#REACTOME_细胞外基质降解通路


###single_gene####
#data <- readRDS("/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/data_pseudotime.rds")
##clone13
diff_data <- read.csv(
  "/public3/DSC/single_cell/DEG_ALL.csv",
  header = TRUE,
  row.names = 1         # 第一列为基因ID（如ENSG00000183853.18_12）
)
rownames(diff_data) <- sub("\\..*", "", rownames(diff_data))

# 转换ENSEMBL ID为Symbol
ensembl_ids <- rownames(diff_data)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")  # 处理多映射问题


diff_data$symbol <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)


dup_symbols <- diff_data$symbol[duplicated(diff_data$symbol)]
print(paste("发现", length(dup_symbols), "个重复Symbol，例如:", head(dup_symbols)))

# 可选：将Symbol列移动到第二列
diff_data <- diff_data[, c(ncol(diff_data), 1:(ncol(diff_data)-1))]

#火山图
library(ggplot2)
library(ggrepel)

# 计算-log10(padj)
diff_data$log10_padj <- -log10(diff_data$padj)
diff_data <- subset(diff_data, !is.na(sig) & sig != "NA")

highlight_genes <- c(
  "CYP27A1", "FGL2", "ID3", "TNFAIP8L2", "PTAFR", "FZD1",
  "SLC11A1", "MMP19", "MMP9", "CCRL2", "CD82", "IL4I1",
  "CCL20", "CXCL8", "IL1R1", "IL1R2", "MGLL", "CEBPB",
  "NLRP1", "WIPI1", "CXCL3"
)
# 在数据框中添加一个标记列
diff_data$highlight <- ifelse(
  diff_data$symbol %in% highlight_genes,
  "Target",  # 标记为目标基因
  "Not Target"  # 非目标基因
)
# 绘制火山图
p <- ggplot(diff_data, aes(x = log2FoldChange, y = log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.8, size = 3) +
   geom_point(
    data = subset(diff_data, highlight == "Target"),
    color = "yellow", size = 4, shape = 1, stroke = 1  # 黄色空心圆，边框加粗
  ) +
  scale_color_manual(values = c("up" = "#E64B35", "down" = "#4DBBD5", "none" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
  theme_bw() +
  geom_text_repel(
    data = subset(diff_data, highlight == "Target"),
    aes(label = symbol),
    size = 4,
    color = "black",
    nudge_x = 0.5,
    direction = "y",
    segment.size = 0,
    max.overlaps = Inf,
    box.padding = 0.8,
    seed = 123
  )+
  theme(legend.position = "right")+
theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 11)
  )
print(p)
# 设置PDF输出参数
pdf(
  file = "/public3/DSC/single_cell/GSE159677/GSVA/matched_markers_clone13.pdf",
  width = 8,
  height = 6,
  pointsize = 12  #
)

print(p)  # 假设p是已生成的ggplot对象

# 关闭图形设备并保存文件
dev.off()


diff_data <- diff_data[diff_data$sig != "none", ]


matched_markers <- markers[rownames(markers) %in% diff_data$symbol, ]
matched_markers$gene <- rownames(matched_markers)
diff_data$gene <- diff_data$symbol
de_cols <- c("gene", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "sig")

matched_markers <- merge(matched_markers,
                     diff_data[, de_cols],
                     by = "gene",
                     all.x = TRUE)

# 恢复行名
#rownames(matched_markers) <- matched_markers$gene
#matched_markers$gene <- NULL
matched_markers$sign_mismatch <- sign(matched_markers$avg_log2FC) != sign(matched_markers$log2FoldChange)

matched_markers <- matched_markers[matched_markers$sign_mismatch, ]
matched_markers$sign_mismatch <- NULL
matched_markers <- matched_markers[order(matched_markers$padj), ]
rownames(head(matched_markers, 200))
write.csv(matched_markers,"/public3/DSC/single_cell/GSE159677/GSVA/matched_markers_clone13.csv")
#MMP9、CSF1R、FOLR2、CTSL、CCL2、PPARG、KLF4、S100A8/S100A9、JUN、THBS1、IL1R1、GPR183、TSC22D3
#JUN,FOSB:  MAPK/NF-kappaB/AP-1 信号通路，TNFAIP8L2：AP-1 下游蛋白，FOS :AP-1 复合物的关键组分,促进巨噬细胞向促炎（M1）型极化
#KLF4:SENP1-KLF4 信号调节 LPS 诱导的巨噬细胞 M1 极化，S100A12：单核标记物
#NLRP1：抑制 NLRP1 炎性小体激活来减轻 OX-LDL 诱导的巨噬细胞炎症。
#CCL2:促进M2激化，TREM1，IL1R1 ：促炎，KLF6，ID3 ,FADS1,CCL20 ：M2激化
#MMP9：细胞外基质降解,FCN1:敲除后显著下调
# MARCO：巨噬细胞受体，KO后显著低表达，抑制泡沫细胞的形成
#CSF1R：巨噬细胞发育存活所必需的,CEBPD：促炎转录因子
#PPARG：上下调不对，CSF2RA ：在脆弱斑块上调 FABP4：上升下降了5倍， 加剧巨噬细胞 M1 型极化和脂代谢紊乱
#PLPP3:基因变异与动脉粥样硬化易感性相关
#TFPI 存在于早期泡沫细胞形成中,M2 的 TFPI mRNA 表达升高.


##clone37##

diff_data <- read.csv(
  "/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/ko/3A_total.csv",
  header = TRUE,
  row.names = 1         # 第一列为基因ID（如ENSG00000183853.18_12）
)
rownames(diff_data) <- sub("\\..*", "", rownames(diff_data))

# 转换ENSEMBL ID为Symbol
ensembl_ids <- rownames(diff_data)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")  # 处理多映射问题


diff_data$symbol <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)


dup_symbols <- diff_data$symbol[duplicated(diff_data$symbol)]
print(paste("发现", length(dup_symbols), "个重复Symbol，例如:", head(dup_symbols)))

# 可选：将Symbol列移动到第二列
diff_data <- diff_data[, c(ncol(diff_data), 1:(ncol(diff_data)-1))]
diff_data <- diff_data[diff_data$sig != "none", ]


matched_markers <- markers[rownames(markers) %in% diff_data$symbol, ]
matched_markers$gene <- rownames(matched_markers)
diff_data$gene <- diff_data$symbol
de_cols <- c("gene", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "sig")

matched_markers <- merge(matched_markers,
                         diff_data[, de_cols],
                         by = "gene",
                         all.x = TRUE)

# 恢复行名
#rownames(matched_markers) <- matched_markers$gene
#matched_markers$gene <- NULL
matched_markers$sign_mismatch <- sign(matched_markers$avg_log2FC) != sign(matched_markers$log2FoldChange)
matched_markers <- matched_markers[matched_markers$sign_mismatch, ]
matched_markers$sign_mismatch <- NULL
matched_markers1 <- matched_markers[order(matched_markers$padj), ]
head(matched_markers1$gene, 200)
write.csv(matched_markers1,"/public3/DSC/single_cell/GSE159677/GSVA/matched_markers_clone37.csv")
#敲除后下降：CCL20、CXCL8、IL1R1、IL1RN、NFKBIA、MGLL、CEBPB、SLC11A1、THBS1、MMP9、MMP19、
#FN1、IL6R、LCP1
#敲除后上升：KLF4、UCP2、PDK4：促炎

common_genes <- intersect(matched_markers1$gene, matched_markers$gene)
print(common_genes)
# "CCL20"       "CXCL8"       "EREG"        "FCAR"        "IL1R1"       "MGLL"        "MMP19"
# "MMP9"        "PLK2"        "C5AR1"       "IGSF6"       "SPATA13"     "GPR35"       "SLC6A6"
# "CD82"        "KLF4"        "CXCL3"       "TUBA1A"      "LUCAT1"      "MRAS"        "RHOB"
# "NRIP3"       "C3AR1"       "UPP1"        "ITGB8"       "G0S2"        "TFPI"        "CEMIP2"
# "THBS1"       "SERPINA1"    "CEBPB"       "ANPEP"       "PRKACB"      "SASH1"       "CYSLTR1"
#"IL4I1"       "CHST15"      "PLA2G7"      "QPCT"        "FZD1"        "BCL2A1"      "ID3"
# "SLC11A1"     "MIR4435-2HG" "GCH1"        "TMCO3"       "EMB"         "RNF125"      "ATF5"
# "CTSL"        "CCL2"        "ADRB2"       "PTGER2"      "CYP27A1"     "RGCC"        "SESN3"
#"TSC22D3"     "RPS6KA4"     "PTAFR"       "MS4A4A"      "AQP9"        "CCRL2"       "FGL2"
# "PLD4"        "EGLN3"       "IL1R2"       "FN1"         "CLMN"        "ZFP36L2"     "WIPI1"
# "CD69"        "OTULINL"     "MT2A"        "TGFBR1"      "KIF13A"      "CLEC4E"      "CST3"
#"FGR"         "KLF2"        "TEX30"       "CCNG2"       "NLRP1"       "OTOA"        "GFRA2"
#"PCOLCE2"     "CPED1"       "FOLR2"       "CALHM6"      "RELL1"       "ZNF33A"      "ZNF124"
# "FILIP1L"     "TNFAIP8L2"   "MS4A6A"      "S1PR4"       "TMEM205"     "FAM13B"      "CD302"
#CYP27A1:胆固醇代谢
#FGL2: M1极化
#ID3：上调，M2极化
#TNFAIP8L2：AP-1 下游蛋白，M1极化
#PTAFR：预后蛋白
#	FZD1： wnt信号通路，M1
#SLC11A1：脂质过氧化
#	MMP19、MMP9
#CCRL2，CD82，IL4I1，CCL20, CXCL8，  m2,敲除后下调
#IL1R1,IL1R2
#MGLL，CEBPB M2,脂质分解酶
#NLRP1：抑制 NLRP1 炎性小体激活来减轻 OX-LDL 诱导的巨噬细胞炎症。
#WIPI1 自噬 敲除后下调
#CXCL3：CXCL3–CXCR2信号轴