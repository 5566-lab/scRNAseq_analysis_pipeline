#####macSpectrum####
library(ggplot2)
library(macSpectrum)
library(Seurat)
library(ggpubr)
#?macspec

expression_matrix <- GetAssayData(data, slot = "data") 
mac_mtx <- as.data.frame(as.matrix(expression_matrix)) 
mac_mtx <- tibble::rownames_to_column(mac_mtx, var = "geneid") 
library(org.Hs.eg.db)  # 或 org.Mm.eg.db
ensembl_ids <- mapIds(org.Hs.eg.db,
                      keys = mac_mtx$geneid,
                      column = "ENSEMBL",
                      keytype = "SYMBOL")
mac_mtx$geneid <- ensembl_ids  # 替换原基因符号为Ensembl ID
mac_mtx <- na.omit(mac_mtx) 
data@meta.data$Celltype_raw1 <- as.character(data@meta.data$Celltype_raw1)
data@meta.data$Celltype_raw1[data@meta.data$Celltype_raw1 == "FOAM_cells2"] <- "Foam cells2"
data@meta.data$Celltype_raw1[data@meta.data$Celltype_raw1 == "FOAM_cells1"] <- "Foam cells1"
data@meta.data$Celltype_raw1 <- as.factor(data@meta.data$Celltype_raw1)
feature<- data@meta.data$Celltype_raw1  #分组信息存储在Celltype_raw1

result <- macspec(mac_mtx, feature, select_hu_mo = "hum")

data@meta.data <- merge(data@meta.data, 
                        result[, c("MPI", "AMDI")], 
                        by.x = "row.names", 
                        by.y = "row.names",
                        all.x = TRUE) %>%
  column_to_rownames("Row.names")

saveRDS(data,"/dell_1/dsc/single_cell/GSE159677/hdWGCNA/Mo_Ma/data.rds")


FeaturePlot(data, features = "MPI", #巨噬细胞极化指数 (MPI) 
            pt.size = 0.8,                     
            max.cutoff = 'q98')+
  scale_colour_gradientn(
    colours = c("lightgrey", "lightgrey", "#FF0000"),  
    values = scales::rescale(c(0, 0.5, 1)),       
    breaks = seq(0, 10, 2)
  )              
FeaturePlot(data, features = "AMDI",  #巨噬细胞分化指数 (AMDI) 
            pt.size = 0.8,                     
            max.cutoff = 'q98') +
  scale_colour_gradientn(
    colours = c("lightgrey", "lightgrey", "#FF0000"),  
    values = scales::rescale(c(0, 0.8, 1)),       
    breaks = seq(0, 10, 2)
  )  
VlnPlot(data, features = "AMDI", group.by = "Celltype_raw", pt.size = 0)
VlnPlot(data, features = "MPI", group.by = "Celltype_raw", pt.size = 0)

VlnPlot(data, features = "AMDI", group.by = "Celltype_raw1", pt.size = 0)
VlnPlot(data, features = "MPI", group.by = "Celltype_raw1", pt.size = 0)

VlnPlot(data, features = "AMDI", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(data, features = "MPI", group.by = "seurat_clusters", pt.size = 0)

# 提取metadata数据
meta_data <- data@meta.data

# 基础散点图
ggplot(meta_data, aes(x = MPI, y = AMDI)) + 
  geom_point(size = 1.5, alpha = 0.6, color = "#1E90FF") +  # 调整点的大小和透明度
  theme_classic() +  # 简洁主题
  labs(x = "MPI", y = "AMDI", title = "MPI vs AMDI Distribution")
ggplot(meta_data, aes(x = MPI, y = AMDI, color = Celltype_raw1)) + 
  geom_point(size = 1.5, alpha = 0.6) 

ggplot(meta_data, aes(x = MPI, y = AMDI)) + 
  geom_jitter(width = 0.2, height = 0) +  # 水平抖动
  facet_wrap(~Celltype_raw1)

sub_data <- subset(data, subset = Celltype_raw == "Macrophage")

FeaturePlot(sub_data, features = "MPI", 
            pt.size = 0.8,  
            max.cutoff = 'q98')+
  scale_colour_gradientn(
    colours = c("lightgrey", "lightgrey", "#FF0000"),  
    values = scales::rescale(c(0, 0.5, 1)),       
    breaks = seq(0, 10, 2)
  )            
FeaturePlot(sub_data, features = "AMDI", 
            pt.size = 0.8,                     
            max.cutoff = 'q98') +
  scale_colour_gradientn(
    colours = c("lightgrey", "lightgrey", "#FF0000"),  
    values = scales::rescale(c(0, 0.8, 1)),       
    breaks = seq(0, 10, 2)
  )  


VlnPlot(sub_data, features = "AMDI", group.by = "Celltype_raw1", pt.size = 0)
VlnPlot(sub_data, features = "MPI", group.by = "Celltype_raw1", pt.size = 0)

VlnPlot(sub_data, features = "AMDI", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(sub_data, features = "MPI", group.by = "seurat_clusters", pt.size = 0)


meta_data <- sub_data@meta.data

# 散点图
ggplot(meta_data %>% filter(seurat_clusters %in% c(0,1,2,5)), 
       aes(x = MPI, y = AMDI)) + 
  geom_density_2d(
    aes(color = seurat_clusters),  # 按组着色等高线
    size = 0.8,                  # 调整线宽
    alpha = 0.7                  # 设置透明度
  ) +
  scale_color_brewer(
    palette = "Set1",            # 使用高对比度配色
    name = "Cell Type"           # 修改图例标题
  ) +
  theme_bw() +                   # 保留坐标轴和网格线
  labs(
    x = "MPI", 
    y = "AMDI", 
    title = "2D Density by Cell Type"
  ) +
  # 强制显示 x=0 和 y=0 轴线
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +  # y=0 轴线
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +  # x=0 轴线
  coord_cartesian(ylim = c(-30, 30))

ggplot(
  meta_data %>% filter(Celltype_raw1 %in% c("Macrophage", "Foam cells1", "Foam cells2")), 
  aes(x = MPI, y = AMDI)
) + 
  geom_density_2d(
    aes(color = Celltype_raw1),
    size = 0.8,
    alpha = 0.7
  ) +
  scale_color_manual(
    name = " ",  # 图例标题
    values = c(
      "Macrophage" = "#1B9E77", 
      "Foam cells1" = "#D95F02", 
      "Foam cells2" = "#7570B3"
    )
  ) +
  # 强制显示 x=0 和 y=0 轴线（与主题协调）
  geom_hline(
    yintercept = 0, 
    color = "black", 
    linewidth = 0.5, 
    linetype = "solid"  # 确保与 combined_plot 的轴线风格一致
  ) +  
  geom_vline(
    xintercept = 0, 
    color = "black", 
    linewidth = 0.5, 
    linetype = "solid"
  ) +
  # 坐标轴范围控制
  coord_cartesian(ylim = c(-30, 30)) +
  # 应用统一主题
  theme_custom +
  # 标签设置（与 combined_plot 一致）
  labs(
    x = "MPI", 
    y = "AMDI"
  )


ggplot(meta_data , 
       aes(x = MPI, y = AMDI)) + 
  geom_density_2d(
    aes(color = Celltype_raw),  # 按组着色等高线
    size = 0.8,                  # 调整线宽
    alpha = 0.7                  # 设置透明度
  ) +
  scale_color_brewer(
    palette = "Set1",            # 使用高对比度配色
    name = "Cell Type"           # 修改图例标题
  ) +
  theme_bw() +                   # 保留坐标轴和网格线
  labs(
    x = "MPI", 
    y = "AMDI", 
    title = "2D Density by Cell Type"
  ) +
  # 强制显示 x=0 和 y=0 轴线
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +  # y=0 轴线
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +  # x=0 轴线
  coord_cartesian(ylim = c(-30, 30))

ggplot(meta_data , 
       aes(x = MPI, y = AMDI)) + 
  geom_density_2d(
    aes(color = AC_PA),  # 按组着色等高线
    size = 0.8,                  # 调整线宽
    alpha = 0.7                  # 设置透明度
  ) +
  scale_color_brewer(
    palette = "Set1",            # 使用高对比度配色
    name = "Cell Type"           # 修改图例标题
  ) +
  theme_bw() +                   # 保留坐标轴和网格线
  labs(
    x = "MPI", 
    y = "AMDI", 
    title = "2D Density by Cell Type"
  ) +
  # 强制显示 x=0 和 y=0 轴线
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +  # y=0 轴线
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +  # x=0 轴线
  coord_cartesian(ylim = c(-30, 30))

ggplot(meta_data %>% filter(Celltype_raw1 %in% c("Macrophage", "Foam cells1", "Foam cells2")), 
       aes(x = MPI, y = AMDI)) + 
  geom_density_2d(
    aes(color = AC_PA),  # 按组着色等高线
    size = 0.8,                  # 调整线宽
    alpha = 0.7                  # 设置透明度
  ) +
  scale_color_brewer(
    palette = "Set1",            # 使用高对比度配色
    name = "Cell Type"           # 修改图例标题
  ) +
  theme_bw() +                   # 保留坐标轴和网格线
  labs(
    x = "MPI", 
    y = "AMDI", 
    title = "2D Density by Cell Type"
  ) +
  # 强制显示 x=0 和 y=0 轴线
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +  # y=0 轴线
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +  # x=0 轴线
  coord_cartesian(ylim = c(-30, 30))



library(clusterProfiler)
library(org.Hs.eg.db)
# 查看原始数据结构
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
head(count_data)


entrez_ids <- rownames(count_data)  # 确保行名为字符型（非数值型）

# 执行转换（自动过滤无效ID）
symbol_df <- bitr(
  entrez_ids,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

# 合并结果到原数据框
count_data$SYMBOL <- symbol_df$SYMBOL[match(rownames(count_data), symbol_df$ENTREZID)]
count_data <- na.omit(count_data)  # 移除未匹配的基因
rownames(count_data) <- make.unique(count_data$SYMBOL)  # 处理重复SYMBOL（如TP53.1）
count_data <- subset(count_data, select = -SYMBOL) 
count_data$geneid <- rownames(count_data)
count_data <- count_data[, c("geneid", setdiff(colnames(count_data), "geneid"))]
sample_cols <- colnames(count_data)[-1]  # geneid 是第一列

# 生成分组标签
sample_groups <- ifelse(grepl("^KO", sample_cols), "KO", "WT")
feature1 <- factor(sample_groups)
result <- macspec(count_data, feature1, select_hu_mo = "hum")  #result==Nan



