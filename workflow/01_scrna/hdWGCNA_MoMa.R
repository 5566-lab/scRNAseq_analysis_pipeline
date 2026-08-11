library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
#devtools::install_github("NightingaleHealth/ggforestplot")
#devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)
align_umap_plots <- function(plot_AC, plot_PA) {
  # 提取两组数据的UMAP坐标范围
  ac_range <- layer_scales(plot_AC)$x$range$range  # 获取AC组的x轴范围
  pa_range <- layer_scales(plot_PA)$x$range$range  # 获取PA组的x轴范围
  x_min <- min(ac_range[1], pa_range[1])          # 计算x轴最小值
  x_max <- max(ac_range[2], pa_range[2])          # 计算x轴最大值

  # 同理获取y轴范围
  y_range <- range(layer_scales(plot_AC)$y$range$range,
                   layer_scales(plot_PA)$y$range$range)
  y_min <- y_range[1]
  y_max <- y_range[2]

  # 应用统一坐标范围和比例
  plot_AC <- plot_AC +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme(
      aspect.ratio = 1,           # 确保画布为正方形
      plot.margin = margin(5,5,5,5)  # 统一边距
    )

  plot_PA <- plot_PA +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5)
    )

  # 使用patchwork精确对齐
  p <- (plot_AC | plot_PA) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom')

  # 返回对齐后的图
  return(p)
}
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(123)

# Workflow起点：sub_integrated_data
sub_integrated_data <- readRDS('/public3/DSC/single_cell/Result/figer_new/group_result.rds')
data <- sub_integrated_data

enableWGCNAThreads(nThreads = 8)
# set up seurat object for WGCNA
# Use assay = 'RNA' and layer = 'data' for SeuratObject 5.0.0+ compatibility
data <- SetupForWGCNA(
  data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "wgcna" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = "AC_PA", # specify the columns in seurat_obj@meta.data to group by
  k = 20, # nearest-neighbors parameter 通常是20-75，10万个细胞可以用50
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'AC_PA' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
data <- NormalizeMetacells(data)

# set up the expression matrix
data <- SetDatExpr(
  data,
  group_name = c("atherosclerotic core","proximal adjacent"), # the name of the group of interest in the group.by column
  group.by='AC_PA', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data (replaces deprecated 'slot' parameter in SeuratObject 5.0.0+)
)

# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

# assemble with patchwork
p <- wrap_plots(plot_list, ncol=2)
# Output paths
output_root <- "/public3/DSC/single_cell/Result/openclaw"
output_hdWGCNA <- file.path(output_root, "hdWGCNA")

ggsave(file.path(output_hdWGCNA, "data_softpowers.pdf"), p, width=10 ,height=8)

power_table <- GetPowerTable(data)
head(power_table)

# construct co-expression network:
data <- ConstructNetwork(
  data,
  soft_power = 10,
  overwrite_tom = TRUE,
  tom_name = 'CAD_MM' # name of the topoligical overlap matrix written to disk
)

pdf(file = file.path(output_hdWGCNA, "data_Dendrogram.pdf"),width=10, height=8)
PlotDendrogram(data, main='hdWGCNA Dendrogram')
dev.off()

modules <- data@misc$wgcna$wgcna_modules
table(modules$module)
write.csv(modules,file.path(output_hdWGCNA,"data_modules.csv"))

#每个细胞对每个模块的特征值
data <- ModuleEigengenes(
  data
)

#鉴定模块内hub基因
data <- ModuleConnectivity(
  data
)

# 计算每个细胞对于每个模块hub基因的表达活性(module score),可使用seurat包或者Ucell包
data <- ModuleExprScore(
  data,
  n_genes = 25,
  method='Seurat'
)
#将细胞对于模块的特征值，整合到seurat的meta.data中
MEs <- GetMEs(data)
mods <- colnames(MEs)
data@meta.data <- cbind(data@meta.data, MEs)


saveRDS(data,file.path(output_hdWGCNA,"data.rds"))

data <- readRDS(file.path(output_hdWGCNA,"data.rds"))
hub_df <- GetHubGenes(data, n_hubs = 10)
write.csv(hub_df,file.path(output_hdWGCNA,"data_genes.csv"))

# 绘制树状图与热图
pdf(file = file.path(output_hdWGCNA, "data_Dendrogram1.pdf"),width=10, height=8)
plotEigengeneNetworks(
  MEs,
  "Eigengene dendrogram",
  plotHeatmaps = FALSE,
  marDendro = c(0, 4, 2, 0)
)

dev.off()

pdf(file = file.path(output_hdWGCNA, "data_Heatmaps.pdf"),width=10, height=8)
par(mar = c(5, 5, 4, 2))
plotEigengeneNetworks(
  MEs,
  "Eigengene adjacency heatmap",
  plotDendrograms = FALSE,
  xLabelsAngle = 90,
  marHeatmap = c(6, 6, 4, 2)  # 热图边距：下、左、上、右（原下边距10→6，上边距1→4）
)

dev.off()

#每个细胞对于每个模块的特征值
pdf(file = file.path(output_hdWGCNA, "data_ModuleFeaturePlot.pdf"), width = 15, height = 12)  # 根据模块数量调整宽高

# 绘制模块特征图并调整边距
plot_list <- ModuleFeaturePlot(
  data,
  features = 'MEs',  # 绘制模块特征基因（hMEs）
  order = TRUE,       # 按 hMEs 值从高到低排序点

)

# 组合图形并调整布局
combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "Module Feature Plots",
    theme = theme(
      plot.title = element_text(
        size = 25,
        hjust = 0.5,
        margin = margin(b = 15)
      )
    )
  ) &
  theme( # 使用&运算符统一设置所有子图主题
    plot.title = element_text(
      size = 18,  # 子标题字体大小（建议值）
      face = "bold", # 加粗
      hjust = 0.5 # 保持标题居中
    )
  )

# 输出图形
print(combined_plot)
dev.off()
#文件过于巨大，无法常规模式保存，png都有93mb,用浏览器储存:ModuleFeaturePlot.png
desired_order <- c("blue", "brown", "green","yellow",'turquoise','purple','grey','black',"magenta",'greenyellow','pink','tan','red')  # 按实际组名修改
mods_ordered <- mods[match(desired_order, mods)]
p <- DotPlot(
  data,
  features = mods_ordered,
  group.by = 'metacell_grouping',
  scale = FALSE,
  dot.scale = 6,
  cols = c('blue', 'red')
) +
  # 添加标题并设置样式
  labs(title = "Macrophage Module Expression Heatmap") +
  # 交换XY轴
  coord_flip() +
  # 统一主题设置
  theme(
    plot.title = element_text(
      face = "bold",          # 加粗
      size = 16,              # 字号
      hjust = 0.5,            # 水平居中[2,3,6](@ref)
      color = "black",        # 颜色
      margin = margin(b = 10) # 下边距调整
    ),
    axis.title.x = element_blank(),   # 隐藏X轴标题[3](@ref)
    axis.title.y = element_blank(),   # 隐藏Y轴标题[3](@ref)
    axis.text.x = element_text(
      angle = 45,            # X轴标签旋转45度
      hjust = 1,             # 右对齐[3](@ref)
      vjust = 1              # 垂直对齐微调
    ),
    axis.text.y = element_text(
      face = "italic",        # Y轴基因名斜体
      margin = margin(r = 5) # 右侧增加间距
    )
  )

# 输出图形
print(p)

# 保存高清PDF(推荐矢量图格式)
ggsave(
  file.path(output_hdWGCNA, "Mono_module.pdf"),
  plot = p,
  width = 8,   # 加宽以适应Y轴长标签[3](@ref)
  height = 7,  # 增加高度显示完整模块名
  device = "pdf",
  limitsize = FALSE
)

# 提取元数据
meta_data <- data@meta.data


# 检查分组类别
table(meta_data$AC_PA)

# 如果存在其他类别或NA值，需要清理数据：
meta_data <- meta_data[meta_data$AC_PA %in% c("atherosclerotic core", "proximal adjacent"), ]

# 定义模块列表
modules <- c("magenta", "red", "turquoise", "greenyellow", "tan", "pink", "yellow"
             ,"brown", "green", "blue",  "purple", "black")

results <- data.frame(
  Module = character(),
  p_value = numeric(),
  mean_AC = numeric(),
  mean_PA = numeric(),
  stringsAsFactors = FALSE
)
# 循环计算每个模块的统计检验
for (module in modules) {
  # 提取模块表达值并移除NA值
  expr_AC <- na.omit(meta_data[meta_data$AC_PA == "atherosclerotic core", module])
  expr_PA <- na.omit(meta_data[meta_data$AC_PA == "proximal adjacent", module])

  # 检查有效样本量（每组至少1个样本）
  if (length(expr_AC) < 1 | length(expr_PA) < 1) {
    warning(paste0("Skipping module ", module, ": PA组或AC组有效样本不足（AC=", length(expr_AC), ", PA=", length(expr_PA), "）"))
    next  # 跳过当前模块
  }

  # 执行Wilcoxon秩和检验
  wilcox_test <- wilcox.test(expr_AC, expr_PA)

  # 存储结果
  results <- rbind(results, data.frame(
    Module = module,
    p_value = wilcox_test$p.value,
    mean_AC = mean(expr_AC),
    mean_PA = mean(expr_PA)
  ))
}

# 输出结果
print(results)
results$p_adj <- p.adjust(results$p_value, method = "BH")
print(results[order(results$p_adj), ])

Module_Colors <- c(
  "magenta" = "magenta",
  "red" = "red",
  "turquoise" = "turquoise",
  "greenyellow" = "greenyellow",
  "tan" = "tan",
  "pink" = "pink",
  "yellow" = "yellow",
  "brown" = "brown",
  "green" = "green",
  "blue" = "blue",
  "purple" = "purple",
  "black" = "black"
)

# 绘制火山图
ggplot(results, aes(x = mean_AC - mean_PA, y = log_p, color = Module)) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Module Differential Expression",
       x = "Effect Size (AC - PA)",
       y = "-log10(Adjusted p-value)") +
  theme_classic() +
  scale_color_manual(values = Module_Colors)



#### 每个module的GO富集 ####
file_paths <- file.path(output_hdWGCNA, 'data_modules.csv')
modules_list <- lapply(file_paths, read.csv)

# 给列表元素命名（可选）
names(modules_list) <- 'MM'
output_root_GO <- file.path(output_hdWGCNA, "GO")

for (i in seq_along(modules_list)) {
  # 获取当前模块名称（假设已命名）
  module_name <- names(modules_list)[i]
  # 若未命名，自动生成名称（如"Module1"）
  if(is.null(module_name)) module_name <- paste0("Module", i)

  # 创建模块专属文件夹
  module_dir <- file.path(output_root, module_name)
  if(!dir.exists(module_dir)) dir.create(module_dir, recursive = TRUE)

  # 提取当前模块数据
  module_df <- modules_list[[i]]

  # 遍历模块中的颜色分组
  co <- as.data.frame(unique(module_df$color))

  for (j in 1:nrow(co)) {
    mc <- co[j, 1]
    ge <- data.frame(Name = module_df[module_df$color == mc, 2])

    # GO富集分析
    GO <- enrichGO(
      gene = ge$Name,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = 'BH',
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    )

    if (nrow(GO@result) > 0) {
      # 构建文件路径前缀
      file_prefix <- file.path(module_dir, mc)

      # 保存CSV结果
      write.csv(
        data.frame(ID = row.names(GO@result), GO@result),
        file = paste0(file_prefix, "_GO.csv"),
        row.names = FALSE
      )

      # 生成并保存可视化
      p <- barplot(GO, drop = TRUE, showCategory = 30) +
        ggtitle(paste0("GO Enrichment for : ", mc)) +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 8)
        ) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

      ggsave(
        filename = paste0(file_prefix, ".pdf"),
        plot = p,
        width = 10,
        height = 8
      )
    }

    # 清理临时对象
    rm(mc, ge, GO, p)
    gc()
  }
}

turquoise <- read.csv(file.path(output_root_GO, 'Mono/turquoise_GO.csv'))

turquoise$Score <- -log10(turquoise$p.adjust)
turquoise <- turquoise[c(1,2,4,7,8,10,11,17,27,29),]
turquoise <- turquoise %>% arrange(desc(Score))
turquoise$Description <- factor(turquoise$Description, levels = rev(turquoise$Description))
p <- ggplot(data = turquoise, aes(x = Description, y = Score)) +
  geom_col(fill = "#ec8574") +
  theme_bw() +
  scale_y_continuous(limits = c(0,max(turquoise$Score))) +
  labs(x = 'Go Term', y = '-Log10(adj.p.value)', title = 'Monocyte Module turquoise : GO - BP', color = '') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  geom_text(data = turquoise, aes(x = Description, y = max(Score)/30, label = Description),
            hjust = 0, size = 4.5,  color = "black")

ggsave(file.path(output_root_GO, 'Mono/turquoise_BAR10.pdf'),p,width = 5, height = 5)
red <- read.csv(file.path(output_root_GO, 'Mono/red_GO.csv'))














# 定义函数
align_umap_plots <- function(plot_AC, plot_PA) {
  # 提取两组数据的UMAP坐标范围
  ac_range <- layer_scales(plot_AC)$x$range$range  # 获取AC组的x轴范围
  pa_range <- layer_scales(plot_PA)$x$range$range  # 获取PA组的x轴范围
  x_min <- min(ac_range[1], pa_range[1])          # 计算x轴最小值
  x_max <- max(ac_range[2], pa_range[2])          # 计算x轴最大值

  # 同理获取y轴范围
  y_range <- range(layer_scales(plot_AC)$y$range$range,
                   layer_scales(plot_PA)$y$range$range)
  y_min <- y_range[1]
  y_max <- y_range[2]

  # 应用统一坐标范围和比例
  plot_AC <- plot_AC +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme(
      aspect.ratio = 1,           # 确保画布为正方形
      plot.margin = margin(5,5,5,5)  # 统一边距
    )

  plot_PA <- plot_PA +
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5)
    )

  # 使用patchwork精确对齐
  p <- (plot_AC | plot_PA) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom')

  # 返回对齐后的图
  return(p)
}




data <- readRDS("/public3/DSC/single_cell/Result/openclaw/hdWGCNA/data.rds")
expression_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")
cell_metadata <- data@meta.data

groups <- unique(cell_metadata$AC_PA) # 获取分组类别
# 创建gene_metadata（确保行名与表达矩阵一致）
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  module = "Default_Module",  # 替换为实际的模块信息
  row.names = rownames(expression_matrix)  # 关键步骤！
)




cell_ids <- rownames(cell_metadata)

sub_matrix <- expression_matrix

cds_MM<- new_cell_data_set(
  sub_matrix,
  cell_metadata = cell_metadata[cell_ids, ],
  gene_metadata = gene_metadata
)

# 导入对应的UMAP坐标（需与当前分组匹配）
reducedDims(cds_MM)$UMAP <- data@reductions$umap@cell.embeddings[cell_ids, ]
# 归一化与特征选择
#cds_MM <- cds_Mono



cds_MM <- preprocess_cds(cds_MM, num_dim = 50)  # 建议根据主成分解释方差图调整维度数[1,2](@ref)
plot_pc_variance_explained(cds_MM)              # 选择拐点前的维度
preprocess_cds(cds_MM, num_dim=15)
# 降维（UMAP/PCA/t-SNE）
cds_MM <- reduce_dimension(cds_MM, reduction_method = "UMAP")

cds_MM <- cluster_cells(cds_MM, resolution = 1e-5)
cds_MM <- learn_graph(
  cds_MM,
  close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =5 ,  # 增加该值以修剪短分支（数值需根据数据调整）
    prune_graph = TRUE         # 启用自动修剪分支
  )
)

cds_partition <- plot_cells(cds_MM, color_cells_by="partition",group_label_size = 4,)

ggsave(cds_partition,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/S1.5_monocle_partition.png',width = 18,height = 8)
ggsave(cds_partition,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/S1.5_monocle_partition.pdf',width = 18,height = 8)






partitions <- cds_MM@clusters$UMAP$partitions

######拆分partition####
#######################

partition1_indices <- cds_MM@clusters$UMAP$partitions == "1"
cds_MM_foam  <- cds_MM[, partition1_indices]
cds_MM_foam  <- preprocess_cds(cds_MM_foam , num_dim = 50)  # 建议根据主成分解释方差图调整维度数[1,2](@ref)
plot_pc_variance_explained(cds_MM_foam )              # 选择拐点前的维度
preprocess_cds(cds_MM_foam , num_dim=20)
# 降维（UMAP/PCA/t-SNE）
cds_MM_foam  <- reduce_dimension(cds_MM_foam , reduction_method = "UMAP")


#cds_MM_foam <- cds_MM


cds_MM_foam <- cluster_cells(cds_MM_foam, resolution = 1e-5)
cds_MM_foam <- learn_graph(
  cds_MM_foam,
  close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =3 #,  # 增加该值以修剪短分支（数值需根据数据调整）
    #prune_graph = TRUE         # 启用自动修剪分支
  )
)
plot_cells(cds_MM_foam,
           color_cells_by ="seurat_clusters", #'AC_PA',#"Celltype_raw1", #"pseudotime",
           #genes = "APOBEC3A",
           cell_size = 1.5,
           label_groups_by_cluster = FALSE,
           group_label_size = 4,
           show_trajectory_graph = T)

p <- plot_cells(
  cds_MM_foam,
  color_cells_by = "Celltype_raw1",  # 按分组着色
  cell_size = 1.5,
  label_groups_by_cluster = FALSE,
  group_label_size = 4,
  show_trajectory_graph = TRUE
) +
  facet_wrap(~AC_PA, nrow = 1) +  # 横向分面（左右布局）
  theme(
    strip.text = element_text(size = 12),  # 分面标题字体
    strip.background = element_blank()     # 分面标题背景透明
  )
print(p)
ggsave(p,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/F2.5_monocle_celltype_cluster.png',width = 18,height = 8)
ggsave(p,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/F2.5_monocle_celltype_cluster.pdf',width = 18,height = 8)




p <- plot_cells(
  cds_MM_foam,
  genes = "APOBEC3A",
  label_groups_by_cluster = FALSE,
  cell_size = 1.5,
  group_label_size = 4,
  show_trajectory_graph = TRUE
) +
  facet_wrap(~AC_PA, nrow = 1) +  # 横向分面（左右布局）
  theme(
    strip.text = element_text(size = 12),  # 分面标题字体
    strip.background = element_blank()     # 分面标题背景透明
  )+
  scale_color_gradient(
    low = "gray90",            # 低表达设为浅灰色
    high = "red3",             # 高表达设为深红色
    #breaks = c(0, 2, 5),       # 自定义颜色断点（根据实际表达范围调整）
    #limits = c(0, 5)           # 限制颜色映射范围（避免极端值影响对比度）
  ) +
  labs(color = "APOBEC3A Expression")+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

print(p)
ggsave(p,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/F2.5_monocle_APOBEC3A.png',width = 18,height = 8)
ggsave(p,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/F2.5_monocle_APOBEC3A.pdf',width = 18,height = 8)





cds_MM_foam <- order_cells(cds_MM_foam)

plot_AC <- plot_cells(
  cds_MM_foam,
  color_cells_by = "pseudotime",    # 按伪时间着色
  label_cell_groups = FALSE,       # 关闭细胞群标签
  label_leaves = TRUE,             # 开启轨迹叶节点标签[1](@ref)
  label_branch_points = TRUE,      # 开启轨迹分支点标签[1](@ref)
  trajectory_graph_color = "black",
  # trajectory_graph_label_size = 5,  # 增大轨迹标签字体[1](@ref)
  cell_size = 1,
  alpha = 0.8
) +
  facet_wrap(~AC_PA, nrow = 1) +   # 分面展示AC和PA组
  theme(
    strip.text = element_text(size = 14, face = "bold", color = "black"),  # 分面标题加粗[2](@ref)
    strip.background = element_rect(fill = "#F5F5F5", color = NA),        # 分面标题背景色
    plot.margin = margin(20, 20, 20, 20)                                   # 调整边距防止标签被裁剪
  ) +
  scale_color_gradientn(
    colours = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red'),       # 自定义伪时间色阶
    name = "Pseudotime",
    guide = guide_colorbar(barwidth = 1.5, title.position = "top")         # 调整图例位置和样式[2](@ref)
  ) +
  theme_dr() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),     # 坐标轴标题
    axis.text = element_text(size = 10),       # 坐标轴刻度
    legend.title = element_text(size = 12),    # 图例标题
    legend.text = element_text(size = 10)      # 图例数值
  )

plot_AC

plot_cells(cds_MM_foam, color_cells_by="partition",group_label_size = 4,)
###添加伪时序信息到 Seurat 中 #####
data$pseudotime <- pseudotime(cds_MM_foam)
summary(data$pseudotime)


###选择不同分支的细胞 #####
cdsMM_AC_subset1 <- choose_graph_segments(cds_MM_foam)#健康_marcophage
cdsMM_AC_subset2 <- choose_graph_segments(cds_MM_foam)#疾病_foam


subset_list <- list(
  subset1 = cdsMM_AC_subset1,
  subset2 = cdsMM_AC_subset2

)

# 定义通用标记函数
mark_subset_cells <- function(data, subset, prefix) {
  subset_cells <- colnames(subset)
  data@meta.data[[prefix]] <- ifelse(colnames(data) %in% subset_cells, "Yes", "No")
  return(data)
}
# 循环处理每个子集
for (i in seq_along(subset_list)) {
  data <- mark_subset_cells(
    data = data,
    subset = subset_list[[i]],
    prefix = paste0("cdsMM_sub", i)  # 生成动态列名（如cdsMM_sub1）
  )
}


data@meta.data[, grep("cdsMM_sub", colnames(data@meta.data))]

# 创建合并绘图的数据框架
combined_data <- data.frame()

# 定义颜色调色板（根据子集数量扩展）
subset_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")[1:length(subset_list)]

# 循环处理每个子集
for(i in seq_along(subset_list)){
  # 提取当前子集名称对应的元数据列
  subset_col <- paste0("cdsMM_sub", i)

  # 筛选属于当前子集的细胞
  subset_cells <- data@meta.data %>%
    filter(!!sym(subset_col) == "Yes") %>%
    rownames()

  # 提取基因表达与伪时间数据
  subset_df <- data.frame(
    pseudotime = data@meta.data[subset_cells, "pseudotime"],
    APOBEC3A = GetAssayData(data, assay = "RNA", slot = "data")["APOBEC3A", subset_cells],
    subset = names(subset_list)[i]
  )

  combined_data <- rbind(combined_data, subset_df)
}

# 绘制联合轨迹图
ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = subset)) +
  #geom_point(alpha = 0.6, size = 1.2) +  # 散点显示细胞分布
  geom_smooth(
    method = "gam",   # 采用GAM模型适应复杂轨迹[3](@ref)
    formula = y ~ s(x, bs = "tp"),
    se = FALSE,
    linewidth = 1.5,
    alpha = 0.8
  ) +
  # geom_smooth(method = "loess", se = FALSE, linewidth = 1.5) +  # 趋势线
  scale_color_manual(values = subset_colors,labels = c( "AC_Mono", "AC_Foam1","AC_Foam2") ) +
  labs(x = "Pseudotime", y = "APOBEC3A Expression (log-normalized)",
       title = "APOBEC3A Expression Dynamics Across Trajectory Subsets") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "grey90"))










#################
#################
#################



# 创建分组cds列表
cds_list <- lapply(groups, function(group) {
  # 筛选当前组的细胞ID
  cell_ids <- rownames(cell_metadata)[cell_metadata$AC_PA == group]

  # 子集化表达矩阵（注意行列对应关系）
  sub_matrix <- expression_matrix[, colnames(expression_matrix) %in% cell_ids]

  # 构建分组cds对象
  group_cds <- new_cell_data_set(
    sub_matrix,
    cell_metadata = cell_metadata[cell_ids, ],
    gene_metadata = gene_metadata
  )

  # 导入对应的UMAP坐标（需与当前分组匹配）
  reducedDims(group_cds)$UMAP <- data@reductions$umap@cell.embeddings[cell_ids, ]

  # 执行必须的聚类步骤
  cluster_cells(group_cds, reduction_method = "UMAP")
})

# 命名列表方便后续调用
names(cds_list) <- groups

# 示例调用AC组的数据
cdsMM_AC <- cds_list[["atherosclerotic core"]]
cdsMM_PA <- cds_list[["proximal adjacent"]]

cdsMM_AC <- learn_graph(
  cdsMM_AC,
  close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =6 ,  # 增加该值以修剪短分支（数值需根据数据调整）
    prune_graph = TRUE         # 启用自动修剪分支
  )
)

cdsMM_PA <- learn_graph(
  cdsMM_PA,
  #close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =3 ,  # 增加该值以修剪短分支（数值需根据数据调整）
    prune_graph = TRUE         # 启用自动修剪分支
  )
)

plot_AC <- plot_cells(
  cdsMM_AC,
  color_cells_by = "seurat_clusters",  # 或 "Celltype_raw"
  #label_groups_by_cluster = FALSE,
  cell_size = 1,
  show_trajectory_graph = TRUE  # 确保轨迹图显示（默认已开启）
)+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

plot_PA<- plot_cells(
  cdsMM_PA,
  color_cells_by = "seurat_clusters",  # 或 "Celltype_raw"
  #label_groups_by_cluster = FALSE,
  cell_size = 1,
  show_trajectory_graph = TRUE  # 确保轨迹图显示（默认已开启）
)+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

seurat_clusters_plot <- align_umap_plots(plot_AC, plot_PA)
seurat_clusters_plot
ggsave(seurat_clusters_plot,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/F2.5_MM_celltype_cluster.png',width = 18,height = 8)
ggsave(seurat_clusters_plot,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/F2.5_MM_celltype_cluster.pdf',width = 18,height = 8)





plot_AC <- plot_cells(
  cdsMM_AC,
  genes = "APOBEC3A",          # 指定目标基因
  label_cell_groups = FALSE,   # 关闭聚类标签
  cell_size = 1.5,             # 调整点大小
  show_trajectory_graph = TRUE # 显示轨迹骨架
) +
  scale_color_gradient(
    low = "gray90",            # 低表达设为浅灰色
    high = "red3",             # 高表达设为深红色
    #breaks = c(0, 2, 5),       # 自定义颜色断点（根据实际表达范围调整）
    #limits = c(0, 5)           # 限制颜色映射范围（避免极端值影响对比度）
  ) +
  labs(color = "APOBEC3A Expression")+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

plot_PA <- plot_cells(
  cdsMM_PA,
  genes = "APOBEC3A",          # 指定目标基因
  label_cell_groups = FALSE,   # 关闭聚类标签
  cell_size = 1.5,             # 调整点大小
  show_trajectory_graph = TRUE # 显示轨迹骨架
) +
  scale_color_gradient(
    low = "gray90",            # 低表达设为浅灰色
    high = "red3",             # 高表达设为深红色
    #breaks = c(0, 2, 5),       # 自定义颜色断点（根据实际表达范围调整）
    #limits = c(0, 5)           # 限制颜色映射范围（避免极端值影响对比度）
  ) +
  labs(color = "APOBEC3A Expression")+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

plot_APOBEC3A <- align_umap_plots(plot_AC, plot_PA)
ggsave(plot_APOBEC3A,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/F2.5_MM_APOBEC3A.png',width = 18,height = 8)
ggsave(plot_APOBEC3A,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/F2.5_MM_APOBEC3A.pdf',width = 18,height = 8)






cdsMM_AC <- order_cells(cdsMM_AC)
cdsMM_PA <- order_cells(cdsMM_PA)

plot_AC <- plot_cells(cdsMM_AC,
                      color_cells_by = "pseudotime",
                      label_cell_groups = FALSE,
                      label_leaves = FALSE,
                      label_branch_points = FALSE,
                      trajectory_graph_color = "black",
                      cell_size = 1)+
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red'),name = "Pseudotime",
  )+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

plot_PA <- plot_cells(cdsMM_PA,
                      color_cells_by = "pseudotime",
                      label_cell_groups = FALSE,
                      label_leaves = FALSE,
                      label_branch_points = FALSE,
                      trajectory_graph_color = "black",
                      cell_size = 1)+
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red'),name = "Pseudotime",
    limits = c(0, 40)
  )+
  theme_dr() + theme(panel.grid=element_blank(),
                     plot.title = element_blank()) + NoLegend()

pseudotime_plot <- align_umap_plots(plot_AC, plot_PA)
p<- plot_AC|plot_PA

ggsave(pseudotime_plot,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/F2.5_MM_Pseudotim.png',width = 18,height = 8)
ggsave(pseudotime_plot,filename = '/public3/DSC/single_cell/GSE159677/monocle3/monocle/F2.5_MM_Pseudotim.pdf',width = 18,height = 8)




ciliated_AC_test_res <- graph_test(cdsMM_AC, neighbor_graph = "principal_graph", cores = 4)
ciliated_PA_test_res <- graph_test(cdsMM_AC, neighbor_graph = "principal_graph", cores = 4)
write.csv(ciliated_AC_test_res,"/public3/DSC/single_cell/GSE159677/monocle3/monocle/ciliated_AC_test_res.csv")
write.csv(ciliated_PA_test_res,"/public3/DSC/single_cell/GSE159677/monocle3/monocle/ciliated_PA_test_res.csv")

plot_cells(cdsMM_AC, color_cells_by="partition")
plot_cells(cdsMM_PA, color_cells_by="partition")

###添加伪时序信息到 Seurat 中 #####
data$pseudotime <- pseudotime(cdsMM_PA)
data$pseudotime <- pseudotime(cdsMM_AC)
summary(data$pseudotime)


###选择不同分支的细胞 #####
cdsMM_AC_subset1 <- choose_graph_segments(cdsMM_AC)#疾病_泡沫细胞
cdsMM_AC_subset2 <- choose_graph_segments(cdsMM_AC)#疾病_健康细胞1
cdsMM_AC_subset3 <- choose_graph_segments(cdsMM_AC)#疾病_健康细胞2
cdsMM_PA_subset1 <- choose_graph_segments(cdsMM_PA)#健康_健康细胞
cdsMM_PA_subset2 <- choose_graph_segments(cdsMM_PA)#健康_健康细胞1
cdsMM_PA_subset3 <- choose_graph_segments(cdsMM_PA)#健康_健康细胞2


subset_list <- list(
  subset1 = cdsMM_AC_subset1,
  subset2 = cdsMM_AC_subset2,
  subset3 = cdsMM_AC_subset3,
  subset4 = cdsMM_PA_subset1,
  subset5 = cdsMM_PA_subset2,
  subset6 = cdsMM_PA_subset3
)


subset_list <- list(
  subset1 = cdsMM_AC_subset1,#疾病_泡沫细胞
  subset4 = cdsMM_PA_subset1,#健康_健康细胞1
  subset5 = cdsMM_PA_subset2#健康_泡沫细胞1
)

# 定义通用标记函数
mark_subset_cells <- function(data, subset, prefix) {
  subset_cells <- colnames(subset)
  data@meta.data[[prefix]] <- ifelse(colnames(data) %in% subset_cells, "Yes", "No")
  return(data)
}
# 循环处理每个子集
for (i in seq_along(subset_list)) {
  data <- mark_subset_cells(
    data = data,
    subset = subset_list[[i]],
    prefix = paste0("cdsMM_sub", i)  # 生成动态列名（如cdsMM_sub1）
  )
}


data@meta.data[, grep("cdsMM_sub", colnames(data@meta.data))]

# 创建合并绘图的数据框架
combined_data <- data.frame()

subset_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                   "#8C564B", "#E377C2")[1:length(subset_list)]

# 循环处理每个子集
for(i in seq_along(subset_list)){
  # 提取当前子集名称对应的元数据列
  subset_col <- paste0("cdsMM_sub", i)

  # 筛选属于当前子集的细胞
  subset_cells <- data@meta.data %>%
    filter(!!sym(subset_col) == "Yes") %>%
    rownames()

  # 提取基因表达与伪时间数据
  subset_df <- data.frame(
    pseudotime = data@meta.data[subset_cells, "pseudotime"],
    APOBEC3A = GetAssayData(data, assay = "RNA", slot = "data")["APOBEC3A", subset_cells],
    subset = names(subset_list)[i]
  )

  combined_data <- rbind(combined_data, subset_df)
}

# 绘制联合轨迹图
ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = subset)) +
  #geom_point(alpha = 0.6, size = 1.2) +  # 散点显示细胞分布
  geom_smooth(
    method = "gam",   # 采用GAM模型适应复杂轨迹[3](@ref)
    formula = y ~ s(x, bs = "tp"),
    se = FALSE,
    linewidth = 1.5,
    alpha = 0.8
  ) +
 # geom_smooth(method = "loess", se = FALSE, linewidth = 1.5) +  # 趋势线
  scale_color_manual(values = subset_colors,labels = c( "AC_Foam", "PA_mono", "PA_foam") ) +
  labs(x = "Pseudotime", y = "APOBEC3A Expression (log-normalized)",
       title = "APOBEC3A Expression Dynamics Across Trajectory Subsets") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "grey90"))





# 提取元数据
meta_data <- data@meta.data


# 检查分组类别
table(meta_data$AC_PA)

# 如果存在其他类别或NA值，需要清理数据：
meta_data <- meta_data[meta_data$AC_PA %in% c("atherosclerotic core", "proximal adjacent"), ]

# 定义模块列表
modules <- c("magenta", "red", "turquoise", "greenyellow", "tan", "pink", "yellow"
             ,"brown", "green", "blue",  "purple", "black")

results <- data.frame(
  Module = character(),
  p_value = numeric(),
  mean_AC = numeric(),
  mean_PA = numeric(),
  stringsAsFactors = FALSE
)
# 循环计算每个模块的统计检验
for (module in modules) {
  # 提取模块表达值并移除NA值
  expr_AC <- na.omit(meta_data[meta_data$AC_PA == "atherosclerotic core", module])
  expr_PA <- na.omit(meta_data[meta_data$AC_PA == "proximal adjacent", module])

  # 检查有效样本量（每组至少1个样本）
  if (length(expr_AC) < 1 | length(expr_PA) < 1) {
    warning(paste0("Skipping module ", module, ": PA组或AC组有效样本不足（AC=", length(expr_AC), ", PA=", length(expr_PA), "）"))
    next  # 跳过当前模块
  }

  # 执行Wilcoxon秩和检验
  wilcox_test <- wilcox.test(expr_AC, expr_PA)

  # 存储结果
  results <- rbind(results, data.frame(
    Module = module,
    p_value = wilcox_test$p.value,
    mean_AC = mean(expr_AC),
    mean_PA = mean(expr_PA)
  ))
}

# 输出结果
print(results)
results$p_adj <- p.adjust(results$p_value, method = "BH")
print(results[order(results$p_adj), ])

Module_Colors <- c(
  "magenta" = "magenta",
  "red" = "red",
  "turquoise" = "turquoise",
  "greenyellow" = "greenyellow",
  "tan" = "tan",
  "pink" = "pink",
  "yellow" = "yellow",
  "brown" = "brown",
  "green" = "green",
  "blue" = "blue",
  "purple" = "purple",
  "black" = "black"
)

# 绘制火山图
ggplot(results, aes(x = mean_AC - mean_PA, y = log_p, color = Module)) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Module Differential Expression",
       x = "Effect Size (AC - PA)",
       y = "-log10(Adjusted p-value)") +
  theme_classic() +
  scale_color_manual(values = Module_Colors)




######双起点
selected_genes <- unique(c(module_gene_list$brown, module_gene_list$red,
                           module_gene_list$tan,module_gene_list$magenta,module_gene_list$blue
))

selected_genes <- unique(c(module_gene_list$brown, module_gene_list$red,
                           module_gene_list$tan,module_gene_list$blue
))




###选择不同分支的细胞 #####
main_umap <- reducedDims(cds_MM_foam)$UMAP  # 主图UMAP坐标
main_pseudotime <- pseudotime(cds_MM_foam)  # 主图伪时间

cdsMM_subset1 <- choose_graph_segments(cds_MM_foam)#疾病_foam
subset1_cells <- colnames(cdsMM_subset1)
reducedDims(cdsMM_subset1)$UMAP <- main_umap[subset1_cells, ]
colData(cdsMM_subset1)$pseudotime <- main_pseudotime[subset1_cells]


cdsMM_subset2 <- choose_graph_segments(cds_MM_foam)#健康_marcophage
subset2_cells <- colnames(cdsMM_subset2)
reducedDims(cdsMM_subset2)$UMAP <- main_umap[subset2_cells, ]
cdsMM_subset2 <- cluster_cells(cdsMM_subset2, resolution = 1e-5)
cdsMM_subset2 <- learn_graph(
  cdsMM_subset2,
  close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =3 ,  # 增加该值以修剪短分支（数值需根据数据调整）
    prune_graph = TRUE         # 启用自动修剪分支
  )
)

p <- plot_cells(
  cdsMM_subset2,
  color_cells_by = "Celltype_raw1",  # 按分组着色
  cell_size = 1.5,
  label_groups_by_cluster = FALSE,
  group_label_size = 4,
  show_trajectory_graph = TRUE
) +
  facet_wrap(~AC_PA, nrow = 1) +  # 横向分面（左右布局）
  theme(
    strip.text = element_text(size = 12),  # 分面标题字体
    strip.background = element_blank()     # 分面标题背景透明
  )
print(p)
colData(cdsMM_subset2)$pseudotime <- main_pseudotime[subset2_cells]



trace_genes2 <- graph_test(cdsMM_subset2,
                           neighbor_graph = "principal_graph",
                           cores = 4)
write.csv(trace_genes,"/public3/DSC/single_cell/GSE159677/monocle3/monocle_MM/trace_genes2.csv")
