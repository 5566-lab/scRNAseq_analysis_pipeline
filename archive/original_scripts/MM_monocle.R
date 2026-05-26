library(hdWGCNA)
library(monocle3)
library(tidydr)

# Output paths
output_root <- "/public3/DSC/single_cell/Result/openclaw"
output_monocle <- file.path(output_root, "monocle3")
output_hdWGCNA <- file.path(output_root, "hdWGCNA")

data <- readRDS(file.path(output_hdWGCNA, "data.rds"))
expression_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")
cell_metadata <- data@meta.data

gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  module = "Default_Module",  
  row.names = rownames(expression_matrix)  
)




cell_ids <- rownames(cell_metadata)

sub_matrix <- expression_matrix

cds_MM<- new_cell_data_set(
  sub_matrix,
  cell_metadata = cell_metadata[cell_ids, ],
  gene_metadata = gene_metadata
)

cds_MM <- preprocess_cds(cds_MM, num_dim=20)
cds_MM <- reduce_dimension(
  cds_MM,
  reduction_method = "UMAP",
  preprocess_method = "PCA",        
  umap.n_neighbors = 30,            
  umap.min_dist = 0.3              
)
#提取 Seurat 的 UMAP 坐标（确保细胞顺序一致）
seurat_umap <- Embeddings(data@reductions$umap)

# 替换 Monocle3 的 UMAP 坐标
cds_MM@int_colData$reducedDims$UMAP <- seurat_umap[colnames(cds_MM), ]
cds_MM <- cluster_cells(cds_MM, resolution = 1e-5)
cds_MM <- learn_graph(
  cds_MM,
  close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =5,  # 增加该值以修剪短分支（数值需根据数据调整）
    prune_graph = TRUE         # 启用自动修剪分支
  )
)



modules <- GetModules(data)

# 创建颜色列表（排除未分组的grey模块）
module_colors <- setdiff(unique(modules$color), "grey")

# 生成模块基因列表
module_gene_list <- lapply(module_colors, function(col){
  modules %>% 
    dplyr::filter(color == col) %>% 
    .$gene_name
})
names(module_gene_list) <- module_colors



#####单起点

selected_genes <- unique(c(module_gene_list$red,module_gene_list$magenta,module_gene_list$yellow,
                           module_gene_list$pink,module_gene_list$brown,module_gene_list$black
                           ,module_gene_list$turquoise #APOBEC3A还行
                           
))

selected_genes <- unique(c(module_gene_list$red,module_gene_list$magenta,module_gene_list$yellow,
                           module_gene_list$green,module_gene_list$brown,module_gene_list$black
                           ,module_gene_list$turquoise#不行
))

selected_genes <- unique(c(module_gene_list$red,module_gene_list$magenta,module_gene_list$yellow,
                           module_gene_list$black,module_gene_list$turquoise,module_gene_list$brown
                           #不行
))

selected_genes <- unique(c(module_gene_list$red,module_gene_list$magenta,module_gene_list$yellow,
                           module_gene_list$turquoise,module_gene_list$brown
))#APOBEC3A不好看


                        
#####双起点
selected_genes <- unique(c(module_gene_list$magenta,module_gene_list$yellow,module_gene_list$black,
                           module_gene_list$turquoise,module_gene_list$brown,module_gene_list$purple
                           
))#APOBEC3A可以

selected_genes <- unique(c(module_gene_list$magenta,module_gene_list$yellow,module_gene_list$black,
                           module_gene_list$turquoise,module_gene_list$brown,module_gene_list$green
                           
))


# 检查 selected_genes 是否在 fData(cds_MM)$gene_short_name 中
selected_genes_in_data <- selected_genes[selected_genes %in% fData(cds_MM)$gene_short_name]
print(paste("Number of selected genes in data:", length(selected_genes_in_data)))

# 如果 selected_genes_in_data 为空，停止分析
if (length(selected_genes_in_data) == 0) {
  stop("No selected genes found in the dataset.")
}


cds_MM <-preprocess_cds(cds_MM, num_dim = 50) 
plot_pc_variance_explained(cds_MM)              # 选择拐点前的维度
cds_MM <- preprocess_cds(cds_MM, num_dim=20, use_genes = selected_genes_in_data)
#cds_MM <- preprocess_cds(cds_MM, num_dim=20)
cds_MM <- reduce_dimension(cds_MM, reduction_method="PCA")
harmony_emb <- HarmonyMatrix(
  data_mat = reducedDims(cds_MM)$PCA,
  meta_data = colData(cds_MM),
  vars_use = "orig.ident",  # 指定批次字段
  do_pca = FALSE
)

# 替换原有PCA坐标
reducedDims(cds_MM)$PCA <- harmony_emb
# 降维（UMAP/PCA/t-SNE）
cds_MM <- reduce_dimension(cds_MM, reduction_method = "UMAP",
                           preprocess_method = "PCA",       
                           umap.n_neighbors = 30,     
                           umap.min_dist = 0.3     )
plot_cells(cds_MM, color_cells_by="orig.ident", group_label_size=5)
cds_MM <- cluster_cells(cds_MM, resolution = 1e-5)
cds_MM <- learn_graph(
  cds_MM,
  close_loop = FALSE,          # 禁用闭环，防止成环
  learn_graph_control = list(
    minimal_branch_len =2,  # 增加该值以修剪短分支（数值需根据数据调整）
    prune_graph = TRUE         # 启用自动修剪分支
  )
)
p <- plot_cells(
  cds_MM,
  color_cells_by = "Celltype_raw1",  # 按分组着色
  cell_size = 1.0,
  label_groups_by_cluster = FALSE,
  group_label_size = 4,
  show_trajectory_graph = TRUE
) + 
  # 上下翻转（Y轴反向）
  scale_y_reverse() +
  # 左右翻转（X轴反向）
  scale_x_reverse() +
  facet_wrap(~AC_PA, nrow = 1) +  # 横向分面（左右布局）
  theme(
    strip.text = element_text(size = 12),  # 分面标题字体
    strip.background = element_blank()     # 分面标题背景透明
  )
print(p)
ggsave(p,filename = file.path(output_monocle, "F2.5_monocle_celltype_cluster.png"),width = 18,height = 8)
ggsave(p,filename = file.path(output_monocle, "F2.5_monocle_celltype_cluster.pdf"),width = 18,height = 8)



######拆分partition####
#######################
partition1_indices <- cds_MM@clusters$UMAP$partitions == "1"
cds_MM_foam  <- cds_MM[, partition1_indices]
cds_MM_foam <- order_cells(cds_MM_foam)

plot_PT <- plot_cells(
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

plot_PT
ggsave(plot_PT,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.5_MM_Pseudotim.png',width = 18,height = 8)
ggsave(plot_PT,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.5_MM_Pseudotim.pdf',width = 18,height = 8)


###添加伪时序信息到 Seurat 中 #####
data$pseudotime <- pseudotime(cds_MM_foam)
summary(data$pseudotime)
saveRDS(data,"/public3/DSC/single_cell/Result/openclaw/monocle3/data_pseudotime.rds")

monocle_Pseudotim  <- PlotModuleTrajectory(
  data,
  pseudotime_col = 'pseudotime'
)

monocle_Pseudotim
ggsave(monocle_Pseudotim,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.6_MM_monocle_Pseudotim.png',width = 18,height = 8)
ggsave(monocle_Pseudotim,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.6_MM_monocle_Pseudotim.pdf',width = 18,height = 8)


#trace_genes<- graph_test(cds_MM_foam, 
#              neighbor_graph = "principal_graph", 
#              cores = 8)
#sorted_res <- trace_genes %>% 
#  arrange(desc(morans_I))


pseudotime_values <- pseudotime(cds_MM_foam)
selected_cells <- colnames(cds_MM_foam)[pseudotime_values >= 0 & pseudotime_values <= 9]
cds_subset <- cds_MM_foam[, selected_cells]

selected_cells <- colnames(cds_subset)
colData(cds_MM_foam)$selected <- colnames(cds_MM_foam) %in% selected_cells
cds_MM_foam@int_metadata$subset_info <- list(
  selected_cells = selected_cells,
  selection_time = Sys.time()
)
p_selected_cells <- plot_cells(cds_MM_foam,
           color_cells_by = "selected",  # 引用colData中的逻辑标记列
           cell_size = 1.0,             # 放大选中细胞点
           group_label_size = 4,
           alpha = ifelse(colData(cds_MM_foam)$selected, 0.8, 0.3),         # 未选中细胞透明度0.3，选中0.8
           trajectory_graph_segment_size = 0.5) +
  scale_color_manual(
    name = "Selected",
    values = c("gray80", "red"),        # 自定义颜色
    labels = c("Unselected", "Selected")
  ) +
  theme(legend.position = "right")
p_selected_cells
ggsave(p_selected_cells,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.7_p_selected_cells.png',width = 18,height = 8)
ggsave(p_selected_cells,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.7_p_selected_cells.pdf',width = 18,height = 8)



subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
sorted_res <- subset_pr_test_res %>% 
  arrange(desc(morans_I))
apobec_row <- sorted_res %>% 
  filter(gene_short_name == "APOBEC3A") 
print(paste("APOBEC3A 的行数为：", which(sorted_res$gene_short_name == "APOBEC3A")))
subset_pr_test_res <-  subset(subset_pr_test_res, 
                              q_value < 0.01)
write.csv(subset_pr_test_res,"/public3/DSC/single_cell/Result/openclaw/monocle3/subset_genes42.csv")



subset_pr_test_res <- read.csv("/public3/DSC/single_cell/Result/openclaw/monocle3/subset_genes40.csv", row.names = 1)
###选择不同分支的细胞 #####

cdsMM_subset1 <- choose_graph_segments(cds_MM_foam)#疾病_foam
cdsMM_subset2 <- choose_graph_segments(cds_MM_foam)#健康_marcophage
saveRDS(cdsMM_subset1, file = "/public3/DSC/single_cell/Result/openclaw/monocle3/cdsMM_subset1_disease_foam.rds")
saveRDS(cdsMM_subset2, file = "/public3/DSC/single_cell/Result/openclaw/monocle3/cdsMM_subset2_healthy_macrophage.rds")

#cdsMM_subset1 <- readRDS("/public3/DSC/single_cell/Result/openclaw/monocle3/cdsMM_subset1_disease_foam.rds")
#cdsMM_subset2 <- readRDS("/public3/DSC/single_cell/Result/openclaw/monocle3/cdsMM_subset2_healthy_macrophage.rds")

subset_list <- list(
  subset1 = cdsMM_subset1,
  subset2 = cdsMM_subset2
  
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
saveRDS(data,"/public3/DSC/single_cell/Result/openclaw/monocle3/data_pseudotime.rds")


#####AC_PA在不同细胞分支上随时间的分布变化#####

p <- ggplot(data@meta.data, aes(x = pseudotime, fill = AC_PA)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Pseudotime by AC_PA") +
  theme_classic()

p <- ggplot(data@meta.data, aes(x = pseudotime, fill = cdsMM_sub1)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Pseudotime by cdsMM_sub1") +
  theme_classic()
p <- ggplot(data@meta.data, aes(x = pseudotime, fill = cdsMM_sub2)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Pseudotime by cdsMM_sub1") +
  theme_classic()
ggsave(p,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/S2.1_Pseudotime_AC_PA.png',width = 18,height = 8)
ggsave(p,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/S2.1_Pseudotime_AC_PA.pdf',width = 18,height = 8)

p1 <- ggplot() +
  # 绘制 cdsMM_sub1 的密度
  geom_density(
    data = data@meta.data[data@meta.data$cdsMM_sub1 == "Yes", ],
    aes(x = pseudotime, fill = "cdsMM_sub1"),
    alpha = 0.5
  ) +
  # 绘制 cdsMM_sub2 的密度
  geom_density(
    data = data@meta.data[data@meta.data$cdsMM_sub2 == "Yes", ],
    aes(x = pseudotime, fill = "cdsMM_sub2"),
    alpha = 0.5
  ) +
  scale_fill_manual(
    name = "Subset",
    values = c("cdsMM_sub1" = "red", "cdsMM_sub2" = "blue"),
    labels = c("cdsMM_sub1 (Yes)", "cdsMM_sub2 (Yes)")
  ) +
  labs(title = "Pseudotime Density Overlay",
       x = "Pseudotime",
       y = "Density") +
  theme_classic()



p1 <- ggplot(
  data = subset(data@meta.data, cdsMM_sub1 == "Yes"),  # 筛选subset1子集
  aes(x = pseudotime, fill = AC_PA)
) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density of Pseudotime by AC_PA (Cell_fate1 only)",  # 更新标题
    x = "Pseudotime",
    y = "Density"
  ) +
  theme_classic()
p2 <- ggplot(
  data = subset(data@meta.data, cdsMM_sub2 == "Yes"),   # 筛选subset2子集
  aes(x = pseudotime, fill = AC_PA)
) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density of Pseudotime by AC_PA (Cell_fate2 only)",  # 更新标题
    x = "Pseudotime",
    y = "Density"
  ) +
  theme_classic()

####cell_fate在不同细胞分支上随时间的分布变化#####

p <- ggplot(combined_data, aes(x = pseudotime, fill = subset)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Pseudotime by cell_fate") +
  theme_classic()

combined_data$cell_type[combined_data$AC_PA == "atherosclerotic core" & 
                          combined_data$subset == "subset1"] <- "cell_type1"
combined_data$cell_type[combined_data$AC_PA == "atherosclerotic core" & 
                          combined_data$subset == "subset2"] <- "cell_type2"
combined_data$cell_type[combined_data$AC_PA == "proximal adjacent" & 
                          combined_data$subset == "subset1"] <- "cell_type3"
combined_data$cell_type[combined_data$AC_PA == "proximal adjacent" & 
                          combined_data$subset == "subset2"] <- "cell_type4"
plot_data1 <- combined_data %>%
  filter(cell_type %in% c("cell_type1", "cell_type2")) %>%
  mutate(cell_type = factor(cell_type)) 

p1 <- ggplot(plot_data1, aes(x = pseudotime, fill = cell_type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Pseudotime by cell_fate") +
  theme_classic()
plot_data2 <- combined_data %>%
  filter(cell_type %in% c("cell_type4")) %>%
  mutate(cell_type = factor(cell_type)) 
p2 <- ggplot(plot_data2, aes(x = pseudotime, fill = cell_type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Pseudotime by cell_fate") +
  theme_classic()


####subset genes#####
# 定义基因列表
gene_list <- c("CD163", 'IL10',"APOBEC3A",'CCL18','CCL22',
               'TNF', "IL6",'IL1B','CXCL10','CD80'
)#clone13敲除后变化趋势相同(上升)：TNF、IL1B
#clone37敲除后变化趋势应该相同
gene_list <- c("NFKB1", 'RELA',"MSR1",
               "LIPA",  "CD36", "APOC1", "CD9", "TREM2",
               "OLR1", "PLIN2", "MARCO", "IL1RN","CCL2","FABP5", "CTSB", "SPP1"
)#敲除后变化趋势应该相同的gene：CD36、OLR1、PLIN2、MARCO、IL1RN、CCL2、FABP5、SPP1
#clone37敲除后变化趋势应该相同的gene：
#clone13敲除后变化趋势应该相同的gene：OLR1,IL1RN,CCL2

gene_list <- c("FCN1", 'VCAN',"APOBEC3A",'OLR1','TIMP1','FN1',
               'CD52', "EREG",'PLIN2','MARCO','C1QA',
               "PLTP", 'SELENOP','IGSF21','FOLR2',
              "SLCO2B1", "FILIP1L", "HSP90AA1","IGF1","BHLHE41")
#clone13敲除后变化趋势应该相同的gene：FCN1、VCAN、FN1、EREG、MARCO、FOLR2、FILIP1L
#clone37敲除后变化趋势应该相同的gene：OLR1、TIMP1、FN1、EREG、SELENOP、FOLR2、FILIP1L
#相同的gene:FOLR2、FILIP1L、FN1、EREG
# 绘制分面图
p <- plot_cells(
  cds_MM_foam,
  genes = gene_list,
  label_groups_by_cluster = FALSE,
  cell_size = 1.0,
  group_label_size = 4,
  show_trajectory_graph = TRUE
) + 
  facet_wrap(~feature_label, nrow = 4) +  
  theme(
    strip.text = element_text(size = 12),  
    strip.background = element_blank()     
  ) + 
  scale_color_gradient(
    low = "gray90",            
    high = "red3",             
    #breaks = c(0, 2, 5),       
    #limits = c(0, 5)           
  ) +
  labs(color = "Expression") + 
  theme_dr() + 
  theme(
    panel.grid = element_blank(), 
    plot.title = element_blank()
  ) + 
  NoLegend()
print(p)
ggsave(p,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/S2.2_Pseudotime_plot_cells.png',width = 20,height = 16)
ggsave(p,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/S2.2_Pseudotime_plot_cells.pdf',width = 20,height = 16)



#####Pseudotime_Expression###### 
combined_data <- data.frame()
# 循环处理每个子集
for (i in seq_along(subset_list)) {
  # 提取当前子集名称对应的元数据列
  subset_col <- paste0("cdsMM_sub", i)
  
  # 筛选属于当前子集的细胞
  subset_cells <- data@meta.data %>% 
    filter(!!sym(subset_col) == "Yes") %>% 
    rownames()
  
  # 初始化一个数据框，存储当前子集的伪时间和其他信息
  subset_df <- data.frame(
    pseudotime = data@meta.data[subset_cells, "pseudotime"],
    subset = names(subset_list)[i],
    AC_PA = data@meta.data[subset_cells, "AC_PA"]
  )
  
  # 循环处理 gene_list 中的每个基因，提取表达数据
  for (gene in gene_list) {
    subset_df[[gene]] <- GetAssayData(data, assay = "RNA", slot = "data")[gene, subset_cells]
  }
  
  # 将当前子集的数据合并到 combined_data 中
  combined_data <- rbind(combined_data, subset_df)
}

combined_data_long <- combined_data %>%
  pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression")
combined_data_long$gene <- factor(combined_data_long$gene, levels = gene_list)

# 绘制联合轨迹图
# 定义颜色调色板（根据子集数量扩展）
subset_colors <- c( "#FF7F0E","#1F77B4", "#2CA02C", "#D62728", "#9467BD")[1:length(subset_list)]
gene_expression_plot <- ggplot(combined_data_long, aes(x = pseudotime, y = expression, color = subset)) +
  geom_smooth(
    method = "gam",   # 采用GAM模型适应复杂轨迹
    formula = y ~ s(x, bs = "tp"), 
    se = FALSE, 
    linewidth = 1.5,
    alpha = 0.8
  ) +
  scale_color_manual(values = subset_colors, labels = c("FOAM_cells", "Macrophage")) +
  labs(x = "Pseudotime", y = "Expression (log-normalized)", 
       title = "Gene Expression Dynamics Across Trajectory Subsets") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "grey90")) +
  facet_wrap(~ gene, scales = "free_y")  # 按基因分面，y轴独立

# 显示图形
gene_expression_plot
ggsave(gene_expression_plot,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/S2.1_Pseudotime_Expression.png',width = 20,height = 16)
ggsave(gene_expression_plot,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/S2.1_Pseudotime_Expression.pdf',width = 20,height = 16)






#####APOBEC3A####

p <- plot_cells(
  cds_MM,
  genes = "APOBEC3A",
  label_groups_by_cluster = FALSE,
  cell_size = 1.0,
  group_label_size = 4,
  show_trajectory_graph = TRUE
)  + 
  # 上下翻转（Y轴反向）
  scale_y_reverse() +
  # 左右翻转（X轴反向）
  scale_x_reverse() +
  facet_wrap(~AC_PA, ncol = 1) +  # 横向分面（左右布局）
  theme(
    strip.text = element_text(size = 12),  # 分面标题字体
    strip.background = element_blank()     # 分面标题背景透明
  )+ 
  scale_color_gradient(
    low = "gray90",            # 低表达设为浅灰色
    high = "red3",             # 高表达设为深红色
    breaks = c(0, 0.1, 1),       # 自定义颜色断点（根据实际表达范围调整）
    #limits = c(0, 5)           # 限制颜色映射范围（避免极端值影响对比度）
  ) +
  labs(color = "APOBEC3A Expression")+ 
  theme_dr() + theme(panel.grid=element_blank(), 
                     plot.title = element_blank()) + NoLegend() 
print(p)
ggsave(p,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.5_monocle_APOBEC3A.png',width = 18,height = 8)
ggsave(p,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F2.5_monocle_APOBEC3A.pdf',width = 18,height = 8)


# 创建合并数据框架
combined_data <- data.frame()


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
    subset = names(subset_list)[i],
    AC_PA = data@meta.data[subset_cells, "AC_PA"]
  )
  
  combined_data <- rbind(combined_data, subset_df)
}





ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = subset)) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~AC_PA, ncol = 2) +  # 修正分面变量为AC_PA
  scale_color_manual(
    values = c("#FF7F0E","#1F77B4",  "#2CA02C")[1:length(subset_list)],  # 添加逗号
    labels = names(subset_list)  # 标签与subset对应
  ) +
  labs(title = "APOBEC3A Expression Dynamics by AC/PA Status") +
  theme_minimal()

# 绘制联合轨迹图
APOBEC3A_EXP_TIME <- ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = subset)) +
  #geom_point(alpha = 0.6, size = 1.2) +  # 散点显示细胞分布
  geom_smooth(
    method = "gam",   # 采用GAM模型适应复杂轨迹
    formula = y ~ s(x, bs = "tp"), 
    se = FALSE, 
    linewidth = 1.5,
    alpha = 0.8
  ) +
  # geom_smooth(method = "loess", se = FALSE, linewidth = 1.5) +  # 趋势线
  scale_color_manual(values = subset_colors,labels = c( "cell_fate1", "cell_fate2") ) +
  labs(x = "Pseudotime", y = "APOBEC3A Expression (log-normalized)", 
       title = "APOBEC3A Expression Dynamics Across Trajectory Subsets") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "grey90"))

APOBEC3A_EXP_TIME
ggsave(APOBEC3A_EXP_TIME,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F3.1_APOBEC3A_EXP_TIME.png',width = 18,height = 8)
ggsave(APOBEC3A_EXP_TIME,filename = '/public3/DSC/single_cell/Result/openclaw/monocle3/F3.1_APOBEC3A_EXP_TIME.pdf',width = 18,height = 8)

trace_genes <- graph_test(cds_MM_foam, 
                          neighbor_graph = "principal_graph", 
                          cores = 4)
write.csv(trace_genes,"/public3/DSC/single_cell/Result/openclaw/monocle3/trace_genes.csv")
