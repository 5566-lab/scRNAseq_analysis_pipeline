library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(harmony)
library(plyr)
library(tidydr)
library(SingleR)
library(pheatmap)
library(celldex)
library(SeuratObject)
library(SeuratData)
#library(monocle)
library(monocle3)
library(MAST)
library(ggbump)
packageVersion("Seurat") #5.0.1

packageVersion("dbplyr") #需要为2.3.4
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
library(ggpubr)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
#devtools::install_github("NightingaleHealth/ggforestplot")
#devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)

#GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)
library(GOplot)
library(stringr)
library(showtext)


font_add(
  family = "Arial",
  regular = "/home/dengsc/.fonts/truetype/msttcorefonts/ARIAL.TTF",
  bold    = "/home/dengsc/.fonts/truetype/msttcorefonts/ARIALBD.TTF", 
  italic  = "/home/dengsc/.fonts/truetype/msttcorefonts/ARIALI.TTF",
  bolditalic = "/home/dengsc/.fonts/truetype/msttcorefonts/ARIALBI.TTF"
)

# 启动字体渲染引擎
showtext_auto()


# ==============================================================================
# INITIALIZATION: Initialize the master list for all Seurat objects
# ==============================================================================
cat("\n>>> Initializing master Seurat object list...\n")
all_individual_seurat_objects <- list()

# ==============================================================================
# PART 1: 处理 GSE260657 (Carotid - Multiple Files)
# ==============================================================================
cat("\n>>> PART 1: Starting processing for GSE260657 (Multiple Files)...\n")

data_dir_260 <- "/public3/DSC/single_cell/GSE260657_Carotid"
file_list_260 <- list.files(path = data_dir_260, pattern = "\\.txt\\.gz$", full.names = TRUE)

# 1.1 初始化临时列表
temp_list_260 <- list()

# 1.2 循环读取
for (file_path in file_list_260) {
  tryCatch({
    # 解析文件名
    file_name <- basename(file_path) # 循环变量：file_name
    patient_num_str <- sub(".*human_([0-9]+)\\.txt\\.gz", "\\1", file_name)
    patient_num <- as.numeric(patient_num_str)
    project_name <- paste0("GSE260657_Human_", patient_num)
    
    cat(paste0("  Loading ", project_name, "...\n"))
    
    # 读取数据
    raw_data <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
    
    # 定义分组 (Human 1-7: Stable, 8-15: Unstable)
    if (!is.na(patient_num) && patient_num <= 7) {
      plaque <- "Carotid Stable Plaque"
    } else {
      plaque <- "Carotid Unstable Plaque"
    }
    
    # 创建对象
    seurat_obj <- CreateSeuratObject(counts = raw_data, project = project_name, min.cells = 3, min.features = 200)
    seurat_obj$Patient_ID <- paste0("Human_", patient_num)
    seurat_obj$AC_PA <- plaque
    seurat_obj$Source_GSE <- "GSE260657"
    
    # 加入临时列表
    temp_list_260[[project_name]] <- seurat_obj
    
  }, error = function(e) {
    cat(paste("  Error processing file", file_path, ":", e$message, "\n"))
  })
}

# 1.3 合并 GSE260657 为单个对象并加入总列表
if (length(temp_list_260) > 0) {
  cat("  Merging GSE260657 samples into one object...\n")
  GSE260657_combined <- merge(
    x = temp_list_260[[1]],
    y = temp_list_260[-1],
    add.cell.ids = names(temp_list_260), # 防止Barcode冲突
    project = "GSE260657_Combined"
  )
  all_individual_seurat_objects[["GSE260657"]] <- GSE260657_combined
  rm(temp_list_260); gc() # 清理内存
} else {
  cat("Warning: No samples loaded for GSE260657.\n")
}


# ==============================================================================
# PART 2: 处理 GSE247238 (Carotid - Matrix Files + GEO Metadata)
# ==============================================================================
cat("\n>>> PART 2: Starting processing for GSE247238 (Matrix Files).\n")

# 2.1 定义临床信息表
patient_db <- list(
  "1" = list(type="Carotid Stable Plaque", age=70, sex="Male"),
  "2" = list(type="Carotid Stable Plaque", age=72, sex="Female"),
  "3" = list(type="Carotid Unstable Plaque", age=74, sex="Male"),
  "4" = list(type="Carotid Unstable Plaque", age=81, sex="Male"),
  "5" = list(type="Carotid Stable Plaque", age=71, sex="Male"),
  "6" = list(type="Carotid Unstable Plaque", age=80, sex="Male"),
  "7" = list(type="Carotid Stable Plaque", age=78, sex="Male"),
  "8" = list(type="Carotid Stable Plaque", age=69, sex="Male"),
  "9" = list(type="Carotid Unstable Plaque", age=82, sex="Male"),
  "10" = list(type="Carotid Stable Plaque", age=71, sex="Male")
)

# 2.2 获取 GEO 元数据
cat("  Fetching GEO metadata for GSE247238...\n")
gse_info <- tryCatch(
  getGEO("GSE247238", GSEMatrix = TRUE, getGPL = FALSE),
  error = function(e) { cat("  Error fetching GEO data:", e$message, "\n"); NULL }
)

gsm_to_title <- c()
if (!is.null(gse_info)) {
  gsm_to_title <- Biobase::pData(gse_info[[1]])$title
  names(gsm_to_title) <- Biobase::pData(gse_info[[1]])$geo_accession
}

# 2.3 初始化临时列表
temp_list_247 <- list()
data_dir_247 <- "/public3/DSC/single_cell/GSE247238_Carotid"
matrix_files <- list.files(data_dir_247, pattern = "_matrix.mtx.gz$", full.names = FALSE)
sample_ids <- str_remove(matrix_files, "_matrix.mtx.gz")

# 2.4 循环读取
for (sample in sample_ids) {
  tryCatch({
    gsm_id <- unlist(strsplit(sample, "_"))[1]
    cat(paste0("  Loading ", sample, " (GSM: ", gsm_id, ")...\n"))
    
    # 读取数据
    mtx_path   <- file.path(data_dir_247, paste0(sample, "_matrix.mtx.gz"))
    barcodes_path <- file.path(data_dir_247, paste0(sample, "_barcodes.tsv.gz"))
    features_path <- file.path(data_dir_247, paste0(sample, "_features.tsv.gz"))
    
    counts <- ReadMtx(mtx = mtx_path, cells = barcodes_path, features = features_path, feature.column = 1)
    sobj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
    
    # 匹配临床信息
    sobj$GSM_ID <- gsm_id
    sobj$Source_GSE <- "GSE247238"
    
    # 默认值
    car_id <- "Unknown"; plaque <- "Unknown"; age <- NA; sex <- "Unknown"
    
    if (gsm_id %in% names(gsm_to_title)) {
      title_str <- gsm_to_title[[gsm_id]]
      car_match <- str_extract(title_str, "CAR[0-9]+")
      if (!is.na(car_match)) {
  car_id <- str_remove(car_match, "CAR")
  if (car_id %in% names(patient_db)) {
    info <- patient_db[[car_id]]
    plaque <- info$type
    age <- info$age
    sex <- info$sex
  }
      }
    }
    
    sobj$Patient_ID <- paste0("Patient_", car_id)
    sobj$AC_PA <- plaque
    sobj$Age <- age
    sobj$Sex <- sex
    
    temp_list_247[[sample]] <- sobj
    
  }, error = function(e) {
    cat(paste("  Error processing", sample, ":", e$message, "\n"))
  })
}

# 2.5 合并 GSE247238 并加入总列表
if (length(temp_list_247) > 0) {
  cat("  Merging GSE247238 samples into one object...\n")
  GSE247238_combined <- merge(
    x = temp_list_247[[1]],
    y = temp_list_247[-1],
    add.cell.ids = names(temp_list_247),
    project = "GSE247238_Combined"
  )
  all_individual_seurat_objects[["GSE247238"]] <- GSE247238_combined
  rm(temp_list_247); gc()
} else {
  cat("Warning: No samples loaded for GSE247238.\n")
}


# ==============================================================================
# PART 3: 处理 GSE131778 (Coronary - Single File)
# ==============================================================================
cat("\n>>> PART 3: Processing GSE131778_AC (Single file)...\n")

GSE131778_combined <- NULL # 初始化变量

tryCatch({
  file_path_131 <- "/public3/DSC/single_cell/GSE131778_Coronary_AC/GSE131778_human_coronary_scRNAseq.txt"
  
  cat("  Reading GSE131778 data table...\n")
  raw_data_131 <- read.table(file_path_131, header = TRUE, row.names = 1, sep = "\t")
  
  # 创建 Seurat 对象
  seurat_obj_131 <- CreateSeuratObject(
    counts = raw_data_131,
    project = "GSE131778",
    min.cells = 3,
    min.features = 200
  )
  
  # 添加元数据 (Metadata)
  seurat_obj_131$Source_GSE <- "GSE131778"
  seurat_obj_131$AC_PA <- "Coronary Atherosclerotic Core" # 指定元数据
  seurat_obj_131$Patient_ID <- "GSE131778_Sample1"    # 只有一个样本
  
  # 因为是单文件，它本身就是 Combined 对象
  GSE131778_combined <- seurat_obj_131
  all_individual_seurat_objects[["GSE131778"]] <- GSE131778_combined
  
  cat("  Successfully processed GSE131778_AC. Cells:", ncol(GSE131778_combined), "\n")
  
}, error = function(e) {
  cat(paste("Warning: Error processing GSE131778_AC:", e$message, "\n"))
})


# ==============================================================================
# PART 4: 处理 GSE210152 (RDS File)
# ==============================================================================
cat("\n>>> PART 4: Processing GSE210152_AC (RDS file)...\n")

GSE210152_combined <- NULL

tryCatch({
  file_path_210 <- "/public3/DSC/single_cell/GSE210152_Carotid_AC/GSE210152_raw.RDS"
  cat("  Reading RDS file...\n")
  seurat_obj_210 <- readRDS(file_path_210)
  
  # 设置元数据
  seurat_obj_210$AC_PA <- "Carotid Atherosclerotic Core"
  seurat_obj_210$Source_GSE <- "GSE210152"
  # 检查是否已有 Patient_ID，否则创建默认ID
  if (!"Patient_ID" %in% colnames(seurat_obj_210@meta.data)) {
    seurat_obj_210$Patient_ID <- "GSE210152_Sample1"
  }
  
  # 设置项目名称
  Project(seurat_obj_210) <- "GSE210152"
  
  GSE210152_combined <- seurat_obj_210
  all_individual_seurat_objects[["GSE210152"]] <- GSE210152_combined
  
  cat("  Successfully processed GSE210152_AC. Cells:", ncol(GSE210152_combined), "\n")
  
}, error = function(e) {
  cat(paste("Warning: Error processing GSE210152_AC:", e$message, "\n"))
})


# ==============================================================================
# PART 5: 处理 GSE155468 (AscAorta Control - File-based)
# ==============================================================================
cat("\n>>> PART 5: Processing GSE155468_AscAorta_Control_PA (File-based)...\n")
data_dir_155 <- "/public3/DSC/single_cell/GSE155468_AscAorta_Control_PA/"
file_names_155 <- c("GSM4704931_Con4.txt.gz", "GSM4704932_Con6.txt.gz", "GSM4704933_Con9.txt.gz")

if (dir.exists(data_dir_155)) {
  for (file_name in file_names_155) { # 循环变量：file_name
    file_path <- file.path(data_dir_155, file_name)
    # 提取 GSM 和 ConID 作为样本 ID
    sample_id <- gsub(".txt.gz", "", file_name) 
    unique_sample_name <- paste0("GSE155468_", sample_id) 
    
    cat(paste("  Reading file:", file_name, "for sample:", unique_sample_name, "\n"))
    
    if (file.exists(file_path)) {
      tryCatch({
  raw_data <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
  # 创建 Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = raw_data,
    project = "GSE155468", 
    min.cells = 3,
    min.features = 200
  )
  # 添加元数据
  seurat_obj$Source_GSE <- "GSE155468"
  # 设置为升主动脉对照
  seurat_obj$AC_PA <- "Control AscAorta" 
  seurat_obj$Patient_ID <- sample_id
  
  # 将对象加入总列表
  all_individual_seurat_objects[[unique_sample_name]] <- seurat_obj
  cat(paste("  Successfully processed file:", file_name, ". Cells:", ncol(seurat_obj), "\n"))
  
      }, error = function(e) {
  cat(paste("  Error processing file:", file_name, ":", e$message, "\n"))
      })
    } else {
      cat(paste("  Warning: File not found, skipping:", file_path, "\n"))
    }
  }
} else {
  cat(paste("Warning: Directory not found, skipping:", data_dir_155, "\n"))
}


# ==============================================================================
# PART 6: 处理 GSE159677 (Carotid AC/PA - 10X Format)
# ==============================================================================
cat("\n>>> PART 6: Starting processing for GSE159677 (Carotid AC/PA)...\n")

gse_root_dir_159 <- '/public3/DSC/single_cell/GSE159677_Carotid_MainProject/'
patient_ids <- 1:3
sample_types <- c("AC", "PA")

for (patient in patient_ids) {
  for (type in sample_types) {
    # 构造目录名，例如 Patient_1_AC
    sample_dir_name <- paste0("Patient_", patient, "_", type)
    dir_path <- file.path(gse_root_dir_159, sample_dir_name)
    
    desired_orig_ident <- paste0("patient", patient, type)
    
    cat(paste0("  Loading ", desired_orig_ident, " from ", dir_path, "...\n"))
    
    if (dir.exists(dir_path)) {
      tryCatch({
  # 1. Read 10X data
  # FIXED: Was dir.path, changed to dir_path
  counts_data <- Read10X(dir_path)
  
  # 2. Define metadata
  if (type == "AC") {
    plaque_type <- "Carotid Atherosclerotic Core"
  } else {
    plaque_type <- "Carotid Proximal Adjacent"
  }
  
  # 3. Create Seurat Object
  sobj <- CreateSeuratObject(
    counts = counts_data, 
    project = desired_orig_ident, 
    min.cells = 3, 
    min.features = 200
  )
  
  # 4. Add Metadata
  sobj$Source_GSE <- "GSE159677"
  sobj$Patient_ID <- paste0("Patient_", patient)
  sobj$Sample_Type <- type
  sobj$AC_PA <- plaque_type
  
  # 5. Add individual object to master list
  all_individual_seurat_objects[[desired_orig_ident]] <- sobj
  cat(paste("  Successfully processed sample:", desired_orig_ident, ". Cells:", ncol(sobj), "\n"))
  
      }, error = function(e) {
  cat(paste("  Error processing sample", desired_orig_ident, ":", e$message, "\n"))
      })
    } else {
      cat(paste("  Warning: Directory not found, skipping:", dir_path, "\n"))
    }
  }
}


# ==============================================================================
# PART 7: 处理 Directory-based GSEs (10X Genomics Format)
# ==============================================================================
cat("\n>>> PART 7: Starting processing for Directory-based GSEs (10X Format)...\n")

# --- Function to process Directory-based GSEs ---
process_directory_gse <- function(root_dir, individual_objects_list) {
  
  # 1. 识别GSE名称和设置AC_PA值
  gse_name <- basename(root_dir) 
  
  ac_pa_value <- "Unknown"
  
  # 匹配 Control AscAorta
  if (grepl("GSE213740|GSE216860|GSE155468", gse_name)) {
    ac_pa_value <- "Control AscAorta"
    # 匹配 Coronary Proximal Adjacent
  } else if (grepl("GSE234077|GSE224273|GSE253903", gse_name)) {
    ac_pa_value <- "Coronary Proximal Adjacent"
  } else if (grepl("Control", gse_name)) {
    ac_pa_value <- "Control"
  } else if (grepl("AC|Plaque", gse_name)) {
    ac_pa_value <- "Atherosclerotic Core"
  }
  
  cat(paste("  Processing Directory-based GSE:", gse_name, "with AC_PA:", ac_pa_value, "\n"))
  
  # 2. 查找样本目录
  sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
  
  if (length(sample_dirs) == 0) {
    cat(paste("  Warning: No subdirectories found in", root_dir, ". Assuming root contains the 10X files.\n"))
    sample_dirs <- root_dir
  }
  
  # 3. 循环读取样本
  for (sample_dir in sample_dirs) {
    sample_id <- basename(sample_dir)
    unique_sample_name <- paste0(gse_name, "_", sample_id)
    
    cat(paste("    Reading data for sample:", unique_sample_name, " from", sample_dir, "\n"))
    
    tryCatch({
      data <- Read10X(data.dir = sample_dir)
      
      seurat_obj <- CreateSeuratObject(
  counts = data,
  project = gse_name, 
  min.cells = 3,
  min.features = 200
      )
      seurat_obj$Source_GSE <- gse_name
      seurat_obj$AC_PA <- ac_pa_value
      seurat_obj$Patient_ID <- sample_id 
      
      individual_objects_list[[unique_sample_name]] <- seurat_obj
      cat(paste("    Successfully processed sample:", unique_sample_name, "\n"))
      
    }, error = function(e) {
      cat(paste("    Error processing sample:", unique_sample_name, ":", e$message, "\n"))
    })
  }
  return(individual_objects_list)
}

# --- Process all Directory-based GSEs using the function ---
dir_based_gses_root_dirs <- c(
  "/public3/DSC/single_cell/GSE213740_AscAorta_Control_PA", 
  "/public3/DSC/single_cell/GSE234077_Carotid_AC", 
  "/public3/DSC/single_cell/GSE224273_Carotid_AC", 
  "/public3/DSC/single_cell/GSE253903_Carotid_AC", 
  "/public3/DSC/single_cell/GSE216860_Carotid_PA"     
)

for (root_dir in dir_based_gses_root_dirs) {
  if (dir.exists(root_dir)) {
    all_individual_seurat_objects <- process_directory_gse(root_dir, all_individual_seurat_objects)
  } else {
    cat(paste("Warning: Directory not found, skipping:", root_dir, "\n"))
  }
}


# ==============================================================================
# PART 8: 检查结果
# ==============================================================================
cat("\n=== Processing Complete ===\n")
cat("List 'all_individual_seurat_objects' now contains:", length(all_individual_seurat_objects), "datasets.\n")
print(names(all_individual_seurat_objects))


# ==============================================================================
# PART 9: 完整流程 - 内存优化版 (Memory Optimized)
# ==============================================================================
cat("\n>>> PART 9: 开始最终整合流程 (内存优化模式)...\n")
cat("  检查并修复列表中的特定对象...\n")

if ("GSE210152" %in% names(all_individual_seurat_objects)) {
  cat("  [Fix] 正在修复 GSE210152 格式问题...\n")
  obj_210 <- all_individual_seurat_objects[["GSE210152"]]
  
  # 1. 获取当前 Assay 名字
  current_assay_names <- names(obj_210@assays)
  # 假设第一个就是主 Assay (比如 originalexp 或 spatial)
  old_name <- current_assay_names[1] 
  cat(paste("  当前 Assay 名:", old_name, "\n"))
  
  # 2. 如果不是 RNA，手动创建一个新的 RNA Assay
  if (old_name != "RNA") {
    cat("  [重命名] 检测到 Assay 名不一致，正在手动迁移数据到 'RNA' Assay...\n")
    
    # --- 手动迁移三部曲 ---
    
    # A. 提取原始矩阵 (Counts)
    # 尝试多种方式获取 counts，确保万无一失
    raw_counts <- tryCatch({
      LayerData(obj_210, assay = old_name, layer = "counts")
    }, error = function(e) {
      # 如果 LayerData 失败 (旧版本 Seurat)，尝试 GetAssayData
      GetAssayData(obj_210, assay = old_name, slot = "counts")
    })
    
    # B. 创建全新的 RNA Assay
    # 这会自动设置 Key 为 "rna_", 避免旧 Key (如 "spatial_") 的干扰
    new_rna_assay <- CreateAssayObject(counts = raw_counts)
    
    # C. 将新 Assay 加入对象
    obj_210[["RNA"]] <- new_rna_assay
    
    # D. 设为默认
    DefaultAssay(obj_210) <- "RNA"
    
    # E. 删除旧 Assay 以释放内存
    obj_210[[old_name]] <- NULL
    
    cat("  [成功] 数据已迁移至 'RNA' Assay，旧 Assay 已删除。\n")
  } else {
    DefaultAssay(obj_210) <- "RNA"
  }
  
  # 3. 清除可能冲突的元数据 (Metadata)
  cat("  [清理] 删除旧的 QC 元数据...\n")
  cols_to_remove <- c("nFeature_RNA", "nCount_RNA", "percent.mt", 
    "nFeature_Spatial", "nCount_Spatial", 
    paste0("nFeature_", old_name), paste0("nCount_", old_name))
  
  meta_cols <- colnames(obj_210@meta.data)
  for (col in cols_to_remove) {
    if (col %in% meta_cols) {
      obj_210[[col]] <- NULL
    }
  }
  
  # 4. 更新回列表
  all_individual_seurat_objects[["GSE210152"]] <- obj_210
  cat("  GSE210152 修复完成。\n")
  
  # 清理临时变量 (内存优化)
  rm(obj_210, raw_counts, new_rna_assay)
  gc()
  
} else {
  cat("  列表中未找到 GSE210152，跳过修复。\n")
}

# ==============================================================================
# [补丁修复版] 修复 GSE247238 基因名 (Seurat v5 兼容)
# ==============================================================================
cat("\n>>> 开始修复 GSE247238 的基因命名问题...\n")

if ("GSE247238" %in% names(all_individual_seurat_objects)) {
  
  # 1. 提取对象
  obj_247 <- all_individual_seurat_objects[["GSE247238"]]
  
  # [关键修复] Seurat v5 必须先合并图层，才能提取完整的 counts 矩阵
  cat("  合并图层 (JoinLayers) 以提取完整矩阵...\n")
  obj_247[["RNA"]] <- JoinLayers(obj_247[["RNA"]])
  
  # 2. 找到原始特征文件
  data_dir_247 <- "/public3/DSC/single_cell/GSE247238_Carotid"
  feature_files <- list.files(data_dir_247, pattern = "_features.tsv.gz", full.names = TRUE)
  
  if (length(feature_files) > 0) {
    cat(paste("  读取特征文件用于映射:", basename(feature_files[1]), "\n"))
    gene_table <- read.table(feature_files[1], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # 假设: V1=ID, V2=Symbol
    id_to_symbol <- gene_table$V2
    names(id_to_symbol) <- gene_table$V1
    
    # 3. 准备新的基因名
    # [关键修复] 使用 LayerData 提取合并后的矩阵
    # 此时对象只有一层 counts，可以直接提取
    raw_counts <- LayerData(obj_247, assay = "RNA", layer = "counts")
    
    current_ids <- rownames(raw_counts)
    
    # 匹配
    new_gene_names <- id_to_symbol[current_ids]
    
    # 处理匹配失败的情况 (保留原ID)
    na_idx <- is.na(new_gene_names)
    new_gene_names[na_idx] <- current_ids[na_idx]
    
    # 处理重复 Symbol (添加 .1, .2)
    if (any(duplicated(new_gene_names))) {
      cat(paste("  [注意] 发现", sum(duplicated(new_gene_names)), "个重复的 Gene Symbol，已自动添加后缀。\n"))
      new_gene_names <- make.unique(new_gene_names)
    }
    
    cat(paste("  转换示例: ", current_ids[1], " -> ", new_gene_names[1], "\n"))
    
    # 4. 修改矩阵行名
    rownames(raw_counts) <- new_gene_names
    
    # 5. 重建 Seurat 对象
    # 这是最安全的做法，防止旧的 feature metadata 冲突
    new_obj_247 <- CreateSeuratObject(
      counts = raw_counts,
      project = "GSE247238",
      meta.data = obj_247@meta.data # 保留所有临床信息
    )
    
    # 恢复关键标识
    new_obj_247$Source_GSE <- "GSE247238"
    new_obj_247$AC_PA <- obj_247$AC_PA
    
    # 6. 更新回总列表
    all_individual_seurat_objects[["GSE247238"]] <- new_obj_247
    
    cat("  GSE247238 基因名已成功转换为 Gene Symbol！\n")
    
    # 清理内存
    rm(obj_247, new_obj_247, raw_counts, gene_table, id_to_symbol)
    gc()
    
  } else {
    cat("  [错误] 在目录下找不到 _features.tsv.gz 文件，无法自动转换基因名！\n")
  }
} else {
  cat("  列表中未找到 GSE247238，跳过。\n")
}


cat("  >>> 修复步骤完成，请继续运行 9.1 合并步骤。\n")
# ------------------------------------------------------------------------------
# 9.1 初步合并所有数据
# ------------------------------------------------------------------------------
if (length(all_individual_seurat_objects) < 2) {
  stop("无法合并: 列表中至少需要两个 Seurat 对象。")
}

cat("  合并所有单独的 Seurat 对象...\n")
merge_data <- merge(
  x = all_individual_seurat_objects[[1]],
  y = all_individual_seurat_objects[-1],
  add.cell.ids = names(all_individual_seurat_objects),
  project = "Athero_Integrated"
)
merge_data$Original_orig.ident <- merge_data$orig.ident

# [内存优化 1] 合并完成后，立即删除原始的列表
cat("  [Memory] Deleting 'all_individual_seurat_objects' to free RAM...\n")
rm(all_individual_seurat_objects)
gc() # 强制垃圾回收

# ------------------------------------------------------------------------------
# 9.2 QC 指标计算与初始统计
# ------------------------------------------------------------------------------
cat("\n  计算线粒体比例...\n")
merge_data[["percent.mt"]] <- PercentageFeatureSet(merge_data, pattern = "^MT-|^mt-", assay = "RNA")

# --- 统计步骤 1: QC前的细胞数 ---
cat("  Generating Pre-QC cell count statistics...\n")
meta_df <- merge_data@meta.data
if (!"GSM_ID" %in% colnames(meta_df)) { meta_df$GSM_ID <- NA }
gsm_indices <- grepl("^GSM", meta_df$Patient_ID)
meta_df$GSM_ID[gsm_indices] <- meta_df$Patient_ID[gsm_indices]

stats_before_qc <- meta_df %>%
  group_by(Source_GSE, Patient_ID, GSM_ID, Original_orig.ident) %>%
  summarise(Count_Before_QC = n(), .groups = 'drop')

# ------------------------------------------------------------------------------
# 9.3 可视化 QC 指标 (过滤前)
# ------------------------------------------------------------------------------
save_path <- "/public3/DSC/single_cell/Result/"
if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

cat("  Visualizing QC metrics...\n")
tryCatch({
  p_vln <- VlnPlot(merge_data, 
       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
       pt.size = 0.1, 
       group.by = "Source_GSE") + labs(group = "Source GSE")
  ggsave(filename = paste0(save_path, "QC_Violin_Source_GSE.pdf"), 
   plot = p_vln, width = 12, height = 6)
  
  plot1 <- FeatureScatter(merge_data, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Source_GSE")
  plot2 <- FeatureScatter(merge_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Source_GSE")
  p_scatter <- plot1 + plot2
  print(p_scatter)
  ggsave(filename = paste0(save_path, "QC_FeatureScatter_Source_GSE.pdf"), 
   plot = p_scatter, width = 12, height = 6)
  
  # 清理绘图变量
  rm(p_vln, plot1, plot2, p_scatter)
  gc(verbose = FALSE)
  
}, error = function(e) {
  cat(paste("  Warning: Visualization failed:", e$message, "\n"))
})


# ------------------------------------------------------------------------------
# 9.4 执行自定义过滤
# ------------------------------------------------------------------------------
cat("\n  Filtering cells based on CUSTOM thresholds...\n")
cells_before <- ncol(merge_data)

merge_data <- subset(
  merge_data,
  subset = 
    nFeature_RNA > 200 & 
    percent.mt < 15 & 
    (
      (Source_GSE == "GSE260657" & nFeature_RNA < 12000) | 
  (Source_GSE != "GSE260657" & nFeature_RNA < 5000)
    )
)

cells_after <- ncol(merge_data)
cat(paste("  过滤前:", cells_before, " -> 过滤后:", cells_after, "\n"))

# --- 统计步骤 2: QC后的细胞数 ---
cat("  Generating Post-QC cell count statistics...\n")
meta_df_post <- merge_data@meta.data
if (!"GSM_ID" %in% colnames(meta_df_post)) { meta_df_post$GSM_ID <- NA }
gsm_indices_post <- grepl("^GSM", meta_df_post$Patient_ID)
meta_df_post$GSM_ID[gsm_indices_post] <- meta_df_post$Patient_ID[gsm_indices_post]

stats_after_qc <- meta_df_post %>%
  group_by(Source_GSE, Patient_ID, GSM_ID, Original_orig.ident) %>%
  summarise(Count_After_QC = n(), .groups = 'drop')

cat("  Merging and saving statistics table...\n")
final_stats <- full_join(stats_before_qc, stats_after_qc, 
       by = c("Source_GSE", "Patient_ID", "GSM_ID", "Original_orig.ident"))
final_stats$Count_After_QC[is.na(final_stats$Count_After_QC)] <- 0
write.csv(final_stats, file = paste0(save_path, "Sample_Cell_Counts_QC_Summary.csv"), row.names = FALSE)

# [内存优化 2] 删除不需要的元数据统计表
rm(meta_df, meta_df_post, stats_before_qc, stats_after_qc, final_stats)
gc()

# ==============================================================================
# 步骤 1: 提取并处理参考数据集 (GSE159677)
# ==============================================================================
cat("\n>>> 步骤 1: 提取并处理参考数据集 (GSE159677)...\n")
ref_name <- "GSE159677"

# 拆分
seurat_ref <- subset(merge_data, subset = Source_GSE == ref_name)
seurat_others <- subset(merge_data, subset = Source_GSE != ref_name)

# [内存优化 3] 拆分完成后，立即删除巨大的 merge_data
cat("  [Memory] Deleting 'merge_data' object...\n")
rm(merge_data)
gc()

# 处理 Reference
cat("  GSE159677: Pre-processing & Harmony...\n")
seurat_ref[["RNA"]] <- JoinLayers(seurat_ref[["RNA"]])
DefaultAssay(seurat_ref) <- "RNA"
seurat_ref <- NormalizeData(seurat_ref, verbose = FALSE)
seurat_ref <- FindVariableFeatures(seurat_ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat_ref <- ScaleData(seurat_ref, verbose = FALSE)
seurat_ref <- RunPCA(seurat_ref, verbose = FALSE)
seurat_ref <- RunHarmony(seurat_ref, group.by.vars = "Patient_ID", dims.use = 1:30, verbose = FALSE)

# 1.3 可视化 GSE159677 的内部结构
cat("  GSE159677: 运行 UMAP (基于 Harmony) 并绘图...\n")
seurat_ref <- RunUMAP(seurat_ref, reduction = "harmony", dims = 1:20)

p_ref_1 <- DimPlot(seurat_ref, reduction = "umap", group.by = "Patient_ID") + ggtitle("GSE159677 Internal: By Patient")
p_ref_2 <- DimPlot(seurat_ref, reduction = "umap", group.by = "Sample_Type") + ggtitle("GSE159677 Internal: By Type")
p_ref_combined <- p_ref_1 + p_ref_2
print(p_ref_combined)
ggsave(filename = paste0(save_path, "GSE159677_Internal_Harmony_UMAP.pdf"), 
 plot = p_ref_combined, width = 12, height = 5)
cat("  已保存 GSE159677 内部 Harmony UMAP 图。\n")

seurat_ref$Study_Type <- "Reference"

# ==============================================================================
# 步骤 2: 预处理其他数据集 (查询集)
# ==============================================================================
cat("\n>>> 步骤 2: 预处理其他数据集...\n")
other_studies_list <- SplitObject(seurat_others, split.by = "Source_GSE")

# [内存优化 4] 拆分为列表后，删除 seurat_others 对象
cat("  [Memory] Deleting 'seurat_others'...\n")
rm(seurat_others)
gc()

processed_others_list <- list()
for (study_name in names(other_studies_list)) {
  cat(paste("  Processing:", study_name, "...\n"))
  obj <- other_studies_list[[study_name]]
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]]) # 修复 Layers 问题
  
  if (ncol(obj) < 50) next
  
  obj <- NormalizeData(obj, verbose = FALSE)
  tryCatch({
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj$Study_Type <- "Query"
    processed_others_list[[study_name]] <- obj
  }, error = function(e) { cat(paste0(" Error: ", e$message, "\n")) })
}

# [内存优化 5] 删除原始拆分列表
rm(other_studies_list)
gc()

# ==============================================================================
# 步骤 3: 准备整合列表
# ==============================================================================
cat("\n>>> 步骤 3: 构建基础整合列表...\n")
integration_list <- c(list(seurat_ref), processed_others_list)
names(integration_list)[1] <- ref_name

# [内存优化 6] 列表构建完成后，删除独立的 seurat_ref 和 processed_others_list
# 注意：此时 integration_list 持有这些对象的引用/拷贝
rm(seurat_ref, processed_others_list)
gc()

# ------------------------------------------------------------------------------
# 定义整合特征 (Features)
# ------------------------------------------------------------------------------
cat("  定义整合基因特征...\n")
# 简单使用 SelectIntegrationFeatures (内存占用较小)
features <- SelectIntegrationFeatures(object.list = integration_list, nfeatures = 3000)

# ==============================================================================
# [流程优化] 串行处理：先 RPCA，清理内存，再 CCA
# ==============================================================================

# ------------------------------------------------------------------------------
# BLOCK A: RPCA 流程
# ------------------------------------------------------------------------------
cat("\n>>> BLOCK A: 开始 RPCA 整合流程...\n")

# A.1 准备 RPCA 数据
cat("  [Memory] Creating temporary list for RPCA...\n")
integration_list_for_rpca <- list()
for (name in names(integration_list)) {
  obj <- integration_list[[name]]
  DefaultAssay(obj) <- "RNA"
  actual_feats <- intersect(features, rownames(obj))
  if(length(actual_feats) < 50) next
  
  obj <- ScaleData(obj, features = actual_feats, verbose = FALSE)
  obj <- RunPCA(obj, features = actual_feats, verbose = FALSE, npcs = 30)
  integration_list_for_rpca[[name]] <- obj
}
# 此时 integration_list_for_rpca 会占用额外内存，但我们还没做CCA，所以还行

# A.2 执行 RPCA 整合
cat("  Running FindIntegrationAnchors (RPCA)...\n")
# 自动寻找 Reference 索引
ref_idx <- which(names(integration_list_for_rpca) == ref_name)
if(length(ref_idx)==0) ref_idx <- 1

rpca_anchors <- FindIntegrationAnchors(
  object.list = integration_list_for_rpca,
  reference = ref_idx,
  anchor.features = features,
  dims = 1:30, reduction = "rpca", k.anchor = 5
)

# [内存优化 7] 找到 Anchor 后，立即删除用于准备的列表
cat("  [Memory] Deleting 'integration_list_for_rpca'...\n")
rm(integration_list_for_rpca) 
gc() 

cat("  Running IntegrateData (RPCA)...\n")
final_integrated_rpca <- IntegrateData(anchorset = rpca_anchors, dims = 1:30, new.assay.name = "integrated_rpca")

# [内存优化 8] 整合完成后，立即删除 Anchors
cat("  [Memory] Deleting 'rpca_anchors'...\n")
rm(rpca_anchors)
gc()

# A.3 RPCA 后处理 (Scale, PCA, UMAP)
cat("  Post-processing RPCA result...\n")
DefaultAssay(final_integrated_rpca) <- "integrated_rpca"
final_integrated_rpca <- ScaleData(final_integrated_rpca, verbose = FALSE)
final_integrated_rpca <- RunPCA(final_integrated_rpca, verbose = FALSE, npcs = 30)
final_integrated_rpca <- RunUMAP(final_integrated_rpca, dims = 1:30, reduction = "pca")
final_integrated_rpca <- FindNeighbors(final_integrated_rpca, dims = 1:30, reduction = "pca")
final_integrated_rpca <- FindClusters(final_integrated_rpca, resolution = 0.6)
final_integrated_rpca$Integration_Method <- "RPCA"

# A.4 立即保存 RPCA 结果 (防止后面 CCA 爆内存导致前功尽弃)
cat("  [Save] Saving RPCA object to disk...\n")
saveRDS(final_integrated_rpca, file = paste0(save_path, "RPCA_Integrated_Seurat_Object.rds"))

# ------------------------------------------------------------------------------
# BLOCK B: CCA 流程
# ------------------------------------------------------------------------------
# ==============================================================================
# RESCUE SCRIPT: 从 RPCA 存档恢复并重试 CCA
# ==============================================================================

# 1. 设置路径
save_path <- "/public3/DSC/single_cell/Result/"
rpca_file <- paste0(save_path, "RPCA_Integrated_Seurat_Object.rds")

# 2. 内存与并行设置 (关键修复)
# 调大 future 全局变量限制 (设置为 200GB 或更大，取决于你的服务器)
options(future.globals.maxSize = 200 * 1024^3) 
# 【强制串行】关闭并行计算，这是解决 "FutureInterruptError" 最稳妥的方法
plan("sequential") 

cat("\n>>> [RESCUE] 正在读取已保存的 RPCA 对象...\n")
if (!file.exists(rpca_file)) {
  stop("错误: 找不到 RPCA_Integrated_Seurat_Object.rds 文件，无法恢复。请检查路径。")
}

# 3. 读取数据

final_integrated_rpca <- readRDS(rpca_file)
cat("  成功读取 RPCA 对象。细胞数:", ncol(final_integrated_rpca), "\n")

# ==============================================================================
# 准备 CCA 数据 (从 RPCA 对象中拆分)
# ==============================================================================
cat("\n>>> [RESCUE] 正在从 RPCA 对象中提取数据用于 CCA...\n")

# 切换回 RNA Assay
DefaultAssay(final_integrated_rpca) <- "RNA"

# 确保 Layers 是合并的 (Seurat v5)
final_integrated_rpca[["RNA"]] <- JoinLayers(final_integrated_rpca[["RNA"]])

# 重新拆分为列表 (基于 Source_GSE)
# 注意：这里我们不需要重新做 Normalize/FindVariableFeatures，因为 RPCA 对象里的 RNA assay 通常已经包含这些信息
# 但为了保险起见，我们在拆分后快速检查一下
integration_list_for_cca <- SplitObject(final_integrated_rpca, split.by = "Source_GSE")

# 释放掉巨大的 RPCA 对象以节省内存给 CCA 用
rm(final_integrated_rpca)
gc()

cat("  数据已拆分为", length(integration_list_for_cca), "个部分。\n")

# 重新选择整合特征 (保证与之前一致)
cat("  重新选择整合特征...\n")
features <- SelectIntegrationFeatures(object.list = integration_list_for_cca, nfeatures = 3000)

# 标准化与缩放检查 (CCA需要数据经过ScaleData)
cat("  正在为 CCA 准备数据 (ScaleData)...\n")
for (i in names(integration_list_for_cca)) {
  # 只对整合特征进行 Scale，节省内存
  integration_list_for_cca[[i]] <- ScaleData(integration_list_for_cca[[i]], features = features, verbose = FALSE)
}

# ==============================================================================
# 执行 CCA 整合 (内存优化模式)
# ==============================================================================
cat("\n>>> [RESCUE] 开始执行 CCA 整合 (串行模式，防止中断)...\n")

ref_name <- "GSE159677" # 你的参考数据集名称
ref_idx <- which(names(integration_list_for_cca) == ref_name)
if(length(ref_idx)==0) ref_idx <- 1

# 1. FindIntegrationAnchors (最容易报错的步骤)
# 这里的 reduction="cca" 是内存杀手
tryCatch({
  cca_anchors <- FindIntegrationAnchors(
    object.list = integration_list_for_cca,
    reference = ref_idx,
    anchor.features = features,
    dims = 1:30, 
    reduction = "cca", 
    k.anchor = 5, # 如果依然爆内存，可以将此值降低为 20 或更低 (默认是5，不用动)
    verbose = TRUE
  )
  cat("  CCA Anchors 寻找成功！\n")
}, error = function(e) {
  stop("CCA FindAnchors 依然失败。原因: ", e$message)
})

# 清理列表
rm(integration_list_for_cca)
gc()

# 2. IntegrateData
cat("  Running IntegrateData (CCA)...\n")
final_integrated_cca <- IntegrateData(anchorset = cca_anchors, dims = 1:30, new.assay.name = "integrated_cca")

# 清理 Anchors
rm(cca_anchors)
gc()
final_integrated_cca <- readRDS('/public3/DSC/single_cell/Result/CCA_Integrated_Seurat_Object_1.rds')
# 3. 后处理
cat("  Post-processing CCA result...\n")
DefaultAssay(final_integrated_cca) <- "integrated_cca"
final_integrated_cca <- ScaleData(final_integrated_cca, verbose = FALSE)
final_integrated_cca <- RunPCA(final_integrated_cca, npcs = 15, verbose = FALSE)

final_integrated_cca <- RunUMAP(final_integrated_cca, dims = 1:15, reduction = "pca",
        min.dist = 0.1,
        spread = 1.5, 
        n.neighbors = 20,   # 减少邻居数，简化结构
        verbose = FALSE)

final_integrated_cca <- FindNeighbors(final_integrated_cca,
        dims = 1:15,
        reduction = "pca",
        verbose = FALSE)
final_integrated_cca <- FindClusters(
  final_integrated_cca,
  resolution = 1.0,  # 降低分辨率，数值越小聚类越少
  verbose = FALSE
)

ref <- HumanPrimaryCellAtlasData()
rds_for_SingleR <- GetAssayData(object = final_integrated_cca, layer = 'data')
clusters <- final_integrated_cca@meta.data$seurat_clusters
rds.hesc <- SingleR(test = rds_for_SingleR, ref = ref, labels = ref$label.main,
        clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
final_integrated_cca$SingleR_labels <- rds.hesc$labels[match(clusters, rownames(rds.hesc))]

p1 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "Source_GSE", raster = TRUE) + 
  ggtitle("CCA: By Study (All Sampled)")

p2 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "AC_PA", raster = TRUE) + 
  ggtitle("CCA: By Condition (All Sampled)")

p3 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = TRUE) + 
  ggtitle("CCA: By Cluster (All Sampled)")
p4 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "SingleR_labels", label = TRUE, repel = TRUE, raster = TRUE) + 
  ggtitle("SingleR: Cell Type Annotation") +
  theme(legend.position = "bottom")

p1
p2
p3
p4


# 4. 保存
cat("  [Save] Saving CCA object to disk...\n")
saveRDS(final_integrated_cca, file = paste0(save_path, "CCA_Integrated_Seurat_Object.rds"))

cat("\n=== CCA 流程修复并完成 ===\n")

# ==============================================================================
# 步骤 4: 绘图 (现在内存中只有两个最终对象，压力小很多)
# ==============================================================================
cat("\n>>> 步骤 4: 可视化...\n")
save_path <- "/public3/DSC/single_cell/Result/"
rpca_file <- paste0(save_path, "RPCA_Integrated_Seurat_Object.rds")
final_integrated_rpca <- readRDS(rpca_file)
file = paste0(save_path, "CCA_Integrated_Seurat_Object.rds")
final_integrated_cca <- readRDS(file)
# 4.1 RPCA 绘图
p_rpca_1 <- DimPlot(final_integrated_rpca, reduction = "umap", group.by = "Source_GSE", raster = TRUE) + ggtitle("RPCA: By Study")
p_rpca_2 <- DimPlot(final_integrated_rpca, reduction = "umap", group.by = "AC_PA", raster = TRUE) + ggtitle("RPCA: By Condition")
p_rpca_3 <- DimPlot(final_integrated_rpca, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = TRUE) + ggtitle("RPCA: By Cluster")
ggsave(filename = paste0(save_path, "RPCA_Final_UMAP.pdf"), plot = (p_rpca_1 + p_rpca_2) / p_rpca_3, width = 12, height = 12)

# 4.2 CCA 绘图
p_cca_1 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "Source_GSE", raster = TRUE) + ggtitle("CCA: By Study")
p_cca_2 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "AC_PA", raster = TRUE) + ggtitle("CCA: By Condition")
p_cca_3 <- DimPlot(final_integrated_cca, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = TRUE) + ggtitle("CCA: By Cluster")
ggsave(filename = paste0(save_path, "CCA_Final_UMAP.pdf"), plot = (p_cca_1 + p_cca_2) / p_cca_3, width = 12, height = 12)

# 4.3 对比绘图
cat("  Saving comparison plot...\n")
combined_plot <- (p_rpca_1 + p_cca_1)
ggsave(filename = paste0(save_path, "RPCA_vs_CCA_Study_Comparison.pdf"), plot = combined_plot, width = 16, height = 6)
rm(GSE131778_combined, GSE210152_combined, GSE247238_combined, GSE260657_combined, obj_check, seurat_obj, seurat_obj_131, seurat_obj_210, sobj); gc()
cat("\n=== 全部完成，内存已清理 ===\n")




#### MACRO_MONO ####
# ==============================================================================
# 0. 准备工作
# ==============================================================================
cat(">>> 正在加载 SingleR 参考数据集 (HumanPrimaryCellAtlas)...\n")
ref <- HumanPrimaryCellAtlasData()

# 1. 路径设置
base_dir <- "/public3/DSC/single_cell/Result/"
output_dir <- paste0(base_dir, "Lightweight_Plots/")
if(!dir.exists(output_dir)) dir.create(output_dir)

big_file_rpca <- paste0(base_dir, "RPCA_Integrated_Seurat_Object.rds")
big_file_cca  <- paste0(base_dir, "CCA_Integrated_Seurat_Object.rds") 

mini_file_rpca <- paste0(output_dir, "RPCA_Mini_umap_only.rds")
mini_file_cca  <- paste0(output_dir, "CCA_Mini_umap_only.rds")

# ==============================================================================
# 终极函数：全量保存 + 迷你抽取
# ==============================================================================
process_and_save_all <- function(input_path, mini_output_path, ref) {
  
  cat(paste0("\n>>> [读取] 正在加载大文件: ", basename(input_path), " ...\n"))
  if (!file.exists(input_path)) { stop("找不到源文件！") }
  
  # 1. 读取完整对象
  full_obj <- readRDS(input_path)
  cat(paste0("    [状态] 原始细胞数: ", ncol(full_obj), "\n"))
  
  # ==========================================================================
  # 步骤 A: 修复 "Plate" 问题 (在完整对象上操作)
  # ==========================================================================
  if (any(full_obj$orig.ident == "Plate")) {
    plate_count <- sum(full_obj$orig.ident == "Plate")
    full_obj$orig.ident[full_obj$orig.ident == "Plate"] <- "GSE260657"
    cat(paste0("    [修复] 已将 ", plate_count, " 个 'Plate' 细胞更名为 'GSE260657'。\n"))
  } else {
    cat("    [检查] 未发现 'Plate' 命名问题，跳过修复。\n")
  }
  
  # ==========================================================================
  # 步骤 B: 高效 SingleR 注释 (Cluster-level) - 避免内存爆炸
  # ==========================================================================
  #### SingleR ####
  
  rds_for_SingleR <- GetAssayData(object = full_obj, layer = 'data')
  
  #注释:基于cluster
  clusters <- full_obj@meta.data$seurat_clusters
  rds.hesc <- SingleR(test = rds_for_SingleR, ref = ref, labels = ref$label.main,
    clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  celltype = data.frame(ClusterID=rownames(rds.hesc), celltype=rds.hesc$labels, stringsAsFactors = FALSE)
  
  #整合注释结果
  full_obj@meta.data$Celltype_raw <- celltype[match(clusters, celltype$ClusterID), 'celltype']
  
  # ==========================================================================
  # 步骤 C: 保存完整大文件 (覆盖或另存)
  # ==========================================================================
  # 为了安全起见，我们保存为 "_Annotated.rds"，如果不想要新文件，可以改回 input_path 覆盖
  save_full_path <- sub(".rds$", "_Annotated.rds", input_path)
  
  cat(paste0("    [保存] 正在保存处理后的完整大文件 (含修复+注释)...\n"))
  cat(paste0("     路径: ", save_full_path, "\n"))
  
  saveRDS(full_obj, save_full_path)
  cat("    [成功] 大文件保存完毕！\n")
  
  # ==========================================================================
  # 步骤 D: 生成 Mini 对象 (用于画图)
  # ==========================================================================
  if (file.exists(mini_output_path)) {
    cat("    [跳过] Mini 文件已存在，不再重新生成。\n")
  } else {
    cat("    [抽样] 正在生成 10% Mini 对象用于绘图...\n")
    
    # 抽样
    set.seed(123)
    cells_to_keep <- sample(Cells(full_obj), size = round(0.1 * ncol(full_obj)))
    small_obj <- subset(full_obj, cells = cells_to_keep)
    
    # 提取 Metadata 和 UMAP
    meta_data <- small_obj@meta.data
    umap_coords <- small_obj@reductions$umap
    
    # 销毁 full_obj 释放内存
    rm(full_obj, small_obj)
    gc()
    
    # 重建轻量对象
    n_cells <- nrow(meta_data)
    dummy_counts <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), dims = c(3, n_cells),
      dimnames = list(c("GeneA", "GeneB", "GeneC"), rownames(meta_data))
    )
    
    current_option <- getOption("Seurat.object.assay.version")
    options(Seurat.object.assay.version = "v3") 
    
    tryCatch({
      mini_obj <- CreateSeuratObject(counts = dummy_counts, meta.data = meta_data)
      if (!is.null(umap_coords)) mini_obj@reductions$umap <- umap_coords
      
      saveRDS(mini_obj, mini_output_path)
      cat(paste0("    [成功] Mini 对象已保存: ", mini_output_path, "\n"))
      
    }, finally = {
      options(Seurat.object.assay.version = current_option)
    })
    
    rm(mini_obj, dummy_counts, meta_data, umap_coords)
    gc()
  }
}

# ==============================================================================
# 执行任务
# ==============================================================================

# 1. 处理 RPCA (大文件会被另存为 *_Annotated.rds)
process_and_save_all(big_file_rpca, mini_file_rpca, ref)

cat("\n----------------------------------------------------\n")
cat("等待内存回收...\n")
Sys.sleep(3)
cat("----------------------------------------------------\n")

# 2. 处理 CCA
process_and_save_all(big_file_cca, mini_file_cca, ref)

cat("\n=== 全部处理完成 ===\n")


mini_output_path <- "/public3/DSC/single_cell/Result/Lightweight_Plots/RPCA_Mini_umap_only.rds"
cat(">>> [绘图] 生成 RPCA 最终图表...\n")

data <- readRDS(mini_output_path)
metadata <- data@meta.data

# 定义统一的因子顺序（与 F1.3 保持绝对一致）
target_levels <- c(
  "Epithelial Cell",  
  "Endothelial Cell",
  "VSMC",
  "Macrophage",
  "Monocyte",
  "Neutrophils",
  "T Cell",
  "NK Cell",
  "B Cell",
  "cDC1",             
  "pDC"              
)

metadata <- metadata %>% 
  mutate(
    Celltype_raw = case_when(
      Celltype_raw %in% c("Smooth_muscle_cells", "Tissue_stem_cells",'Chondrocytes','MSC') ~ "VSMC",
      Celltype_raw == "B_cell" ~ "B Cell",
      Celltype_raw == "T_cells" ~ "T Cell",
      Celltype_raw == "NK_cell" ~ "NK Cell",
      # 【修改点1】：彻底统一 Endothelial 和 Epithelial 的命名格式
      Celltype_raw %in% c("Endothelial_cells", "Endothelial Cell") ~ "Endothelial Cell",
      Celltype_raw %in% c("Epithelial_cells", "Epithelial cells", "Epithelial Cell") ~ "Epithelial Cell",
      TRUE ~ as.character(Celltype_raw)
    )
  ) %>%
  mutate(
    # 【修改点2】：应用全局统一的 levels 顺序
    Celltype_raw = factor(
      Celltype_raw,
      levels = intersect(target_levels, unique(Celltype_raw))
    )
  )
data@meta.data <- metadata

# 图1: SingleR 结果 (Cluster)
p_cluster <- DimPlot(
  data, 
  group.by = "seurat_clusters", 
  label = TRUE, 
  label.size = 8,      # 恢复为您最初的 8
  reduction = 'umap',
  pt.size = 0.5        # 【统一强制设定点大小】
) + 
  theme_dr() +         # 保持原汁原味的箭头
  theme(
    panel.grid = element_blank(), 
    plot.title = element_blank()
  ) + 
  NoLegend() 

# 图2: Celltype 分布
p_celltype <- DimPlot(
  data, 
  group.by = "Celltype_raw", 
  label = TRUE, 
  label.size = 8,      # 恢复为您最初的 8
  repel = TRUE,        # 【修复】：开启防重叠，B-Cell和Epithelial会自动弹开
  reduction = 'umap',
  pt.size = 0.5        # 【统一强制设定点大小】
) +
  theme_dr() +         # 保持原汁原味的箭头
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(size = 24, family = "Arial"), 
    legend.key.height = unit(1.25, "cm"), 
    text = element_text()
  )

p1 <- p_cluster | p_celltype
print(p1)

ggplot2::ggsave(p1,filename = '/public3/DSC/single_cell/Result/figer_new/F1.1_UMAP_celltype_cluster.png',width = 18,height = 8)
ggplot2::ggsave(p1,filename = '/public3/DSC/single_cell/Result/figer_new/F1.1_UMAP_celltype_cluster.pdf',width = 20,height = 8, 
                device = cairo_pdf)


# 图2: 检查清洗后的 AC_PA
p2 <- DimPlot(data, group.by = "AC_PA", cols = c("#E41A1C", "#377EB8"), reduction='umap') + 
  theme_dr() + ggtitle("Plaque vs Control")

cat(">>> 正在深度清洗 orig.ident 并生成纯净的 GSE_ID...\n")

# 1. 提取原始 ID
clean_ids <- data$orig.ident

# 2. 手动收敛特定格式的样本到它们对应的 GSE 项目
# 把所有的 GSMxxxxxx 都归为 GSE247238
clean_ids[grepl("^GSM", clean_ids)] <- "GSE247238"

# 把 patient1AC, patient2PA 等归为 GSE159677
clean_ids[grepl("^patient", clean_ids)] <- "GSE159677"

# 3. 剥离冗长后缀，只保留纯净的 "GSE+数字" (例如 GSE213740_AscAorta_Control_PA -> GSE213740)
clean_ids <- str_extract(clean_ids, "GSE[0-9]+")

# 4. 容错兜底：如果有既不是 GSM/patient，也没有 GSE 前缀的未知数据，保留原名防 NA
clean_ids[is.na(clean_ids)] <- data$orig.ident[is.na(clean_ids)]

# 5. 将洗净的 ID 存入 Seurat 对象
data$GSE_Clean_ID <- clean_ids
data <- subset(data, subset = GSE_Clean_ID != "GSE213740")
cat(">>> 清洗完成！现在的分类包含：\n")
print(table(data$GSE_Clean_ID))

# ==============================================================================
# 重新绘制去批次效应评估图 (UMAP)
# ==============================================================================
p3_batch_clean <- DimPlot(
  data, 
  group.by = "GSE_Clean_ID", 
  shuffle = TRUE,      
  reduction = 'umap',
  pt.size = 0.5        # 【统一强制设定点大小，和 p1 绝对一致】
) + 
  theme_dr() +         # 同样使用带有箭头的坐标系
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),                             
    legend.text = element_text(size = 30, family = "Arial"),  
    legend.title = element_blank(),                           
    legend.key.height = unit(1.25, "cm"),                     
    text = element_text(family = "Arial")                     
  )

print(p3_batch_clean)

# 保存 (此时图例很长，设置 width=11 来保证左侧 UMAP 的视觉比例与 p1 中单张图一致)
ggplot2::ggsave(p3_batch_clean, 
                filename = '/public3/DSC/single_cell/Result/figer_new/F1.2_UMAP_Batch_GSE_Clean.png', 
                width = 11,   
                height = 8)

ggplot2::ggsave(p3_batch_clean, 
                filename = '/public3/DSC/single_cell/Result/figer_new/F1.2_UMAP_Batch_GSE_Clean.pdf', 
                width = 11, 
                height = 8, 
                device = cairo_pdf)
# ==============================================================================
# 在大文件 (Big File) 上运行 FindAllMarkers
# ==============================================================================
# 在 R 脚本开头添加这一行
Sys.setenv(SLURM_CONF = "/home/dsc/my_slurm_conf/slurm.conf")
cat("\n>>> [分析] 准备在 RPCA 大文件上寻找差异基因...\n")

# 1. 设置大文件路径
big_annotated_path <- "/public3/DSC/single_cell/Result/RPCA_Integrated_Seurat_Object_Annotated.rds"

if (!file.exists(big_annotated_path)) {
  stop("错误：找不到已注释的大文件: ", big_annotated_path)
}

# 2. 读取大文件
cat("    [读取] 正在加载完整 RPCA 对象 (这可能需要一点时间)...\n")
big_data <- readRDS(big_annotated_path)

# 3. [关键] 同步修改细胞类型注释
cat("    [处理] 正在应用细胞类型重命名规则...\n")

meta_big <- big_data@meta.data
meta_big <- meta_big %>% 
  mutate(
    Celltype_raw = case_when(
      Celltype_raw %in% c("Smooth_muscle_cells", "Tissue_stem_cells", 'Chondrocytes','MSC') ~ "VSMC",
      Celltype_raw == "B_cell" ~ "B Cell",
      Celltype_raw == "T_cells" ~ "T Cell",
      Celltype_raw == "NK_cell" ~ "NK Cell",
      # 【大文件同步修改】：统一命名格式
      Celltype_raw %in% c("Endothelial_cells", "Endothelial Cell") ~ "Endothelial Cell",
      Celltype_raw %in% c("Epithelial_cells", "Epithelial cells", "Epithelial Cell") ~ "Epithelial Cell",
      TRUE ~ as.character(Celltype_raw)
    )
  ) %>%
  mutate(
    # 【大文件同步修改】：应用全局统一的 levels 顺序
    Celltype_raw = factor(
      Celltype_raw, 
      levels = intersect(target_levels, unique(Celltype_raw))
    )
  )

# 将修改后的元数据写回对象
big_data@meta.data <- meta_big

# 5. 设置 Idents
Idents(big_data) <- "Celltype_raw"

# 6. 开启并行计算 (强烈建议！大文件跑 Marker 很慢)
cat("    [计算] 正在配置并行计算 (多核加速)...\n")
library(future)
plan("multisession", workers = 8) # 根据您的服务器情况调整核数，建议 4-8 核
options(future.globals.maxSize = 100 * 1024^3) # 设置最大内存为 100GB，防止报错

# 7. 运行 FindAllMarkers
cat("    [计算] 开始寻找差异基因 (FindAllMarkers)... 请耐心等待...\n")

# 确保使用 RNA assay 的 data 层
DefaultAssay(big_data) <- "RNA"
# 如果是 Seurat v5 且层是分离的，合并一下
if (packageVersion("Seurat") >= "5.0.0") {
  big_data[["RNA"]] <- JoinLayers(big_data[["RNA"]])
}

all_markers <- FindAllMarkers(
  big_data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)

# 8. 保存结果
output_csv_path <- "/public3/DSC/single_cell/Result/RPCA_All_Cell_Markers.csv"
cat(paste0("    [保存] 正在将 Marker 结果保存至: ", output_csv_path, "\n"))
write.csv(all_markers, output_csv_path, row.names = FALSE)
# 关闭并行策略，释放资源
plan("sequential")

all_markers <- read.csv("/public3/DSC/single_cell/Result/RPCA_All_Cell_Markers.csv")

# 1. 确保对象的 Idents 是你的细胞类型，并且已经设置好了正确的因子顺序
Idents(big_data) <- "Celltype_raw"

# 获取当前对象中严格的从左到右的细胞类型顺序
current_cluster_order <- levels(big_data@active.ident)

# 2. 提取每个细胞类型前5个差异基因 (关键修复点)
top20_markers <- all_markers %>%
  # 过滤掉可能存在的 NA
  filter(!is.na(cluster)) %>%
  # 【核心修改】：将表格中的 cluster 转换为 factor，并应用与 UMAP/Heatmap 完全一致的顺序
  mutate(cluster = factor(cluster, levels = current_cluster_order)) %>%
  # 先按细胞类型顺序排，再按 logFC 降序排
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  # 使用 slice_head 提取排在最前面的 5 个，确保顺序不乱
  slice_head(n = 5) %>%
  ungroup()

# 3. 提取基因列表。此时提取出来的基因顺序，完美对应了从左到右的细胞类型
heatmap_genes <- unique(top20_markers$gene)

cat(">>> 正在为热图基因补充 ScaleData...\n")
DefaultAssay(big_data) <- "RNA"

if (packageVersion("Seurat") >= "5.0.0") {
  big_data[["RNA"]] <- JoinLayers(big_data[["RNA"]]) 
}

# 补充 ScaleData
big_data <- ScaleData(big_data, features = heatmap_genes, verbose = FALSE)

p <- DoHeatmap(subset(big_data, downsample = 200),
               features = heatmap_genes,
               size = 3,
               assay = 'RNA',
               slot = "scale.data") +
  scale_fill_gradientn(colors = c("#94C4E1", "white", "red")) +
  theme(
    axis.text.y = element_text(size = 10, face = "italic", family = "Arial"),  
    legend.text = element_text(size = 12, family = "Arial"),      
    legend.title = element_text(size = 12, family = "Arial"),     
    legend.key.size = unit(0.5, "cm")    
  )

print(p)
ggplot2::ggsave(p, 
                filename = "/public3/DSC/single_cell/Result/figer_new/celltype_heatmap_10.pdf",
                width = 10, 
                height = 8,
                device = cairo_pdf)


diff_mono_macro <- FindMarkers(big_data, 
                               ident.1 = "Monocyte", 
                               ident.2 = "Macrophage",
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# 整理 Monocyte 特异基因 (Top 100)
top_mono <- diff_mono_macro %>% 
  filter(avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(100) %>%
  mutate(Type = "High_in_Monocyte")

# 整理 Macrophage 特异基因 (Top 100)
top_macro <- diff_mono_macro %>% 
  filter(avg_log2FC < 0) %>% 
  arrange(avg_log2FC) %>%  
  head(100) %>%
  mutate(Type = "High_in_Macrophage")

top_msc_from_all <- all_markers %>%
  filter(cluster == "MSC") %>%
  filter(avg_log2FC > 0) %>%     
  arrange(desc(avg_log2FC)) %>%  
  head(100)


cat(">>> [1/3] 正在提取表达量数据...\n")

features <- c(
  'GKN2', 'PGC',  # Epithelial
  'VWF', 'CDH5',       # Endothelial
  'ACTA2', 'MYH11',    # VSMC
  'SPP1', 'APOE',      # Macrophage
  'FCN1', 'VCAN',    # Monocyte
  'FCGR3B', 'S100A8',  # Neutrophils
  'CD3D', 'TRAC',      # T Cell
  'KLRD1', 'KLRB1',    # NK Cell
  'CD79A', 'MS4A1'    # B Cell
)

if (!exists("plot_data")) {
  if(exists("big_data")) {
    plot_data <- subset(big_data, downsample = 1000)
  } else {
    stop("找不到 plot_data 或 big_data，请先加载数据！")
  }
}

mat <- FetchData(plot_data, vars = c("Celltype_raw", features))

gene_data <- mat %>%
  as.data.frame() %>%
  pivot_longer(
    cols = all_of(features), 
    names_to = "Gene", 
    values_to = "Expression" 
  )
# ==============================================================================
# 2. 设置细胞顺序 (Factor Levels)
# ==============================================================================
cat(">>> [2/3] 设置细胞类型顺序...\n")

# 定义您的理想顺序
my_levels <- c(
  "Epithelial Cell", "Endothelial Cell", "VSMC", 
  "Macrophage", "Monocyte", "Neutrophils",   
  "T Cell", "NK Cell", "B Cell",        
  "cDC1", "pDC"        
)

# 取交集：只保留数据里真正有的
exist_levels <- intersect(my_levels, unique(gene_data$Celltype_raw))

# 设置因子 (注意：在 coord_flip 下，因子顺序通常要 rev 反转才能从上到下显示)
gene_data$Celltype_raw <- factor(gene_data$Celltype_raw, levels = rev(exist_levels))

# 同样，设置 Gene 的因子水平，保证画图时基因按 features 列表排序，而不是字母顺序
gene_data$Gene <- factor(gene_data$Gene, levels = features)

# ==============================================================================
# 3. 手动绘制 ggplot (Manual Plotting)
# ==============================================================================
cat(">>> [3/3] 开始绘制...\n")

# 建议：在这里定义字体大小变量，方便统一调整
base_font_size <- 14 

vlnplt <- ggplot(gene_data, aes(x = Celltype_raw, y = Expression, fill = Celltype_raw)) +
  
  # 1. 绘制小提琴
  geom_violin(scale = "width", trim = TRUE, size = 0.2, alpha = 0.8) +
  
  # 2. 翻转坐标轴
  coord_flip() +
  
  # 3. 分面
  facet_grid(~Gene, scales = "free_x") +
  
  # 4. 主题美化
  # 修改点A: 设置 base_size = 14 (默认为11)，这会整体放大所有线条和文字
  theme_classic(base_size = base_font_size) + 
  scale_fill_hue() + 
  labs(x = "", y = "Expression Level") +
  
  theme(
    legend.position = "none",
    
    # 调整分面标题 (基因名)
    strip.background = element_blank(),
    # 修改点B: 调大 size，并微调 vjust 确保文字不压线
    strip.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 0.2, vjust = 0.3), 
    
    # 调整 Y 轴 (细胞类型)
    # 修改点C: 调大 size
    axis.text.y = element_text(size = 18, color = "black", face = "italic"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # 调整 X 轴
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title = element_text(size = 25, face = "bold", vjust = 1),
    
    # 调整图之间的间距
    panel.spacing = unit(1.5, "lines"),
    
    # === 关键修改点 D: 增加页边距 ===
    # 顺序为: 上, 右, 下, 左。
    # 因为基因名在上方且向右倾斜，所以必须显著增加 Top 和 Right 的边距
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm") 
  )

# ==============================================================================
# 4. 输出与保存
# ==============================================================================
vlnplt

# 注意：如果字体变大了，可能需要适当增加 height 和 width，以防内容显得拥挤
ggplot2::ggsave(filename = "/public3/DSC/single_cell/Result/figer_new/F1.3_All_cell_markers_violin.pdf", 
    height = 10, width = 26, plot = vlnplt, device = cairo_pdf) # 高度和宽度略微增加
ggplot2::ggsave(filename = "/public3/DSC/single_cell/Result/figer_new/F1.3_All_cell_markers_violin.png", 
    height = 10, width = 26, plot = vlnplt)


#### MACRO_MONO ####
#### MACRO_MONO ####

cat(">>> [0/5] 数据准备与元数据清洗...\n")

# 1. 提取子集 (单核/巨噬/cDC1)
sub_integrated_data <- subset(x = big_data, subset = Celltype_raw %in% c("Monocyte", "Macrophage", "cDC1"))

# 2. 定义 Sample_Type (AC/PA/Control)
sub_integrated_data@meta.data <- sub_integrated_data@meta.data %>%
  mutate(Sample_Type = case_when(
    grepl("Core|Plaque", AC_PA, ignore.case = TRUE) ~ "AC",
    grepl("Adjacent", AC_PA, ignore.case = TRUE) ~ "PA",
    grepl("Control", AC_PA, ignore.case = TRUE) ~ "Control",
    TRUE ~ Sample_Type
  ))

# 检查 Sample_Type
print(table(sub_integrated_data$Sample_Type, useNA = "always"))

# 3. [关键] 合并 GSE247238 为单一批次
# 先备份旧ID
sub_integrated_data$original_GSM_ID <- sub_integrated_data$orig.ident
# 修改 orig.ident
sub_integrated_data$orig.ident <- ifelse(
  sub_integrated_data$Source_GSE == "GSE247238", 
  "GSE247238", 
  sub_integrated_data$orig.ident
)


cat(">>> 更新后的批次分布 (orig.ident):\n")
print(table(sub_integrated_data$orig.ident))

# ==============================================================================
# 步骤 1: Seurat V5 图层重置与预处理
# ==============================================================================
# 1. 提取核心数据
# 确保此时 JoinLayers 拿到最完整的矩阵
final_counts <- JoinLayers(sub_integrated_data[["RNA"]])$counts
final_metadata <- sub_integrated_data@meta.data

# 2. 彻底销毁旧对象，创建一个全新的
# 这样能保证 meta.features, scale.data, layers 全是空白且对齐的
sub_integrated_data <- CreateSeuratObject(counts = final_counts, meta.data = final_metadata)

# 3. 重新拆分 Layers (Seurat v5)
# 此时对象已经是 Assay5 格式，split 不会再报转换警告
sub_integrated_data[["RNA"]] <- split(sub_integrated_data[["RNA"]], f = sub_integrated_data$orig.ident)

# 4. 分步执行（不要用 %>%，以便定位具体哪一步报错）
sub_integrated_data <- NormalizeData(sub_integrated_data, verbose = FALSE)
sub_integrated_data <- FindVariableFeatures(sub_integrated_data, nfeatures = 3000, verbose = FALSE)

# 关键：先看 Normalize 后对象是否合法
# if(!validObject(sub_integrated_data)) stop("Object is invalid after Norm!")

sub_integrated_data <- ScaleData(sub_integrated_data, vars.to.regress = c("percent.mt", "nFeature_RNA"), verbose = FALSE)

# 5. 重新运行标准化与降维
# 此时对象是绝对干净的，不会触发 validObject 报错
sub_integrated_data <- sub_integrated_data %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt", "nFeature_RNA"), verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE)

cat(">>> PCA 运行成功！\n")

# 查看 ElbowPlot 决定维度 (可选)
# print(ElbowPlot(sub_integrated_data, ndims = 50))

# ==============================================================================
# 步骤 2: 运行三种整合方法
# ==============================================================================

# --- 方法 A: Harmony ---
cat(">>> [2/5] 正在运行 Harmony...\n")
sub_integrated_data <- RunHarmony(
  object = sub_integrated_data,
  group.by.vars = "orig.ident",
  reduction = "pca",
  reduction.save = "harmony",
  dims.use = 1:30,     # 建议与后续整合保持一致 (如30)
  theta = 2,     # 默认是2，10可能太强了，除非批次效应极难去除
  max.iter.harmony = 20,     # [修正] 参数名为 max.iter.harmony
  plot_convergence = FALSE, 
  verbose = FALSE
)
options(future.globals.maxSize = 20 * 1024^3)
# --- 方法 B: RPCA ---
cat(">>> [3/5] 正在运行 Seurat RPCA...\n")
sub_integrated_data <- IntegrateLayers(
  object = sub_integrated_data,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  k.anchor = 20,       # 针对合并后的大数据集，稍微调大 anchor 搜索
  verbose = FALSE
)

# --- 方法 C: CCA ---
cat(">>> [4/5] 正在运行 Seurat CCA...\n")
sub_integrated_data <- IntegrateLayers(
  object = sub_integrated_data,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,
  k.anchor = 20,
  verbose = FALSE
)

# ==============================================================================
# 步骤 3: 降维与 LISI 评估
# ==============================================================================
cat(">>> [5/5] 生成 UMAP 和 LISI 评分...\n")
library(lisi)
# 3.1 运行 UMAP
dims_umap <- 1:30
sub_integrated_data <- RunUMAP(sub_integrated_data, reduction = "harmony", dims = dims_umap, reduction.name = "umap.harmony", verbose = FALSE)
sub_integrated_data <- RunUMAP(sub_integrated_data, reduction = "integrated.rpca", dims = dims_umap, reduction.name = "umap.rpca", verbose = FALSE)
sub_integrated_data <- RunUMAP(sub_integrated_data, reduction = "integrated.cca", dims = dims_umap, reduction.name = "umap.cca", verbose = FALSE)

# 3.2 计算 LISI (Batch Mixing Quality)
# 定义函数
calc_lisi_mean <- function(obj, reduction_name) {
  # 提取 Embeddings
  emb <- Embeddings(obj, reduction_name)
  # 确保维度不过大，取前30维
  if(ncol(emb) > 30) emb <- emb[, 1:30]
  
  # 计算 LISI
  lisi_res <- compute_lisi(X = emb, meta_data = obj@meta.data, label_colnames = "orig.ident")
  return(lisi_res$orig.ident)
}

cat("  正在计算 LISI Score (可能需要几分钟)...\n")
sub_integrated_data$LISI_RawPCA  <- calc_lisi_mean(sub_integrated_data, "pca")
sub_integrated_data$LISI_Harmony <- calc_lisi_mean(sub_integrated_data, "harmony")
sub_integrated_data$LISI_RPCA    <- calc_lisi_mean(sub_integrated_data, "integrated.rpca")
sub_integrated_data$LISI_CCA     <- calc_lisi_mean(sub_integrated_data, "integrated.cca")

# 打印平均分 (越高越好，满分是批次总数)
cat("\n=== 平均 LISI 分数 (Batch Mixing, Higher is Better) ===\n")
cat("Raw PCA (Baseline):", round(mean(sub_integrated_data$LISI_RawPCA), 2), "\n")
cat("Harmony:     ", round(mean(sub_integrated_data$LISI_Harmony), 2), "\n")
cat("RPCA:  ", round(mean(sub_integrated_data$LISI_RPCA), 2), "\n")
cat("CCA:   ", round(mean(sub_integrated_data$LISI_CCA), 2), "\n")

# ==============================================================================
# 步骤 4: 可视化输出
# ==============================================================================
library(patchwork)
# 4.1 LISI 小提琴图
lisi_data <- sub_integrated_data@meta.data %>%
  select(starts_with("LISI_")) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Score") %>%
  mutate(Method = factor(Method, levels = c("LISI_RawPCA", "LISI_Harmony", "LISI_RPCA", "LISI_CCA")))

p_lisi <- ggplot(lisi_data, aes(x = Method, y = Score, fill = Method)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = c("gray", "#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(title = "Integration Quality (LISI)", y = "LISI Score", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# 4.2 UMAP 对比图
# 使用 raster=TRUE 加快渲染速度，shuffle=TRUE 避免特定批次遮挡
p1 <- DimPlot(sub_integrated_data, reduction = "umap.harmony", group.by = "orig.ident", pt.size=0.1, raster=TRUE, shuffle=TRUE) + ggtitle("Harmony") + NoLegend() + theme(axis.title = element_blank())
p2 <- DimPlot(sub_integrated_data, reduction = "umap.rpca", group.by = "orig.ident", pt.size=0.1, raster=TRUE, shuffle=TRUE) + ggtitle("RPCA") + NoLegend() + theme(axis.title = element_blank())
p3 <- DimPlot(sub_integrated_data, reduction = "umap.cca", group.by = "orig.ident", pt.size=0.1, raster=TRUE, shuffle=TRUE) + ggtitle("CCA") + NoLegend() + theme(axis.title = element_blank())

# 拼图
final_plot <- (p_lisi | (p1 / p2 / p3)) + plot_layout(widths = c(1, 1))

# 输出
print(final_plot)

# 注意 reduction 名字要对应
sub_integrated_data <- FindNeighbors(
  sub_integrated_data, 
  reduction = "integrated.rpca", 
  dims = 1:20, 
  k.param = 40,
  graph.name = c("integrated_rpca_nn", "integrated_rpca_snn"),
  force.recalc = TRUE # 强制重新计算
)
sub_integrated_data <- RunUMAP(
  sub_integrated_data, 
  reduction = "integrated.rpca", 
  reduction.name = "umap.rpca",
  dims = 1:20,   # <--- 减少维度，排除噪音干扰
  n.neighbors = 30
  
)
sub_integrated_data <- FindClusters(sub_integrated_data, graph.name = "integrated_rpca_snn", resolution = 0.5, cluster.name = "rpca_clusters")

DimPlot(sub_integrated_data, reduction = "umap.rpca", group.by = "rpca_clusters", label = TRUE)

# 1. 建图：必须指定 reduction = "harmony"
sub_integrated_data <- FindNeighbors(sub_integrated_data, reduction = "harmony", dims = 1:20, graph.name = "harmony_snn")

# 2. 聚类：基于上一步建好的图 (默认图名称通常是 RNA_snn 或 harmony_snn，Seurat会自动识别)
# 为了防止覆盖，建议指定 graph.name
sub_integrated_data <- FindClusters(sub_integrated_data, graph.name = "harmony_snn", resolution = 0.5, cluster.name = "harmony_clusters")

# 3. 绘图
#DimPlot(sub_integrated_data, reduction = "umap.harmony", group.by = "harmony_clusters", label = TRUE)

sub_integrated_data <- FindNeighbors(sub_integrated_data, reduction = "integrated.cca", dims = 1:20, graph.name = "cca_snn")
sub_integrated_data <- FindClusters(sub_integrated_data, graph.name = "cca_snn", resolution = 0.5, cluster.name = "cca_clusters")

#DimPlot(sub_integrated_data, reduction = "umap.cca", group.by = "cca_clusters", label = TRUE)
Layers(sub_integrated_data)

# 执行合并
sub_integrated_data <- JoinLayers(sub_integrated_data)

# 再次检查，应该只剩下 counts 和 data
Layers(sub_integrated_data)
StackedVlnPlot <- function(obj, 
         features, 
         group.by = "Celltype_raw", 
         my_levels = NULL, 
         title_size = 14,      # 字体大小直接在这里调！
         spacing = 2) {  # 图片间隙
  
  cat(">>> [1/3] 提取数据...\n")
  mat <- FetchData(obj, vars = c(group.by, features), layer = "data")
  
  if (is.null(my_levels)) my_levels <- unique(mat[[group.by]])
  mat[[group.by]] <- factor(mat[[group.by]], levels = rev(intersect(my_levels, unique(mat[[group.by]]))))
  
  cat(">>> [2/3] 绘图 (Clip Off 模式)...\n")
  
  plot_list <- list()
  
  for (i in seq_along(features)) {
    gene <- features[i]
    
    # 基础绘图
    p <- ggplot(mat, aes(x = .data[[group.by]], y = .data[[gene]], fill = .data[[group.by]])) +
      geom_violin(scale = "width", trim = TRUE, size = 0.3, alpha = 0.8) +
      
      # === 关键点 1: 关闭裁切 ===
      # clip = "off" 允许文字延伸到绘图框之外而不消失
      coord_flip(clip = "off") + 
      
      scale_fill_hue() +
      
      # === 关键点 2: 直接使用原生标题 ===
      ggtitle(gene) +
      
      theme_classic(base_size = 14) + # 基础线条粗细
      theme(
  legend.position = "none",
  
  # --- 标题 (基因名) 设置 ---
  # 这种方式您可以随心所欲控制字体大小、颜色、粗细
  plot.title = element_text(
    size = title_size,   # <--- 您的字体大小参数
    face = "bold.italic", 
    angle = 45,    # 旋转角度
    hjust = 0,     # 左对齐 (根据旋转角度微调)
    vjust = 0,     # 底部对齐
    color = "black"
  ),
  
  # --- 坐标轴清理 ---
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.x = element_blank(),
  
  # === 关键点 3: 边距控制 ===
  # t = 10: 顶部留出 10pt 的巨大空间给标题，防止文字“撞”到上边界
  # r = spacing: 右侧留出缝隙
  plot.margin = margin(t = 10, r = spacing, b = 10, l = 0) 
      )
    
    # --- 第一张图保留 Y 轴，其他去掉 ---
    if (i == 1) {
      p <- p + theme(
  axis.text.y = element_text(size = 12, color = "black"),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank()
      )
    } else {
      p <- p + theme(
  axis.text.y = element_blank(),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  # 只有右边距，左边距为0
  plot.margin = margin(t = 10, r = spacing, b = 10, l = 0)
      )
    }
    
    plot_list[[i]] <- p
  }
  
  cat(">>> [3/3] 拼接...\n")
  # 直接拼成一行，不需要再搞什么高度比例了
  final_plot <- wrap_plots(plot_list, nrow = 1)
  
  return(final_plot)
}


features_combined <- c(
  # --- [原有] 单核细胞 (Monocytes) ---
  'CD14', 'FCGR3A', 'CCR2', 'CX3CR1', 'LYZ', 'S100A8', 'S100A9',
  
  # --- [原有] 巨噬细胞 (Macrophages) ---
  'CD68', 'CD163', 'CD80', 'CD86', 'F13A1', 'APOE', 'SPP1',
  'NOS2', 'IL1B', 'ARG1', 'MRC1', 'TIMD4', 'VSIG4',
  
  # --- [新增] T 细胞 (T cells) ---
  'CD3D', 'CD3E', 'CD3G', 'TRAC', 
  
  # --- [新增] 浆细胞 (Plasma cells) ---
  'MZB1', 'IGKC', 'JCHAIN', 'SDC1', # SDC1即CD138
  
  # --- [新增] 成纤维细胞 (Fibroblasts) ---
  'COL1A1', 'DCN', 'LUM', 'FAP', 'THY1', 
  
  # --- [新增] 间充质基质细胞 (MSC) ---
  'NT5E', 'ENG', 'PDGFRB', 'MCAM' # NT5E=CD73, ENG=CD105
)

# 去除重复基因并检查对象中是否存在这些基因
features_combined <- unique(features_combined)
features_combined <- features_combined[features_combined %in% rownames(sub_integrated_data)]

vlnplt_monocytes_macrophages <- StackedVlnPlot(obj = sub_integrated_data,
           features = features_combined,
           group.by = "rpca_clusters"
)
vlnplt_monocytes_macrophages <- StackedVlnPlot(obj = sub_integrated_data,
           features = features_combined,
           group.by = "Celltype_raw"
)

vlnplt_monocytes_macrophages
#rpca_clusters == 14  间充质基质细胞
#rpca_clusters == 13  #浆细胞
#rpca_clusters == 12  #不知道是啥
#rpca_clusters == 10  成纤维细胞
#rpca_clusters == 11  CD8+ T细胞
saveRDS(sub_integrated_data, "/public3/DSC/single_cell/Result/sub_integrated_data")
#sub_integrated_data<-readRDS('/public3/DSC/single_cell/Result/sub_integrated_data')
rm(big_data)

source("/public3/DSC/single_cell/Result/singleR.R")
res <- run_cell_annotation(
  seurat_obj = sub_integrated_data, 
  return_plots = TRUE
)
sub_integrated_data <- res$obj
View(res$summary)
final_p <-res$plots$comparison
print(final_p)
ggplot2::ggsave(final_p, 
    filename = "/public3/DSC/single_cell/Result/figer_new/Mono_celltype_1.pdf",
    width = 15, 
    height = 10,
    device = cairo_pdf)



#rpca_clusters == 14  内皮细胞
#rpca_clusters == 13  #浆细胞
#rpca_clusters == 12  #不知道是啥
#rpca_clusters == 10  成纤维细胞，血小板？
#rpca_clusters == 11  CD8+ T细胞

saveRDS(sub_integrated_data, "/public3/DSC/single_cell/Result/sub_integrated_data")

#第一次清洗
sub_integrated_data <- subset(
  sub_integrated_data,
  subset = !rpca_clusters %in% c( 10, 11, 12, 13, 14)
)
source("/public3/DSC/single_cell/Result/singleR.R")
res <- run_cell_annotation(
  seurat_obj = sub_integrated_data, 
  return_plots = TRUE
)

sub_integrated_data <- res$obj
sub_integrated_data <- FindNeighbors(sub_integrated_data, reduction = "integrated.rpca",dims = 1:20)
sub_integrated_data <- FindClusters(sub_integrated_data, graph.name = "integrated_rpca_snn", resolution = 0.5, cluster.name = "rpca_clusters")
sub_integrated_data <- RunUMAP(sub_integrated_data, reduction = "integrated.rpca",n.neighbors = 30, dims = 1:15)
DimPlot(sub_integrated_data, reduction = "umap", group.by = "rpca_clusters", label = TRUE)

####去除cDC细胞####
#第二次清洗
sub_integrated_data <- subset(
  sub_integrated_data,
  subset = !rpca_clusters %in% c( 9, 7, 12, 14)
)
res <- run_cell_annotation(
  seurat_obj = sub_integrated_data, 
  return_plots = TRUE
)
sub_integrated_data <- res$obj
View(res$summary)
final_p <-res$plots$comparison
print(final_p)
sub_integrated_data <- FindNeighbors(sub_integrated_data, reduction = "integrated.rpca",dims = 1:30)
sub_integrated_data <- FindClusters(sub_integrated_data, graph.name = "integrated_rpca_snn", resolution = 0.8, cluster.name = "rpca_clusters")
sub_integrated_data <- RunUMAP(sub_integrated_data, reduction = "integrated.rpca",n.neighbors = 20, min.dist = 0.5,dims = 1:30)
DimPlot(sub_integrated_data, reduction = "umap", group.by = "rpca_clusters", label = TRUE)
#saveRDS(sub_integrated_data, "/public3/DSC/single_cell/Result/sub_integrated_data1")
#第三次清洗
res <- run_cell_annotation(
  seurat_obj = sub_integrated_data, 
  return_plots = TRUE
)
sub_integrated_data <- res$obj
final_p <-res$plots$comparison
print(final_p)
rm(res)
gc()
table(sub_integrated_data@meta.data$rpca_clusters,sub_integrated_data@meta.data$predicted.celltype.l2)
# 1. 自动识别现有的注释列
anno_cols <- c("SingleR_Monaco", "predicted.celltype.l2", "ScType_Label", "CellTypist_Label")
existing_cols <- intersect(anno_cols, colnames(sub_integrated_data@meta.data))

# 2. 提取数据并计算每个 Cluster 投票最高的结果
cluster_comparison <- sub_integrated_data@meta.data %>%
  group_by(rpca_clusters) %>%
  summarise(across(all_of(existing_cols), 
       ~ names(which.max(table(.x))), 
       .names = "{.col}"),
      Cell_Count = n()) %>%
  arrange(as.numeric(as.character(rpca_clusters)))

# 3. 打印精简表格
cat("\n>>> 各聚类在不同算法下的注释一致性对比：\n")
print(as.data.frame(cluster_comparison))

contamination_labels <- c( "CD4 TCM", "CD4 TEM", 
         "CD8 TCM", "CD8 TEM", 
)

# 2. 识别并标记需要剔除的细胞
# 逻辑：如果细胞属于 Cluster 1 或 14，且其标签在污染列表里，则标记为剔除
all_cells <- Cells(sub_integrated_data)
metadata <- sub_integrated_data@meta.data

cells_to_remove <- rownames(metadata[
  (metadata$rpca_clusters %in% c("1", "14")) & 
    (metadata$predicted.celltype.l2 %in% contamination_labels), 
])

# 3. 执行过滤
sub_integrated_data <- subset(sub_integrated_data, cells = setdiff(all_cells, cells_to_remove))

forbidden_labels <- c("Eryth", "Platelet", "B intermediate", "B memory", 
    "Plasmablast", "NK", "ILC", "HSPC")

# 2. 执行精准过滤
# 逻辑：保留那些【不在 forbidden_labels 中】且【不是特定冲突身份】的细胞
sub_integrated_data <- subset(sub_integrated_data, 
      subset = predicted.celltype.l2 %in% forbidden_labels, 
      invert = TRUE)

sub_integrated_data <- FindNeighbors(sub_integrated_data, reduction = "integrated.rpca",dims = 1:30)
sub_integrated_data <- FindClusters(sub_integrated_data, graph.name = "integrated_rpca_snn", resolution = 0.8, cluster.name = "rpca_clusters")
sub_integrated_data <- RunUMAP(sub_integrated_data, reduction = "integrated.rpca",n.neighbors = 20, min.dist = 0.5,dims = 1:30)
DimPlot(sub_integrated_data, reduction = "umap", group.by = "rpca_clusters", label = TRUE)

res <- run_cell_annotation(
  seurat_obj = sub_integrated_data, 
  return_plots = TRUE
)
sub_integrated_data <- res$obj
final_p <-res$plots$comparison
print(final_p)
rm(res)
gc()
table(sub_integrated_data@meta.data$rpca_clusters,sub_integrated_data@meta.data$predicted.celltype.l2)
# 1. 自动识别现有的注释列
anno_cols <- c("SingleR_Monaco", "predicted.celltype.l2", "ScType_Label", "CellTypist_Label")
existing_cols <- intersect(anno_cols, colnames(sub_integrated_data@meta.data))

# 2. 提取数据并计算每个 Cluster 投票最高的结果
cluster_comparison <- sub_integrated_data@meta.data %>%
  group_by(rpca_clusters) %>%
  summarise(across(all_of(existing_cols), 
       ~ names(which.max(table(.x))), 
       .names = "{.col}"),
      Cell_Count = n()) %>%
  arrange(as.numeric(as.character(rpca_clusters)))

# 3. 打印精简表格
cat("\n>>> 各聚类在不同算法下的注释一致性对比：\n")
print(as.data.frame(cluster_comparison))
saveRDS(sub_integrated_data, "/public3/DSC/single_cell/Result/sub_integrated_data_Final.rds")


#####读取CCA.R处理好的sub_integrated_data_Final 文件####
sub_integrated_data <- readRDS( "/public3/DSC/single_cell/Result/sub_integrated_data_Final.rds")
##GSE213740的三个样品就是GSE216860六个样品中的三个，是重复的
sub_integrated_data <- subset(
  x = sub_integrated_data, 
  subset = orig.ident != "GSE213740_AscAorta_Control_PA"
)
# 确保当前的 Identity 是 rpca_clusters
Idents(sub_integrated_data) <- "rpca_clusters"
sub_integrated_data <- subset(sub_integrated_data, subset = rpca_clusters != "14")
all_markers <- FindAllMarkers(sub_integrated_data, 
      only.pos = TRUE, 
      min.pct = 0.25, 
      logfc.threshold = 0.25)

# 提取每个 Cluster 前 10 个显著的 Marker 基因进行查看
top20_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

# 2. 完整打印结果 (确保能看到所有 Cluster)
# 使用 as.data.frame 可以避免 tibble 的缩略显示问题
print(as.data.frame(top20_markers))

# 保存完整的 Marker 列表到 Excel，方便后续筛选
write.csv(all_markers, "/public3/DSC/single_cell/Result/rpca_clusters_all_markers.csv", row.names = FALSE)
all_markers <- read.csv("/public3/DSC/single_cell/Result/rpca_clusters_all_markers.csv")
p_umap <- DimPlot(sub_integrated_data, reduction = "umap", label = TRUE, pt.size = 0.5)
print(p_umap)



# 将元数据转换为数据框方便操作
metadata <- sub_integrated_data@meta.data

# 创建新列并修改
# 第一步：创建Celltype_raw列（基于原始rpca_clusters的分类）
metadata <- metadata %>% 
  mutate(
    Celltype_raw = case_when(
      rpca_clusters %in% c(2, 5, 10, 11, 12, 13) ~ "Monocyte",
      rpca_clusters %in% c(1, 3, 4, 7) ~ "LAM",
      rpca_clusters %in% c(0, 6, 8, 9) ~ "Macrophage", 
      TRUE ~ as.character(Celltype_raw) 
    )
  )

# 第二步：基于原始rpca_clusters创建Celltype_raw1的细分类别
metadata <- metadata %>% 
  mutate(
    Celltype_raw1 = case_when(
      rpca_clusters == 10 ~ "Classical Mono",
      rpca_clusters == 2  ~ "Non-classical Mono",
      rpca_clusters %in% c(5,12)  ~ "Inflammatory Mono",
      rpca_clusters == 11 ~ "ISG+ Mono",
      rpca_clusters == 13 ~ "ISG+ Mono",
      rpca_clusters == 7  ~ "Foam cells1",
      rpca_clusters == 1  ~ "Foam cells2",
      rpca_clusters == 6  ~ "LYVE1+ TRM",
      rpca_clusters == 0  ~ "CX3CR1+ TRM",
      rpca_clusters %in% c(8,9) ~ "Transitional Mac",
      rpca_clusters == 4  ~ "LAM",
      rpca_clusters == 3  ~ "LAM",
      TRUE ~ Celltype_raw  # 其他情况使用Celltype_raw的值
    )
  )

# 将修改后的metadata写回Seurat对象
sub_integrated_data@meta.data <- metadata

# 提取rpca_clusters 5和12的细胞
sub_cluster5_12 <- subset(sub_integrated_data, subset = rpca_clusters %in% c(5, 12))

# 可以选择重新进行UMAP（基于原有的integrated.rpca降维）
# 但需要重新找邻居和聚类，因为原来的邻居是基于所有细胞的
sub_cluster5_12 <- FindNeighbors(sub_cluster5_12, reduction = "integrated.rpca", dims = 1:15)
sub_cluster5_12 <- FindClusters(sub_cluster5_12, graph.name = "integrated_rpca_snn", resolution = 0.8, cluster.name = "sub_rpca_clusters")
sub_cluster5_12 <- RunUMAP(sub_cluster5_12, reduction = "integrated.rpca", n.neighbors = 20, min.dist = 0.5, dims = 1:15)

# 绘制UMAP
DimPlot(sub_cluster5_12, reduction = "umap", group.by = "sub_rpca_clusters", label = TRUE)
# 确保 sub_cluster5_12 存在且包含新的聚类列
if (!exists("sub_cluster5_12")) {
  stop("请先运行子集提取和重新聚类的代码。")
}

# 设置身份为 sub_rpca_clusters
Idents(sub_cluster5_12) <- "sub_rpca_clusters"

# 运行 FindAllMarkers（只在这些新聚类之间比较）
cat("\n正在计算 sub_cluster5_12 中新聚类的标记基因...\n")
all_markers_sub <- FindAllMarkers(
  sub_cluster5_12,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 筛选显著基因（p_val_adj < 0.05）
sig_markers_sub <- all_markers_sub %>%
  filter(p_val_adj < 0.05)

# 提取每个聚类的前20个标记基因（按 avg_log2FC 降序）
top20_sub <- sig_markers_sub %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup() %>%
  arrange(cluster, dplyr::desc(avg_log2FC))

# 打印结果
cat("\n========== 每个 sub_rpca_clusters 的前50个标记基因（含 fold change）==========\n")
if (nrow(top20_sub) > 0) {
  clusters <- unique(top20_sub$cluster)
  for (cl in clusters) {
    cat("\n--- 聚类", cl, "---\n")
    # 表格形式打印基因名和 avg_log2FC
    top20_sub %>%
      filter(cluster == cl) %>%
      dplyr::select(gene, avg_log2FC) %>%
      print(row.names = FALSE)
  }
} else {
  cat("\n没有找到显著的标记基因。\n")
}


# 提取子集细胞的 barcode 和新聚类标签
sub_cells <- Cells(sub_cluster5_12)
sub_clusters <- sub_cluster5_12$sub_rpca_clusters
names(sub_clusters) <- sub_cells

# 在大对象中创建新列，默认值为 "Other"
sub_integrated_data$sub_highlight <- "Other"

# 将子集细胞的标签填入（转换为字符，便于着色）
sub_integrated_data$sub_highlight[sub_cells] <- as.character(sub_clusters[sub_cells])


library(ggplot2)
library(Seurat)
library(tidydr)  # 如果您使用 theme_dr()

# 获取所有出现的子集聚类（可能有 0-8）
highlight_clusters <- sort(unique(na.omit(as.character(sub_clusters))))
n_clusters <- length(highlight_clusters)

# 为每个聚类生成颜色（例如使用 RColorBrewer 或自定义）
# 这里用彩虹色，可根据需要调整
cluster_colors <- setNames(
  c(RColorBrewer::brewer.pal(n_clusters, "Set1"), "grey70"),
  c(highlight_clusters, "Other")
)

# 绘制 UMAP
p_highlight <- DimPlot(
  sub_integrated_data,
  group.by = "sub_highlight",
  reduction = "umap",
  label = TRUE,
  repel = TRUE,     # 避免标签重叠
  label.size = 6,
  pt.size = 1,
  cols = cluster_colors   # 手动指定颜色
) +
  theme_dr() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    plot.title = element_blank()
  ) +
  labs(color = "Subcluster")   # 图例标题

# 可选：如果您需要翻转坐标（通常不建议）
p_highlight <- p_highlight + coord_flip()

print(p_highlight)

cluster0_markers <- all_markers_sub %>%
  filter(cluster == 0, p_val_adj < 0.05) %>%
  arrange(dplyr::desc(avg_log2FC))

# 查看前50个基因（含 avg_log2FC）
cat("\n========== cluster0 的前50个标记基因（按 avg_log2FC 降序）==========\n")
print(cluster0_markers %>% 
  dplyr::select(gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
  head(50), 
      row.names = FALSE)



# 从子集对象中提取细胞条形码和新聚类标签
sub_cells <- Cells(sub_cluster5_12)
sub_clusters <- sub_cluster5_12$sub_rpca_clusters
names(sub_clusters) <- sub_cells

# 在大对象中创建新列，初始化为 NA
sub_integrated_data$sub_rpca_clusters <- NA

# 将子集细胞的标签填入（转换为字符以保持一致性）
sub_integrated_data$sub_rpca_clusters[sub_cells] <- as.character(sub_clusters[sub_cells])

# 检查是否成功
table(sub_integrated_data$sub_rpca_clusters, useNA = "ifany")
# 确保 sub_rpca_clusters 是字符型（方便映射）
sub_integrated_data$sub_rpca_clusters <- as.character(sub_integrated_data$sub_rpca_clusters)

# 定义映射：sub_rpca_clusters -> 新的 rpca_clusters 标签
cluster_to_celltype_raw1 <- c(
  "0" = "Non-classical Mono",
  "1" = "Inflammatory Mono",
  "2" = "ISG+ Mono",
  "3" = "Transitional Mac",
  "4" = "LAM",
  "5" = "Non-classical Mono",
  "6" = "Inflammatory Mono",
  "7" = "Classical Mono",
  "8" = "Classical Mono"
)

# 步骤3：获取子集细胞索引（sub_rpca_clusters 非NA）
subset_idx <- !is.na(sub_integrated_data$sub_rpca_clusters)

# 步骤4：更新 Celltype_raw1（先转换为字符，避免因子问题）
sub_integrated_data$Celltype_raw1 <- as.character(sub_integrated_data$Celltype_raw1)
sub_integrated_data$Celltype_raw1[subset_idx] <- 
  cluster_to_celltype_raw1[sub_integrated_data$sub_rpca_clusters[subset_idx]]

# 步骤5：验证更新结果
cat("更新后，子集细胞的 Celltype_raw1 分布：\n")
print(table(sub_integrated_data$Celltype_raw1[subset_idx]))

cat("\n更新后，整体 rpca_clusters vs Celltype_raw1 交叉表：\n")
print(table(sub_integrated_data$rpca_clusters, sub_integrated_data$Celltype_raw1))





all_markers <- FindAllMarkers(
  object = SetIdent(sub_integrated_data, value = "Celltype_raw1"), # 关键在这里
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 选择排序列（根据您的输出是 avg_log2FC）
sort_col <- "avg_log2FC"

# 提取每个细胞类型的前100个marker基因
top100_markers <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = !!rlang::sym(sort_col), n = 100) %>%
  dplyr::ungroup()

# 计算每个基因在多少个细胞类型的top100中出现
gene_cluster_count <- top100_markers %>%
  dplyr::select(gene, cluster) %>%
  dplyr::distinct() %>%    # 每个基因在每个cluster中只保留一条
  dplyr::count(gene, name = "n_clusters")   # 统计每个基因出现的cluster数

# 合并统计信息
top100_markers <- top100_markers %>%
  dplyr::left_join(gene_cluster_count, by = "gene")

# 提取每个细胞类型的独特 marker 基因（n_clusters == 1）
unique_markers_list <- top100_markers %>%
  dplyr::filter(n_clusters == 1) %>%
  dplyr::select(cluster, gene) %>%
  dplyr::arrange(cluster)

# 打印每个细胞类型的独特基因
cat("\n========== 每个细胞类型的独特 marker 基因 ==========\n")
if (nrow(unique_markers_list) > 0) {
  clusters_with_unique <- unique(unique_markers_list$cluster)
  for (cl in clusters_with_unique) {
    genes <- unique_markers_list %>% dplyr::filter(cluster == cl) %>% dplyr::pull(gene)
    cat("\n---", cl, "---\n")
    cat(paste(genes, collapse = ", "), "\n")
  }
  
  all_clusters <- unique(top100_markers$cluster)
  clusters_without_unique <- setdiff(all_clusters, clusters_with_unique)
  if (length(clusters_without_unique) > 0) {
    cat("\n以下细胞类型没有独特 marker 基因:\n")
    cat(paste(clusters_without_unique, collapse = ", "), "\n")
  }
} else {
  cat("\n没有任何细胞类型拥有独特 marker 基因。\n")
}


# 保存完整的 Marker 列表到 Excel，方便后续筛选
write.csv(all_markers, "/public3/DSC/single_cell/Result/Celltype_raw1_all_markers.csv", row.names = FALSE)


# 1. 确保 Idents 设置为 Celltype_raw1
Idents(sub_integrated_data) <- "Celltype_raw1"

# 2. 定义要比较的三个群
groups <- c("Foam cells1", "Foam cells2", "LAM")

# 3. 提取这三个群构建专属子集
# 这一步极其重要，只有 subset 之后，后续的比较才是单纯在这三者之间进行的
lipid_macs_subset <- subset(sub_integrated_data, idents = groups)

# 4. 寻找差异特征基因 (FindAllMarkers 默认会遍历 Idents，将每一个与其剩余的所有进行比较)
# 设置 only.pos = TRUE 因为我们通常只对能作为阳性标记（上调表达）的基因感兴趣
lipid_markers <- FindAllMarkers(lipid_macs_subset, 
        only.pos = TRUE,   
        min.pct = 0.25,    
        logfc.threshold = 0.25)  

# 5. 提取每个细胞群排名前 50 的特征基因（按 avg_log2FC 表达差异倍数降序排列）
# 推荐使用 slice_max，它是 dplyr 较新版本中替代 top_n 的更稳定函数
top50_lipid_markers <- lipid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

# 6. 打印结果 (n = 150 确保三组的 50 个基因都能完整显示在控制台)
print(top50_lipid_markers, n = 150)



# 1. 定义要画的基因列表
features_to_plot <- c(
  # --- Monocytes ---
  "S100A12", "SELL", # Classical Mono (S100A8 比 A12 更干净，SELL 是金标准)
  "MTRNR2L8",  "CCL5", # Inflammatory Mono (PTGS2 是炎症风暴核心)
  "IFIT1", "ISG15",      # ISG+ Mono (RSAD2 背景极低)
  "CLEC10A",   "THBS1",      # Repair Mono (生长因子)
  # --- Metabolic / Foam ---
  "PHLDA1",'CCL2',#"CCL7", # Foam1 (趋化/早期)
  "CD36", "MGLL", # Foam2 (成熟/脂滴)
  "TREM2",   "PLAU", # LAM (溶酶体活跃/分解脂质)
  # --- Resident ---
  'F13A1','SELENOP',     
  "CX3CR1", "C3",   # CX3CR1+ TRM
  "LYVE1", "PDGFC"      # LYVE1+ TRM (IGF1 特异性很高)
  
  
)



my_levels <- c(
  # --- 1. 单核细胞起点 & 过渡 ---
  "Classical Mono",     # Cluster 10: 始祖
  
  # --- 2. 功能性单核亚群 (炎症/抗病毒/修复) ---
  "Inflammatory Mono",  # Cluster 12: 炎症风暴
  "ISG+ Mono",   # Cluster 11/13: 干扰素反应
  "Non-classical Mono",   # Cluster 2:  修复/M2样前体
  
  
  # --- 4. 代谢/病理巨噬细胞 (泡沫化路线) ---
  "Foam cells1",  # Cluster 7:  早期/单核样泡沫
  "Foam cells2",  # Cluster 1:  成熟泡沫
  "LAM",   # Cluster 4:  脂质相关/TREM2+
  
  # --- 5. 组织驻留巨噬细胞 (维稳卫士) ---
  "Transitional Mac" ,  # Cluster 9:  最典型的M2驻留
  "CX3CR1+ TRM",       # Cluster 0:  巡逻/抗原呈递
  "LYVE1+ TRM"       # Cluster 6:  血管旁驻留

)

#泡沫marker1:LIPA+ APOC1+ CD36 APOC1 CD9 TREM2
#泡沫marker2:OLR1 PLIN2+ MARCO IL1RN CCL2
#泡沫marker3:FABP5 CTSB SPP1 CD36

#'OTOA','CLDN1','FCGBP','TREM2','MMP14'
#'TREM2','CD28','PLAU','ITGB5','KCNMA','SYNE1', "CCND1", "BEST1"
#IL2R1, PDGFC EGFL7 CR1 "SGMS1" "RGL1" CPM
#单核细胞marker "FCN1", "S100A8"
#FOLR2: F13A1,SELENOP，SLC40A1，GPR34
#APOE 巨噬细胞marker
#泡沫细胞marker SPP1,APOC1,LPL,FABP5,FABP4,CYP27A1,MMP9,SLC2A1

# 2. 将 Celltype_raw1 转换为指定顺序的因子 (Factor)
vlnplt_monocytes_macrophages <- StackedVlnPlot(obj = sub_integrated_data, 
           features = features_to_plot,
           my_levels = my_levels,
           group.by = "Celltype_raw1")
vlnplt_monocytes_macrophages

ggplot2::ggsave(vlnplt_monocytes_macrophages,filename = '/public3/DSC/single_cell/Result/figer_new/F2.0_vlnplt_celltype1.png',width = 24,height = 8)
ggplot2::ggsave(vlnplt_monocytes_macrophages,filename = '/public3/DSC/single_cell/Result/figer_new/F2.0_vlnplt_celltype1.png.pdf',width = 24,height = 8, device = cairo_pdf)




features_to_plot <- c(
  # --- Monocytes ---
  "FCN1", "S100A8", # Classical Mono (S100A8 比 A12 更干净，SELL 是金标准)
  # --- Foam ---
  "SPP1",'GPNMB',#"CCL7", # Foam1 (趋化/早期)
  # --- marcophage ---
  'SELENOP','F13A1'     
  
)

my_levels <- c(
  "Monocyte", "LAM", "Macrophage"
  
)
# 2. 将 Celltype_raw1 转换为指定顺序的因子 (Factor)
vlnplt_monocytes_macrophages <- StackedVlnPlot(obj = sub_integrated_data, 
           features = features_to_plot,
           my_levels = my_levels,
           group.by = "Celltype_raw")
vlnplt_monocytes_macrophages
ggplot2::ggsave(vlnplt_monocytes_macrophages,filename = '/public3/DSC/single_cell/Result/figer_new/F2.0_vlnplt_celltype.pdf',width = 9,height = 8, device = cairo_pdf)

# 1. 提取并修正详细分类 AC_PA
# 重点：将 GSE216860 明确归类为 Control
metadata <- sub_integrated_data@meta.data %>%
  mutate(orig.ident = as.character(orig.ident)) %>%
  mutate(AC_PA_Fixed = case_when(
    # --- 核心对照组 (Control) ---
    # GSE216860 是正常人升主动脉
    # GSE155468 和 GSE213740 是升主动脉瘤的健康对照
    grepl("GSE216860|GSE155468", orig.ident) ~ "Control AscAorta",
    
    # --- 斑块核心组 (AC) ---
    grepl("GSE131778", orig.ident) ~ "Coronary Atherosclerotic Core",
    grepl("GSE224273|GSE234077|GSE253903|AC$", orig.ident, ignore.case = TRUE) ~ "Carotid Atherosclerotic Core",
    
    # --- 稳定/不稳定斑块组 ---
    grepl("GSE247238|GSE260657", orig.ident) & grepl("Stable", AC_PA) ~ "Carotid Stable Plaque",
    grepl("GSE247238|GSE260657", orig.ident) & grepl("Unstable", AC_PA) ~ "Carotid Unstable Plaque",
    
    # --- 邻近区域组 (PA) ---
    # 仅保留真实的邻近样本（如本地 patientPA）
    grepl("PA$", orig.ident) & !grepl("GSE216860|GSE213740", orig.ident) ~ "Carotid Proximal Adjacent",
    
    TRUE ~ as.character(AC_PA)
  ))

# 将修正结果写回
sub_integrated_data$AC_PA <- metadata$AC_PA_Fixed

# 2. 基于修正后的 AC_PA，统一 Sample_Type 为三大类
sub_integrated_data@meta.data <- sub_integrated_data@meta.data %>%
  mutate(Sample_Type = case_when(
    grepl("Control", AC_PA, ignore.case = TRUE) ~ "Proximal Adjacent",
    grepl("Adjacent", AC_PA, ignore.case = TRUE) ~ "Proximal Adjacent",
    grepl("Core|Plaque", AC_PA, ignore.case = TRUE) ~ "Atherosclerotic Core",
    TRUE ~ "Unknown"
  ))

# 3. 设置因子顺序：按照生物学进展排列 (Control -> PA -> AC)
# 这样在绘图时，Control 会作为基准出现在最左侧
sub_integrated_data$Sample_Type <- factor(sub_integrated_data$Sample_Type, 
      levels = c("Proximal Adjacent", "Atherosclerotic Core"))

# 4. 最终验证
print("修正后的详细分类与数据集对应关系：")
table(sub_integrated_data$AC_PA, sub_integrated_data$orig.ident)

print("修正后的宏观分类统计：")
table(sub_integrated_data$Sample_Type)
p_AC_PA <- DimPlot(sub_integrated_data, 
       group.by = "orig.ident", 
       split.by = "AC_PA",  # 核心修改：按 AC_PA 分面
       ncol = 3,      # 设置分面列数，让 AC/PA/Control 横向排开
       repel = TRUE,  # 启用 ggrepel 算法防止标签重叠
       # label = TRUE, 
       label = FALSE,       # 彻底关闭图上的文字标签
       label.size = 6,      # 分面后建议缩小 label，否则会显得拥挤
       reduction = 'umap', 
       pt.size = 1) + 
  theme_dr() + 
  theme(panel.grid = element_blank(),
  text = element_text(family = "Arial"),
  plot.title = element_blank(),
  strip.text = element_text(size = 20, face = "bold"), # 修改分面标题的大小
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 15)) +       # 分面后建议缩小图例字体
  coord_flip()

# 打印图片
print(p_AC_PA)
# 1. 专门提取 "Control AscAorta" 的细胞子集
control_asc_data <- subset(x = sub_integrated_data, subset = AC_PA == "Control AscAorta")

# 2. 【关键步骤】清理冗余的因子层级 (Factor levels)
# 这一步确保在绘图时，只有 Control 组实际包含的 orig.ident 才会出现，避免生成一堆空白的分面。
control_asc_data$orig.ident <- droplevels(as.factor(control_asc_data$orig.ident))

# 3. 绘制分面 UMAP 图
p_control <- DimPlot(control_asc_data, 
   group.by = "orig.ident", 
   split.by = "orig.ident",     # 核心：按 orig.ident 分面
   ncol = 2,        # 根据您目前剩下的 Control 样本数量（估计是2-3个），设置合适的列数
   label = FALSE,   # 建议关闭标签，直接看顶部标题和图例即可
   reduction = 'umap', 
   pt.size = 1) + 
  theme_dr() + 
  theme(panel.grid = element_blank(),
  text = element_text(family = "Arial"),
  plot.title = element_blank(),       # 关闭默认的整个大图标题
  strip.text = element_text(size = 16, face = "bold"), # 分面框顶部的样本名大小
  legend.position = "right",
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12)) + 
  coord_flip()      # 保留您之前的坐标轴翻转习惯

# 打印图片
print(p_control)
rm(control_asc_data)
gc()
# 重新绘制UMAP
p_clusters <- DimPlot(sub_integrated_data, group.by = "rpca_clusters", label=T, label.size=10, reduction='umap', pt.size = 1) + # 调整 pt.size 的值
  theme_dr() + theme(panel.grid=element_blank(),text = element_text(family = "Arial"),
   plot.title = element_blank()) + NoLegend()+coord_flip()



p_celltype <- DimPlot(sub_integrated_data, 
     group.by = "Celltype_raw", 
     label = TRUE, 
     repel = TRUE,    # 关键修改：将 repel = TRUE 移到 DimPlot() 内部
     label.size = 10, 
     reduction = 'umap', 
     pt.size = 1) + 
  theme_dr() + 
  theme(panel.grid = element_blank(),
  text = element_text(family = "Arial"),
  plot.title = element_blank()) + 
  NoLegend() + 
  coord_flip()




p_celltype1 <- DimPlot(sub_integrated_data, 
     group.by = "Celltype_raw1", 
     label = TRUE, 
     repel = TRUE,    # 关键修改：将 repel = TRUE 移到 DimPlot() 内部
     label.size = 10, 
     reduction = 'umap', 
     pt.size = 1) + 
  theme_dr() + 
  theme(panel.grid = element_blank(),
  text = element_text(family = "Arial"),
  plot.title = element_blank()) + 
  NoLegend() + 
  coord_flip()
sub_integrated_data$Sample_Type <- factor(
  sub_integrated_data$Sample_Type, 
  levels = c(
    "Atherosclerotic Core",   # 放在第一位（图例最上方/第一个颜色）
    "Proximal Adjacent"      # 第二位
  )
)

p_AC_PA <- DimPlot(sub_integrated_data, 
       group.by = "Sample_Type", 
       split.by = "Sample_Type", # 添加这一行实现分面
       label = TRUE, 
       label.size = 10, 
       reduction = 'umap', 
       pt.size = 1) + 
  theme_dr() + 
  theme(panel.grid = element_blank(),
  text = element_text(family = "Arial"),
  plot.title = element_blank(),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 30),
  strip.text = element_text(size = 20, family = "Arial")) + # 可选：调整分面标题的大小和字体
  coord_flip()

p_AC_PA
p_AC_PA <- DimPlot(sub_integrated_data, 
       group.by = "Sample_Type", 
       label=T, 
       label.size=10, 
       reduction='umap', 
       pt.size = 1) + # 调整 pt.size 的值
  theme_dr() + 
  theme(panel.grid=element_blank(),
  text = element_text(family = "Arial"),
  plot.title = element_blank(),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 30))+coord_flip()

p1 <-  p_celltype1 | p_AC_PA

print(p1)

ggplot2::ggsave(p1,filename = '/public3/DSC/single_cell/Result/figer_new/F2.1_MM_UMAP_celltype_cluster.png',width = 21,height = 8)
ggplot2::ggsave(p1,filename = '/public3/DSC/single_cell/Result/figer_new/F2.1_MM_UMAP_celltype_cluster.pdf',width = 21,height = 8, device = cairo_pdf)

p_cluster <- DimPlot(sub_integrated_data, group.by = "seurat_clusters", label=T, label.size=10, pt.size = 1, reduction='umap') + 
  theme_dr() + theme(panel.grid=element_blank(), 
   plot.title = element_blank()) + NoLegend() 

p_cluster <-  p_clusters | p_celltype
p_cluster
ggplot2::ggsave(p_cluster,filename = '/public3/DSC/single_cell/Result/figer_new/S1.1_MM_UMAP_celltype_cluster.png',width = 18,height = 8)
ggplot2::ggsave(p_cluster,filename = '/public3/DSC/single_cell/Result/figer_new/S1.1_MM_UMAP_celltype_cluster.pdf',width = 18,height = 8, device = cairo_pdf)

avg_exp_celltype <- AverageExpression(
  sub_integrated_data,
  group.by = "Celltype_raw1",  # Key change: group by Celltype_raw1
  assays = "RNA",
  slot = "data"  # Use normalized data
)$RNA

# Check dimensions (genes x number of unique Celltype_raw1)
print("Dimensions of average expression matrix (genes x Celltype_raw1):")
print(dim(avg_exp_celltype))
print("Column names (Celltype_raw1 types):")
print(colnames(avg_exp_celltype))

# Ensure there are at least 2 cell types for correlation
if (ncol(avg_exp_celltype) < 2) {
  stop("Need at least two groups in Celltype_raw1 to calculate correlations. Found: ", ncol(avg_exp_celltype))
}

# 2. Calculate the full correlation matrix for Celltype_raw1
# Hmisc::rcorr requires a matrix with at least 2 columns.
# The result 'r' will have rownames and colnames derived from colnames(avg_exp_celltype)
cor_matrix_celltype <- Hmisc::rcorr(as.matrix(avg_exp_celltype), type = "spearman")$r

# Check dimensions (number of Celltype_raw1 x number of Celltype_raw1)
print("Dimensions of correlation matrix (Celltype_raw1 x Celltype_raw1):")
print(dim(cor_matrix_celltype))
print("Row names of correlation matrix:")
print(rownames(cor_matrix_celltype))
print("Column names of correlation matrix:")
print(colnames(cor_matrix_celltype))


# 3. (Optional) Filtering or Reordering based on Celltype_raw1 if needed
# If you need to remove specific Celltype_raw1 or set a specific order:
# Example:
# celltypes_to_remove <- c("Unknown_Celltype", "Doublets")
# celltypes_to_keep <- setdiff(colnames(cor_matrix_celltype), celltypes_to_remove)
# cor_matrix_subset_celltype <- cor_matrix_celltype[celltypes_to_keep, celltypes_to_keep]

# Example of specific ordering:
# desired_celltype_order <- c("Astrocyte", "Neuron", "Microglia", "Oligodendrocyte")
# # Ensure all names in desired_celltype_order are present in colnames(cor_matrix_celltype)
# desired_celltype_order_filtered <- intersect(desired_celltype_order, colnames(cor_matrix_celltype))
# if (length(desired_celltype_order_filtered) > 1) {
#   ordered_matrix_celltype <- cor_matrix_celltype[desired_celltype_order_filtered, desired_celltype_order_filtered]
#   matrix_to_plot <- ordered_matrix_celltype
#   cluster_r <- FALSE # If providing specific order, don't cluster
#   cluster_c <- FALSE
# } else {
#   print("Desired order not fully applicable or results in too few cell types, using original matrix and clustering.")
#   matrix_to_plot <- cor_matrix_celltype
#   cluster_r <- TRUE
#   cluster_c <- TRUE
# }

# For this example, we'll use the full matrix and let pheatmap cluster
matrix_to_plot <- cor_matrix_celltype
cluster_r <- TRUE # Allow pheatmap to cluster rows
cluster_c <- TRUE # Allow pheatmap to cluster columns

# 4. Generate the heatmap
my_colors <- colorRampPalette(c("#4178a7", "white", "#e2201c"))(100)

# Dynamically set breaks based on the data range in the correlation matrix
min_val <- min(matrix_to_plot, na.rm = TRUE)
max_val <- max(matrix_to_plot, na.rm = TRUE)
print(paste("Min correlation:", round(min_val, 3), "Max correlation:", round(max_val, 3)))

# If all values are very similar and close to 1, your original breaks might be fine.
# Otherwise, adjust breaks to span the actual range of correlation values or a relevant sub-range.
# breaks_list <- seq(min_val, max_val, length.out = 100)
# Or, if you want to focus on a specific range, e.g., positive correlations:
if (max_val > 0.5) { # Ensure there's a reasonable positive correlation range
  breaks_list <- seq(max(min_val, 0), max_val, length.out = 100) # Focus on 0 to max, or min_val to max_val
  # If your previous range 0.88 to 1 is still meaningful for this data:
  # breaks_list <- seq(0.88, 1, length.out = 100)
} else { # If max correlation is low, adjust accordingly
  breaks_list <- seq(min_val, max_val, length.out = 100)
}
# If the matrix contains negative correlations and you want to show them distinctly:
# breaks_list <- c(seq(min_val, -0.01, length.out = 49), seq(-0.01, 0.01, length.out=2), seq(0.01, max_val, length.out = 49))


heatmap_plot_celltype <- pheatmap(
  matrix_to_plot,
  color = my_colors,
  breaks = breaks_list,
  cluster_rows = cluster_r,
  cluster_cols = cluster_c,
  display_numbers = TRUE,
  number_color = "black",
  angle_col = 45,   # Keep column labels tilted
  fontfamily = "Arial",
  fontsize_row = 18,      # Adjust as needed, 18 might be too large if many cell types
  fontsize_col = 18,      # Adjust as needed
  fontsize = 18,    # General fontsize for legend, etc.
  fontsize_number = 18,    # Fontsize for numbers in cells, adjust as needed
  legend = FALSE,
  main = "" # Correlation Heatmap of Avg. Expression by Celltype_raw1 (Spearman)
)
ggplot2::ggsave(heatmap_plot_celltype,filename = '/public3/DSC/single_cell/Result/figer_new/F2.4_heatmap.pdf',width = 9,height = 8, device = cairo_pdf)

broad_markers <- FindAllMarkers(
  sub_integrated_data,
  group.by = "Celltype_raw",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 提取并查看每个大类下最具代表性的前 5 个 Marker 基因
top5_broad_markers <- broad_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC) 

print(as.data.frame(top5_broad_markers))

p2 <- FeaturePlot(
  sub_integrated_data,
  features = c(    # --- Monocytes ---
    "FCN1", "S100A8", # Classical Mono (S100A8 比 A12 更干净，SELL 是金标准)
    # --- Foam ---
    "SPP1",'GPNMB',#"CCL7", # Foam1 (趋化/早期)
    # --- marcophage ---
    'SELENOP','F13A1'    
  ),
  reduction = "umap",    # 🌟 核心修复：明确指定使用标准的 umap
  pt.size = 1,
  ncol = 3               # 🌟 核心修改：指定排版为 3 列
) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,         # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,      # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 45, 5), "pt")
  ) &
  NoLegend() &
  coord_flip()

print(p2)

#单核，巨噬细胞marker

p3 <- VlnPlot(
  sub_integrated_data, 
  group.by = "Celltype_raw",
  features = c(    # --- Monocytes ---
    "FCN1", "S100A8", # Classical Mono (S100A8 比 A12 更干净，SELL 是金标准)
    # --- Foam ---
    "SPP1",'GPNMB',#"CCL7", # Foam1 (趋化/早期)
    # --- marcophage ---
    'SELENOP','F13A1'    
  ),
  pt.size = 0,
  ncol = 2
) +
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )

print(p3)

# 🌟 核心修改：调整 3 列排版对应的最佳长宽比例
ggplot2::ggsave(p2, filename = '/public3/DSC/single_cell/Result/figer_new/F2.2_marker_feature.png', width = 24, height = 16)
ggplot2::ggsave(p2, filename = '/public3/DSC/single_cell/Result/figer_new/F2.2_marker_feature.pdf', width = 24, height = 16, device = cairo_pdf)

# p3 保持不变
ggplot2::ggsave(p3, filename = '/public3/DSC/single_cell/Result/figer_new/F2.2_marker_VlnPlot.png', width = 8, height = 8)
ggplot2::ggsave(p3, filename = '/public3/DSC/single_cell/Result/figer_new/F2.2_marker_VlnPlot.pdf', width = 8, height = 8, device = cairo_pdf)

#胚胎来源巨噬细胞与单细胞来源巨噬细胞：'TREM2','FOLR2'！
#cluster5 MHEM marker 基因为 MMP9,NR1H3(LXRα)
##cluster3,1 M1 marker 基因为 IL1B
#CD14+ CD16- CD62L(SELL)+ CM, cluster6


p3 <- FeaturePlot(sub_integrated_data,
      features = c("LYVE1", "PDGFC"),
      reduction = "umap", 
      pt.size = 1,ncol=3) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,   # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,       # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 15, 5), "pt")
  ) &
  NoLegend()&
  coord_flip()#M2a 免疫抑制巨噬细胞
p4 <- VlnPlot(sub_integrated_data,features = c("LYVE1", "PDGFC"),pt.size = 0) +
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )#M2a 免疫抑制巨噬细胞
p3
p4

ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.2_M2_marker_feature.png',width = 18,height = 8)
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.2_M2_marker_feature.pdf',width = 18,height = 8, device = cairo_pdf)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.2_M2_marker_VlnPlot.png',width = 12,height = 4)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.2_M2_marker_VlnPlot.pdf',width = 12,height = 4, device = cairo_pdf)

p3 <- FeaturePlot(sub_integrated_data,features = c("S100A12", "SELL"),  reduction = "umap", 
      pt.size = 1,ncol=3) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,   # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,       # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 15, 5), "pt")
  ) &
  NoLegend()&
  coord_flip()#CM
p4 <- VlnPlot(sub_integrated_data,features = c("S100A12", "SELL"),pt.size = 0)+
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )#CM
print(p3)
print(p4)
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.3_MO_marker_feature.png',width = 18,height = 8)
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.3_MO_marker_feature.pdf',width = 18,height = 8, device = cairo_pdf)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.3_MO_marker_VlnPlot.png',width = 12,height = 4)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.3_MO_marker_VlnPlot.pdf',width = 12,height = 4, device = cairo_pdf)

#泡沫marker1:LIPA+ APOC1+ CD36 APOC1 CD9 TREM2
#泡沫marker2:OLR1 PLIN2+ MARCO IL1RN CCL2
#泡沫marker3:FABP5 CTSB SPP1 CD36

p3 <- FeaturePlot(sub_integrated_data,features = c("PHLDA1",'CCL2'),
      reduction = "umap", 
      pt.size = 1,ncol=3) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,   # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,       # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 15, 5), "pt")
  ) &
  NoLegend()&
  coord_flip()#FC1
p4 <- VlnPlot(sub_integrated_data,
  group.by = "Celltype_raw1",
  features = c("PHLDA1",'CCL2'),
  pt.size = 0)+
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )#FC1
print(p3)
print(p4)

ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker2_feature.png',width = 18,height = 8)
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker2_feature.pdf',width = 18,height = 8, device = cairo_pdf)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker2_VlnPlot.png',width = 12,height = 4)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker2_VlnPlot.pdf',width = 12,height = 4, device = cairo_pdf)

p3 <- FeaturePlot(sub_integrated_data,features = c("CD36","MGLL"),
      reduction = "umap", 
      pt.size = 1,ncol=3) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,   # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,       # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 15, 5), "pt")
  ) &
  NoLegend()&
  coord_flip()#FC1
p4 <- VlnPlot(sub_integrated_data,group.by = "Celltype_raw1",features = c("CD36","MGLL"),pt.size = 0)+
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )#FC2

p3
p4
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker1_feature.png',width = 27,height = 8)
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker1_feature.pdf',width = 27,height = 8, device = cairo_pdf)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker1_VlnPlot.png',width = 12,height = 4)
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker1_VlnPlot.pdf',width = 12,height = 4, device = cairo_pdf)

p3 <- FeaturePlot(sub_integrated_data,features = c( "TREM2","PLAU"),
      reduction = "umap", 
      pt.size = 1,ncol=3) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,   # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,       # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 15, 5), "pt")
  ) &
  NoLegend()&
  coord_flip()#FC1
p4 <- VlnPlot(sub_integrated_data,group.by = "Celltype_raw1",features = c( "TREM2","PLAU"),pt.size = 0)+
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )#FC3
print(p3)
print(p4)

ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_feature.png',width = 18,height = 8,units = "in")
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_feature.pdf',width = 18,height = 8, device = cairo_pdf,units = "in")
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_VlnPlot.png',width = 12,height = 4,units = "in")
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_VlnPlot.pdf',width = 12,height = 4, device = cairo_pdf,units = "in")


p3 <- FeaturePlot(sub_integrated_data,features = c( "IFIT1", "ISG15"),
      reduction = "umap", 
      pt.size = 1,ncol=3) &  # 使用 & 使修改应用于所有子图
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(  # 修改标题样式
      family = "Arial",  # 字体
      face = "bold",     # 加粗
      size = 40,   # 字号
      hjust = 0.5,       # 水平居中
      vjust = 0.5 ,       # 垂直居中
      margin = margin(b = 5) 
    ),
    text = element_text(family = "Arial"),  # 全局字体
    plot.margin = unit(c(15, 5, 15, 5), "pt")
  ) &
  NoLegend()&
  coord_flip()#FC1
p4 <- VlnPlot(sub_integrated_data,group.by = "Celltype_raw1",features = c( "IFIT1", "ISG15"),pt.size = 0)+
  # 3. 统一设置Arial字体
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    axis.title = element_blank(),  # 坐标轴标题
    axis.text = element_blank(),  # 坐标轴刻度
    legend.text = element_blank(),  # 图例文本
    strip.text = element_blank(),  # 移除子图标题文本
    strip.background = element_blank() # 移除子图标题背景
    
  )#ISG+
print(p3)
print(p4)

ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_feature.png',width = 18,height = 8,units = "in")
ggplot2::ggsave(p3,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_feature.pdf',width = 18,height = 8, device = cairo_pdf,units = "in")
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_VlnPlot.png',width = 12,height = 4,units = "in")
ggplot2::ggsave(p4,filename = '/public3/DSC/single_cell/Result/figer_new/S1.4_FC_marker3_VlnPlot.pdf',width = 12,height = 4, device = cairo_pdf,units = "in")


cat(">>> 正在生成所有核心 Marker 的全局 FeaturePlot 大图...\n")

# 1. 汇总所有需要展示的基因 (共 12 个)
# all_marker_genes <- c(
#   "S100A12", "SELL",   # CM / 经典单核
#   "IFIT1", "ISG15",     # ISG+ / 干扰素响应单核
#   "LYVE1", "PDGFC",    # M2a / 免疫抑制
#   "TREM2", "PLAU",     # LAM / 脂质相关巨噬细胞
#   "PHLDA1", "CCL2",    # Foam 1 / 早期趋化泡沫
#   "CD36", "MGLL",      # Foam 2 / 成熟脂滴泡沫
#   "TREM2", "PLAU"     # LAM / 脂质相关巨噬细胞
# 
# )

# # 1. 汇总所有需要展示的基因 (共 20 个)
all_marker_genes <- c(
  # --- Monocytes ---
  "S100A12", "SELL", # Classical Mono (S100A8 比 A12 更干净，SELL 是金标准)
  "MTRNR2L8", "CCL5", # Inflammatory Mono (PTGS2 是炎症风暴核心)
  "IFIT1", "ISG15", # ISG+ Mono (RSAD2 背景极低)
  "CLEC10A", "THBS1",   # Repair Mono (生长因子)

  # --- Metabolic / Foam ---
  "PHLDA1",'CCL2',#"CCL7", # Foam1 (趋化/早期)
  "CD36", "MGLL", # Foam2 (成熟/脂滴)
  "TREM2",   "PLAU", # LAM (溶酶体活跃/分解脂质)

  # --- Resident ---
  'F13A1','SELENOP',
  "CX3CR1", "C3",   # CX3CR1+ TRM
  "LYVE1", "PDGFC"  # LYVE1+ TRM (IGF1 特异性很高)
)
# 2. 一次性绘制 12 个基因
p_all_features <- FeaturePlot(
  sub_integrated_data,
  features = all_marker_genes,
  reduction = "umap", 
  pt.size = 1,
  order = TRUE,
  ncol = 7             # 20个基因，设置为 5列 × 4行 排版
) &  
  theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      family = "Arial",  
      face = "bold",     
      size = 35,         # 考虑到图变多了，标题字号设为35保持协调
      hjust = 0.5,       
      vjust = 0.5,       
      margin = margin(b = 10) 
    ),
    text = element_text(family = "Arial"), 
    plot.margin = unit(c(15, 10, 15, 10), "pt") # 增加各子图之间的呼吸空间
  ) &
  NoLegend() &
  coord_flip()

# 打印预览
print(p_all_features)

# 3. 保存为超大尺寸的高清图片
# 因为包含 5列x4行 共20张子图，我们将画布等比例放大到 24 x 16
out_dir <- '/public3/DSC/single_cell/Result/figer_new/'

ggplot2::ggsave(p_all_features, 
                filename = paste0(out_dir, 'S1_All_Markers_FeaturePlot_Combined.png'), 
                width = 34, 
                height = 12)

ggplot2::ggsave(p_all_features, 
                filename = paste0(out_dir, 'S1_All_Markers_FeaturePlot_Combined.pdf'), 
                width = 34, 
                height = 12, 
                device = cairo_pdf)

# ==========================================
# 堆叠柱状图

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)
library(scales)
library(ggpubr)
library(gtable)

# ==========================================
# 1. 全局配置 (确保两张图风格绝对统一的核心)
# ==========================================
# 统一颜色字典 
cell_colors <- c(
  "Classical Mono"     = "#8c564b", 
  "Inflammatory Mono"  = "#b15928", 
  "ISG+ Mono"          = "#bcbd22", 
  "Non-classical Mono" = "#c7c7c7", 
  "Foam cells1"        = "#d62728", 
  "Foam cells2"        = "#ff7f0e", 
  "LAM"                = "#9467bd", 
  "Transitional Mac"   = "#2ca02c", 
  "CX3CR1+ TRM"        = "#e377c2", 
  "LYVE1+ TRM"         = "#f7b6d2",
  "TrMs"               = "#e377c2", 
  "CM"                 = "#2ca02c", 
  "Macrophage"         = "#9467bd", 
  "Monocyte"           = "#8c564b", 
  "cDC1"               = "#1f77b4"
)

# 统一表格样式 (居中对齐，黑底白字表头)
my_table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontfamily = "Arial", fontsize = 10, hjust = 0.5, x = 0.5),
    bg_params = list(fill = c("white", "#f7f7f7"))
  ),
  colhead = list(
    bg_params = list(fill = "#404040"),
    fg_params = list(col = "white", fontface = "bold", fontfamily = "Arial", fontsize = 10, hjust = 0.5, x = 0.5)
  )
)

# ==========================================
# 2. 绘制图 F (PA vs AC)
# ==========================================
table_counts_F <- table(sub_integrated_data@meta.data$Sample_Type, sub_integrated_data@meta.data$Celltype_raw1)
Df_F <- as.data.frame(table_counts_F[, colSums(table_counts_F) > 0])
colnames(Df_F) <- c("Status", "CellType", "Count")

# 设定顺序：PA在左(上)，AC在右(下)
Df_F$Status <- factor(Df_F$Status, levels = c("Proximal Adjacent", "Atherosclerotic Core"))

Df_F <- Df_F %>% 
  group_by(Status) %>% 
  mutate(Proportion = Count / sum(Count))

p_f <- ggplot(Df_F, aes(x = Status, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) + 
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors, breaks = names(cell_colors)) +
  labs( x = "Group", y = "Proportion") +#title = "MPS Subtype Composition (PA vs AC)",
  theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.margin = margin(t=20, r=20, b=10, l=20, unit="pt"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

table_data_F <- Df_F %>%
  group_by(Status, CellType) %>%
  dplyr::summarise(Proportion = paste0(round(mean(Proportion) * 100, 1), "%"), .groups = "drop") %>%
  pivot_wider(names_from = CellType, values_from = Proportion)

n_cols_F <- ncol(table_data_F)
split_idx_F <- ceiling((n_cols_F - 1) / 2) + 1
grob_f1 <- tableGrob(table_data_F[, 1:split_idx_F], rows = NULL, theme = my_table_theme)
grob_f2 <- tableGrob(table_data_F[, c(1, (split_idx_F+1):n_cols_F)], rows = NULL, theme = my_table_theme)
table_grob_f <- gtable_rbind(grob_f1, grob_f2)

final_f <- grid.arrange(p_f, table_grob_f, nrow = 2, heights = c(0.7, 0.3))


# ==========================================
# 3. 绘制图 G (Stable Plaque vs Vulnerable Plaque)
# ==========================================
symptom_subset <- subset(
  sub_integrated_data,
  subset = Source_GSE %in% c("GSE260657", "GSE247238") & grepl("Stable|Unstable", AC_PA)
)

# 核心修改：更新标签名称
symptom_subset$Status <- ifelse(
  grepl("Stable", symptom_subset$AC_PA), 
  "Stable Plaque",       
  "Vulnerable Plaque"
)
# 设定因子顺序：Stable 在左，Vulnerable 在右
symptom_subset$Status <- factor(symptom_subset$Status, levels = c("Stable Plaque", "Vulnerable Plaque"))

table_counts_G <- table(symptom_subset$Status, symptom_subset$Celltype_raw1)
Df_G <- as.data.frame(table_counts_G[, colSums(table_counts_G) > 0])
colnames(Df_G) <- c("Status", "CellType", "Count")

Df_G <- Df_G %>% 
  group_by(Status) %>% 
  mutate(Proportion = Count / sum(Count))

p_g <- ggplot(Df_G, aes(x = Status, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) + 
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors, breaks = names(cell_colors)) +
  labs( x = "Group", y = "Proportion") + # 更新了标题
  theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.margin = margin(t=20, r=20, b=10, l=20, unit="pt"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

table_data_G <- Df_G %>%
  group_by(Status, CellType) %>%
  dplyr::summarise(Proportion = paste0(round(mean(Proportion) * 100, 1), "%"), .groups = "drop") %>%
  pivot_wider(names_from = CellType, values_from = Proportion)

n_cols_G <- ncol(table_data_G)
if (n_cols_G > 6) {
  split_idx_G <- ceiling((n_cols_G - 1) / 2) + 1
  grob_g1 <- tableGrob(table_data_G[, 1:split_idx_G], rows = NULL, theme = my_table_theme)
  grob_g2 <- tableGrob(table_data_G[, c(1, (split_idx_G+1):n_cols_G)], rows = NULL, theme = my_table_theme)
  table_grob_g <- gtable_rbind(grob_g1, grob_g2)
} else {
  table_grob_g <- tableGrob(table_data_G, rows = NULL, theme = my_table_theme)
}

final_g <- grid.arrange(p_g, table_grob_g, nrow = 2, heights = c(0.7, 0.3))


# ==========================================
# 4. 统一尺寸输出
# ==========================================
save_path <- "/public3/DSC/single_cell/Result/figer_new/"
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# 图 F 保存
ggsave(paste0(save_path, "Fig_F_PA_vs_AC.pdf"), final_f, width = 10, height = 8, device = cairo_pdf)
ggsave(paste0(save_path, "Fig_F_PA_vs_AC.png"), final_f, width = 10, height = 8, dpi = 300)

# 图 G 保存
ggsave(paste0(save_path, "Fig_G_Plaque_Status.pdf"), final_g, width = 10, height = 8, device = cairo_pdf)
ggsave(paste0(save_path, "Fig_G_Plaque_Status.png"), final_g, width = 10, height = 8, dpi = 300)

message("已更新图 G 标签为 Stable/Vulnerable Plaque 并完成保存。")

# ==========================================


Mono_markers <- FindAllMarkers(sub_integrated_data, 
       only.pos = TRUE, 
       min.pct = 0.25,group.by = "Celltype_raw1",
       logfc.threshold = 0.25)
Mono_cluster_markers <- FindAllMarkers(sub_integrated_data, 
         only.pos = TRUE, 
         min.pct = 0.25,group.by = "seurat_clusters",
         logfc.threshold = 0.25)
write.csv(Mono_markers,"/public3/DSC/single_cell/Result/figer_new/ST2_Ma_Mo_Marker.csv")
write.csv(Mono_cluster_markers,"/public3/DSC/single_cell/Result/figer_new/ST2_Mono_cluster_markers.csv")
Mono_markers <- read.csv("/public3/DSC/single_cell/Result/figer_new/ST2_Ma_Mo_Marker.csv")
# 排序顺序
target_order <- c(
  # --- 1. 单核细胞起点 & 过渡 ---
  "Classical Mono",     # Cluster 10: 始祖
  
  # --- 2. 功能性单核亚群 (炎症/抗病毒/修复) ---
  "Inflammatory Mono",  # Cluster 12: 炎症风暴
  "ISG+ Mono",    # Cluster 11/13: 干扰素反应
  "Non-classical Mono", # Cluster 2:  修复/M2样前体
  
  # --- 4. 代谢/病理巨噬细胞 (泡沫化路线) ---
  "Foam cells1",  # Cluster 7:  早期/单核样泡沫
  "Foam cells2",  # Cluster 1:  成熟泡沫
  "LAM",    # Cluster 4:  脂质相关/TREM2+
  
  # --- 5. 组织驻留巨噬细胞 (维稳卫士) ---
  "Transitional Mac" ,  # Cluster 9:  最典型的M2驻留
  "CX3CR1+ TRM",  # Cluster 0:  巡逻/抗原呈递
  "LYVE1+ TRM"    # Cluster 6:  血管旁驻留
)


# 将Celltype_raw1的因子水平重置为目标顺序 (由于 Seurat 默认从下往上画 Y 轴，所以使用 rev)
sub_integrated_data@meta.data$Celltype_raw1 <- factor(
  sub_integrated_data@meta.data$Celltype_raw1,
  levels = rev(target_order) 
)

# === 修改区域开始 ===
# 1. 强制将 marker 数据框按照 target_order 排序
Top5_MARKERS <- Mono_markers %>%
  mutate(cluster = factor(cluster, levels = target_order)) %>%
  # 关键修复：先按细胞群排序，再按 avg_log2FC 降序严格排序！
  arrange(cluster, desc(avg_log2FC)) %>% 
  group_by(cluster) %>%
  # 弃用 top_n，改用 slice_head 确保提取的基因顺序与 arrange 后的完全一致
  slice_head(n = 5) %>% 
  ungroup()

# 2. 此时提取的基因就会严格遵循从 Classical Mono 到 LYVE1+ TRM 的顺序
DotPlot_genes <- unique(Top5_MARKERS$gene)

p6 <- DotPlot(sub_integrated_data,
  group.by = "Celltype_raw1",
  features = DotPlot_genes,
  cluster.idents = FALSE,
  dot.scale = 10) +
  RotatedAxis() +
  theme(
    text = element_text(size = 20, family = "Arial"),  # 设置所有文本的默认大小和字体
    axis.title.x = element_text(size = 0), # 设置X轴标签的大小和字体
    axis.title.y = element_text(size = 0), # 设置Y轴标签的大小和字体
    axis.text.x = element_text(size = 18),  # 设置X轴刻度标签的大小和字体
    axis.text.y = element_text(size = 18)   # 设置Y轴刻度标签的大小和字体
  )

print(p6)

ggplot2::ggsave(p6,filename = '/public3/DSC/single_cell/Result/figer_new/F2.4_DotPlot.png',width = 12,height = 8)
ggplot2::ggsave(p6,filename = '/public3/DSC/single_cell/Result/figer_new/F2.4_DotPlot.pdf',width = 18,height = 8,device = cairo_pdf)
saveRDS(sub_integrated_data,file = '/public3/DSC/single_cell/Result/figer_new/group_result.rds')

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
#sub_integrated_data <- readRDS("/public3/DSC/single_cell/Result/figer_new/group_result.rds")
align_umap_plots <- function(plot_AC, plot_PA) {
  # 提取两组数据的UMAP坐标范围
  ac_range <- layer_scales(plot_AC)$x$range$range  # 获取AC组的x轴范围
  pa_range <- layer_scales(plot_PA)$x$range$range  # 获取PA组的x轴范围
  x_min <- min(ac_range[1], pa_range[1])    # 计算x轴最小值
  x_max <- max(ac_range[2], pa_range[2])    # 计算x轴最大值
  
  # 同理获取y轴范围
  y_range <- range(layer_scales(plot_AC)$y$range$range,
       layer_scales(plot_PA)$y$range$range)
  y_min <- y_range[1]
  y_max <- y_range[2]
  
  # 应用统一坐标范围和比例
  plot_AC <- plot_AC + 
    coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme(
      aspect.ratio = 1,     # 确保画布为正方形
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

data <- readRDS("/public3/DSC/single_cell/Result/figer_new/group_result.rds")
#data <- sub_integrated_data

enableWGCNAThreads(nThreads = 8)
# set up seurat object for WGCNA
data <- SetupForWGCNA(
  data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "wgcna" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
data <- MetacellsByGroups(
  seurat_obj = data,
  group.by = "Sample_Type", # specify the columns in seurat_obj@meta.data to group by
  k = 20, # nearest-neighbors parameter 通常是20-75，10万个细胞可以用50
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Sample_Type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
data <- NormalizeMetacells(data)

# set up the expression matrix
data <- SetDatExpr(
  data,
  group_name = c("Atherosclerotic Core","Proximal Adjacent"), # the name of the group of interest in the group.by column
  group.by='Sample_Type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
data <- TestSoftPowers(
  data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(data)

my_theme <- theme(
  text = element_text(size = 14, family = "Arial"), # 设置整体字体大小和字体
  axis.title = element_text(size = 20),       # 调整轴标题大小
  axis.text = element_text(size = 15),  # 调整轴标签大小
  legend.title = element_text(size = 20),     # 调整图例标题大小
  legend.text = element_text(size = 15),      # 调整图例标签大小
  title = element_text(size = 18)       # 调整图标题大小
)

# 遍历 plot_list 中的每个 ggplot 对象并添加主题
modified_plot_list <- lapply(plot_list, function(plot) {
  return(plot + my_theme)
})
# assemble with patchwork
library(patchwork)
p <- wrap_plots(plot_list, ncol=2)
print(p)
ggplot2::ggsave("/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_softpowers.pdf", p, width=7.5 ,height=6,device = cairo_pdf)

power_table <- GetPowerTable(data)
head(power_table)

# construct co-expression network:
data <- ConstructNetwork(
  data,
  soft_power = 9,
  overwrite_tom = TRUE,
  tom_name = 'CAD_MM' # name of the topoligical overlap matrix written to disk
)

pdf(file = "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_Dendrogram.pdf", width=10, height=8)

# 设置全局字体为 Arial (可能不适用于所有绘图元素)
#par(family = "Arial")

# 同时增大字体大小 (结合之前的建议)
par(cex = 2, cex.main = 3, cex.lab = 2, cex.axis = 2)

PlotDendrogram(data,  setLabels = "",)

dev.off()

modules <- data@misc$wgcna$wgcna_modules
table(modules$module)
write.csv(modules,"/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_modules.csv")

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


saveRDS(data,"/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data.rds")


data <- readRDS("/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data.rds")
hub_df <- GetHubGenes(data, n_hubs = 10)
write.csv(hub_df,"/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_genes.csv")

# 绘制树状图与热图
pdf(file = "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_Dendrogram1.pdf",width=10, height=8)
plotEigengeneNetworks(
  MEs, 
  setLabels = "",
  plotHeatmaps = FALSE, 
  marDendro = c(0, 4, 2, 0)
)

dev.off()
MEs <- MEs[, colnames(MEs) != "grey"]
# 设置PDF输出
pdf(file = "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_Heatmaps.pdf", 
    width=12.5, height=10,)  # 设置字体为Arial

# 设置绘图参数
# 设置全局图形参数，将字体设置为 Arial，标签颜色设置为黑色
par(mar = c(10, 10, 10, 2),  col.lab = "black", col.axis = "black")

# 绘制热图并调整字体大小
plotEigengeneNetworks(
  MEs,
  setLabels = "",
  plotDendrograms = FALSE,
  xLabelsAngle = 90,
  marHeatmap = c(15, 15, 4, 2),
  cex.main = 3,    # 主标题字体大小
  cex.lab = 2.5,     # 坐标轴标签字体大小 (可能会影响轴的文本，不一定是行/列标签)
  cex.axis = 2,    # 坐标轴刻度字体大小
  cexRow = 2,      # 行标签字体大小 
  cexCol = 2 # 列标签字体大小 
)


# 关闭图形设备
dev.off()


desired_order <- c("green",'turquoise',"brown",'purple','greenyellow',
       'pink','red',"black","yellow","magenta" ,"blue")
mods_ordered <- mods[match(desired_order, mods)]

p <- DotPlot(
  data, 
  features = mods_ordered, 
  group.by = 'metacell_grouping',
  scale = FALSE,
  dot.scale = 12,  # 增大点的大小
  cols = c('blue', 'red')
) + 
  labs(title = " ") + 
  coord_flip() + 
  theme_minimal(base_family = "Arial") +  
  scale_y_discrete(expand = expansion(mult = c(1, 1))) +  # 大幅减少y轴扩展空间
  theme(
    text = element_text(family = "Arial", size = 15),  # 全局字体设置
    axis.title = element_blank(),   
    axis.text.x = element_text(
      size = 20,  # 增大X轴标签字号
      angle = 45,      
      hjust = 1,
      color = "black"  ,
      vjust = 1  
    ),
    axis.text.y = element_text(
      size = 20,  # 增大Y轴标签字号
      face = "italic",  
      margin = margin(r = 5) # 增加右侧间距
    ),
    legend.text = element_text(size = 15,color = "black"  ),  # 图例文字大小
    legend.title = element_text(size = 18,color = "black"  )  # 图例标题大小
  )

# 输出图形
print(p)

# 保存高清PDF（确保字体嵌入）
ggplot2::ggsave(
  "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/Mono_module.pdf", 
  plot = p,
  width = 7.5,   # 增加画布宽度
  height = 7.5,   # 增加画布高度
  device = cairo_pdf,  # 使用cairo设备确保字体嵌入
  family = "Arial",    # 指定PDF字体
  limitsize = FALSE)


# 提取元数据
meta_data <- data@meta.data


# 检查分组类别
table(meta_data$Sample_Type)

cat("正在重绘 ModuleFeaturePlot (原默认顺序，4列布局 + XY轴反转)...\n")

# 步骤 A: 正常调用函数生成图表列表 (传入 'MEs' 避免报错)
plot_list <- ModuleFeaturePlot(
  data,
  features = 'MEs',
  order = TRUE,
  label = TRUE
)

# 设置 PDF（增加高度以完美容纳 4 列）
pdf(file = "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_ModuleFeaturePlot_4cols_flipped.pdf", 
    width = 24, height = 20)

# 步骤 B: 直接使用 patchwork 绘制 4列 布局
combined_feature_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(
    theme = theme(plot.margin = margin(20, 20, 20, 20))
  ) & 
  coord_flip() &  # 应用 XY 轴反转
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12)
  )

print(combined_feature_plot)
dev.off()


# =========================================================================
# 修改 2：S2.5_Module_Genes_kME_Plot (原默认顺序 + 4列 布局)
# =========================================================================
cat("正在重绘 kME Plot (原默认顺序，4列布局)...\n")

generate_kme_plots_4cols <- function(module_kme_df) {
  
  # 直接获取所有非 grey 的有效颜色，保持原有默认顺序
  module_colors <- unique(module_kme_df$color) %>%
    grep("grey", ., invert = TRUE, value = TRUE)
  
  kme_plots <- list()
  
  for (module in module_colors) {
    kme_column <- paste0("kME_", module)
    
    if (!kme_column %in% colnames(module_kme_df)) next
    
    # 数据预处理
    module_genes_kme <- module_kme_df %>%
      filter(color == module) %>%
      mutate(
        kME = .data[[kme_column]],
        rank = row_number(-.data[[kme_column]]) 
      ) %>%
      filter(kME > 0.25) %>% 
      arrange(desc(kME)) %>%
      dplyr::select(gene_name, kME, rank)
    
    if (nrow(module_genes_kme) == 0) next
    
    top_genes <- module_genes_kme %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      head(20) 
    
    # 创建主图
    main_plot <- ggplot(module_genes_kme, aes(x = rank, y = kME)) +
      geom_bar(stat = "identity", fill = module, width = 0.8) +
      scale_x_reverse() +
      labs(
        title = paste(module, "Module (kME > 0.25)"),
        subtitle = paste("Top", nrow(module_genes_kme), "Genes"),
        x = "Gene Rank",
        y = "kME"
      ) +
      theme_minimal(base_size = 15, base_family = "Arial") +
      theme(
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", color = "black", size = 20, hjust = 0.5),
        plot.subtitle = element_text(color = "gray40", hjust = 0.5),
        axis.title = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial"),
        plot.margin = margin(1, -10, 1, 1)
      )
    
    # 创建TOP基因列表表格图
    table_plot <- ggplot(top_genes, aes(x = 1, y = rev(seq_along(gene_name)), label = gene_name)) +
      geom_text(size = 4, family = "Arial", hjust = 0) +
      scale_y_continuous(limits = c(0.5, nrow(top_genes) + 0.5)) +
      theme_void() +
      theme(
        plot.margin = margin(1, 1, 1, -10),
        text = element_text(family = "Arial")
      ) +
      labs(title = "") +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    # 将主图和表格图组合
    combined_plot <- main_plot + plot_spacer() + table_plot + 
      plot_layout(widths = c(4, -0.75, 1.5)) + 
      theme(plot.margin = margin(l = 10, r = 10, unit = "pt"))
    
    kme_plots[[module]] <- combined_plot
  }
  
  # 统一 4 列布局
  if (length(kme_plots) > 0) {
    wrap_plots(kme_plots, ncol = 4) + 
      plot_annotation(
        title = "",
        theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
      )
  } else {
    message("No modules found with genes having kME > 0.25.")
    return(NULL)
  }
}

# 调用函数生成 4列布局的 kME 图
combined_kme_plot <- generate_kme_plots_4cols(data@misc$wgcna$wgcna_modules)

# 保存 4x4 kME Plot
if (!is.null(combined_kme_plot)) {
  ggplot2::ggsave(
    filename = "/public3/DSC/single_cell/Result/figer_new/S2.5_Module_Genes_kME_Plot_4cols.pdf",
    plot = combined_kme_plot,
    device = cairo_pdf, 
    width = 25,     
    height = 22     
  )
}

cat("全套分析出图已完成！\n")

# 提取元数据
meta_data <- data@meta.data

# 完整定义 16 个模块
modules_to_test <- c("turquoise", "yellow", "green", "black", "blue", 
                     "magenta", "brown", "lightcyan", "greenyellow", 
                     "salmon", "cyan", "midnightblue", "pink", 
                     "purple", "red", "tan")

results <- data.frame(
  Module = character(),
  mean_AC = numeric(),
  mean_PA = numeric(),
  stringsAsFactors = FALSE
)

# 循环计算每个模块的均值
for (module in modules_to_test) {
  if(module %in% colnames(meta_data)){
    expr_AC <- na.omit(meta_data[meta_data$Sample_Type == "Atherosclerotic Core", module])
    expr_PA <- na.omit(meta_data[meta_data$Sample_Type == "Proximal Adjacent", module])
    
    results <- rbind(results, data.frame(
      Module = module,
      mean_AC = mean(expr_AC),
      mean_PA = mean(expr_PA)
    ))
  }
}

# 计算差值：PA - AC (正值代表 PA 较高，负值代表 AC 较高)
results$diff_PA_AC <- results$mean_PA - results$mean_AC

# 按差值降序排列（确保纵向模块按 Proximal Adjacent 最高到 Atherosclerotic Core 最高排列）
ordered_results <- results[order(results$diff_PA_AC, decreasing = FALSE), ]
mods_ordered <- ordered_results$Module

# 强制设置横轴因子顺序，使横向按 "Proximal Adjacent", "Atherosclerotic Core" 排列
data@meta.data$metacell_grouping <- factor(
  data@meta.data$metacell_grouping, 
  levels = c("Proximal Adjacent", "Atherosclerotic Core")
)

# 使用排序后的 16 个模块重新绘图
p <- DotPlot(
  data, 
  features = mods_ordered, 
  group.by = 'metacell_grouping',
  scale = FALSE,
  dot.scale = 12, 
  cols = c('blue', 'red')
) + 
  labs(title = " ") + 
  coord_flip() + 
  theme_minimal(base_family = "Arial") +  
  scale_y_discrete(expand = expansion(mult = c(0.5, 0.5))) + 
  theme(
    text = element_text(family = "Arial", size = 15), 
    axis.title = element_blank(),   
    axis.text.x = element_text(
      size = 20, 
      angle = 45,      
      hjust = 1,
      color = "black",
      vjust = 1  
    ),
    axis.text.y = element_text(
      size = 20, 
      face = "italic",  
      margin = margin(r = 5)
    ),
    legend.text = element_text(size = 15, color = "black"), 
    legend.title = element_text(size = 18, color = "black")
  )

# 输出图形
print(p)

ggplot2::ggsave(
  "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/Mono_module.pdf", 
  plot = p,
  width = 7.5,   # 增加画布宽度
  height = 7.5,   # 增加画布高度
  device = cairo_pdf,  # 使用cairo设备确保字体嵌入
  family = "Arial",    # 指定PDF字体
  limitsize = FALSE)



# 1. 完整定义 16 个待展示的模块
modules_to_plot <- c("turquoise", "yellow", "green", "black", "blue", 
                     "magenta", "brown", "lightcyan", "greenyellow", 
                     "salmon", "cyan", "midnightblue", "pink", 
                     "purple", "red", "tan")

# 2. 提取并计算每个细胞类型 (Celltype_raw1) 中模块的平均表达值
# 使用 dplyr:: 强行指定命名空间，完美避开 plyr 等包的冲突
avg_exp <- data@meta.data %>%
  dplyr::group_by(Celltype_raw1) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(modules_to_plot), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame() %>%  # 转为基础数据框
  tibble::column_to_rownames("Celltype_raw1")

# 3. 对模块（Y轴）进行数据标准化和数学聚类，自动寻找最佳上下顺序
avg_exp_t <- t(avg_exp)
# 定义 Z-score 标准化函数
cal_z_score <- function(x){ (x - mean(x)) / sd(x) }
avg_exp_scaled <- t(apply(avg_exp_t, 1, cal_z_score))

# 执行层次聚类并提取最佳模块顺序
hc_rows <- hclust(dist(avg_exp_scaled))
optimal_mod_order <- hc_rows$labels[hc_rows$order]

# 4. 手动定义完美的对角线细胞类型顺序（X轴）
ideal_cell_order <- c(
  "Foam cells1", 
  "Foam cells2", 
  "LAM", 
  "LYVE1+ TRM", 
  "Transitional Mac", 
  "CX3CR1+ TRM", 
  "Inflammatory Mono", 
  "Classical Mono", 
  "Non-classical Mono", 
  "ISG+ Mono"
)

# 5. 强制修改 Seurat 对象元数据中的细胞类型因子水平
data@meta.data$Celltype_raw1 <- factor(
  data@meta.data$Celltype_raw1,
  levels = ideal_cell_order
)

# 6. 使用最优的 Y 轴模块顺序和完美的 X 轴细胞顺序绘制 DotPlot
p_diagonal_perfect <- DotPlot(
  data, 
  features = optimal_mod_order, 
  group.by = "Celltype_raw1",   
  scale = FALSE,
  dot.scale = 14, # 放大点，使图像饱满
  cols = c('blue', 'red')
) + 
  #labs(title = "Clustered Module Expression") + 
  coord_flip() + 
  theme_minimal(base_family = "Arial") +  
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) + # 极致压缩 Y 轴留白
  theme(
    text = element_text(family = "Arial", size = 15), 
    axis.title = element_blank(),   
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black", vjust = 1),
    axis.text.y = element_text(size = 18, face = "italic", margin = margin(r = 5)),
    legend.text = element_text(size = 15, color = "black"), 
    legend.title = element_text(size = 18, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5) # 添加网格线增强视觉紧凑感
  )

# 7. 在 RStudio 中预览图形
print(p_diagonal_perfect)

# 8. 保存高清 PDF 结果文件（开启 limitsize = FALSE 防止超幅报错）
ggplot2::ggsave(
  filename = "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/Module_CelltypeRaw1_Diagonal_Perfect.pdf", 
  plot = p_diagonal_perfect, 
  width = 10, 
  height = 8, 
  device = cairo_pdf,
  family = "Arial",
  limitsize = FALSE
)


#### 每个module的GO富集 ####
file_paths <- '/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_modules.csv'
modules_list <- lapply(file_paths, read.csv)

# 给列表元素命名（可选）
names(modules_list) <- 'MM'
output_root <- "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/GO/"

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
      
      ggplot2::ggsave(
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

modules <- list(
  # TURQUOISE - 保留自噬、液泡、糖/脂代谢、囊泡运输和DNA修复
  turquoise = c(1, 2, 6, 7, 8, 11, 12, 14, 25, 26),
  # YELLOW - 保留细胞凋亡、细菌响应、T细胞激活/分化、趋化与NF-kB
  yellow = c(1, 3, 4, 9, 10, 17, 19, 27, 28, 29),
  # GREEN - 保留肌动蛋白、PI3K/AKT、补体激活、内吞、神经发育与细胞增殖
  green = c(1, 2, 6, 8, 9, 10, 11, 16, 20, 23),
  # BLACK - 保留膜定位、细胞质分裂、线粒体组装、TLR4、泛素化及ATP
  black = c(1, 3, 4, 5, 6, 8, 15, 20, 25, 30),
  # BLUE - 移除大量“病毒响应”同义词，保留病毒防御、干扰素信号、受体信号
  blue = c(1, 10, 14, 16, 20, 22, 24, 25, 28, 30),
  # MAGENTA - 保留蛋白折叠/应激、骨髓细胞/红细胞分化、包涵体与RNA剪接
  magenta = c(1, 3, 5, 9, 10, 11, 12, 13, 22, 23),
  # BROWN - 保留受体信号、mRNA代谢/剪接、表观遗传、核转运、染色单体及微管
  brown = c(1, 2, 4, 5, 6, 10, 11, 12, 13, 19),
  # LIGHTCYAN - 重点挑选糖酵解、嘌呤/吡啶代谢、丙酮酸及ATP代谢
  lightcyan = c(1, 5, 6, 7, 14, 22, 26, 27, 28, 29),
  # GREENYELLOW - 保留细胞外基质、金属离子解毒、TGF-beta、淀粉样纤维与成纤维细胞
  greenyellow = c(1, 4, 5, 7, 15, 17, 19, 22, 24, 29),
  # SALMON - 保留脂质/碳水化合物代谢、脂蛋白清除、MHC II、巨噬细胞激活与呼吸爆发
  salmon = c(1, 2, 5, 7, 9, 11, 13, 19, 23, 25),
  # CYAN - 保留T细胞细胞毒性、溶酶体酸化、胞饮作用、神经节苷脂代谢与EGFR负调控
  cyan = c(1, 2, 8, 11, 16, 17, 19, 23, 25, 28),
  # MIDNIGHTBLUE - 保留伤口愈合、单核细胞趋化、内皮细胞凋亡、血管生成与基质粘附
  midnightblue= c(1, 3, 4, 6, 12, 13, 25, 27, 28, 30),
  # PINK - 严格控制核糖体冗余，保留翻译起始与调控
  pink = c(1, 2, 4, 12, 16, 26, 35, 37, 40, 46),
  # PURPLE - ，保留核心能量代谢，深挖膜电位与剪接体组装
  purple = c(1, 6, 7, 11, 15, 27, 38, 40, 49, 50),
  # RED - 融合炎症小体(IL-1β)、趋化性与特定免疫细胞(NK/T细胞)
  red = c(1, 2, 3, 10, 12, 15, 24, 36, 37, 47),
  # TAN - 剔除冗余呈递词，融合ADCC效应、DC/巨噬细胞激活与免疫记忆
  tan = c(1, 3, 12, 14, 21, 24, 33, 36, 40, 44)
)
##
# Create output directory if it doesn't exist
output_dir <- "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/GO/MM/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each module
for (module in names(modules)) {
  # Read the data
  file_path <- paste0(output_dir, module, "_GO.csv")
  data <- read.csv(file_path)
  
  # Process the data
  data$Score <- -log10(data$p.adjust)
  data <- data[modules[[module]], ]  # Select specific rows
  data <- data %>% arrange(desc(Score))
  data$Description <- factor(data$Description, levels = rev(data$Description))
  
  # Create the plot
  bar_plot <- ggplot(data, aes(x = Description, y = Score, fill = Score)) +
    geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.3) +
    scale_fill_gradient(low = "#6BAED6", high = "#FDAE6B") +
    labs(
      title = paste("Top Enriched GO Terms (", module, " Module)", sep = ""),
      x = "GO Biological Process",
      y = "-log10(Adjusted p-value)"
    ) +
    coord_flip(clip = "off") +
    theme_classic(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 15, color = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    geom_text(
      aes(y = 0, label = Description),
      hjust = 0,
      size = 5,
      color = "black",
      nudge_x = 0.1
    )
  
  # Save the plot
  output_file <- paste0(output_dir, module, "_BAR10.pdf")
  ggplot2::ggsave(output_file, bar_plot, width = 8, height = 6)
  
  # Print progress
  message("Created plot for ", module, " module")
}

#####macSpectrum####
# ------------------------------------------------------------------------------
# 0. 环境准备与包加载
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)
library(scales)
library(ggpubr)
library(gtable)
library(AUCell)
library(hdWGCNA)
library(monocle3)
library(tidydr)
library(harmony)
library(patchwork)

# 全局绘图主题设置
theme_custom <- theme_bw(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

out_dir <- "/public3/DSC/single_cell/Result/figer_new/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 模块 1：AUCell 极化指数(MPI)与成熟度指数(AMDI)打分 
# ==============================================================================
# 初始数据加载 (请修改为您的实际初始无打分对象路径)
# data <- readRDS("/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data.rds")

# 1. 读取并提取 Classical Mono 特征基因集 (Top 100)
markers_df <- read.csv("/public3/DSC/single_cell/Result/figer_new/ST2_Ma_Mo_Marker.csv")
top100_mono <- markers_df %>%
  filter(cluster == "Classical Mono") %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%  # 按显著性和 log2FC 排序
  slice_head(n = 100) %>%                   # 提取前 100 行
  pull(gene)                                # 取出基因名向量

# 2. 定义 M1 极化特征基因集 (145 个特征)(合并 Proteomics, Core set, Azizi 后去重，共 145 个特征)
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

# 3. 定义 M2 极化特征基因集 (165 个特征)
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

# 4. 构建 AUCell 输入 List
geneSets <- list(
  M1_Score = m1_features, 
  M2_Score = m2_features,
  Mono_Immaturity_Score = top100_mono
)

# 5. 执行 AUCell 打分
# 取消注释以运行完整的打分流程
# exprMatrix <- GetAssayData(data, layer = "counts")
# message("正在构建细胞排名树，请稍候...")
# cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 16, plotStats = FALSE)
# message("正在计算 M1, M2 以及 单核/AMDI 的 AUC 分数...")
# cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
# auc_scores <- as.data.frame(t(getAUC(cells_AUC)))

# 6. 核心计算逻辑：AMDI_Index 取相反数, Polarization 取差值
# auc_scores <- auc_scores %>%
#   mutate(
#     AMDI_Index = -Mono_Immaturity_Score,  
#     Polarization_Index = M1_Score - M2_Score 
#   )
# data <- AddMetaData(data, metadata = auc_scores)

# 7. 保存带有最终得分的 Seurat 对象
# saveRDS(data, "/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_scored.rds")


# ==============================================================================
# 模块 2：读取打分数据 & 鲁棒性降采样 & 堆叠柱状复合图表
# ==============================================================================
message("\n--- 开始加载已打分的 Seurat 对象 ---")
data <- readRDS("/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_scored.rds")
meta <- data@meta.data

# 设置输出目录
out_dir <- "/public3/DSC/single_cell/Result/figer_new/"
dir.create(out_dir, showWarnings = FALSE)


# --- 2. 绘制 Macrophage_Polarization_Direction (极化方向柱状图 p2) ---
plot_data_polar <- data@meta.data %>%
  dplyr::group_by(Celltype_raw1) %>%
  dplyr::summarise(Mean_Polar = mean(Polarization_Index, na.rm = TRUE)) %>%
  dplyr::arrange(desc(Mean_Polar)) %>%
  dplyr::mutate(Celltype_raw1 = factor(Celltype_raw1, levels = Celltype_raw1))

p2 <- ggplot(plot_data_polar, aes(x = Celltype_raw1, y = Mean_Polar, fill = Mean_Polar > 0)) +
  geom_bar(stat = "identity", color = "white", width = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#BC3C29FF", "FALSE" = "#0072B5FF"), 
                    labels = c("TRUE" = "M1 Dominant", "FALSE" = "M2 Dominant")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  theme_classic() +
  labs(#title = "Macrophage Polarization Direction (AUCell)",
       #subtitle = "Net Score: M1_AUC - M2_AUC",
       x = "Cell Subtypes",
       y = "Macrophage Polarization Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        legend.position = "top",
        legend.title = element_blank())

ggsave(paste0(out_dir, "Macrophage_Polarization_Direction.pdf"), plot = p2, device = cairo_pdf, width = 12, height = 6)


# --- 3. 绘制 Macrophage_Maturation_Index (成熟度阶梯图 p3) ---
plot_amdi <- data@meta.data %>%
  dplyr::group_by(Celltype_raw1) %>%
  dplyr::summarise(Mean_AMDI = mean(AMDI_Index, na.rm = TRUE)) %>%
  dplyr::arrange(Mean_AMDI) %>%
  dplyr::mutate(Celltype_raw1 = factor(Celltype_raw1, levels = Celltype_raw1))

p3 <- ggplot(plot_amdi, aes(x = Celltype_raw1, y = Mean_AMDI, fill = Mean_AMDI)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_gradient(low = "#E5F5E0", high = "#31A354") +
  theme_classic() +
  labs(#title = "Macrophage Maturation Index (AUCell)", 
       #subtitle = "MMI Score: 1 - Monocyte_Immaturity_AUC",
       x = "Cell Subtypes", 
       y = "Macrophage Maturation Index (MMI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        legend.position = "none")

ggsave(paste0(out_dir, "Macrophage_Maturation_Index.pdf"), plot = p3, device = cairo_pdf, width = 12, height = 6)


# --- 4. 绘制 FeaturePlot (空间映射图 p_mmi / p_mpi) ---
# 修复：加入 max.cutoff = 'q98' 控制极值饱和，明确 name 标签
p_mpi <- FeaturePlot(data, features = "Polarization_Index", pt.size = 0.8, 
                     reduction = "umap", max.cutoff = 'q98') +
  scale_colour_gradient2(low = "#0072B5FF", mid = "lightgrey", high = "#BC3C29FF", 
                         midpoint = 0, name = "MPI") +
  ggtitle("Macrophage Polarization Index (MPI)") + 
  coord_flip()

ggsave(paste0(out_dir, "FeaturePlot_MPI_coord_flipped.pdf"), plot = p_mpi, device = cairo_pdf, width = 10, height = 8)

p_mmi <- FeaturePlot(data, features = "AMDI_Index", pt.size = 0.8, 
                     reduction = "umap", max.cutoff = 'q98') +
  scale_colour_gradientn(
    colours = c("#0072B5FF", "lightgrey", "#FF0000"), 
    values = scales::rescale(c(0.4, 0.8, 1)), 
    name = "MMI"
  ) +
  ggtitle("Macrophage Maturation Index (MMI)") + 
  coord_flip()

p_mmi
ggsave(paste0(out_dir, "FeaturePlot_MMI_coord_flipped.pdf"), plot = p_mmi, device = cairo_pdf, width = 10, height = 8)


# --- 5. 绘制 Macrophage_AUCell_M1_vs_M2_Scores (背靠背柱状图 p1_reversed) ---
plot_data_m1_m2 <- data@meta.data %>%
  dplyr::group_by(Celltype_raw1) %>%
  dplyr::summarise(Mean_M1 = mean(M1_Score, na.rm = TRUE),
                   Mean_M2 = mean(M2_Score, na.rm = TRUE))

plot_long <- plot_data_m1_m2 %>% 
  tidyr::pivot_longer(cols = c(Mean_M1, Mean_M2), names_to = "Metric", values_to = "Score") %>%
  dplyr::mutate(Score = ifelse(Metric == "Mean_M2", -Score, Score))

p1_reversed <- ggplot(plot_long, aes(x = Celltype_raw1, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("Mean_M1" = "#E64B35", "Mean_M2" = "#4DBBD5"),
                    labels = c("Mean_M1" = "M1 AUC Score", "Mean_M2" = "M2 AUC Score")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_y_continuous(labels = abs) +
  theme_classic() +
  labs( x = "Cell Subtypes", y = "AUC Score") + #title = "AUCell M1 vs M2 Scoring",
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        legend.position = "top", 
        legend.title = element_blank())

ggsave(paste0(out_dir, "Macrophage_AUCell_M1_vs_M2_Scores.pdf"), plot = p1_reversed, device = cairo_pdf, width = 12, height = 6)

# --- 2.1 降采样前：提取原始比例 ---
table_counts_ori <- table(data@meta.data$Sample_Type, data@meta.data$Celltype_raw1)
df_ori <- as.data.frame(table_counts_ori)
colnames(df_ori) <- c("Status", "CellType", "Count")
df_ori <- df_ori %>%
  dplyr::group_by(Status) %>%
  dplyr::mutate(Proportion_Original = Count / sum(Count)) %>%
  dplyr::ungroup()

# --- 2.2 执行高鲁棒性降采样 (针对 Control AscAorta) ---
set.seed(123)
# 修复: 从 AC_PA 列提取主动脉细胞，使用 which 防 NA 报错
cells_asc <- rownames(data@meta.data)[which(data@meta.data$AC_PA == "Control AscAorta")]
cat("\n成功抓取到的 AscAorta 细胞数量: ", length(cells_asc), "\n")

cells_asc_downsampled <- sample(cells_asc, size = 1000)
cells_other <- setdiff(rownames(meta), cells_asc)
cells_keep <- c(cells_asc_downsampled, cells_other)

meta_ds <- meta[cells_keep, ]
data_ds <- subset(data, cells = cells_keep) # 生成降采样后的对象，供后续 DotPlot 使用

# --- 2.3 降采样后：计算新的绘图数据 ---
table_counts_ds <- table(meta_ds$Sample_Type, meta_ds$Celltype_raw1)
table_counts_ds <- table_counts_ds[, colSums(table_counts_ds) > 0]

df <- as.data.frame(table_counts_ds)
colnames(df) <- c("Status", "CellType", "Count")
df$Status <- factor(df$Status, levels = c("Proximal Adjacent", "Atherosclerotic Core"))

df <- df %>%
  dplyr::group_by(Status) %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::ungroup()

# --- 2.4 生成降采样前后细胞比例差异对比表 ---
df_compare <- df %>%
  dplyr::rename(Proportion_Downsampled = Proportion) %>%
  dplyr::left_join(df_ori[, c("Status", "CellType", "Proportion_Original")], by = c("Status", "CellType")) %>%
  dplyr::filter(Status == "Proximal Adjacent") %>% 
  dplyr::mutate(
    Original_Pct = paste0(round(Proportion_Original * 100, 2), "%"),
    Downsampled_Pct = paste0(round(Proportion_Downsampled * 100, 2), "%"),
    Difference = paste0(ifelse(Proportion_Downsampled > Proportion_Original, "+", ""),
                        round((Proportion_Downsampled - Proportion_Original) * 100, 2), "%")
  ) %>%
  dplyr::arrange(desc(abs(Proportion_Downsampled - Proportion_Original))) %>%
  dplyr::select(CellType, Original_Pct, Downsampled_Pct, Difference)

cat("\n[分析结果] 降采样前后 Proximal Adjacent 组细胞比例变化对比：\n")
print(as.data.frame(df_compare))

# 统一颜色字典 (全集)
cell_colors <- c(
  "Classical Mono"     = "#8c564b", 
  "Inflammatory Mono"  = "#b15928", 
  "ISG+ Mono"          = "#bcbd22", 
  "Non-classical Mono" = "#c7c7c7", 
  "Foam cells1"        = "#d62728", 
  "Foam cells2"        = "#ff7f0e", 
  "LAM"                = "#9467bd", 
  "Transitional Mac"   = "#2ca02c", 
  "CX3CR1+ TRM"        = "#e377c2", 
  "LYVE1+ TRM"         = "#f7b6d2",
  "TrMs"               = "#e377c2", 
  "CM"                 = "#2ca02c", 
  "Macrophage"         = "#9467bd", 
  "Monocyte"           = "#8c564b", 
  "cDC1"               = "#1f77b4"
)

# 确保数据因子的 level 顺序正确
df$CellType <- factor(df$CellType, levels = names(cell_colors))

# 统一表格样式
my_table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontfamily = "Arial", fontsize = 10, hjust = 0.5, x = 0.5),
    bg_params = list(fill = c("white", "#f7f7f7"))
  ),
  colhead = list(
    bg_params = list(fill = "#404040"),
    fg_params = list(col = "white", fontface = "bold", fontfamily = "Arial", fontsize = 10, hjust = 0.5, x = 0.5)
  )
)

# ================= 数据汇总处理 (剔除全 0 列) =================
table_data <- df %>%
  dplyr::group_by(Status, CellType) %>%
  dplyr::summarise(Proportion = mean(Proportion), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = CellType, values_from = Proportion, values_fill = 0) %>%
  dplyr::select(Status, dplyr::any_of(names(cell_colors))) %>% 
  # 【剔除逻辑】：保留字符列以及总和大于 0 的细胞列
  dplyr::select(where(~ is.character(.) || is.factor(.) || sum(., na.rm = TRUE) > 0)) %>%
  # 转为百分比字符串
  dplyr::mutate(dplyr::across(where(is.numeric), ~ paste0(round(. * 100, 1), "%")))

# ================= 1. 绘制堆叠柱状图 =================
p_stack <- ggplot(df, aes(x = Status, y = Proportion, fill = CellType)) +
  # 【反转逻辑】：position_fill(reverse = TRUE) 实现堆叠反转，使之与参考图一致
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.6) + 
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors, breaks = names(cell_colors)) +
  labs(title = " ", x = "Group", y = "Proportion") +
  theme_minimal(base_family = "Arial", base_size = 12) + 
  theme(
    plot.margin = margin(t=20, r=20, b=10, l=20, unit="pt"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12, color = "black"), 
    legend.position = "right", 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()  
  )

# ================= 2. 动态拆分表格 (防报错版) =================
n_cols <- ncol(table_data)

if (n_cols > 6) {
  split_idx <- ceiling((n_cols - 1) / 2) + 1
  
  part1 <- table_data[, 1:split_idx, drop = FALSE]
  part2 <- table_data[, c(1, (split_idx+1):n_cols), drop = FALSE]
  
  # 【安全补齐逻辑】：如果 part2 列数少于 part1，用隐形空格列补齐
  diff_cols <- ncol(part1) - ncol(part2)
  if (diff_cols > 0) {
    for (i in 1:diff_cols) {
      empty_col_name <- paste0(rep(" ", i), collapse = "") # 制造不同长度的空格防重名
      part2[[empty_col_name]] <- ""
    }
  }
  
  grob_1 <- tableGrob(part1, rows = NULL, theme = my_table_theme)
  grob_2 <- tableGrob(part2, rows = NULL, theme = my_table_theme)
  # 现在由于列数绝对相等，rbind 再也不会报错了
  table_grob <- gtable_rbind(grob_1, grob_2)
  
} else {
  table_grob <- tableGrob(table_data, rows = NULL, theme = my_table_theme)
}

# ================= 3. 组合排版 =================
final_stack_plot <- grid.arrange(p_stack, table_grob, nrow = 2, heights = c(0.7, 0.3))
ggplot2::ggsave(paste0(out_dir, "/Cell_Proportion_StackedBar_Table.pdf"), plot = final_stack_plot, device = cairo_pdf, width = 12, height = 10)


# ==============================================================================
# 模块 3：模块表达 DotPlot (降采样前后相对百分比修正验证)
# ==============================================================================
message("\n--- 开始计算降采样前后模块表达变化 (DotPlot) ---")
# 注意：假设您的全集模块名称保存在变量 mods 中
mods <- unique(GetModules(data)$color) 
mods <- setdiff(mods, "grey")

dp_data_ori <- DotPlot(data, features = mods, group.by = 'metacell_grouping', scale = FALSE)$data
dp_data_ds <- DotPlot(data_ds, features = mods, group.by = 'metacell_grouping', scale = FALSE)$data

PA_ori <- subset(dp_data_ori, id == "Proximal Adjacent")[, c("features.plot", "pct.exp")]
AC_ori <- subset(dp_data_ori, id == "Atherosclerotic Core")[, c("features.plot", "pct.exp")]
PA_ds <- subset(dp_data_ds, id == "Proximal Adjacent")[, c("features.plot", "pct.exp")]
AC_ds <- subset(dp_data_ds, id == "Atherosclerotic Core")[, c("features.plot", "pct.exp")]

core_vs_adj <- data.frame(
  Module = PA_ori$features.plot,
  Gap_Original = ifelse(PA_ori$pct.exp == 0, NA, (AC_ori$pct.exp - PA_ori$pct.exp) / PA_ori$pct.exp * 100),
  Gap_Downsampled = ifelse(PA_ds$pct.exp == 0, NA, (AC_ds$pct.exp - PA_ds$pct.exp) / PA_ds$pct.exp * 100)
)

core_vs_adj_formatted <- core_vs_adj %>%
  dplyr::mutate(Correction_Magnitude = Gap_Downsampled - Gap_Original) %>%
  dplyr::arrange(desc(abs(Correction_Magnitude))) %>%
  dplyr::mutate(
    Gap_Original_Pct = paste0(ifelse(Gap_Original > 0, "+", ""), round(Gap_Original, 2), "%"),
    Gap_Downsampled_Pct = paste0(ifelse(Gap_Downsampled > 0, "+", ""), round(Gap_Downsampled, 2), "%"),
    Correction_Magnitude_Pct = paste0(ifelse(Correction_Magnitude > 0, "+", ""), round(Correction_Magnitude, 2), "%")
  ) %>%
  dplyr::select(Module, Gap_Original_Pct, Gap_Downsampled_Pct, Correction_Magnitude_Pct)

cat("\n[核心对比] 降采样前后：斑块核心相对于对照的表达 [相对百分比变化]：\n")
print(core_vs_adj_formatted)

# --- 绘制真实分布的 DotPlot ---
proximal_ds <- subset(dp_data_ds, id == "Proximal Adjacent")
proximal_data_sorted <- proximal_ds[order(proximal_ds$pct.exp, decreasing = FALSE), ]
ordered_features_by_expr <- as.character(proximal_data_sorted$features.plot)

p_dot <- DotPlot(
  data_ds, features = ordered_features_by_expr, group.by = 'metacell_grouping',
  scale = FALSE, dot.scale = 12, cols = c('blue', 'red')
) + 
  labs(title = " ") + coord_flip() + theme_minimal(base_family = "Arial") +  
  scale_y_discrete(expand = expansion(mult = c(1, 1))) +
  theme(
    text = element_text(family = "Arial", size = 15), axis.title = element_blank(),   
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1, color = "black", vjust = 1),
    axis.text.y = element_text(size = 20, face = "italic", margin = margin(r = 5)),
    legend.text = element_text(size = 15, color = "black"), legend.title = element_text(size = 18, color = "black")
  )
ggplot2::ggsave(paste0(out_dir, "Mono_module_DotPlot_Corrected.pdf"), plot = p_dot, device = cairo_pdf, width = 7.5, height = 7.5)


# ==============================================================================
# 模块 4：极化与分化指数展示 (VlnPlot & Density)
# ==============================================================================
message("\n--- 开始绘制极化(MPI)与分化(AMDI)指数 ---")
Sample_Type_comparisons <- list(c("Atherosclerotic Core", "Proximal Adjacent"))

# --- VlnPlot: MPI (降采样后) ---
mpi_plot_after <- VlnPlot(data, features = "Polarization_Index", group.by = "Sample_Type", pt.size = 0, cols = c("#D95F02" ,"#1B9E77")) + 
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  labs(x = "Macrophage", y = "MPI") +
  stat_compare_means(comparisons = Sample_Type_comparisons, method = "wilcox.test", label = "p.signif", size = 5, bracket.size = 0.6) +
  ylim(NA, max(data_ds$Polarization_Index, na.rm = TRUE) * 1.5) + ggtitle("MPI: After Downsampling")

# --- VlnPlot: AMDI (降采样后) ---
amdi_plot_after <- VlnPlot(data, features = "AMDI_Index", group.by = "Sample_Type", pt.size = 0, cols = c("#D95F02" ,"#1B9E77")) + 
  theme_custom +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  labs(x = "Macrophage", y = "AMDI") +
  stat_compare_means(comparisons = Sample_Type_comparisons, method = "wilcox.test", label = "p.signif", size = 5, bracket.size = 0.6) +
  ylim(min(data_ds$AMDI_Index, na.rm = TRUE) * 1.1, max(data_ds$AMDI_Index, na.rm = TRUE) * 1.1) + ggtitle("AMDI: After Downsampling")

ggplot2::ggsave(paste0(out_dir, "MPI_AMDI_VlnPlot.pdf"), plot = (mpi_plot_after | amdi_plot_after), device = cairo_pdf, width = 10, height = 6)

# --- 2D Density: 状态分布 ---
p_density <- ggplot(
  data@meta.data %>% filter(Celltype_raw1 %in% c("LAM", "Foam cells1", "Foam cells2")), 
  aes(x = Polarization_Index, y = AMDI_Index) 
) + 
  geom_density_2d(aes(color = Celltype_raw1), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("LAM" = "#1B9E77", "Foam cells1" = "#D95F02", "Foam cells2" = "#7570B3")) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  theme_custom + 
  theme(legend.position = "right", legend.title = element_text(face = "bold"), legend.key = element_blank()) +
  labs(x = "Polarization Index (MPI)", y = "Maturation Index (AMDI)")
ggplot2::ggsave(paste0(out_dir, "Trajectory_Density.pdf"), plot = p_density, device = cairo_pdf, width = 8, height = 6)

# 严选的 10 个完美契合生物学故事的 Top Marker
target_genes <- c("APOE", "PLIN2", "TREM2", "IL1B", "CXCL8", "CTSS", "TIMP1", "C1QB", "HLA-DRA", "LDHA")# 批量绘制
plots_foam <- lapply(target_genes, function(gene) {
  VlnPlot(
    sub_data_foam,
    features = gene,
    group.by = "Celltype_raw1",
    pt.size = 0,
    cols = c("#1B9E77", "#D95F02")  
  ) +
    theme_custom +
    labs(y = "Expression Level") +
    stat_compare_means(
      comparisons = list(c("Foam cells1", "Foam cells2")),
      method = "wilcox.test",
      label = "p.signif",
      bracket.size = 0.6,
      tip.length = 0.02,
      size = 5,
      vjust = 0.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) 
})# 拼图展示 (2行5列)
combined_violin_plots <- wrap_plots(plots_foam, ncol = 5) 
print(combined_violin_plots)# 导出高清 PDF
ggplot2::ggsave(
  filename = "/public3/DSC/single_cell/Result/figer_new/Top10_Foam_Markers_VlnPlot.pdf",
  plot = combined_violin_plots,
  device = cairo_pdf,
  width = 20,  
  height = 8
)

#####monocle3
library(hdWGCNA)
library(monocle3)
library(tidydr)

data <- readRDS("/public3/DSC/single_cell/Result/figer_new/group_result.rds")
expression_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
cell_metadata <- data@meta.data

gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  module = "Default_Module",  
  row.names = rownames(expression_matrix)  
)

cell_ids <- rownames(cell_metadata)
sub_matrix <- expression_matrix

cds_MM <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

seurat_umap <- Embeddings(data, reduction = "umap")
rownames(seurat_umap) <- colnames(cds_MM)
reducedDims(cds_MM)[["UMAP"]] <- seurat_umap

##跳过 Monocle3 降维，直接导入 Seurat 的 UMAP 与 PCA 坐标

if ("integrated.rpca" %in% names(data@reductions)) {
  reducedDims(cds_MM)[["PCA"]] <- Embeddings(data, reduction = "integrated.rpca")
} else if ("rpca" %in% names(data@reductions)) {
  reducedDims(cds_MM)[["PCA"]] <- Embeddings(data, reduction = "rpca")
} else {
  reducedDims(cds_MM)[["PCA"]] <- Embeddings(data, reduction = "pca")
}

# 1. 聚类 (直接基于已有的 UMAP)
cds_MM <- cluster_cells(cds_MM, resolution = 5e-05)

# 2. 拟时序建树
cds_MM <- learn_graph(
  cds_MM,
  close_loop = TRUE,#FALSE,
  learn_graph_control = list(
    minimal_branch_len = 5, 
    prune_graph = TRUE
  )
)

p <- plot_cells(
  cds_MM,
  color_cells_by = "Celltype_raw1",
  cell_size = 1.0,
  label_groups_by_cluster = FALSE,
  group_label_size = 8,
  show_trajectory_graph = TRUE
) + 
  # 【核心修改】：交换 X 轴和 Y 轴
  coord_flip() + 
  
  # 注意：使用了 coord_flip 后，scale_y_reverse 实际上反转的是视觉上的水平轴，
  # scale_x_reverse 反转的是视觉上的垂直轴。如果方向不对，可以尝试删掉这两行或只留一行。
  #scale_y_reverse() +
  #scale_x_reverse() +
  
  facet_wrap(~Sample_Type, nrow = 1, scales = "free_y") +
  theme_dr() +
  theme(
    strip.text = element_text(size = 20),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    legend.key.height = unit(1.25, "cm"),
    text = element_text(family = "Arial"),
    # 由于坐标轴交换，原来的 y 轴现在变成了底部的 x 轴。
    # 如果你想隐藏右侧（或顶部/底部）的刻度，可能需要改为 axis.text.x.top / axis.text.x.bottom 等
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank()
  ) +
  annotate("label",  
           x = -5.5, y = 2.5,  # 提示：翻转后，这里的 x=11.5 将对应视觉上的垂直方向，y=4.5 对应水平方向
           label = "Cell fate 2",
           hjust = -0.5, vjust = 2,
           size = 6, color = "black",
           fill = "grey90",  
           label.padding = unit(0.15, "lines"),  
           label.r = unit(0.05, "lines")) +  
  annotate("label",
           x = -8, y = -3,    # 同理，这里的坐标可能需要根据翻转后的实际图形重新调整
           label = "Cell fate 1",
           hjust = 1.2, vjust = -1,
           size = 6, color = "black",
           fill = "grey90",
           label.padding = unit(0.15, "lines"),
           label.r = unit(0.05, "lines"))

print(p)


ggplot2::ggsave(p, 
                filename = '/public3/DSC/single_cell/Result/figer_new/F2.5_monocle_celltype_cluster.pdf',
                width = 14, 
                height = 7, 
                device = cairo_pdf)



######拆分partition####
#######################
cds_MM_foam <- order_cells(cds_MM)

plot_PT <- plot_cells(
  cds_MM_foam,
  color_cells_by = "pseudotime",    # 按伪时间着色
  label_cell_groups = FALSE, # 关闭细胞群标签
  label_leaves = TRUE,       # 开启轨迹叶节点标签[1](@ref)
  label_branch_points = TRUE,      # 开启轨迹分支点标签[1](@ref)
  trajectory_graph_color = "black", 
  # trajectory_graph_label_size = 5,  # 增大轨迹标签字体[1](@ref)
  cell_size = 1,
  alpha = 0.8
)   + 
  # 【核心修改】：交换 X 轴和 Y 轴
  coord_flip() + 
  
  # 注意：使用了 coord_flip 后，scale_y_reverse 实际上反转的是视觉上的水平轴，
  # scale_x_reverse 反转的是视觉上的垂直轴。如果方向不对，可以尝试删掉这两行或只留一行。
  #scale_y_reverse() +
  facet_wrap(~Sample_Type, nrow = 1,scales = "free_y") +   # 分面展示AC和PA组
  scale_color_gradientn(
    colours = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red'), # 自定义伪时间色阶
    name = "Pseudotime",
    guide = guide_colorbar(barwidth = 1.5, title.position = "top")   # 调整图例位置和样式[2](@ref)
  ) + 
  theme_dr() + 
  theme(
    strip.text = element_text(size = 20),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(
      size = 15,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    text = element_text(family = "Arial")
    
  ) +
  annotate("label",  
           x = -5.5, y = 2.5,  # 提示：翻转后，这里的 x=11.5 将对应视觉上的垂直方向，y=4.5 对应水平方向
           label = "Cell fate 2",
           hjust = -0.5, vjust = 2,
           size = 6, color = "black",
           fill = "grey90",  
           label.padding = unit(0.15, "lines"),  
           label.r = unit(0.05, "lines")) +  
  annotate("label",
           x = -8, y = -3,    # 同理，这里的坐标可能需要根据翻转后的实际图形重新调整
           label = "Cell fate 1",
           hjust = 1.2, vjust = -1,
           size = 6, color = "black",
           fill = "grey90",
           label.padding = unit(0.15, "lines"),
           label.r = unit(0.05, "lines"))



plot_PT
ggplot2::ggsave(plot_PT,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_MM_Pseudotim.png',width = 18,height = 8)
ggplot2::ggsave(plot_PT, filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_MM_Pseudotim.pdf', width = 14, height = 7, device = cairo_pdf)

p <- plot_cells(
  cds_MM,
  genes = "APOBEC3A",
  label_groups_by_cluster = FALSE,
  cell_size = 1.0,
  group_label_size = 4,
  show_trajectory_graph = TRUE
)  +
  coord_flip() + 
  
  # 注意：使用了 coord_flip 后，scale_y_reverse 实际上反转的是视觉上的水平轴，
  # scale_x_reverse 反转的是视觉上的垂直轴。如果方向不对，可以尝试删掉这两行或只留一行。
  #scale_y_reverse() +
  facet_wrap(~Sample_Type, ncol = 2, scales = "free") + 
  scale_color_gradientn(
    colors = c("gray90", "red2", "red3"),
    values = scales::rescale(c(0, 0.3, 1)),  # 低表达区0-30%映射到灰色到粉色，30-100%快速过渡到红色
    breaks = c(0, 0.5, 1)        # 强制色阶以0.5为中间断点
  )+
  labs(color = "APOBEC3A Expression") + 
  theme_dr() + 
  theme(text = element_text(family = "Arial"), 
        panel.grid = element_blank(),
        strip.text = element_text(size = 20),  
        strip.background = element_blank(),  
        plot.title = element_blank()) +
  annotate("label",  
           x = -5.5, y = 2.5,  # 提示：翻转后，这里的 x=11.5 将对应视觉上的垂直方向，y=4.5 对应水平方向
           label = "Cell fate 2",
           hjust = -0.5, vjust = 2,
           size = 6, color = "black",
           fill = "grey90",  
           label.padding = unit(0.15, "lines"),  
           label.r = unit(0.05, "lines")) +  
  annotate("label",
           x = -8, y = -3,    # 同理，这里的坐标可能需要根据翻转后的实际图形重新调整
           label = "Cell fate 1",
           hjust = 1.2, vjust = -1,
           size = 6, color = "black",
           fill = "grey90",
           label.padding = unit(0.15, "lines"),
           label.r = unit(0.05, "lines"))
p


###添加伪时序信息到 Seurat 中 #####
data<-readRDS("/public3/DSC/single_cell/Result/figer_new/hdWGCNA/Mo_Ma/data_scored.rds")
data$pseudotime <- pseudotime(cds_MM_foam)
summary(data$pseudotime)

monocle_Pseudotim  <- PlotModuleTrajectory(
  data,
  pseudotime_col = 'pseudotime'
) +
  theme(
    text = element_text(family = "Arial"))
monocle_Pseudotim <- monocle_Pseudotim + 
  theme(axis.title.x = element_text(margin = margin(t = 15)))
monocle_Pseudotim
ggplot2::ggsave(monocle_Pseudotim,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.6_MM_monocle_Pseudotim.png',width = 7,height = 4)
ggplot2::ggsave(monocle_Pseudotim,
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.6_MM_monocle_Pseudotim.pdf',
                width = 10,
                height = 5,
                device = cairo_pdf)


#trace_genes<- graph_test(cds_MM_foam, 
#  neighbor_graph = "principal_graph", 
#  cores = 8)
#sorted_res <- trace_genes %>% 
#  arrange(desc(morans_I))
#data<-readRDS("/public3/DSC/single_cell/Result/figer_new/monocle3/data_pseudotime.rds")

# 确保已加载 patchwork 用于拼图 (如未安装请运行 install.packages("patchwork"))
library(patchwork)
# 1. 弹出窗口选第 1 次
cat("请在弹出的窗口中选择第 1 个分支 (Fate 1)，完成后点击 Done...\n")
cds_subset_1 <- choose_graph_segments(cds_MM_foam)
selected_cells_1 <- colnames(cds_subset_1)

cat("请在弹出的窗口中选择第 2 个分支 (Fate 2)，完成后点击 Done...\n")
cds_subset_2 <- choose_graph_segments(cds_MM_foam)
selected_cells_2 <- colnames(cds_subset_2)
# ==========================================
# 1. 准备数据：为两张图分别创建标记列
# ==========================================
# 标记 Subset 1
colData(cds_MM_foam)$subset_1 <- ifelse(
  colnames(cds_MM_foam) %in% selected_cells_1, 
  "Fate 1", 
  "Unselected"
)

# 标记 Subset 2
colData(cds_MM_foam)$subset_2 <- ifelse(
  colnames(cds_MM_foam) %in% selected_cells_2, 
  "Fate 2", 
  "Unselected"
)

# ==========================================
# 2. 绘制左图：突出显示 cds_subset_1 (蓝色)
# ==========================================
p1 <- plot_cells(cds_MM_foam,
                 color_cells_by = "subset_1",  
                 label_cell_groups = FALSE,
                 cell_size = 1,        
                 trajectory_graph_segment_size = 0.5)  + 
  coord_flip() +
  # 移除 facet_wrap，改为独立的颜色映射
  scale_color_manual(
    name = "Selection",
    values = c("Unselected" = "gray80", "Fate 1" = "blue") # subset_1 为蓝色
  ) + 
  theme_dr() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5), # 增加居中标题
    legend.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    text = element_text(family = "Arial"),
    legend.position = "bottom" # 图例放到底部，避免左右拼图时挤占空间
  ) +
  ggtitle("Cell fate 1") + # 为左图添加标题
  # 只保留 Fate 1 的标签
  annotate("label",
           x = -8, y = -3,    
           label = "Cell fate 1",
           hjust = 1.2, vjust = -1,
           size = 6, color = "black",
           fill = "grey90",
           label.padding = unit(0.15, "lines"),
           label.r = unit(0.05, "lines"))

# ==========================================
# 3. 绘制右图：突出显示 cds_subset_2 (红色)
# ==========================================
p2 <- plot_cells(cds_MM_foam,
                 color_cells_by = "subset_2",  
                 label_cell_groups = FALSE,
                 cell_size = 1,        
                 trajectory_graph_segment_size = 0.5)  + 
  coord_flip() +
  scale_color_manual(
    name = "Selection",
    values = c("Unselected" = "gray80", "Fate 2" = "red") # subset_2 为红色
  ) + 
  theme_dr() + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    text = element_text(family = "Arial"),
    legend.position = "bottom" 
  ) +
  ggtitle("Cell fate 2") + 
  # 只保留 Fate 2 的标签
  annotate("label",  
           x = -5.5, y = 2.5,  
           label = "Cell fate 2",
           hjust = -0.5, vjust = 2,
           size = 6, color = "black",
           fill = "grey90",  
           label.padding = unit(0.15, "lines"),  
           label.r = unit(0.05, "lines"))

# ==========================================
# 4. 拼接并查看最终图像
# ==========================================
# 使用 patchwork 的 + 号语法将左右两张图拼在一起
p_selected_cells <- p1 + p2

p_selected_cells

ggplot2::ggsave(p_selected_cells,
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.7_MM_selected_cells.pdf', # 你可以根据需要修改文件名编号
                width = 14,  # 左右拼图，宽度设为 10 或 12 比较合适
                height = 7,
                device = cairo_pdf)
# 将对象保存到指定路径
# 1. 保存完整的整体对象 (会生成一个名为 cds_MM_foam 的文件夹)
save_monocle_objects(cds_MM_foam, 
                     directory_path = "/public3/DSC/single_cell/Result/figer_new/monocle3/cds_MM_foam")

# 2. 保存 Subset 2 (Disease Foam) (会生成一个名为 cdsMM_subset2_disease_foam 的文件夹)
save_monocle_objects(cds_subset_2, 
                     directory_path = "/public3/DSC/single_cell/Result/figer_new/monocle3/cdsMM_subset2_disease_foam")

# 3. 保存 Subset 1 (Healthy Macrophage) (会生成一个名为 cdsMM_subset1_healthy_macrophage 的文件夹)
save_monocle_objects(cds_subset_1, 
                     directory_path = "/public3/DSC/single_cell/Result/figer_new/monocle3/cdsMM_subset1_healthy_macrophage")

# 读取 Subset 2 (Disease Foam)
#cds_subset_2 <- load_monocle_objects(directory_path = "/public3/DSC/single_cell/Result/figer_new/monocle3/cdsMM_subset2_disease_foam")

# 读取 Subset 1 (Healthy Macrophage)
#cds_subset_1 <- load_monocle_objects(directory_path = "/public3/DSC/single_cell/Result/figer_new/monocle3/cdsMM_subset1_healthy_macrophage")

subset_list <- list(
  subset1 = cds_subset_1,
  subset2 = cds_subset_2
  
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
saveRDS(data,"/public3/DSC/single_cell/Result/figer_new/monocle3/data_pseudotime.rds")
#data<-readRDS("/public3/DSC/single_cell/Result/figer_new/monocle3/data_pseudotime.rds")
#####Sample_Type在不同细胞分支上随时间的分布变化#####

p <- ggplot( data = data@meta.data, aes(x = pseudotime, fill = Sample_Type)) +
  geom_density(alpha = 0.5) +
  # 1. 在这里输入你想要的标题内容
  labs(title = "Density of Pseudotime by Sample Type", 
       x = "Pseudotime", # 如果需要自定义X轴名字可保留，不需要则删去
       y = "Density") +  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),  
    axis.text = element_text(size = 12),    
    axis.title = element_text(size = 14),   
    legend.title = element_blank(),   
    legend.text = element_text(size = 12),  
    legend.position = "right",
    # 2. 增加这行代码，让新加的标题居中、加粗，保持排版一致性
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
  )

print(p)
ggplot2::ggsave(p,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_monocle_APOBEC3A.png',width = 12,height = 5)
ggplot2::ggsave(p,
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_monocle_APOBEC3A.pdf',
                width = 12,
                height = 5,
                device = cairo_pdf) # 加上这个关键参数

#########
# 提取 Fate 1 的数据
plot_data_fate1 <- data@meta.data %>%
  filter(cdsMM_sub1 == "Yes") %>%
  mutate(Cell_Fate = "Cell Fate 1")

# 提取 Fate 2 的数据
plot_data_fate2 <- data@meta.data %>%
  filter(cdsMM_sub2 == "Yes") %>%
  mutate(Cell_Fate = "Cell Fate 2")

# 将两组数据上下合并
plot_data_combined <- bind_rows(plot_data_fate1, plot_data_fate2)

p2 <- ggplot(plot_data_combined, aes(x = pseudotime, fill = Cell_Fate)) +
  geom_density(alpha = 0.5) +  
  scale_fill_manual(values = c("Cell Fate 1" = "#00BFC4", "Cell Fate 2" = "#F8766D")) + 
  labs(
    title = "Density of Pseudotime by Cell Fates",
    x = "Pseudotime",
    y = "Density"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),  
    axis.text = element_text(size = 12),    
    axis.title = element_text(size = 14),   
    legend.title = element_blank(),         # 隐藏图例标题，与上一张图保持一致
    legend.text = element_text(size = 12),  
    legend.position = "right",              # 图例放在右侧
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 标题居中显示
  )

# 打印图片
print(p2)

# 保存为 PNG
ggplot2::ggsave(p2, 
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_monocle_CellFates.png', 
                width = 12, 
                height = 5)

# 保存为 PDF (使用 cairo_pdf 保证字体正常导出)
ggplot2::ggsave(p2,
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_monocle_CellFates.pdf',
                width = 12,
                height = 5,
                device = cairo_pdf)







###### pr_graph_test_res #####
library(dplyr)
library(DescTools)
library(monocle3)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DescTools)      # 用于计算 AUC
library(splines)       # 用于拟合自然样条
library(clusterProfiler)
library(org.Hs.eg.db)
#####自然样条函数（Natural Splines）的广义线性模型
#主效应模型分析（Additive Model）显示，控制拟时间变量后，APOBEC3A 在 Fate 2 分支中的全轨迹平均表达量极显著高于 Fate 1（AUC 差异 = +182，q-value < 1e-99）。
#同时，交互效应模型（Interaction Model）的检验结果表明，该基因在两条命运分支间的动态变化趋势（形状）并无显著差异（q > 0.05）。
# ==============================================================================
# Step 1: 环境设置与数据导入
# ==============================================================================
out_dir <- "/public3/DSC/single_cell/Result/figer_new/monocle3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[1/7] 正在加载数据...\n")
# 加载 Monocle 对象
cds_MM_foam <- load_monocle_objects(directory_path = file.path(out_dir, "cds_MM_foam"))
cds_subset_1 <- load_monocle_objects(directory_path = file.path(out_dir, "cdsMM_subset1_healthy_macrophage"))
cds_subset_2 <- load_monocle_objects(directory_path = file.path(out_dir, "cdsMM_subset2_disease_foam"))

# 加载额外的拟时间数据（如有需要）
data_pt <- readRDS(file.path(out_dir, "data_pseudotime.rds"))

# ==============================================================================
# Step 2: 全局趋势基因筛选 (Graph Test)
# ==============================================================================
cat("[2/7] 正在执行全局轨迹相关基因筛选...\n")
pr_graph_test_res <- graph_test(cds_MM_foam, 
                                neighbor_graph = "principal_graph", 
                                cores = 4)

# 筛选显著动态基因 (q < 0.01)，按 Moran's I 排序
sig_dynamic_genes <- pr_graph_test_res %>% 
  filter(q_value < 0.01) %>% 
  arrange(desc(morans_I))

if(!"gene_short_name" %in% colnames(sig_dynamic_genes)) {
  sig_dynamic_genes$gene_short_name <- rownames(sig_dynamic_genes)
}

# 提取前 1000 个高变基因用于后续分支比对
target_genes <- head(rownames(sig_dynamic_genes), 1000)
write.csv(sig_dynamic_genes, file.path(out_dir, "T1_Global_Trajectory_Dependent_Genes.csv"), row.names = FALSE)

# ==============================================================================
# Step 3: 分支数据预处理 (命运锚定)
# ==============================================================================
cat("[3/7] 正在构建分支比对模型数据集...\n")

# 固化拟时间变量
colData(cds_MM_foam)$pseudo_t <- pseudotime(cds_MM_foam)

# 识别并分配祖细胞与分支细胞
selected_cells_1 <- colnames(cds_subset_1)
selected_cells_2 <- colnames(cds_subset_2)
overlapping_cells <- intersect(selected_cells_1, selected_cells_2)
pure_fate1_cells  <- setdiff(selected_cells_1, overlapping_cells)
pure_fate2_cells  <- setdiff(selected_cells_2, overlapping_cells)

# 分支 1 数据集
cds_f1 <- cds_MM_foam[, c(overlapping_cells, pure_fate1_cells)]
colData(cds_f1)$Branch_Fate <- "Fate 1"
colnames(cds_f1) <- paste0(colnames(cds_f1), "_f1")

# 分支 2 数据集
cds_f2 <- cds_MM_foam[, c(overlapping_cells, pure_fate2_cells)]
colData(cds_f2)$Branch_Fate <- "Fate 2"
colnames(cds_f2) <- paste0(colnames(cds_f2), "_f2")

# 合并为测试对象
cds_branch_test <- cbind(cds_f1, cds_f2)
colData(cds_branch_test)$Branch_Fate <- factor(colData(cds_branch_test)$Branch_Fate, levels = c("Fate 1", "Fate 2"))

# ==============================================================================
# Step 4: 拟合交互效应模型 (Interaction Model)
# ==============================================================================
cat("[4/7] 正在拟合分支特异性动态模型 (Time x Fate Interaction)...\n")

cds_target_test <- cds_branch_test[rowData(cds_branch_test)$gene_short_name %in% target_genes, ]

# 建立模型：考虑拟时间样条、分支分配、以及技术协变量
gene_fits <- fit_models(cds_target_test, 
                        model_formula_str = "~ splines::ns(pseudo_t, df=3) * Branch_Fate + orig.ident + nCount_RNA", 
                        cores = 8)

fit_coefs <- coefficient_table(gene_fits)

# 提取具有显著交互作用（即在不同分支间趋势不同）的基因
hetero_genes_clean <- fit_coefs %>%
  filter(grepl("Branch_Fate", term) & grepl(":", term)) %>%
  filter(status == "OK") %>%
  group_by(gene_short_name) %>%
  slice_min(order_by = q_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(q_value) %>%
  mutate(rank_q = row_number()) %>%
  arrange(desc(abs(estimate))) %>% 
  mutate(rank_est = row_number()) %>%
  arrange(rank_est)

# ==============================================================================
# Step 5: 计算 AUC 差异并生成最终汇总表
# ==============================================================================
cat("[5/7] 正在计算分支间表达量曲线的面积差 (AUC Diff)...\n")

# 准备预测用的标准化网格
max_pt <- max(colData(cds_branch_test)$pseudo_t, na.rm = TRUE)
pt_grid <- seq(0, max_pt, length.out = 100)
mean_ncount <- mean(colData(cds_branch_test)$nCount_RNA, na.rm = TRUE)
mean_sf <- mean(colData(cds_branch_test)$Size_Factor, na.rm = TRUE)
base_ident <- names(sort(table(colData(cds_branch_test)$orig.ident), decreasing = TRUE))[1]

pred_df_base <- expand.grid(
  pseudo_t = pt_grid,
  Branch_Fate = c("Fate 1", "Fate 2"),
  nCount_RNA = mean_ncount,
  orig.ident = base_ident,
  Size_Factor = mean_sf
)

# 计算显著基因的 AUC
sig_genes_for_auc <- hetero_genes_clean %>% filter(q_value < 0.01)

auc_results <- lapply(sig_genes_for_auc$gene_short_name, function(gene) {
  m_obj <- (gene_fits %>% filter(gene_short_name == gene) %>% pull(model))[[1]]
  if(is.null(m_obj)) return(NULL)
  
  temp_pred <- pred_df_base
  temp_pred$pred_expr <- predict(m_obj, newdata = temp_pred, type = "response")
  
  auc_1 <- AUC(x = temp_pred$pseudo_t[temp_pred$Branch_Fate == "Fate 1"], 
               y = temp_pred$pred_expr[temp_pred$Branch_Fate == "Fate 1"], method = "trapezoid")
  auc_2 <- AUC(x = temp_pred$pseudo_t[temp_pred$Branch_Fate == "Fate 2"], 
               y = temp_pred$pred_expr[temp_pred$Branch_Fate == "Fate 2"], method = "trapezoid")
  
  data.frame(gene_short_name = gene, AUC_Fate1 = auc_1, AUC_Fate2 = auc_2, AUC_Diff = auc_2 - auc_1)
})

auc_df <- do.call(rbind, auc_results)
final_summary_table <- sig_genes_for_auc %>% left_join(auc_df, by = "gene_short_name") %>% arrange(rank_est)

write.csv(final_summary_table, file.path(out_dir, "T2_Model_Based_Hetero_Genes_Final.csv"), row.names = FALSE)
save(gene_fits, final_summary_table, pred_df_base, cds_branch_test, 
     file = "/public3/DSC/single_cell/Result/figer_new/monocle3/Plot_Environment_Backup.RData")

cat("绘图核心环境已保存！\n")
# ==============================================================================
# Step 6: GO/KEGG 功能富集分析
# ==============================================================================
cat("[6/7] 正在进行基因功能富集分析...\n")

# 定义上调 (Fate 2) 和 下调 (Fate 1) 基因集 (基于 AUC 差值)
up_genes <- final_summary_table %>% filter(AUC_Diff > 0) %>% arrange(desc(AUC_Diff)) %>% pull(gene_short_name) %>% head(300)
down_genes <- final_summary_table %>% filter(AUC_Diff < 0) %>% arrange(AUC_Diff) %>% pull(gene_short_name) %>% head(300)

# 定义分析辅助函数
perform_enrichment <- function(genes, type = "GO") {
  gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (type == "GO") {
    res <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
    if(!is.null(res)) res <- simplify(res)
  } else {
    res <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa')
    if(!is.null(res)) res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  return(res)
}

go_up <- perform_enrichment(up_genes, "GO")
go_down <- perform_enrichment(down_genes, "GO")
kegg_up <- perform_enrichment(up_genes, "KEGG")
kegg_down <- perform_enrichment(down_genes, "KEGG")

# 保存富集结果表格
write.csv(as.data.frame(go_up), file.path(out_dir, "GO_Fate2_Up.csv"))
write.csv(as.data.frame(go_down), file.path(out_dir, "GO_Fate1_Up.csv"))
write.csv(as.data.frame(kegg_up), file.path(out_dir, "KEGG_Fate2_Up.csv"))
write.csv(as.data.frame(kegg_down), file.path(out_dir, "KEGG_Fate1_Up.csv"))
# ==============================================================================
# Step 7: 可视化绘制
# ==============================================================================
cat("[7/7] 正在生成分析图表...\n")

# 1. 绘制富集气泡图
p1 <- dotplot(go_up, showCategory = 15) + ggtitle("GO BP: Disease Foam (Fate 2) Up")
p2 <- dotplot(go_down, showCategory = 15) + ggtitle("GO BP: Healthy Mono (Fate 1) Up")
p <- p1 + p2
print(p)
ggsave(file.path(out_dir, "Plot_GO_Combined.pdf"), p1 + p2, width = 14, height = 7)

p3 <- dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG: Disease Foam (Fate 2) Up")
p4 <- dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG: Healthy Mono (Fate 1) Up")
p <- p3 + p4
print(p)

ggsave(file.path(out_dir, "Plot_KEGG_Combined.pdf"), p3 + p4, width = 14, height = 7)



# 2. 一键加载之前保存的拟合模型和环境 (加载后，gene_fits 等对象会自动出现)
load("/public3/DSC/single_cell/Result/figer_new/monocle3/Plot_Environment_Backup.RData")
# 2. 绘制 Top 10 基因趋势预测图
plot_model_prediction <- function(target_gene, model_fits, cds_obj, pred_base) {
  m_obj <- (model_fits %>% filter(gene_short_name == target_gene) %>% pull(model))[[1]]
  if(is.null(m_obj)) return(NULL)
  
  df <- pred_base
  df$Predicted_Exp <- predict(m_obj, newdata = df, type = "response")
  
  ggplot(df, aes(x = pseudo_t, y = Predicted_Exp, color = Branch_Fate)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("Fate 1" = "#1F77B4", "Fate 2" = "#FF7F0E")) +
    labs(title = target_gene, x = "Pseudotime", y = "Fitted Expression") +
    theme_classic() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
}

# ==========================================================
# 1. 精心挑选的 20个最具代表性 上调 & 下调 基因列表
# ==========================================================

# Fate 2 (Disease Foam) 极具代表性上调 20 基因 (炎症、脂质代谢、趋化)
top_20_up <- c(
  "CD36", "CCL2", "TNF", "FABP5", "S100A6", 
  "LGALS3", "STAT1", "APOBEC3A", "PFKP", "PLAU", 
  "ISG15", "GBP1", "SPHK1", "S100A12", "ALDH1A1", 
  "PLA2G7", "FLT1", "FBP1", "MX1", "CSTB"
)

###"S100A6", "FABP5","CSTB","LGALS3","CCL2","APOBEC3A","FBP1","GBP1", "ISG15","CMPK2", "MX1", "IFI44L","RSAD2","GBP4","OAS2","STAT1", "LY6E", "CD36","TNF" 
###CD36,FABP5,PFKP,STAT1,ISG15,APOBEC3A,CCL2,TNF,S100A6,LGALS3，

# Fate 1 (Healthy Macro) 极具代表性下调 20 基因 (稳态、修复、驻留)
top_20_down <- c(
  "PHACTR1", "CPM", "CD302", "CLEC10A", "THBS1", 
  "F13A1", "C3", "AXL", "FPR2", "FGL2", 
  "CD1D", "NR4A3", "SELL", "CCR2", "CLEC4E", 
  "OSM", "SLC40A1", "CH25H", "AOAH", "PTGS2"
)

###"FGL2", "AREG","THBS1","F13A1","C3","CLEC10A","AOAH","RGS18","CH25H", "FCER1A", "PLD4","SDS",
###CLEC10A，CD302，AXL，AOAH，CH25H，PHACTR1，FGL2，THBS1，F13A1，AREG


####对 CH25H 的备注（最新研究显示可能促炎）在传统的认知里，CH25H 产生的 25-羟基胆固醇是抗炎的（抑制炎症小体）；
####但在最新的研究（如某些 Immunity 或 Cell Metabolism 的文章）中，发现在特定的高脂或晚期斑块微环境中，它也可能放大炎症。
top_20_all <- c(
  
  # ----------------------------------------------------------------------------
  # 【上半区：Fate 2 (疾病泡沫细胞) —— 恶化五部曲】
  # 写作逻辑：强调微环境压力下，巨噬细胞如何从脂质过载走向彻底崩溃与破坏。
  # ----------------------------------------------------------------------------
  
  # 模块 1：原罪 (脂质吞噬与胞内超载)
  # -> 文章落脚点：证明分群的病理基础是氧化脂质的过度摄取。
  "CD36",      # 核心清道夫受体，无限制吞噬 oxLDL 的直接推手
  "FABP5",     # 脂肪酸结合蛋白，印证胞内脂质过载与早期脂毒性
  
  # 模块 2：窒息 (缺氧与代谢重编程)
  # -> 文章落脚点：经典“瓦博格效应”，细胞放弃有氧呼吸，转向促炎代谢。
  "PFKP",      # 糖酵解核心限速酶，驱动巨噬细胞向 M1 型/促炎极化的代谢引擎
  
  # 模块 3：错乱 (脂毒性引发的核酸应激与假性感染) —— ★ 文章最大创新点 (Novelty)
  # -> 文章落脚点：结晶刺破溶酶体导致的 DNA 损伤与极端应激反应。
  "STAT1",     # 干扰素信号核心转录因子，全面拉响细胞级“抗病毒”警报
  "ISG15",     # 强效干扰素刺激基因，在无菌微环境中异常放大炎症信号
  "APOBEC3A",  # 关键 RNA/DNA 编辑酶！连接脂毒性、核酸损伤与极端促炎的特色靶点
  
  # 模块 4：摇人 (致命趋化与炎症放大)
  # -> 文章落脚点：病态细胞如何通过旁分泌机制导致局部炎症失控。
  "CCL2",      # 最强单核细胞趋化因子(MCP-1)，呼叫血液单核细胞进入斑块“送死”
  
  # 模块 5：毁灭 (组织毒杀与物理崩塌)
  # -> 文章落脚点：泡沫细胞对正常血管结构的最终破坏，促发斑块破裂风险。
  "TNF",       # 炎症与细胞毒性总开关，诱导平滑肌细胞死亡及坏死核心形成
  "S100A6",    # 极强促炎钙结合蛋白，介导极端氧化应激与斑块失稳
  "LGALS3",    # 半乳糖凝集素-3，驱动晚期斑块的病理纤维化与坏死重塑
  
  
  # ----------------------------------------------------------------------------
  # 【下半区：Fate 1 (健康巨噬细胞) —— 守护者的五项美德】
  # 写作逻辑：镜像对比！强调疾病发生(Fate 2)本质上是这些健康保护机制的丧失。
  # ----------------------------------------------------------------------------
  
  # 模块 6：身份 (和平哨兵的良民证)
  # -> 文章落脚点：证明它们是原生的驻留型管家，而非促炎新兵。
  "CLEC10A",   # (CD301) 经典 M2 型/替代激活标志物，宣告抗炎与修复属性
  "CD302",     # 和平哨兵受体，负责非炎症性内吞与温和的免疫监视
  
  # 模块 7：清理 (无声的微环境打扫)
  # -> 文章落脚点：高效的日常清理机制，防止继发性坏死和无菌性炎症。
  "AXL",       # 胞葬作用(Efferocytosis)核心受体，安静吞噬死细胞，防止坏死核心形成
  "AOAH",      # 内毒素拆弹专家，就地灭活脂多糖，维持微环境免疫耐受
  
  # 模块 8：维稳 (脂质平衡与遗传学抗风险) —— ★ 临床转化价值提升
  # -> 文章落脚点：结合人类遗传学(GWAS)，探讨血管稳态的深层维持机制。
  "CH25H",     # 胆固醇羟化酶，维持早期局部脂质调度平衡与抗炎网络
  "PHACTR1",   # 冠心病顶流 GWAS 风险/保护基因！抑制巨噬细胞 M1 极化，维持血管稳定
  
  # 模块 9：镇静 (压制适应性免疫暴动)
  # -> 文章落脚点：巨噬细胞如何通过细胞间通讯平息 T 细胞的过度反应。
  "FGL2",      # 纤维介素蛋白-2，强效免疫镇静剂，压制周围 T 细胞的促炎活性
  
  # 模块 10：修缮 (细胞外基质与内皮缝合)
  # -> 文章落脚点：巨噬细胞作为组织重塑者，在斑块破裂前进行的最后挽救。
  "THBS1",     # 组织通讯使者，激活 TGF-β 促进微环境深度抗炎与伤口愈合
  "F13A1",     # 血管壁的基质缝合师，交联受损细胞外基质(ECM)，物理稳定血管结构
  "AREG"       # 双调蛋白，终极治愈因子，直接促进内皮/平滑肌生长与创伤修复
  
)

cat("40 个最具代表性的核心明星基因已加载！\n")

# ==========================================================
# 2. 批量生成动态预测图
# (调用之前定义的 plot_model_prediction 函数)
# ==========================================================
plots_up <- lapply(top_20_up, plot_model_prediction, gene_fits, cds_branch_test, pred_df_base)
plots_down <- lapply(top_20_down, plot_model_prediction, gene_fits, cds_branch_test, pred_df_base)

# 去除空图防止报错
plots_up <- Filter(Negate(is.null), plots_up)
plots_down <- Filter(Negate(is.null), plots_down)

# ==========================================================
# 3. 拼接并保存为 PDF (5列 x 4行 布局)
# ==========================================================
library(patchwork)

master_plot_up <- wrap_plots(plots_up, ncol = 5) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")

master_plot_down <- wrap_plots(plots_down, ncol = 5) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")

master_plot_up
master_plot_down
# 1. 批量生成 20 个基因的动态曲线图
# (请确保你的环境里还有 plot_model_prediction、gene_fits 等基础变量)
plots_all <- lapply(top_20_all, plot_model_prediction, gene_fits, cds_branch_test, pred_df_base)

# 2. 过滤可能为空的图（防止报错）
plots_all <- Filter(Negate(is.null), plots_all)

# 3. 拼接为一张宏大的主图 (4 行 x 5 列)
master_figure_all <- wrap_plots(plots_all, ncol = 5) + 
  plot_layout(guides = 'collect') & 
  labs(color = "Cell Fate") &  
  theme(
    legend.position = "bottom",
    # 保持你设定的加粗和大字体样式
    legend.title = element_text(size = 17, face = "bold"), 
    legend.text = element_text(size = 15),
    # 增加图例 key 的尺寸，确保在高分辨率 PDF 中清晰可见
    legend.key.size = unit(1.5, "cm") 
  )

master_figure_all
# 4. 设置输出路径并保存
out_dir <- "/public3/DSC/single_cell/Result/figer_new/monocle3"

# 保存为大尺寸 PDF，保证 20 张小图的高清排版 (推荐 16x12)
ggsave(file.path(out_dir, "MainFigure_Top20_Complete_Story.pdf"), 
       plot = master_figure_all, 
       width = 12, height = 8)

cat("完美包含 20 个基因的终极主图已成功保存至：", out_dir, "\n")

# ==============================================================================
# 使用 top_20_all 绘制轨迹空间上的基因表达分布图 (Feature Plot)
# ==============================================================================

p <- plot_cells(
  cds_MM_foam,
  genes = top_20_all,                   # 直接传入 20 个基因的列表
  label_groups_by_cluster = FALSE,
  cell_size = 1.0,
  group_label_size = 4,
  show_trajectory_graph = TRUE
) +
  coord_flip() +                        # 翻转坐标轴
  facet_wrap(~feature_label, ncol = 5, scales = "free") +  # 按基因名分面，排成 5 列
  scale_color_gradientn(
    colors = c("gray90", "red2", "red3"),
    values = c(0, 0.5, 1)               # 按照你的要求更新颜色映射区间
  ) +
  labs(color = "Expression") + 
  theme_dr() + 
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 20, face = "bold"),  # 基因名字体放大并加粗
    strip.background = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    text = element_text(family = "Arial"),
    legend.position = "right"           # 将图例放在右侧，让出更多空间给基因图
  )

# 屏幕预览
print(p)

# ==============================================================================
# 保存高清图片
# ==============================================================================
out_png <- '/public3/DSC/single_cell/Result/figer_new/monocle3/S2.2_Pseudotime_plot_cells_Top20.png'
out_pdf <- '/public3/DSC/single_cell/Result/figer_new/monocle3/S2.2_Pseudotime_plot_cells_Top20.pdf'

# 保存 PNG (20x16 完美适配 4行5列)
ggplot2::ggsave(p, filename = out_png, width = 20, height = 16)

# 保存 PDF (建议和 PNG 尺寸一致，保证矢量图排版比例完美)
ggplot2::ggsave(
  p, 
  filename = out_pdf, 
  width = 20, 
  height = 16, 
  device = grDevices::cairo_pdf,  # 强制指定使用 Cairo 引擎
  limitsize = FALSE
)


# ====================================================================
# 1. 数据准备 (保持你的筛选逻辑)
# ====================================================================
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

out_dir <- "/public3/DSC/single_cell/Result/figer_new/monocle3"
csv_path <- file.path(out_dir, "T2_Model_Based_Hetero_Genes_Final.csv")
final_summary_table <- read.csv(csv_path, stringsAsFactors = FALSE)

# 筛选驱动基因
fate2_top_genes <- final_summary_table %>% 
  filter(AUC_Diff > 0) %>% arrange(desc(AUC_Diff)) %>% head(25) %>% pull(gene_short_name)
fate1_top_genes <- final_summary_table %>% 
  filter(AUC_Diff < 0) %>% arrange(AUC_Diff) %>% head(25) %>% pull(gene_short_name)
heatmap_genes <- unique(c(fate1_top_genes, fate2_top_genes))

# ====================================================================
# 2. 数据平滑与标准化 (使用 bin 方法)
# ====================================================================
get_unscaled_smooth_mat <- function(cds, global_pt, genes, num_bins = 100) {
  cell_names <- colnames(cds)
  pt <- global_pt[cell_names]
  pt <- pt[!is.na(pt) & is.finite(pt)]
  valid_cells <- names(pt)
  expr_mat <- as.matrix(normalized_counts(cds)[genes, valid_cells])
  pt_ordered <- sort(pt)
  expr_mat <- expr_mat[, names(pt_ordered)]
  bins <- cut(pt_ordered, breaks = num_bins, labels = FALSE)
  bins_factor <- factor(bins, levels = 1:num_bins)
  
  smooth_mat <- apply(expr_mat, 1, function(x) tapply(x, bins_factor, mean, na.rm = TRUE))
  smooth_mat <- t(smooth_mat) 
  smooth_mat <- t(apply(smooth_mat, 1, function(row_data) {
    if (all(is.na(row_data))) return(rep(0, length(row_data)))
    valid_idx <- which(!is.na(row_data))
    valid_vals <- row_data[!is.na(row_data)]
    stats::approx(x = valid_idx, y = valid_vals, xout = seq_along(row_data), rule = 2)$y
  }))
  return(smooth_mat)
}


unscaled_f1 <- get_unscaled_smooth_mat(cds_f1, pseudotime(cds_MM_foam), heatmap_genes, num_bins = 100)
unscaled_f2 <- get_unscaled_smooth_mat(cds_f2, pseudotime(cds_MM_foam), heatmap_genes, num_bins = 100)

# Fate 1 翻转 (从共同起点向左发散)，Fate 2 正常 (从共同起点向右)
unscaled_f1_rev <- unscaled_f1[, ncol(unscaled_f1):1]
combined_unscaled <- cbind(unscaled_f1_rev, unscaled_f2)
combined_scaled <- t(scale(t(combined_unscaled)))
combined_scaled[combined_scaled > 2.5] <- 2.5
combined_scaled[combined_scaled < -2.5] <- -2.5

# ====================================================================
# 3. 聚类与排序 (彻底解决双重 Label 问题)
# ====================================================================
set.seed(123)
km_res <- kmeans(combined_scaled, centers = 4, nstart = 25)
row_clusters <- km_res$cluster

# 按聚类和表达趋势排序
row_order <- order(row_clusters, rowMeans(combined_scaled[, 101:200])) 
combined_scaled_sorted <- combined_scaled[row_order, ]
row_clusters_sorted <- row_clusters[row_order]

# 拆分
fate1_scaled <- combined_scaled_sorted[, 1:100]
fate2_scaled <- combined_scaled_sorted[, 101:200]

# ====================================================================
# 4. 坐标轴标签准备 (不再使用 decorate，直接注入)
# ====================================================================
# 【修正点】：同样使用正确的 cds_f1 和 cds_f2 变量
t1_r <- range(pseudotime(cds_MM_foam)[colnames(cds_f1)], na.rm = TRUE)
t2_r <- range(pseudotime(cds_MM_foam)[colnames(cds_f2)], na.rm = TRUE)

# 创建 100 个点的标签向量，只在特定位置显示数字
make_lab <- function(r, reverse = FALSE) {
  vals <- if(reverse) seq(r[2], r[1], length.out = 100) else seq(r[1], r[2], length.out = 100)
  labs <- rep("", 100)
  at <- c(1, 50, 100)
  labs[at] <- as.character(round(vals[at], 1))
  return(labs)
}

fate1_labs <- make_lab(t1_r, reverse = TRUE)
fate2_labs <- make_lab(t2_r, reverse = FALSE)

# ====================================================================
# 5. 构建热图
# ====================================================================
col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("#4575b4", "white", "#d73027"))
cluster_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3")

# 左侧热图 (Fate 1)
ht1 <- Heatmap(fate1_scaled, 
               name = "Z-score", col = col_fun,
               cluster_columns = FALSE, cluster_rows = FALSE,
               row_split = row_clusters_sorted, 
               row_title = NULL, # 禁用默认的分组标题防止干扰
               column_title = "Cell Fate 1",
               column_labels = fate1_labs,
               column_names_gp = gpar(fontsize = 9),
               column_names_rot = 0,
               show_row_names = FALSE,
               left_annotation = rowAnnotation(
                 Cluster = anno_block(gp = gpar(fill = cluster_colors),
                                      labels = c("C1", "C2", "C3", "C4"),
                                      labels_gp = gpar(col = "white", fontface = "bold"))
               ),
               border = TRUE)

# 中间基因名 (紧凑型)
ht_mid <- rowAnnotation(
  gene = anno_text(rownames(combined_scaled_sorted), 
                   gp = gpar(fontsize = 8, fontface = "italic"),
                   just = "center", location = 0.5),
  width = unit(1.5, "cm")
)

# 右侧热图 (Fate 2)
ht2 <- Heatmap(fate2_scaled, 
               col = col_fun, show_heatmap_legend = FALSE,
               cluster_columns = FALSE, cluster_rows = FALSE,
               row_split = row_clusters_sorted,
               row_title = NULL,
               column_title = "Cell Fate 2",
               column_labels = fate2_labs,
               column_names_gp = gpar(fontsize = 9),
               column_names_rot = 0,
               show_row_names = FALSE,
               border = TRUE)

final_plot <- ht1 + ht_mid + ht2

# ====================================================================
# 6. 保存导出
# ====================================================================
pdf_path <- file.path(out_dir, "F4.1_Final_Clean_Heatmap.pdf")
cairo_pdf(pdf_path, width = 10, height = 8)

draw(final_plot, 
     #column_title = "Divergent Gene Dynamics",
     column_title_gp = gpar(fontsize = 16, fontface = "bold"),
     ht_gap = unit(1, "mm")) # 缩小热图间的物理间隙

dev.off()


# 加载必备包
library(dplyr)
library(ggplot2)
library(DescTools)

cat("正在启动命运抉择窗口 (Decision Window) 扫描分析...\n")

# ==============================================================================
# 1. 精准定位“命运岔路口”的时间坐标 (完全脱离 cds_MM_foam 的终极安全版)
# ==============================================================================
# 提取测试对象的完整元数据
meta_df <- as.data.frame(colData(cds_branch_test))

# 提取 Fate 1 和 Fate 2 的细胞数据
meta_f1 <- meta_df %>% dplyr::filter(Branch_Fate == "Fate 1")
meta_f2 <- meta_df %>% dplyr::filter(Branch_Fate == "Fate 2")

# 还原基础细胞名
base_f1 <- gsub("_f1$", "", rownames(meta_f1))
base_f2 <- gsub("_f2$", "", rownames(meta_f2))

# 取交集，找到共同的祖细胞
overlap_base <- intersect(base_f1, base_f2)

if (length(overlap_base) == 0) {
  stop("错误：没有找到共同的祖细胞！请检查之前的分支选择是否正确。")
}

# 映射回带后缀的细胞名（随便取 _f1 或 _f2 都可以，因为它们在岔路口前的 pseudo_t 是完全一样的）
overlap_f1_names <- paste0(overlap_base, "_f1")

# 【核心修复】：直接从 cds_branch_test 中提取物理坐标，绝不报错
branch_t <- max(meta_df[overlap_f1_names, "pseudo_t"], na.rm = TRUE)
max_t <- max(meta_df$pseudo_t, na.rm = TRUE)

cat(sprintf("-> 发现命运分叉点坐标: Pseudotime = %.2f (总轨迹长度: %.2f)\n", branch_t, max_t))

# 定义“抉择窗口”：岔路口之后 15% 的时间段
window_end <- branch_t + (max_t * 0.15)

# ==============================================================================
# 2. 扫描所有显著基因，提取“波峰坐标”与“早期散度”
# ==============================================================================
# sig_genes_for_auc 是之前生成的显著基因表，如果丢失可替换为 final_summary_table
bifurcation_stats <- lapply(final_summary_table$gene_short_name, function(gene) {
  
  # 提取该基因的拟合模型 (严格使用 dplyr::)
  m_obj <- (gene_fits %>% dplyr::filter(gene_short_name == gene) %>% dplyr::pull(model))[[1]]
  if(is.null(m_obj)) return(NULL)
  
  # 使用标准化网格进行预测
  temp_pred <- pred_df_base
  temp_pred$pred_expr <- predict(m_obj, newdata = temp_pred, type = "response")
  
  # 指标 A：寻找“瞬时波峰” (Peak Time)
  peak_row <- temp_pred %>% dplyr::arrange(desc(pred_expr)) %>% dplyr::slice(1)
  peak_time <- peak_row$pseudo_t
  
  # 判断波峰是否正好落在岔路口附近 (正负 10% 范围内)
  is_transient_driver <- abs(peak_time - branch_t) < (max_t * 0.1)
  
  # 指标 B：寻找“极速开关” (Early Divergence)
  early_pred <- temp_pred %>% dplyr::filter(pseudo_t >= branch_t & pseudo_t <= window_end)
  
  auc_f1_early <- AUC(x = early_pred$pseudo_t[early_pred$Branch_Fate == "Fate 1"], 
                      y = early_pred$pred_expr[early_pred$Branch_Fate == "Fate 1"], method = "trapezoid")
  auc_f2_early <- AUC(x = early_pred$pseudo_t[early_pred$Branch_Fate == "Fate 2"], 
                      y = early_pred$pred_expr[early_pred$Branch_Fate == "Fate 2"], method = "trapezoid")
  
  early_diff <- auc_f2_early - auc_f1_early
  
  # 返回单个基因的统计结果
  data.frame(
    gene_short_name = gene,
    peak_time = peak_time,
    is_transient_driver = is_transient_driver,
    early_auc_diff = early_diff,
    max_expr = peak_row$pred_expr
  )
})

# 合并结果
bifurcation_df <- do.call(rbind, bifurcation_stats)

cat(sprintf("瞬时波峰 (Transient Driver) 的伪时间有效范围是: [%.2f,  %.2f]\n", 
            branch_t - (max_t * 0.1), 
            branch_t + (max_t * 0.1)))
# 统计 TRUE 和 FALSE 的基因数量
table(bifurcation_df$is_transient_driver, useNA = "ifany")

# ==============================================================================
# 3. 筛选并导出两类关键候选分子
# ==============================================================================
# 候选池 1：瞬时波峰型 (Transient Drivers)
transient_drivers <- bifurcation_df %>% 
  dplyr::filter(is_transient_driver == TRUE) %>% 
  dplyr::arrange(desc(max_expr)) 

# 候选池 2：极速开关型 (Early Switchers) - 走向 Fate 2 (疾病泡沫) 瞬间激增的基因
early_switch_fate2 <- bifurcation_df %>% 
  dplyr::arrange(desc(early_auc_diff))

# 候选池 3：极速开关型 (Early Switchers) - 走向 Fate 1 (健康巨噬) 瞬间激增的基因
early_switch_fate1 <- bifurcation_df %>% 
  dplyr::arrange(early_auc_diff)

library(openxlsx)
out_dir <- "/public3/DSC/single_cell/Result/figer_new/monocle3"

# 将三个表格装入一个命名列表 (List) 中
# 列表前面的名字 (如 "Transient_Drivers") 就会自动变成 Excel 里的 Sheet 名称
export_list <- list(
  "Transient_Drivers" = transient_drivers,
  "Early_Switch_Fate2" = early_switch_fate2,
  "Early_Switch_Fate1" = early_switch_fate1
)

# 一键导出为 Excel 文件
openxlsx::write.xlsx(
  export_list, 
  file = file.path(out_dir, "T3_T5_Bifurcation_Key_Genes.xlsx"), 
  rowNames = FALSE
)

cat("岔路口关键基因提取完毕！所有表格已合并保存至 T3_T5_Bifurcation_Key_Genes.xlsx 的不同 Sheet 中。\n")


# ==============================================================================
# 4. 绘制重点基因并标注“分叉警戒线”
# ==============================================================================
plot_bifurcation_dynamics <- function(target_gene, model_fits, pred_base, branch_t) {
  
  # 【安全检查 1】：确认基因是否在模型列表中 (严格使用 dplyr:: 前缀)
  target_data <- model_fits %>% dplyr::filter(gene_short_name == target_gene)
  
  if(nrow(target_data) == 0) {
    message(paste("⚠️ 警告: 基因", target_gene, "不在拟合模型中，已自动跳过。"))
    return(NULL)
  }
  
  # 提取模型 (严格使用 dplyr::)
  m_obj <- (target_data %>% dplyr::pull(model))[[1]]
  
  # 【安全检查 2】：确认模型是否拟合成功
  if(is.null(m_obj)) {
    message(paste("⚠️ 警告: 基因", target_gene, "的模型为空，已自动跳过。"))
    return(NULL)
  }
  
  # 预测表达量
  df <- pred_base
  df$Predicted_Exp <- predict(m_obj, newdata = df, type = "response")
  
  # 绘图
  p <- ggplot(df, aes(x = pseudo_t, y = Predicted_Exp, color = Branch_Fate)) +
    # 添加代表“岔路口”的垂直虚线
    geom_vline(xintercept = branch_t, linetype = "dashed", color = "grey50", size = 1) +
    # 标注岔路口文本
    annotate("text", x = branch_t - 0.9, 
             y = min(df$Predicted_Exp) + (max(df$Predicted_Exp) - min(df$Predicted_Exp)) * 0.2, 
             label = "Bifurcation", angle = 90, color = "grey30", size = 4) +
    
    geom_line(size = 1.5, alpha = 0.8) +
    scale_color_manual(values = c("Fate 1" = "#1F77B4", "Fate 2" = "#FF7F0E")) +
    labs(title = target_gene, x = "Pseudotime", y = "Fitted Expression") +
    theme_classic(base_size = 14) + 
    theme(
      legend.position = "bottom", 
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
  
  return(p)
}

# 从你筛选出的瞬时或早期开关基因中挑几个最有意思的画图 (比如假设挑了这三个)
# target_tfs <- c(early_switch_fate2$gene_short_name[1:2], transient_drivers$gene_short_name[1])
# 这里以之前的明星分子为例演示：
demo_genes <- c("OLR1", "CLEC7A", "PTPRC", "PIM3", "APOBEC3A", 
                "IFI16", "ISG15", "IFITM2", "GBP1", "CARD16", 
                "SRGN", "CCL4", "TIMP1", "FTH1", "TNFSF10")

cat("开始批量绘制 10 个明星基因的分支动力学曲线...\n")

plots_bifurcation <- lapply(demo_genes, plot_bifurcation_dynamics, 
                            gene_fits, pred_df_base, branch_t)

library(patchwork)
p_bifurcation_combined <- wrap_plots(plots_bifurcation, ncol = 5) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")

print(p_bifurcation_combined)

ggsave(file.path(out_dir, "F5_Bifurcation_Key_TFs.pdf"), p_bifurcation_combined, width = 15, height = 10, device = cairo_pdf)


# 为了深入探究巨噬细胞向病理性表型（Fate 2）分化的核心转录驱动机制，我们沿拟时间轴（Pseudotime）重构了关键命运决定因子的动态表达轨迹
# 。基因表达的分叉动力学（Bifurcation dynamics）揭示了一场具有高度时序特异性的病理级联反应。在命运抉择的极早期阶段（分叉点前及附近），
# 微环境压力感受器（OLR1, CLEC7A）与跨膜信号转导枢纽（PTPRC, PIM3）率先响应，构成了细胞应激的前置信号。跨越分叉点后，
# 病理分支（Fate 2）瞬间爆发了剧烈的、具有高度分支特异性的转录突变：核酸编辑酶 APOBEC3A 与胞内核酸感受器 IFI16 呈现典型的瞬时波峰（Transient peaks），
# 精准触发了强烈的 I 型干扰素（Type I IFN）级联放大反应，导致 ISG15、IFITM2 和 GBP1 等干扰素刺激基因如同剪刀般与健康分支（Fate 1）彻底割裂并呈指数级飙升。
# 伴随这种由内源性核酸紊乱引发的“假性感染（Pseudo-infection）”状态，巨噬细胞迅速启动了不可逆的促炎破坏程序，
# 集中表现为炎症小体调节阀（CARD16）、炎症颗粒包装与趋化因子（SRGN, CCL4）的相继激活，
# 以及与组织基质破坏（TIMP1）、代谢崩塌（FTH1）和死亡诱导（TNFSF10）直接相关的终局效应分子被持续性极度上调。
# 这一完美的时序动态轨迹有力地证明，由极早期细胞应激诱发的内源性核酸与干扰素风暴，是驱动巨噬细胞向破坏性表型极化的关键“命运扳机”。

library(hdWGCNA)
library(dplyr)
library(tidyr)
library(ggplot2)

cat("正在重构 F2.6 全局模块双重分化轨迹 (深度复刻 hdWGCNA 原版视觉风格)...\n")

# ==============================================================================
# 1. 定义严谨的目标顺序与提取数据
# ==============================================================================
# 严格按照你提供的顺序（已修正拼写：meganta -> magenta, salmaon -> salmon）
target_order <- c("yellow", "blue", "turquoise", "purple", "brown", "pink", 
                  "magenta", "green", "lightcyan", "red", "tan", "greenyellow", 
                  "salmon", "cyan", "black", "midnightblue")

# 获取原始模块数据列名，防守型编程：只提取 target_order 中真实存在的列
all_mes <- colnames(GetMEs(data, harmonized = FALSE))
valid_modules <- intersect(target_order, all_mes)

# ==============================================================================
# 2. 提取并复制祖细胞 (基于 meta.data，完美规避重名报错)
# ==============================================================================
# 提取 Fate 1 的数据
df_fate1 <- data@meta.data %>%
  dplyr::filter(cdsMM_sub1 == "Yes") %>%
  dplyr::mutate(Cell_Fate = "Cell Fate 1")

# 提取 Fate 2 的数据
df_fate2 <- data@meta.data %>%
  dplyr::filter(cdsMM_sub2 == "Yes") %>%
  dplyr::mutate(Cell_Fate = "Cell Fate 2")

# 上下合并两组数据，在分叉点完美复制共同祖细胞
merged_df <- dplyr::bind_rows(df_fate1, df_fate2)

# ==============================================================================
# 3. 数据结构重塑 (转换为长数据)
# ==============================================================================
long_df <- merged_df %>%
  dplyr::select(pseudotime, Cell_Fate, dplyr::all_of(valid_modules)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(valid_modules), 
    names_to = "Module",
    values_to = "ME_Score"
  ) %>%
  # 强制转换为因子，严格锁定为你指定的 16 个模块顺序
  dplyr::mutate(Module = factor(Module, levels = target_order))

# ==============================================================================
# 4. 绘制深度复刻 hdWGCNA 风格的大图
# ==============================================================================
p_all_modules <- ggplot(long_df, aes(x = pseudotime, y = ME_Score, color = Cell_Fate, fill = Cell_Fate)) +
  
  # 复刻 1：添加 y = 0 的水平基准线
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.8) +
  
  # 保留我们之前设计的 F5 命运岔路口垂直虚线
  geom_vline(xintercept = branch_t, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # 复刻 2：巧妙使用 stat_summary_bin 模拟原图稀疏的散点 (避免全量细胞打点导致画面糊成一团)
  #stat_summary_bin(fun = "mean", geom = "point", bins = 25, size = 1.5, alpha = 0.8) +
  
  # 复刻 3：添加带有“预测范围”的平滑曲线 (se = TRUE)
  geom_smooth(
    method = "gam",    
    formula = y ~ s(x, bs = "tp"), 
    se = TRUE,       # 开启置信区间阴影
    linewidth = 1.2,
    alpha = 0.2      # 调整阴影透明度，显得高级且不遮挡
  ) +
  
  # 按原图设置为 4 列布局
  facet_wrap(~ Module, scales = "free_y", ncol = 4) + 
  
  # 依然保持红蓝配色以区分双重命运（如果强行用原图的模块色，你将无法分辨哪条是Fate1，哪条是Fate2）
  scale_color_manual(values = c("Cell Fate 1" = "#1F77B4", "Cell Fate 2" = "#FF7F0E"), name = "Fate Branch") +
  scale_fill_manual(values = c("Cell Fate 1" = "#1F77B4", "Cell Fate 2" = "#FF7F0E"), name = "Fate Branch") +
  
  labs(
    title = "hdWGCNA Module Trajectories across Cell Fates",
    x = "Pseudotime",
    y = "Module Eigengene"
  ) +
  
  # 复刻 4：使用 theme_bw() 为基础，手工调校回 hdWGCNA 原版的灰色质感
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 20)),
    
    # 恢复原图典型的浅灰色背景板与白色网格线
    panel.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.6),
    panel.grid.minor = element_line(color = "white", linewidth = 0.3),
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
    
    # 恢复原图经典的灰色分面标题框 (Facet strip)
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 12, face = "plain", color = "black"), 
    
    axis.title.x = element_text(size = 14, margin = margin(t = 15)), 
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),  
    axis.text = element_text(color = "black"),
    
    # 图例放底部
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.width = unit(2, "cm")
  )

print(p_all_modules)

# ==============================================================================
# 5. 高清输出保存
# ==============================================================================
ggplot2::ggsave(p_all_modules, 
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.6_MM_monocle_Pseudotim_AllModules_v2.png', 
                width = 16, height = 12)

ggplot2::ggsave(p_all_modules, 
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.6_MM_monocle_Pseudotim_AllModules_v2.pdf', 
                width = 16, height = 12, device = cairo_pdf)

cat("完美！带有预测区间、经典灰底质感且红蓝分离的 F2.6 模块总图已保存！\n")



# ==============================================================================
# 究极进化版：F2.6 拟时序命运分叉密度图 + 底部 Celltype_raw1 主导权轨道
# ==============================================================================
library(ggplot2)
library(dplyr)
library(patchwork) # 必须加载 patchwork 用于无缝拼图

cat("正在计算两条命运分支中，各 Celltype_raw1 的动态主导权...\n")

# ==========================================================
# 0. 严谨的数据检查：确保数据中包含 Celltype_raw1 列
# ==========================================================
if(!"Celltype_raw1" %in% colnames(plot_data_combined)) {
  # 如果 combined 之前没合并这个信息，我们直接从原始 data 的 metadata 抓取
  plot_data_combined$Celltype_raw1 <- data@meta.data[rownames(plot_data_combined), "Celltype_raw1"]
}

# ==========================================================
# 1. 核密度统筹与主导权计算 (核心算法)
# ==========================================================
pt_min <- min(plot_data_combined$pseudotime, na.rm = TRUE)
pt_max <- max(plot_data_combined$pseudotime, na.rm = TRUE)
pt_grid <- seq(pt_min, pt_max, length.out = 1000)

# 定义高阶算法：计算特定 Fate 内部，哪个 Celltype 密度最高
calc_dominant_type <- function(fate_name) {
  df <- plot_data_combined %>% 
    dplyr::filter(Cell_Fate == fate_name) %>%
    dplyr::filter(!is.na(Celltype_raw1) & !is.na(pseudotime))
  
  total_cells <- nrow(df)
  types <- unique(df$Celltype_raw1)
  
  # 计算每种细胞类型的拟时间加权密度
  dens_matrix <- sapply(types, function(t) {
    sub_df <- df %>% dplyr::filter(Celltype_raw1 == t)
    # 如果某种细胞数量极少，密度计算会报错，直接赋 0
    if(nrow(sub_df) < 5) return(rep(0, length(pt_grid))) 
    
    d <- density(sub_df$pseudotime, na.rm = TRUE)
    # 核心：将密度乘以该细胞类型的数量占比，还原真实的绝对主导地位
    weight <- nrow(sub_df) / total_cells
    y_val <- approx(d$x, d$y, xout = pt_grid, rule = 2)$y * weight
    y_val[is.na(y_val)] <- 0
    return(y_val)
  })
  
  # 找出在每个拟时间点上，密度最大（主导）的细胞类型
  max_indices <- apply(dens_matrix, 1, which.max)
  dominant_types <- types[max_indices]
  
  data.frame(
    pseudotime = pt_grid, 
    Dominant_Type = dominant_types, 
    Cell_Fate = fate_name
  )
}

# 分别计算两条轨道的状况，并合并
dom_fate1 <- calc_dominant_type("Cell Fate 1")
dom_fate2 <- calc_dominant_type("Cell Fate 2")
dom_combined <- dplyr::bind_rows(dom_fate1, dom_fate2)

# 计算进度条上文字标签的最佳居中位置 (取各自主导区间的拟时间中位数)
labels_combined <- dom_combined %>%
  dplyr::group_by(Cell_Fate, Dominant_Type) %>%
  dplyr::summarize(mid_x = median(pseudotime), .groups = 'drop')

# ==========================================================
# 2. 构建上半部分：主密度图 (宏观命运演化)
# ==========================================================
# 自动寻找文字锚点
dens_f1 <- density((plot_data_combined %>% dplyr::filter(Cell_Fate == "Cell Fate 1"))$pseudotime, na.rm=TRUE)
dens_f2 <- density((plot_data_combined %>% dplyr::filter(Cell_Fate == "Cell Fate 2"))$pseudotime, na.rm=TRUE)
max_density <- max(c(dens_f1$y, dens_f2$y), na.rm = TRUE)

p_main <- ggplot(plot_data_combined, aes(x = pseudotime, fill = Cell_Fate, color = Cell_Fate)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  geom_vline(xintercept = branch_t, linetype = "dashed", color = "grey50", linewidth = 1) +
  annotate("text", x = branch_t - 0.9, y = max_density * 0.85, 
           label = "Bifurcation", angle = 90, color = "grey30", size = 5, fontface = "italic") +
  scale_fill_manual(values = c("Cell Fate 1" = "#1F77B4", "Cell Fate 2" = "#FF7F0E")) +
  scale_color_manual(values = c("Cell Fate 1" = "#1F77B4", "Cell Fate 2" = "#FF7F0E")) +
  labs(title = "Divergent Fates and Dominant Cell Types", y = "Fate Density") +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(color = "black", size = 12),
    # 去除主图的 X 轴，将其无缝过渡给下方的进度条
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    # 图例放上面，平衡底部长条的视觉比重
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_x_continuous(limits = c(pt_min, pt_max), expand = c(0, 0))

# ==========================================================
# 3. 构建下半部分：双轨道主导权进度条 (微观类型接力)
# ==========================================================
p_bars <- ggplot(dom_combined, aes(x = pseudotime, y = Cell_Fate, fill = Dominant_Type)) +
  # height = 0.8 可以在两条进度条之间留出极其高级的白色空隙
  geom_tile(height = 0.8, alpha = 0.9) + 
  # 在进度条上也贯穿分叉虚线
  geom_vline(xintercept = branch_t, linetype = "dashed", color = "black", linewidth = 1) +
  # 打印细胞类型名字 (居中显示)
  geom_text(data = labels_combined, aes(x = mid_x, y = Cell_Fate, label = Dominant_Type),
            color = "black", fontface = "bold", size = 4, inherit.aes = FALSE) +
  # 翻转 Y 轴顺序，使得 Fate 1 轨道在上方，Fate 2 在下方
  scale_y_discrete(limits = rev(c("Cell Fate 1", "Cell Fate 2"))) +
  labs(x = "Pseudotime", fill = "Dominant Cell Type") +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(color = "black", size = 12),
    # 底部图例用于展示各色块对应的具体细胞亚群名称
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5) 
  ) +
  scale_x_continuous(limits = c(pt_min, pt_max), expand = c(0, 0))

# ==========================================================
# 4. 使用 Patchwork 进行终极无缝拼接
# ==========================================================
# 按照 4:1 的物理高度比例，将宏观密度图与双轨道进度条合并
p_final <- p_main / p_bars + plot_layout(heights = c(4, 1.2))

print(p_final)

# ==========================================================
# 5. 高清导出
# ==========================================================
ggplot2::ggsave(p_final, 
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.6_Pseudotime_Density_Celltype_Tracks.pdf', 
                width = 12, height = 7, device = cairo_pdf)
ggplot2::ggsave(p_final, 
                filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.6_Pseudotime_Density_Celltype_Tracks.png', 
                width = 12, height = 7)

cat("带有双轨道 Celltype_raw1 主导权进度条的顶级图表绘制完成！\n")









##### APOBEC3A ####
p <- plot_cells(
  cds_MM_foam,
  genes = "APOBEC3A",
  label_groups_by_cluster = FALSE,
  cell_size = 1.0,
  group_label_size = 6,           # 调大潜在的分组标签
  show_trajectory_graph = TRUE
) +
  coord_flip() + 
  
  # 按样本类型分面
  facet_wrap(~Sample_Type, ncol = 2, scales = "free") + 
  
  scale_color_gradientn(
    colors = c("gray90", "red2", "red3"),
    values = c(0, 0.25, 1) 
  ) +
  
  labs(color = "APOBEC3A") + 
  theme_dr() + 
  theme( 
    panel.grid = element_blank(),
    # 1. 分面标题字体加大加粗 (Normal vs Plaque 等)
    strip.text = element_text(size = 25),  
    strip.background = element_blank(),
    
    # 2. 坐标轴标题字体加大 (UMAP_1, UMAP_2)
    axis.title = element_text(size = 20),
    # 3. 坐标轴刻度字体加大
    axis.text = element_text(size = 15),
    
    # 4. 图例标题和标签加大
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    
    plot.title = element_blank(),
    text = element_text(family = "Arial")
  ) +
  
  # 5. 调整 annotate 标签的大小 (size 提高到 8)
  # 注意：annotate 中的 size 单位与 theme 不同，8 已经非常醒目
  annotate("label",  
           x = -5.5, y = 1.8, 
           label = "Cell fate 2",
           hjust = -0.5, vjust = 2,
           size = 8,           # 从 6 调到 8
           color = "black",
           fill = "grey90",  
           label.padding = unit(0.2, "lines"),  
           label.r = unit(0.05, "lines")) +  
  annotate("label",
           x = -9, y = -2.5,
           label = "Cell fate 1",
           hjust = 1.2, vjust = -1,
           size = 8,           # 从 6 调到 8
           color = "black",
           fill = "grey90",
           label.padding = unit(0.2, "lines"),
           label.r = unit(0.05, "lines"))

# 打印预览
print(p)

# 保存文件
ggplot2::ggsave(p, filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_monocle_APOBEC3A.png', width = 12, height = 5)
ggplot2::ggsave(p, filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F2.5_monocle_APOBEC3A.pdf', width = 18, height = 7, device = cairo_pdf)


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
    APOBEC3A = GetAssayData(data, assay = "RNA", layer = "data")["APOBEC3A", subset_cells],
    subset = names(subset_list)[i],
    Sample_Type = data@meta.data[subset_cells, "Sample_Type"]
  )
  
  combined_data <- rbind(combined_data, subset_df)
}


fate_colors_named <- c(
  "subset1" = "#1F77B4",  # 蓝线
  "subset2" = "#FF7F0E"   # 橙线
)

# 2. 定义图例显示名称的命名向量 (数据真名 = 图例上想显示的漂亮名字)
fate_labels_named <- c(
  "subset1" = "Cell Fate 1",
  "subset2" = "Cell Fate 2"
)

APOBEC3A_AC_PA <- ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = subset)) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~Sample_Type, ncol = 2) +  # 按斑块位置 (AC vs PA) 进行分面
  
  # 核心修改：应用命名向量，彻底删除 labels
  scale_color_manual(
    values = fate_colors_named,
    labels = fate_labels_named,
    name = "Fate Branch"
  ) +
  
  labs(
    title = "APOBEC3A Expression Dynamics by AC/PA Status",
    x = "Pseudotime",
    y = "APOBEC3A Expression (log-normalized)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    strip.text = element_text(size = 14, face = "bold"), # 放大分面标题使其更醒目
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold") # 建议将 8 稍微调大为 10，方便阅读
  )

# 打印预览
print(APOBEC3A_AC_PA)


fate_colors_named <- c(
  "subset1" = "#1F77B4",  # 蓝线
  "subset2" = "#FF7F0E"   # 橙线
)

# 2. 定义图例显示名称的命名向量 (数据真名 = 图例上想显示的漂亮名字)
fate_labels_named <- c(
  "subset1" = "Cell Fate 1",
  "subset2" = "Cell Fate 2"
)

APOBEC3A_EXP_TIME <- ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = subset)) +
  geom_smooth(
    method = "gam",   
    formula = y ~ s(x, bs = "tp"), 
    se = FALSE, 
    linewidth = 1.5,
    alpha = 0.8
  ) +
  # 核心修改：使用命名向量，删除 labels
  scale_color_manual(
    values = fate_colors_named,
    labels = fate_labels_named,
    name = "Fate Branch"
  ) +
  labs(
    x = "Pseudotime", 
    y = "APOBEC3A Expression (log-normalized)", 
    title = ""
  ) + 
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),  
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),  
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90")
  )


# ==============================================================================
# 2. 绘制联合轨迹图：按斑块空间位置 (Sample_Type)
# ==============================================================================

# 定义空间位置的严格命名向量
sample_colors_named <- c(
  "Atherosclerotic Core" = "#FF7F0E", # 红色对应斑块核心
  "Proximal Adjacent" = "#1F77B4"     # 蓝色对应邻近区域
)

APOBEC3A_EXP_Sample_Type <- ggplot(combined_data, aes(x = pseudotime, y = APOBEC3A, color = Sample_Type)) +
  geom_smooth(
    method = "gam",   
    formula = y ~ s(x, bs = "tp"), 
    se = FALSE, 
    linewidth = 1.5,
    alpha = 0.8
  ) +
  # 核心修改：使用命名向量，删除 labels，不再混用 subset_colors
  scale_color_manual(
    values = sample_colors_named,
    name = "Plaque Location"
  ) +
  labs(
    x = "Pseudotime", 
    y = "APOBEC3A Expression (log-normalized)", 
    title = ""
  ) + 
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),  
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),  
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90")
  )

APOBEC3A_EXP <- APOBEC3A_EXP_Sample_Type | APOBEC3A_EXP_TIME
APOBEC3A_EXP
ggplot2::ggsave(APOBEC3A_EXP,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F3.1_APOBEC3A_EXP_TIME.png',width = 18,height = 8)
ggplot2::ggsave(APOBEC3A_EXP,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F3.1_APOBEC3A_EXP_TIME.pdf',width = 12,height = 5,
                device = cairo_pdf)
ggplot2::ggsave(APOBEC3A_AC_PA,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F3.1_APOBEC3A_AC_PA_TIME.png',width = 18,height = 8)
ggplot2::ggsave(APOBEC3A_AC_PA,filename = '/public3/DSC/single_cell/Result/figer_new/monocle3/F3.1_APOBEC3A_AC_PA_TIME.pdf',width = 15,height = 6,
                device = cairo_pdf)
