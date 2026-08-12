#!/usr/bin/env Rscript

# ============================================================
# Xenium Mo/Ma label transfer: publication-ready complete version
# 使用 sc_ref$Celltype_raw 进行注释转移
# 目标标签：Monocyte / Macrophage / LAM/Foam Cell
#
# 功能：
# 1. 只读取一次 sc_ref
# 2. 使用 sc_ref$Celltype_raw 作为 label transfer 注释列
# 3. 保留原始 xen FOV，不删除 FOV
# 4. 临时 image-free copy 用于 label transfer
# 5. transfer 结果写回 xen@meta.data
# 6. 不保存完整 xen 对象，只保存 metadata / summary / figures
# 7. 画所有 FOV 的 before/after 空间图，副标题显示 severity
# 8. 画不同 severity 进展时期的 Monocyte / Macrophage / LAM/Foam Cell 比例堆叠柱状图
# 9. 所有图片中统一将 LAM 改为 LAM/Foam Cell
# 10. 图片标题、子图标题、legend title 改为可发表风格
# 11. 新增坏死核心样区域识别：局部 LAM/Foam Cell 富集 + permutation + DBSCAN + concave hull
# 12. 新增细胞级 NC_like_region / NC_region_id 标记，并输出 FOV 最大 core 面积柱状图
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(tidyr)
  library(scales)
  library(gridExtra)
  library(grid)
  library(gtable)
  library(ggpubr)
})

# ============================================================
# 0. 基础参数检查
# ============================================================

env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

ROOT <- env_path("AST_ROOT", "/public3/DSC/single_cell/spatial/atherosclerosis_spatial_publication")
SPATIAL_ROOT <- env_path("AST_SPATIAL_ROOT", "/public3/DSC/single_cell/spatial")
SC_ROOT <- env_path("AST_SC_ROOT", "/public3/DSC/single_cell")
OUTPUT_ROOT <- env_path("AST_OUTPUT_ROOT", file.path(ROOT, "results"))
TABLE_DIR <- file.path(OUTPUT_ROOT, "tables")
FIG_DIR <- file.path(OUTPUT_ROOT, "figures")
DOC_DIR <- file.path(ROOT, "docs")
REF_DIR <- file.path(ROOT, "reference", "geomx")
OUT_PREFIX <- env_path("AST_OUT_PREFIX", "myeloid_story")

SAVE_FULL_METADATA <- env_flag("AST_SAVE_FULL_METADATA", FALSE)
SAVE_CELL_LEVEL_TABLES <- env_flag("AST_SAVE_CELL_LEVEL_TABLES", FALSE)
SAVE_TRANSFER_RDS <- env_flag("AST_SAVE_TRANSFER_RDS", FALSE)
RUN_EXPLORATORY_ISG_CORE_RELATION <- env_flag("AST_RUN_EXPLORATORY_ISG_CORE_RELATION", FALSE)

if (!exists("SC_ROOT")) {
  stop("没有找到 SC_ROOT，请先定义 SC_ROOT。")
}

if (!exists("SPATIAL_ROOT")) {
  stop("没有找到 SPATIAL_ROOT，请先定义 SPATIAL_ROOT。")
}

if (!exists("TABLE_DIR")) {
  TABLE_DIR <- file.path(SPATIAL_ROOT, "atherosclerosis_spatial_publication", "tables")
  warning("没有找到 TABLE_DIR，已自动设置为: ", TABLE_DIR)
}

if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(SPATIAL_ROOT, "atherosclerosis_spatial_publication", "figures")
  warning("没有找到 FIG_DIR，已自动设置为: ", FIG_DIR)
}

if (!exists("OUT_PREFIX")) {
  OUT_PREFIX <- "myeloid_story"
  warning("没有找到 OUT_PREFIX，已自动设置为: ", OUT_PREFIX)
}

if (!exists("TARGET_GENE")) {
  TARGET_GENE <- "APOBEC3A"
}

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0.1 统一标签设置
# ============================================================

LAM_DISPLAY_LABEL <- "LAM/Foam Cell"

rename_moma_labels <- function(x) {
  x <- as.character(x)
  x[x == "LAM"] <- LAM_DISPLAY_LABEL
  x
}

# ============================================================
# 0.2 单细胞部分一致的颜色字典
# ============================================================

cell_colors <- c(
  "Classical Mono"     = "#8c564b",
  "Inflammatory Mono"  = "#b15928",
  "ISG+ Mono"          = "#bcbd22",
  "Non-classical Mono" = "#c7c7c7",
  "Foam cells1"        = "#d62728",
  "Foam cells2"        = "#ff7f0e",
  "LAM"                = "#ff7f0e",
  "LAM/Foam Cell"      = "#ff7f0e",
  "Transitional Mac"   = "#2ca02c",
  "CX3CR1+ TRM"        = "#e377c2",
  "LYVE1+ TRM"         = "#f7b6d2",
  "TrMs"               = "#e377c2",
  "CM"                 = "#2ca02c",
  "Macrophage"         = "#9467bd",
  "Monocyte"           = "#2ca02c",
  "cDC1"               = "#1f77b4"
)

my_table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(
      fontfamily = "Arial",
      fontsize = 10,
      hjust = 0.5,
      x = 0.5
    ),
    bg_params = list(fill = c("white", "#f7f7f7"))
  ),
  colhead = list(
    bg_params = list(fill = "#404040"),
    fg_params = list(
      col = "white",
      fontface = "bold",
      fontfamily = "Arial",
      fontsize = 10,
      hjust = 0.5,
      x = 0.5
    )
  )
)

moma_levels <- c("Monocyte", "Macrophage", LAM_DISPLAY_LABEL)

severity_levels <- c(
  "Normal",
  "Mild",
  "Moderate",
  "Severe",
  "not annotated",
  "erosion plaque"
)


# ============================================================
# 0.3 坏死核心样区域识别参数（只改这里）
# ============================================================
# 调参逻辑：
# 1) 先看 *_threshold_diagnostics.csv：确认 local seed 是否存在。
# 2) 再看 *_DBSCAN_polygon_diagnostics.csv：定位 DBSCAN / polygon / cluster filter 哪一步筛没。
# 3) 最后再根据 *_cluster_summary.csv 的 polygon_area 分布调整 NC_CORE_AREA_CUTOFF。

NC_RUN_THRESHOLD_DIAGNOSTICS <- TRUE
NC_RUN_DBSCAN_POLYGON_DIAGNOSTICS <- TRUE
NC_RUN_FINAL_DETECTION <- TRUE
NC_PLOT_ALL_FOV_NC <- TRUE

# A3A / APOBEC3A 额外图：
# 1) 每个 FOV 的 APOBEC3A 空间表达热图
# 2) 每个 FOV 的 APOBEC3A 空间表达热图 + core-like boundary
# 3) 髓系细胞整体 APOBEC3A 表达量 / 阳性率 与疾病进展关系柱状图
A3A_RUN_SPATIAL_HEATMAPS <- TRUE
A3A_RUN_SPATIAL_HEATMAPS_WITH_CORE <- TRUE
A3A_RUN_DISEASE_BARPLOTS <- TRUE

# NULL = 画全部 FOV；也可以只画重点 FOV，例如：c("fov.8", "fov.9", "fov.10", "fov.11")
A3A_HEATMAP_IMAGE_NAMES <- NULL

# NULL = 跑全部 FOV；也可以只跑重点 FOV 加速调参，例如：c("fov.8", "fov.9", "fov.10", "fov.11")
NC_IMAGE_NAMES <- NULL

# 最终判断 FOV 是否存在 core-like region 的面积 cutoff。
# 注意：这个 cutoff 只影响 FOV-level positive/negative 和最终展示，不影响前面的诊断输出。
NC_RANDOM_SEED <- 123

# Step 1：局部 LAM/Foam Cell 富集 seed 参数
# 如果 pass_all_seed = 0：优先调这里。
NC_LOCAL_PARAMS <- list(
  radius = 120,
  n_perm = 200,
  min_total = 20,
  min_foam = 10,
  min_foam_fraction = 0.45,
  min_z = 2,
  max_p = 0.05
)

# Step 2：把 seed 连成不规则区域的 DBSCAN / polygon / cluster filter 参数
# 如果 pass_all_seed > 0 但 cluster_summary 为空：优先调这里。
NC_CLUSTER_PARAMS <- list(
  dbscan_eps = 120,              # 可试 80 / 100 / 120
  dbscan_minPts = 5,             # 可试 10 / 5 / 4
  min_cluster_cells = 12,
  min_cluster_foam = 8,
  min_cluster_foam_fraction = 0.45,
  concavity = 4
)

# polygon 生成后向外扩展的距离；0 = 不扩展。
# 注意：它只改变最终边界、polygon_area、NC_like_region，不改变 local seed / DBSCAN 结果。
NC_POLYGON_BUFFER <- 20
NC_CORE_AREA_CUTOFF <- 20000

# 备用宽松参数示例。正式分析时不要同时打开，复制到上面的 NC_LOCAL_PARAMS / NC_CLUSTER_PARAMS 中使用。
# NC_LOCAL_PARAMS <- list(radius = 80, n_perm = 200, min_total = 20, min_foam = 10, min_foam_fraction = 0.60, min_z = 2, max_p = 0.05)
# NC_CLUSTER_PARAMS <- list(dbscan_eps = 120, dbscan_minPts = 5, min_cluster_cells = 15, min_cluster_foam = 8, min_cluster_foam_fraction = 0.45, concavity = 3)

# ============================================================
# 1. 工具函数
# ============================================================

get_assay_data_safe <- function(object, assay = NULL, layer = "counts") {
  if (is.null(assay)) assay <- DefaultAssay(object)

  tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e1) {
      tryCatch(
        GetAssayData(object, assay = assay, slot = layer),
        error = function(e2) {
          stop(
            "无法从 assay=", assay, " 提取 layer/slot=", layer,
            "\nlayer error: ", e1$message,
            "\nslot error: ", e2$message
          )
        }
      )
    }
  )
}

set_default_assay_safely <- function(obj, preferred = c("Xenium", "RNA", "Spatial", "SCT")) {

  if (!inherits(obj, "Seurat")) {
    stop("输入对象不是 Seurat object。当前 class: ", paste(class(obj), collapse = ", "))
  }

  assay_names <- as.character(names(obj@assays))

  if (length(assay_names) == 0 || all(is.na(assay_names))) {
    stop("这个 Seurat 对象里没有可用 assay。")
  }

  message("Available assays: ", paste(assay_names, collapse = ", "))

  chosen <- intersect(preferred, assay_names)

  if (length(chosen) > 0) {
    DefaultAssay(obj) <- chosen[1]
  } else {
    DefaultAssay(obj) <- assay_names[1]
    warning(
      "没有找到 preferred assay: ",
      paste(preferred, collapse = ", "),
      "；已自动使用第一个 assay: ",
      assay_names[1]
    )
  }

  message("DefaultAssay set to: ", DefaultAssay(obj))
  obj
}

join_layers_safely <- function(obj, assay = NULL, object_name = "object") {

  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }

  assay_obj <- obj[[assay]]

  if (packageVersion("Seurat") >= "5.0.0" && inherits(assay_obj, "Assay5")) {

    obj[[assay]] <- tryCatch(
      JoinLayers(assay_obj),
      error = function(e) {
        message(object_name, " JoinLayers skipped: ", e$message)
        assay_obj
      }
    )

  } else {

    message(
      object_name,
      " JoinLayers skipped: assay class is ",
      paste(class(assay_obj), collapse = "/"),
      "; this is OK for Seurat v4-style Assay."
    )
  }

  return(obj)
}

strip_images_copy <- function(obj) {
  obj2 <- obj
  if ("images" %in% slotNames(obj2)) {
    obj2@images <- list()
  }
  obj2
}

save_plot_both <- function(p, filename_base, width = 8, height = 6) {

  pdf_file <- file.path(FIG_DIR, paste0(filename_base, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename_base, ".png"))

  ggsave(
    pdf_file,
    p,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    png_file,
    p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )

  message("Saved plot: ", pdf_file)
  message("Saved plot: ", png_file)

  invisible(c(pdf_file, png_file))
}

save_grob_both <- function(grob_obj, filename_base, width = 8, height = 6, dpi = 320) {

  pdf_file <- file.path(FIG_DIR, paste0(filename_base, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename_base, ".png"))

  grDevices::cairo_pdf(pdf_file, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(grob_obj)
  grDevices::dev.off()

  grDevices::png(
    png_file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    type = "cairo",
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(grob_obj)
  grDevices::dev.off()

  message("Saved plot: ", pdf_file)
  message("Saved plot: ", png_file)

  invisible(c(pdf_file, png_file))
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-\\.]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

choose_first_meta_col <- function(obj, candidates) {
  existing <- intersect(candidates, colnames(obj@meta.data))
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

choose_first_reduction <- function(obj, candidates = c(
  "umap", "UMAP", "umap.rpca", "umap.harmony", "umap.cca",
  "integrated.umap", "pca"
)) {
  reds <- names(obj@reductions)
  existing <- intersect(candidates, reds)
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

repair_fov_slots_safely <- function(obj, object_name = "xen") {

  if (!inherits(obj, "Seurat")) {
    stop(object_name, " is not a Seurat object.")
  }

  if (!"images" %in% slotNames(obj) || length(obj@images) == 0) {
    message(object_name, ": no images/FOV found.")
    return(obj)
  }

  for (img_name in names(obj@images)) {

    fov <- obj@images[[img_name]]
    fov_class <- class(fov)[1]

    message("Checking FOV: ", img_name, " / class: ", fov_class)

    fov_slots <- tryCatch(
      methods::getSlots(fov_class),
      error = function(e) NULL
    )

    if (is.null(fov_slots)) {
      message("Cannot inspect FOV slots for: ", img_name)
      next
    }

    if ("misc" %in% names(fov_slots)) {
      fov <- tryCatch(
        {
          methods::slot(fov, "misc") <- list()
          fov
        },
        error = function(e) {
          message("misc slot patch skipped for ", img_name, ": ", e$message)
          fov
        }
      )
    }

    if ("coords_x_orientation" %in% names(fov_slots)) {

      slot_class <- unname(fov_slots[["coords_x_orientation"]])

      default_value <- switch(
        slot_class,
        character = "right",
        logical = TRUE,
        numeric = 1,
        integer = 1L,
        list = list(),
        "right"
      )

      fov <- tryCatch(
        {
          methods::slot(fov, "coords_x_orientation") <- default_value
          fov
        },
        error = function(e) {
          message("coords_x_orientation slot patch skipped for ", img_name, ": ", e$message)
          fov
        }
      )
    }

    obj@images[[img_name]] <- fov
  }

  valid_after <- tryCatch(
    {
      methods::validObject(obj)
      TRUE
    },
    error = function(e) {
      message(object_name, " is still not fully valid after FOV patch: ", e$message)
      FALSE
    }
  )

  if (valid_after) {
    message(object_name, " FOV repaired successfully.")
  } else {
    warning(
      object_name,
      " FOV patch did not fully validate the object, but FOV was preserved. ",
      "Label transfer will use a temporary image-free copy."
    )
  }

  return(obj)
}

get_severity_col <- function(meta_df) {
  candidates <- c(
    "severity",
    "disease",
    "Disease",
    "grade",
    "Grade",
    "category",
    "Category",
    "condition",
    "Condition",
    "sample_type",
    "Sample_Type",
    "lesion_type",
    "plaque_type"
  )
  existing <- intersect(candidates, colnames(meta_df))
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

get_fov_severity_label <- function(df) {

  severity_col <- get_severity_col(df)

  if (is.na(severity_col)) {
    return("not available")
  }

  sev_vec <- as.character(df[[severity_col]])
  sev_vec <- sev_vec[!is.na(sev_vec) & sev_vec != ""]

  if (length(sev_vec) == 0) {
    return("not annotated")
  }

  sev_tab <- sort(table(sev_vec), decreasing = TRUE)

  if (length(sev_tab) == 1) {
    return(names(sev_tab)[1])
  } else {
    return(paste0(
      names(sev_tab),
      " n=",
      as.integer(sev_tab),
      collapse = "; "
    ))
  }
}

# ============================================================
# 2. 单细胞参考读取与整理
#    使用 sc_ref$Celltype_raw 作为转移标签
#    并将 LAM 重命名为 LAM/Foam Cell
# ============================================================

prepare_sc_ref_for_transfer <- function(
    sc_ref_path = file.path(SC_ROOT, "Result", "sub_integrated_data_Final.rds"),
    ref_label_source_col = "Celltype_raw"
) {

  if (!file.exists(sc_ref_path)) {
    stop("找不到单细胞参考对象: ", sc_ref_path)
  }

  message("Loading single-cell Mo/Ma reference: ", sc_ref_path)
  sc_ref <- readRDS(sc_ref_path)

  if (!inherits(sc_ref, "Seurat")) {
    stop("sc_ref 不是 Seurat object。当前 class: ", paste(class(sc_ref), collapse = ", "))
  }

  if (!ref_label_source_col %in% colnames(sc_ref@meta.data)) {
    stop(
      "sc_ref@meta.data 中没有 ", ref_label_source_col,
      "，不能进行 Mo/Ma 细胞注释映射。\n当前可用列包括:\n",
      paste(colnames(sc_ref@meta.data), collapse = "\n")
    )
  }

  sc_ref <- set_default_assay_safely(
    sc_ref,
    preferred = c("RNA", "SCT", "integrated")
  )

  sc_ref <- join_layers_safely(
    sc_ref,
    assay = DefaultAssay(sc_ref),
    object_name = "sc_ref"
  )

  sc_ref[[ref_label_source_col]] <- as.character(sc_ref@meta.data[[ref_label_source_col]])

  raw_labels <- as.character(sc_ref@meta.data[[ref_label_source_col]])
  raw_labels <- rename_moma_labels(raw_labels)

  sc_ref$Celltype_transfer <- raw_labels

  # 只保留 Monocyte / Macrophage / LAM/Foam Cell
  sc_ref$Celltype_transfer[!sc_ref$Celltype_transfer %in% moma_levels] <- NA_character_

  keep_ref_cells <- rownames(sc_ref@meta.data)[!is.na(sc_ref$Celltype_transfer)]

  if (length(keep_ref_cells) < 20) {
    stop(
      "使用 ", ref_label_source_col, " 后可用于 transfer 的参考细胞少于 20 个。\n",
      "请检查 unique(sc_ref$", ref_label_source_col, ")。"
    )
  }

  sc_ref <- subset(sc_ref, cells = keep_ref_cells)

  sc_ref$Celltype_transfer <- factor(
    sc_ref$Celltype_transfer,
    levels = moma_levels
  )

  message("Reference source label column: ", ref_label_source_col)
  message("Reference labels for transfer:")
  print(table(sc_ref$Celltype_transfer, useNA = "ifany"))

  return(sc_ref)
}

# ============================================================
# 3. Label transfer：保留 FOV，只用临时 image-free copy 计算
# ============================================================

transfer_moma_labels_to_spatial_keep_fov <- function(
    ref_obj,
    query_obj,
    ref_label_col = "Celltype_transfer",
    query_assay = NULL,
    out_prefix = "Xenium_MoMa",
    dims_use = 1:30,
    restrict_query_cells = NULL,
    exclude_anchor_genes = c("APOBEC3A_B", "CD45")
) {

  message("========== Transfer Mo/Ma labels to ", out_prefix, " while preserving FOV ==========")

  if (is.null(query_assay)) {
    query_assay <- DefaultAssay(query_obj)
  }

  query_transfer <- strip_images_copy(query_obj)
  DefaultAssay(query_transfer) <- query_assay

  query_transfer <- join_layers_safely(
    query_transfer,
    assay = query_assay,
    object_name = paste0(out_prefix, "_temporary_query")
  )

  if (!is.null(restrict_query_cells)) {
    keep_cells <- intersect(restrict_query_cells, Cells(query_transfer))
    query_sub <- subset(query_transfer, cells = keep_cells)
  } else {
    keep_cells <- Cells(query_transfer)
    query_sub <- query_transfer
  }

  if (ncol(query_sub) < 20) {
    warning(out_prefix, ": query_sub cells too few; skipped label transfer.")
    return(list(
      query_obj = query_obj,
      pred = data.frame()
    ))
  }

  ref_assay <- DefaultAssay(ref_obj)

  common_genes <- intersect(rownames(ref_obj), rownames(query_sub))
  common_genes <- setdiff(common_genes, exclude_anchor_genes)

  message(out_prefix, ": common genes = ", length(common_genes))

  if (length(common_genes) < 100) {
    warning(out_prefix, ": common genes < 100, transfer may be unstable.")
  }

  ref_obj <- NormalizeData(ref_obj, verbose = FALSE)
  ref_obj <- FindVariableFeatures(ref_obj, nfeatures = 3000, verbose = FALSE)

  transfer_features <- intersect(VariableFeatures(ref_obj), common_genes)
  transfer_features <- intersect(transfer_features, rownames(query_sub))

  if (length(transfer_features) < 100) {
    transfer_features <- common_genes
  }

  message(out_prefix, ": transfer features = ", length(transfer_features))

  ref_obj <- ScaleData(ref_obj, features = transfer_features, verbose = FALSE)
  ref_obj <- RunPCA(
    ref_obj,
    features = transfer_features,
    npcs = max(dims_use),
    verbose = FALSE
  )

  query_sub <- NormalizeData(query_sub, verbose = FALSE)
  query_sub <- FindVariableFeatures(query_sub, nfeatures = 3000, verbose = FALSE)

  anchors <- FindTransferAnchors(
    reference = ref_obj,
    query = query_sub,
    reference.assay = ref_assay,
    query.assay = query_assay,
    normalization.method = "LogNormalize",
    features = transfer_features,
    dims = dims_use,
    reduction = "pcaproject",
    reference.reduction = "pca",
    verbose = TRUE
  )

  pred <- TransferData(
    anchorset = anchors,
    refdata = ref_obj[[ref_label_col, drop = TRUE]],
    dims = dims_use,
    prediction.assay = FALSE,
    verbose = TRUE
  )

  pred_col <- paste0(out_prefix, "_scPred_Celltype")
  score_col <- paste0(out_prefix, "_scPred_score")

  query_obj@meta.data[[pred_col]] <- NA_character_
  query_obj@meta.data[[score_col]] <- NA_real_

  query_obj@meta.data[rownames(pred), pred_col] <- as.character(pred$predicted.id)
  query_obj@meta.data[rownames(pred), score_col] <- as.numeric(pred$prediction.score.max)

  score_cols <- grep("^prediction.score.", colnames(pred), value = TRUE)

  for (cc in score_cols) {
    new_cc <- paste0(out_prefix, "_", cc)
    query_obj@meta.data[[new_cc]] <- NA_real_
    query_obj@meta.data[rownames(pred), new_cc] <- as.numeric(pred[[cc]])
  }

  message(out_prefix, ": label transfer done and written back to original object with FOV.")

  return(list(
    query_obj = query_obj,
    pred = pred
  ))
}

# ============================================================
# 4. UMAP 绘图函数：publication-ready
#    左边为完整细胞注释，右边为髓系细胞注释
# ============================================================

plot_before_after_dim <- function(
    obj,
    before_col,
    after_col,
    title_prefix = "Xenium",
    cells = NULL,
    filename_base = NULL,
    width = 13,
    height = 6
) {

  plot_obj <- strip_images_copy(obj)

  reduction_use <- choose_first_reduction(plot_obj)

  if (is.na(reduction_use)) {
    warning("没有找到可用 reduction，跳过 DimPlot。")
    return(NULL)
  }

  if (!is.null(cells)) {
    cells <- intersect(cells, Cells(plot_obj))
    if (length(cells) < 5) {
      warning("指定 cells 少于 5 个，跳过绘图。")
      return(NULL)
    }
    plot_obj <- subset(plot_obj, cells = cells)
  }

  if (!before_col %in% colnames(plot_obj@meta.data)) {
    warning("没有找到 before_col: ", before_col)
    return(NULL)
  }

  if (!after_col %in% colnames(plot_obj@meta.data)) {
    warning("没有找到 after_col: ", after_col)
    return(NULL)
  }

  plot_obj@meta.data[[after_col]] <- factor(
    as.character(plot_obj@meta.data[[after_col]]),
    levels = moma_levels
  )

  p_before <- DimPlot(
    plot_obj,
    reduction = reduction_use,
    group.by = before_col,
    raster = TRUE,
    shuffle = TRUE,
    pt.size = 0.25
  ) +
    labs(color = "Cell annotation") +
    ggtitle("Comprehensive cell annotation") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )

  p_after <- DimPlot(
    plot_obj,
    reduction = reduction_use,
    group.by = after_col,
    raster = TRUE,
    shuffle = TRUE,
    pt.size = 0.25,
    cols = cell_colors[moma_levels]
  ) +
    labs(color = "Myeloid annotation") +
    ggtitle("Myeloid cell annotation") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )

  p <- (p_before | p_after) +
    patchwork::plot_annotation(
      title = title_prefix,
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
      )
    )

  print(p)

  if (!is.null(filename_base)) {
    save_plot_both(p, filename_base, width = width, height = height)
  }

  return(p)
}

plot_transfer_score_dim <- function(
    obj,
    score_col = "Xenium_MoMa_scPred_score",
    title = "Myeloid annotation transfer confidence",
    filename_base = NULL,
    width = 7,
    height = 6
) {

  plot_obj <- strip_images_copy(obj)

  reduction_use <- choose_first_reduction(plot_obj)

  if (is.na(reduction_use)) {
    warning("没有找到可用 reduction，跳过 confidence DimPlot。")
    return(NULL)
  }

  if (!score_col %in% colnames(plot_obj@meta.data)) {
    warning("没有找到 score_col: ", score_col)
    return(NULL)
  }

  p <- FeaturePlot(
    plot_obj,
    features = score_col,
    reduction = reduction_use,
    raster = TRUE,
    pt.size = 0.25
  ) +
    ggtitle(title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  print(p)

  if (!is.null(filename_base)) {
    save_plot_both(p, filename_base, width = width, height = height)
  }

  return(p)
}

# ============================================================
# 5. 所有 FOV 空间图：publication-ready
#    左图：完整细胞注释
#    右图：髓系细胞注释
#    标题和 legend 改为适合发表
# ============================================================

get_one_fov_df <- function(
    obj,
    image_use,
    before_col = "predicted.id",
    after_col = "Xenium_MoMa_scPred_Celltype"
) {

  coord <- GetTissueCoordinates(obj, image = image_use)
  coord <- as.data.frame(coord)

  if (!"cell" %in% colnames(coord)) {
    stop(
      "FOV ", image_use, " 的 GetTissueCoordinates() 结果没有 cell 列。当前列名: ",
      paste(colnames(coord), collapse = ", ")
    )
  }

  if (!all(c("x", "y") %in% colnames(coord))) {
    stop(
      "FOV ", image_use, " 的坐标表没有 x/y 列。当前列名: ",
      paste(colnames(coord), collapse = ", ")
    )
  }

  coord$cell <- as.character(coord$cell)

  meta <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  meta$cell <- as.character(meta$cell)

  n_match <- length(intersect(coord$cell, meta$cell))

  if (n_match == 0) {
    stop(
      "FOV ", image_use, " 的 coord$cell 与 metadata rownames 没有匹配。",
      "\ncoord$cell head: ", paste(head(coord$cell), collapse = ", "),
      "\nmetadata head: ", paste(head(meta$cell), collapse = ", ")
    )
  }

  df <- coord %>%
    left_join(meta, by = "cell")

  df
}

summarise_one_fov_simple <- function(
    df,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score"
) {

  after_tab <- table(as.character(df[[after_col]]), useNA = "no")

  subtype_counts <- sapply(moma_levels, function(x) {
    if (x %in% names(after_tab)) as.integer(after_tab[[x]]) else 0L
  })

  n_cells <- nrow(df)
  n_transferred <- sum(!is.na(df[[after_col]]))

  n_high_conf <- if (score_col %in% colnames(df)) {
    sum(!is.na(df[[after_col]]) & as.numeric(df[[score_col]]) >= 0.5, na.rm = TRUE)
  } else {
    NA_integer_
  }

  n_a3a_transferred <- if ("APOBEC3A_detected_transfer" %in% colnames(df)) {
    sum(df$APOBEC3A_detected_transfer & !is.na(df[[after_col]]), na.rm = TRUE)
  } else {
    NA_integer_
  }

  severity_label <- get_fov_severity_label(df)

  tibble(
    FOV = image_use,
    severity = severity_label,
    n_cells_in_FOV = n_cells,
    n_transferred = n_transferred,
    n_high_confidence = n_high_conf,
    n_APOBEC3A_positive_transferred = n_a3a_transferred,
    n_Monocyte = as.integer(subtype_counts[["Monocyte"]]),
    n_Macrophage = as.integer(subtype_counts[["Macrophage"]]),
    n_LAM_Foam_Cell = as.integer(subtype_counts[[LAM_DISPLAY_LABEL]])
  )
}

plot_one_xenium_fov_with_severity_subtitle <- function(
    obj,
    image_use,
    before_col = "predicted.id",
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    only_transferred_after = TRUE,
    max_points = 200000,
    save_plot = TRUE,
    width = 14,
    height = 6
) {

  message("\n========== Plotting FOV: ", image_use, " ==========")

  df <- get_one_fov_df(
    obj = obj,
    image_use = image_use,
    before_col = before_col,
    after_col = after_col
  )

  severity_label <- get_fov_severity_label(df)

  if (nrow(df) > max_points) {
    set.seed(123)
    df <- df[sample(seq_len(nrow(df)), max_points), , drop = FALSE]
  }

  df[[after_col]] <- factor(
    as.character(df[[after_col]]),
    levels = moma_levels
  )

  p_before <- ggplot(
    df,
    aes(x = x, y = y, color = .data[[before_col]])
  ) +
    geom_point(size = 0.10, alpha = 0.80) +
    coord_fixed() +
    theme_void(base_size = 12) +
    labs(
      color = "Cell annotation",
      title = "Comprehensive cell annotation"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  if (only_transferred_after) {

    df_after <- df %>%
      filter(!is.na(.data[[after_col]]))

    after_plot_col <- after_col

  } else {

    df_after <- df %>%
      mutate(
        after_plot_label = ifelse(
          is.na(.data[[after_col]]),
          "Not transferred",
          as.character(.data[[after_col]])
        )
      )

    df_after$after_plot_label <- factor(
      df_after$after_plot_label,
      levels = c(moma_levels, "Not transferred")
    )

    after_plot_col <- "after_plot_label"
  }

  if (nrow(df_after) == 0) {

    p_after <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "No myeloid cells were assigned in this section",
        size = 5
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void() +
      ggtitle("Myeloid cell annotation") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

  } else {

    p_after <- ggplot(
      df_after,
      aes(x = x, y = y, color = .data[[after_plot_col]])
    ) +
      geom_point(size = 0.13, alpha = 0.90) +
      coord_fixed() +
      theme_void(base_size = 12) +
      scale_color_manual(
        values = cell_colors[moma_levels],
        breaks = moma_levels,
        drop = FALSE
      ) +
      labs(
        color = "Myeloid annotation",
        title = "Myeloid cell annotation"
      ) +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )
  }

  p_all <- (p_before | p_after) +
    patchwork::plot_annotation(
      title = "Xenium spatial annotation",
      subtitle = paste0("Section: ", image_use, " | Severity: ", severity_label),
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )

  print(p_all)

  if (save_plot) {

    safe_img <- sanitize_filename(image_use)

    pdf_file <- file.path(
      FIG_DIR,
      paste0(OUT_PREFIX, "_", safe_img, "_FOV_publication_annotation.pdf")
    )

    png_file <- file.path(
      FIG_DIR,
      paste0(OUT_PREFIX, "_", safe_img, "_FOV_publication_annotation.png")
    )

    ggsave(
      pdf_file,
      p_all,
      width = width,
      height = height,
      device = cairo_pdf,
      bg = "white"
    )

    ggsave(
      png_file,
      p_all,
      width = width,
      height = height,
      dpi = 320,
      bg = "white"
    )

    message("Saved FOV plot PDF: ", pdf_file)
    message("Saved FOV plot PNG: ", png_file)
  }

  list(
    plot = p_all,
    summary_row = summarise_one_fov_simple(
      df = df,
      image_use = image_use,
      after_col = after_col,
      score_col = score_col
    )
  )
}

plot_all_xenium_fovs_with_severity_subtitle <- function(
    obj,
    before_col = "predicted.id",
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    only_transferred_after = TRUE,
    max_points = 200000
) {

  if (length(obj@images) == 0) {
    stop("xen@images 为空，没有 FOV 可以画。")
  }

  image_names <- names(obj@images)

  message("Total FOVs to plot: ", length(image_names))
  print(image_names)

  all_rows <- list()

  for (image_use in image_names) {

    res <- tryCatch(
      {
        plot_one_xenium_fov_with_severity_subtitle(
          obj = obj,
          image_use = image_use,
          before_col = before_col,
          after_col = after_col,
          score_col = score_col,
          only_transferred_after = only_transferred_after,
          max_points = max_points,
          save_plot = TRUE,
          width = 14,
          height = 6
        )
      },
      error = function(e) {
        warning("FOV ", image_use, " plotting failed: ", e$message)
        NULL
      }
    )

    if (!is.null(res)) {
      all_rows[[image_use]] <- res$summary_row
    }
  }

  fov_summary <- bind_rows(all_rows)

  fov_summary_out <- file.path(
    TABLE_DIR,
    paste0(OUT_PREFIX, "_all_FOV_publication_annotation_summary.csv")
  )

  write.csv(
    fov_summary,
    fov_summary_out,
    row.names = FALSE
  )

  message("Saved all FOV summary: ", fov_summary_out)

  fov_summary
}

# ============================================================
# 6. 不同进展时期髓系细胞比例堆叠柱状图
#    修正版：不用 count()，避免冲突
#    并统一使用 LAM/Foam Cell
# ============================================================

plot_moma_composition_by_severity <- function(
    obj,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    use_high_confidence_only = FALSE,
    confidence_cutoff = 0.5,
    filename_base = paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity"),
    width = 9,
    height = 6.5
) {

  meta <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  severity_col <- get_severity_col(meta)

  if (is.na(severity_col)) {
    stop(
      "xen@meta.data 中没有 severity/disease/grade/category 等进展时期列。\n当前列名包括:\n",
      paste(colnames(meta), collapse = "\n")
    )
  }

  if (!after_col %in% colnames(meta)) {
    stop(
      "metadata 中没有 after_col: ", after_col,
      "\n当前列名包括:\n",
      paste(colnames(meta), collapse = "\n")
    )
  }

  if (use_high_confidence_only && !score_col %in% colnames(meta)) {
    stop(
      "use_high_confidence_only = TRUE，但 metadata 中没有 score_col: ",
      score_col
    )
  }

  plot_meta <- meta
  plot_meta$Severity <- as.character(plot_meta[[severity_col]])
  plot_meta$CellType <- as.character(plot_meta[[after_col]])

  if (score_col %in% colnames(plot_meta)) {
    plot_meta$TransferScore <- as.numeric(plot_meta[[score_col]])
  } else {
    plot_meta$TransferScore <- NA_real_
  }

  plot_meta <- plot_meta %>%
    dplyr::filter(
      !is.na(Severity),
      Severity != "",
      !is.na(CellType),
      CellType %in% moma_levels
    )

  if (use_high_confidence_only) {
    plot_meta <- plot_meta %>%
      dplyr::filter(
        !is.na(TransferScore),
        TransferScore >= confidence_cutoff
      )
  }

  if (nrow(plot_meta) == 0) {
    stop("没有可用于绘制堆叠柱状图的 transferred myeloid cells。")
  }

  observed_severity <- unique(plot_meta$Severity)

  severity_order <- c(
    severity_levels[severity_levels %in% observed_severity],
    setdiff(observed_severity, severity_levels)
  )

  plot_meta$Status <- factor(
    plot_meta$Severity,
    levels = severity_order
  )

  plot_meta$CellType <- factor(
    plot_meta$CellType,
    levels = moma_levels
  )

  table_counts <- base::table(
    plot_meta$Status,
    plot_meta$CellType,
    useNA = "no"
  )

  Df_comp <- as.data.frame(table_counts, stringsAsFactors = FALSE)
  colnames(Df_comp) <- c("Status", "CellType", "Count")

  Df_comp$Status <- factor(
    as.character(Df_comp$Status),
    levels = severity_order
  )

  Df_comp$CellType <- factor(
    as.character(Df_comp$CellType),
    levels = moma_levels
  )

  Df_comp$Count <- as.numeric(Df_comp$Count)

  Df_comp <- Df_comp %>%
    dplyr::group_by(Status) %>%
    dplyr::mutate(
      Total = sum(Count),
      Proportion = ifelse(Total > 0, Count / Total, 0)
    ) %>%
    dplyr::ungroup()

  composition_out <- file.path(
    TABLE_DIR,
    paste0(filename_base, "_counts_and_proportions.csv")
  )

  write.csv(
    Df_comp,
    composition_out,
    row.names = FALSE
  )

  message("Saved composition table: ", composition_out)

  message("========== Myeloid composition by severity: counts ==========")
  print(
    Df_comp %>%
      dplyr::select(Status, CellType, Count) %>%
      tidyr::pivot_wider(
        names_from = CellType,
        values_from = Count,
        values_fill = 0
      ) %>%
      as.data.frame()
  )

  message("========== Myeloid composition by severity: proportions ==========")
  print(
    Df_comp %>%
      dplyr::mutate(
        Proportion_percent = paste0(round(Proportion * 100, 1), "%")
      ) %>%
      dplyr::select(Status, CellType, Proportion_percent) %>%
      tidyr::pivot_wider(
        names_from = CellType,
        values_from = Proportion_percent,
        values_fill = "0%"
      ) %>%
      as.data.frame()
  )

  p_stack <- ggplot(
    Df_comp,
    aes(x = Status, y = Proportion, fill = CellType)
  ) +
    geom_bar(
      stat = "identity",
      position = "fill",
      width = 0.6
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = cell_colors[moma_levels],
      breaks = moma_levels,
      drop = FALSE
    ) +
    labs(
      x = "Progression stage",
      y = "Proportion",
      fill = "Myeloid annotation"
    ) +
    ggtitle("Myeloid cell composition across lesion progression stages") +
    theme_minimal(base_family = "Arial", base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(t = 20, r = 20, b = 10, l = 20, unit = "pt"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  table_data <- Df_comp %>%
    dplyr::group_by(Status, CellType) %>%
    dplyr::summarise(
      Proportion = paste0(round(mean(Proportion, na.rm = TRUE) * 100, 1), "%"),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = CellType,
      values_from = Proportion,
      values_fill = "0%"
    ) %>%
    dplyr::arrange(Status)

  table_data <- table_data %>%
    dplyr::select(
      Status,
      dplyr::any_of(moma_levels)
    )

  n_cols <- ncol(table_data)

  if (n_cols <= 6) {

    table_grob <- gridExtra::tableGrob(
      table_data,
      rows = NULL,
      theme = my_table_theme
    )

  } else {

    split_idx <- ceiling((n_cols - 1) / 2) + 1

    grob_1 <- gridExtra::tableGrob(
      table_data[, 1:split_idx, drop = FALSE],
      rows = NULL,
      theme = my_table_theme
    )

    grob_2 <- gridExtra::tableGrob(
      table_data[, c(1, (split_idx + 1):n_cols), drop = FALSE],
      rows = NULL,
      theme = my_table_theme
    )

    table_grob <- gridExtra::arrangeGrob(
      grob_1,
      grob_2,
      nrow = 2
    )
  }

  final_plot <- gridExtra::arrangeGrob(
    p_stack,
    table_grob,
    nrow = 2,
    heights = c(0.7, 0.3)
  )

  grid::grid.newpage()
  grid::grid.draw(final_plot)

  save_grob_both(
    final_plot,
    filename_base,
    width = width,
    height = height,
    dpi = 320
  )

  return(list(
    plot = p_stack,
    final_plot = final_plot,
    composition_table = Df_comp,
    table_data = table_data
  ))
}


# ============================================================
# 6.1 Necrotic-core-like region detection utilities
#     局部 LAM/Foam Cell 富集 + permutation + DBSCAN + concave hull
#     注意：已有图片文件名不改；本节只新增坏死核心样区域相关输出
# ============================================================

check_nc_required_packages <- function() {
  required_pkgs <- c("dbscan", "sf")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "缺少坏死核心样区域识别所需 R 包: ", paste(missing_pkgs, collapse = ", "),
      "\n请先安装，例如 install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
    )
  }
  if (!requireNamespace("concaveman", quietly = TRUE)) {
    warning(
      "没有安装 concaveman 包；将使用 convex hull 作为备用边界。",
      "建议安装：install.packages('concaveman')"
    )
  }
  invisible(TRUE)
}

mix_color_with_grey <- function(col, grey = "#CFCFCF", grey_weight = 0.68) {
  col_rgb <- grDevices::col2rgb(col) / 255
  grey_rgb <- grDevices::col2rgb(grey) / 255
  out <- (1 - grey_weight) * col_rgb + grey_weight * grey_rgb
  grDevices::rgb(out[1, ], out[2, ], out[3, ])
}

make_severity_palette <- function(severity_values) {
  severity_values <- unique(as.character(severity_values))
  severity_values <- severity_values[!is.na(severity_values) & severity_values != ""]
  if (length(severity_values) == 0) severity_values <- "not annotated"

  base_palette <- c(
    "Normal" = "#4DAF4A",
    "Mild" = "#377EB8",
    "Moderate" = "#984EA3",
    "Severe" = "#E41A1C",
    "not annotated" = "#8C8C8C",
    "erosion plaque" = "#FF7F00"
  )

  missing_levels <- setdiff(severity_values, names(base_palette))
  if (length(missing_levels) > 0) {
    extra_cols <- scales::hue_pal()(length(missing_levels))
    names(extra_cols) <- missing_levels
    base_palette <- c(base_palette, extra_cols)
  }

  base_palette[severity_values]
}

make_nc_area_fill_palette <- function(severity_values) {
  sev_cols <- make_severity_palette(severity_values)
  bright_cols <- sev_cols
  names(bright_cols) <- paste0(names(sev_cols), "__above")
  grey_cols <- vapply(sev_cols, mix_color_with_grey, character(1))
  names(grey_cols) <- paste0(names(sev_cols), "__below")
  c(bright_cols, grey_cols)
}

nc_empty_cluster_summary <- function(image_use = character()) {
  tibble(
    FOV = character(),
    severity = character(),
    core_id = character(),
    area_pass_cutoff = logical(),
    n_cells_in_polygon = integer(),
    n_myeloid_in_polygon = integer(),
    n_foam_in_polygon = integer(),
    foam_fraction_in_polygon = numeric(),
    polygon_area = numeric(),
    centroid_x = numeric(),
    centroid_y = numeric(),
    mean_local_foam_fraction_seed = numeric(),
    mean_local_foam_z_seed = numeric(),
    min_local_foam_p_seed = numeric()
  )
}

nc_empty_fov_summary <- function(
    image_use,
    severity_label = "not annotated",
    n_total_cells = 0L,
    n_myeloid_cells = 0L,
    area_cutoff = NA_real_
) {
  tibble(
    FOV = image_use,
    severity = severity_label,
    NC_like_positive = FALSE,
    NC_area_cutoff = area_cutoff,
    n_candidate_core_regions = 0L,
    n_core_regions = 0L,
    max_core_area = 0,
    total_core_area = 0,
    n_total_cells = n_total_cells,
    n_myeloid_cells = n_myeloid_cells,
    n_core_like_cells = 0L,
    n_myeloid_core_like_cells = 0L,
    n_foam_core_like_cells = 0L
  )
}

polygon_sf_to_df <- function(poly_sf) {
  if (is.null(poly_sf) || nrow(poly_sf) == 0) {
    return(data.frame())
  }

  out_list <- lapply(seq_len(nrow(poly_sf)), function(i) {
    cc <- as.data.frame(sf::st_coordinates(poly_sf[i, ]))
    if (nrow(cc) == 0) return(NULL)

    cc$core_id <- poly_sf$core_id[i]

    grp_cols <- intersect(c("L1", "L2", "L3"), colnames(cc))
    if (length(grp_cols) == 0) {
      cc$group_path <- poly_sf$core_id[i]
    } else {
      grp_val <- apply(cc[, grp_cols, drop = FALSE], 1, paste, collapse = "_")
      cc$group_path <- paste0(poly_sf$core_id[i], "__", grp_val)
    }
    cc
  })

  bind_rows(out_list)
}

points_to_core_polygon <- function(df_pts, core_id, concavity = 2) {
  if (nrow(df_pts) < 3) return(NULL)

  sf_pts <- sf::st_as_sf(df_pts, coords = c("x", "y"), crs = NA)

  poly <- tryCatch(
    {
      if (requireNamespace("concaveman", quietly = TRUE)) {
        concaveman::concaveman(sf_pts, concavity = concavity, length_threshold = 0)
      } else {
        sf::st_sf(
          core_id = core_id,
          geometry = sf::st_sfc(sf::st_convex_hull(sf::st_union(sf_pts)))
        )
      }
    },
    error = function(e) {
      message("concave hull failed for ", core_id, "; fallback to convex hull. Error: ", e$message)
      sf::st_sf(
        core_id = core_id,
        geometry = sf::st_sfc(sf::st_convex_hull(sf::st_union(sf_pts)))
      )
    }
  )

  poly <- sf::st_make_valid(poly)
  poly$core_id <- core_id
  poly <- poly[, "core_id", drop = FALSE]

  geom_type <- as.character(sf::st_geometry_type(poly, by_geometry = FALSE))
  if (!grepl("POLYGON", geom_type)) {
    return(NULL)
  }

  return(poly)
}

get_one_fov_df_for_nc <- function(
    obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    a3a_count_col = "APOBEC3A_count_transfer",
    a3a_detect_col = "APOBEC3A_detected_transfer"
) {
  df <- get_one_fov_df(
    obj = obj,
    image_use = image_use,
    before_col = choose_first_meta_col(obj, c("predicted.id", "celltype", "Celltype", "seurat_clusters")),
    after_col = after_col
  )

  df$FOV <- image_use
  df$CellType_transfer <- as.character(df[[after_col]])
  df$is_myeloid_transfer <- !is.na(df$CellType_transfer) & df$CellType_transfer %in% moma_levels
  df$is_foam_transfer <- !is.na(df$CellType_transfer) & df$CellType_transfer == LAM_DISPLAY_LABEL

  if (score_col %in% colnames(df)) {
    df$TransferScore <- as.numeric(df[[score_col]])
  } else {
    df$TransferScore <- NA_real_
  }

  if (a3a_count_col %in% colnames(df)) {
    df$APOBEC3A_count_plot <- as.numeric(df[[a3a_count_col]])
  } else {
    df$APOBEC3A_count_plot <- NA_real_
  }

  if (a3a_detect_col %in% colnames(df)) {
    df$APOBEC3A_detect_plot <- as.logical(df[[a3a_detect_col]])
  } else {
    df$APOBEC3A_detect_plot <- NA
  }

  df
}

# ------------------------------------------------------------
# 单个 FOV：识别坏死核心样区域
# 逻辑：
# 1. 在 transferred myeloid cells 内计算局部 LAM/Foam Cell 富集
# 2. permutation 判断局部富集是否超过随机分布
# 3. DBSCAN 把候选 seed 连成多个空间 cluster
# 4. concave hull 重建不规则区域
# 5. 用 min_core_area_cutoff 作为 FOV 是否存在坏死核心样区域的判定条件
# ------------------------------------------------------------
detect_necrotic_core_like_one_fov <- function(
    obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    radius = 80,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.60,
    min_z = 2,
    max_p = 0.05,
    dbscan_eps = NULL,
    dbscan_minPts = 10,
    min_cluster_cells = 25,
    min_cluster_foam = 15,
    min_cluster_foam_fraction = 0.60,
    min_core_area_cutoff = 20000,
    concavity = 2,
    polygon_buffer = 0,
    seed = 123
) {

  check_nc_required_packages()
  if (is.null(dbscan_eps)) dbscan_eps <- radius

  df <- get_one_fov_df_for_nc(
    obj = obj,
    image_use = image_use,
    after_col = after_col,
    score_col = score_col
  )

  severity_label <- get_fov_severity_label(df)

  df$NC_local_total <- NA_integer_
  df$NC_local_foam_count <- NA_integer_
  df$NC_local_foam_fraction <- NA_real_
  df$NC_local_foam_z <- NA_real_
  df$NC_local_foam_p <- NA_real_
  df$NC_candidate_seed <- FALSE
  df$NC_region_id <- NA_character_
  df$NC_like_region <- FALSE

  myeloid_df <- df %>% dplyr::filter(is_myeloid_transfer)

  if (nrow(myeloid_df) < min_total) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = data.frame(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = data.frame(),
      cluster_summary = nc_empty_cluster_summary(image_use),
      fov_summary = nc_empty_fov_summary(
        image_use = image_use,
        severity_label = severity_label,
        n_total_cells = nrow(df),
        n_myeloid_cells = nrow(myeloid_df),
        area_cutoff = min_core_area_cutoff
      )
    ))
  }

  coords <- as.matrix(myeloid_df[, c("x", "y")])
  nn <- dbscan::frNN(coords, eps = radius)

  obs_total <- lengths(nn$id)
  obs_foam <- sapply(nn$id, function(idx) sum(myeloid_df$is_foam_transfer[idx], na.rm = TRUE))
  obs_frac <- ifelse(obs_total > 0, obs_foam / obs_total, NA_real_)

  set.seed(seed)
  perm_mat <- replicate(n_perm, {
    shuffled <- sample(myeloid_df$is_foam_transfer)
    sapply(nn$id, function(idx) sum(shuffled[idx], na.rm = TRUE))
  })

  if (is.null(dim(perm_mat))) {
    perm_mat <- matrix(perm_mat, ncol = 1)
  }

  perm_mean <- rowMeans(perm_mat, na.rm = TRUE)
  perm_sd <- apply(perm_mat, 1, sd, na.rm = TRUE)
  perm_sd[is.na(perm_sd) | perm_sd == 0] <- 1e-8

  obs_z <- (obs_foam - perm_mean) / perm_sd
  obs_p <- (rowSums(perm_mat >= obs_foam) + 1) / (n_perm + 1)

  seed_flag <- (
    obs_total >= min_total &
      obs_foam >= min_foam &
      obs_frac >= min_foam_fraction &
      obs_z >= min_z &
      obs_p <= max_p
  )

  myeloid_df$NC_local_total <- as.integer(obs_total)
  myeloid_df$NC_local_foam_count <- as.integer(obs_foam)
  myeloid_df$NC_local_foam_fraction <- as.numeric(obs_frac)
  myeloid_df$NC_local_foam_z <- as.numeric(obs_z)
  myeloid_df$NC_local_foam_p <- as.numeric(obs_p)
  myeloid_df$NC_candidate_seed <- seed_flag

  idx_match <- match(myeloid_df$cell, df$cell)
  df[idx_match, c(
    "NC_local_total",
    "NC_local_foam_count",
    "NC_local_foam_fraction",
    "NC_local_foam_z",
    "NC_local_foam_p",
    "NC_candidate_seed"
  )] <- myeloid_df[, c(
    "NC_local_total",
    "NC_local_foam_count",
    "NC_local_foam_fraction",
    "NC_local_foam_z",
    "NC_local_foam_p",
    "NC_candidate_seed"
  )]

  candidate_df <- myeloid_df %>% dplyr::filter(NC_candidate_seed)

  if (nrow(candidate_df) < dbscan_minPts) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = data.frame(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = data.frame(),
      cluster_summary = nc_empty_cluster_summary(image_use),
      fov_summary = nc_empty_fov_summary(
        image_use = image_use,
        severity_label = severity_label,
        n_total_cells = nrow(df),
        n_myeloid_cells = nrow(myeloid_df),
        area_cutoff = min_core_area_cutoff
      )
    ))
  }

  db <- dbscan::dbscan(
    as.matrix(candidate_df[, c("x", "y")]),
    eps = dbscan_eps,
    minPts = dbscan_minPts
  )

  candidate_df$raw_cluster <- db$cluster
  candidate_df <- candidate_df %>% dplyr::filter(raw_cluster > 0)

  if (nrow(candidate_df) == 0) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = data.frame(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = data.frame(),
      cluster_summary = nc_empty_cluster_summary(image_use),
      fov_summary = nc_empty_fov_summary(
        image_use = image_use,
        severity_label = severity_label,
        n_total_cells = nrow(df),
        n_myeloid_cells = nrow(myeloid_df),
        area_cutoff = min_core_area_cutoff
      )
    ))
  }

  raw_clusters <- sort(unique(candidate_df$raw_cluster))
  all_pts_sf <- sf::st_as_sf(df, coords = c("x", "y"), crs = NA)

  candidate_poly_list <- list()
  accepted_poly_list <- list()
  summary_list <- list()
  core_counter <- 0L

  for (cl in raw_clusters) {
    cl_pts <- candidate_df %>% dplyr::filter(raw_cluster == cl)
    if (nrow(cl_pts) < 3) next

    core_counter <- core_counter + 1L
    core_id <- paste0(image_use, "_core", core_counter)

    poly_sf <- tryCatch(
      points_to_core_polygon(cl_pts[, c("x", "y")], core_id = core_id, concavity = concavity),
      error = function(e) NULL
    )

    if (is.null(poly_sf) || nrow(poly_sf) == 0) next

    poly_sf <- sf::st_make_valid(poly_sf)

    if (!is.null(polygon_buffer) && is.finite(polygon_buffer) && polygon_buffer > 0) {
      poly_sf <- sf::st_buffer(poly_sf, dist = polygon_buffer)
      poly_sf <- sf::st_make_valid(poly_sf)
    }

    inside_list <- sf::st_within(all_pts_sf, poly_sf, sparse = TRUE)
    inside_flag <- lengths(inside_list) > 0

    sub_all <- df[inside_flag, , drop = FALSE]
    sub_myeloid <- sub_all[sub_all$is_myeloid_transfer, , drop = FALSE]
    sub_foam <- sub_all[sub_all$is_foam_transfer, , drop = FALSE]

    n_cells_in_polygon <- nrow(sub_all)
    n_myeloid_in_polygon <- nrow(sub_myeloid)
    n_foam_in_polygon <- nrow(sub_foam)
    foam_fraction_in_polygon <- ifelse(
      n_myeloid_in_polygon > 0,
      n_foam_in_polygon / n_myeloid_in_polygon,
      0
    )

    basic_pass <- (
      n_myeloid_in_polygon >= min_cluster_cells &
        n_foam_in_polygon >= min_cluster_foam &
        foam_fraction_in_polygon >= min_cluster_foam_fraction
    )

    if (!basic_pass) next

    area_val <- tryCatch(as.numeric(sf::st_area(poly_sf)), error = function(e) NA_real_)
    centroid <- tryCatch(
      sf::st_coordinates(sf::st_centroid(poly_sf)),
      error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1)
    )

    area_pass <- !is.na(area_val) && area_val >= min_core_area_cutoff

    candidate_poly_list[[core_id]] <- poly_sf
    if (area_pass) {
      accepted_poly_list[[core_id]] <- poly_sf
    }

    summary_list[[core_id]] <- tibble(
      FOV = image_use,
      severity = severity_label,
      core_id = core_id,
      area_pass_cutoff = area_pass,
      n_cells_in_polygon = n_cells_in_polygon,
      n_myeloid_in_polygon = n_myeloid_in_polygon,
      n_foam_in_polygon = n_foam_in_polygon,
      foam_fraction_in_polygon = foam_fraction_in_polygon,
      polygon_area = area_val,
      centroid_x = centroid[1, 1],
      centroid_y = centroid[1, 2],
      mean_local_foam_fraction_seed = mean(cl_pts$NC_local_foam_fraction, na.rm = TRUE),
      mean_local_foam_z_seed = mean(cl_pts$NC_local_foam_z, na.rm = TRUE),
      min_local_foam_p_seed = min(cl_pts$NC_local_foam_p, na.rm = TRUE)
    )
  }

  if (length(summary_list) == 0) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = data.frame(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = data.frame(),
      cluster_summary = nc_empty_cluster_summary(image_use),
      fov_summary = nc_empty_fov_summary(
        image_use = image_use,
        severity_label = severity_label,
        n_total_cells = nrow(df),
        n_myeloid_cells = nrow(myeloid_df),
        area_cutoff = min_core_area_cutoff
      )
    ))
  }

  cluster_summary <- bind_rows(summary_list)
  candidate_polygons_sf <- do.call(rbind, candidate_poly_list)
  candidate_polygons_df <- polygon_sf_to_df(candidate_polygons_sf)

  accepted_polygons_sf <- NULL
  accepted_polygons_df <- data.frame()

  if (length(accepted_poly_list) > 0) {
    accepted_polygons_sf <- do.call(rbind, accepted_poly_list)
    accepted_polygons_df <- polygon_sf_to_df(accepted_polygons_sf)

    inside_all <- sf::st_within(all_pts_sf, accepted_polygons_sf, sparse = TRUE)
    region_id_vec <- rep(NA_character_, nrow(df))
    for (i in seq_len(nrow(df))) {
      hit <- inside_all[[i]]
      if (length(hit) > 0) {
        region_id_vec[i] <- accepted_polygons_sf$core_id[hit[1]]
      }
    }

    df$NC_region_id <- region_id_vec
    df$NC_like_region <- !is.na(df$NC_region_id)
  }

  max_core_area <- max(cluster_summary$polygon_area, na.rm = TRUE)
  if (!is.finite(max_core_area)) max_core_area <- 0
  total_core_area <- sum(cluster_summary$polygon_area[cluster_summary$area_pass_cutoff], na.rm = TRUE)

  fov_summary <- tibble(
    FOV = image_use,
    severity = severity_label,
    NC_like_positive = max_core_area >= min_core_area_cutoff,
    NC_area_cutoff = min_core_area_cutoff,
    n_candidate_core_regions = nrow(cluster_summary),
    n_core_regions = sum(cluster_summary$area_pass_cutoff, na.rm = TRUE),
    max_core_area = max_core_area,
    total_core_area = total_core_area,
    n_total_cells = nrow(df),
    n_myeloid_cells = nrow(myeloid_df),
    n_core_like_cells = sum(df$NC_like_region, na.rm = TRUE),
    n_myeloid_core_like_cells = sum(df$NC_like_region & df$is_myeloid_transfer, na.rm = TRUE),
    n_foam_core_like_cells = sum(df$NC_like_region & df$is_foam_transfer, na.rm = TRUE)
  )

  return(list(
    cell_df = df,
    polygons_sf = accepted_polygons_sf,
    polygons_df = accepted_polygons_df,
    candidate_polygons_sf = candidate_polygons_sf,
    candidate_polygons_df = candidate_polygons_df,
    cluster_summary = cluster_summary,
    fov_summary = fov_summary
  ))
}

# ------------------------------------------------------------
# 所有 FOV：识别坏死核心样区域，并写回 xen@meta.data
# ------------------------------------------------------------
detect_necrotic_core_like_all_fovs <- function(
    obj,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    image_names = NULL,
    radius = 80,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.60,
    min_z = 2,
    max_p = 0.05,
    dbscan_eps = NULL,
    dbscan_minPts = 10,
    min_cluster_cells = 25,
    min_cluster_foam = 15,
    min_cluster_foam_fraction = 0.60,
    min_core_area_cutoff = 20000,
    concavity = 2,
    polygon_buffer = 0,
    seed = 123
) {

  check_nc_required_packages()

  if (is.null(image_names)) {
    image_names <- names(obj@images)
  }

  obj@meta.data$NC_local_total <- NA_integer_
  obj@meta.data$NC_local_foam_count <- NA_integer_
  obj@meta.data$NC_local_foam_fraction <- NA_real_
  obj@meta.data$NC_local_foam_z <- NA_real_
  obj@meta.data$NC_local_foam_p <- NA_real_
  obj@meta.data$NC_candidate_seed <- FALSE
  obj@meta.data$NC_region_id <- NA_character_
  obj@meta.data$NC_like_region <- FALSE
  obj@meta.data$NC_FOV_has_core <- FALSE
  obj@meta.data$NC_FOV_max_core_area <- 0
  obj@meta.data$NC_area_cutoff <- min_core_area_cutoff

  per_fov <- list()

  for (image_use in image_names) {
    message("Detecting necrotic-core-like regions in FOV: ", image_use)

    res <- tryCatch(
      detect_necrotic_core_like_one_fov(
        obj = obj,
        image_use = image_use,
        after_col = after_col,
        score_col = score_col,
        radius = radius,
        n_perm = n_perm,
        min_total = min_total,
        min_foam = min_foam,
        min_foam_fraction = min_foam_fraction,
        min_z = min_z,
        max_p = max_p,
        dbscan_eps = dbscan_eps,
        dbscan_minPts = dbscan_minPts,
        min_cluster_cells = min_cluster_cells,
        min_cluster_foam = min_cluster_foam,
        min_cluster_foam_fraction = min_cluster_foam_fraction,
        min_core_area_cutoff = min_core_area_cutoff,
        concavity = concavity,
        polygon_buffer = polygon_buffer,
        seed = seed
      ),
      error = function(e) {
        warning("Necrotic-core-like detection failed for ", image_use, ": ", e$message)
        df <- get_one_fov_df_for_nc(
          obj = obj,
          image_use = image_use,
          after_col = after_col,
          score_col = score_col
        )
        list(
          cell_df = df,
          polygons_sf = NULL,
          polygons_df = data.frame(),
          candidate_polygons_sf = NULL,
          candidate_polygons_df = data.frame(),
          cluster_summary = nc_empty_cluster_summary(image_use),
          fov_summary = nc_empty_fov_summary(
            image_use = image_use,
            severity_label = get_fov_severity_label(df),
            n_total_cells = nrow(df),
            n_myeloid_cells = sum(df$is_myeloid_transfer, na.rm = TRUE),
            area_cutoff = min_core_area_cutoff
          )
        )
      }
    )

    per_fov[[image_use]] <- res

    md <- res$cell_df %>%
      dplyr::select(
        cell,
        NC_local_total,
        NC_local_foam_count,
        NC_local_foam_fraction,
        NC_local_foam_z,
        NC_local_foam_p,
        NC_candidate_seed,
        NC_region_id,
        NC_like_region
      )

    common_cells <- intersect(md$cell, rownames(obj@meta.data))
    md <- md[match(common_cells, md$cell), , drop = FALSE]

    obj@meta.data[common_cells, c(
      "NC_local_total",
      "NC_local_foam_count",
      "NC_local_foam_fraction",
      "NC_local_foam_z",
      "NC_local_foam_p",
      "NC_candidate_seed",
      "NC_region_id",
      "NC_like_region"
    )] <- md[, c(
      "NC_local_total",
      "NC_local_foam_count",
      "NC_local_foam_fraction",
      "NC_local_foam_z",
      "NC_local_foam_p",
      "NC_candidate_seed",
      "NC_region_id",
      "NC_like_region"
    )]

    obj@meta.data[common_cells, "NC_FOV_has_core"] <- res$fov_summary$NC_like_positive[1]
    obj@meta.data[common_cells, "NC_FOV_max_core_area"] <- res$fov_summary$max_core_area[1]
    obj@meta.data[common_cells, "NC_area_cutoff"] <- min_core_area_cutoff
  }

  fov_summary <- bind_rows(lapply(per_fov, function(x) x$fov_summary))
  cluster_summary <- bind_rows(lapply(per_fov, function(x) x$cluster_summary))

  return(list(
    obj = obj,
    per_fov = per_fov,
    fov_summary = fov_summary,
    cluster_summary = cluster_summary
  ))
}

# ------------------------------------------------------------
# 单个 FOV：画坏死核心样区域
# 左图：髓系注释 + 不规则 core 边界
# 右图：所有细胞中 NC_like_region 标记
# ------------------------------------------------------------
plot_one_fov_necrotic_core_like <- function(
    det_obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    save_plot = TRUE,
    filename_base = NULL,
    width = 13,
    height = 6
) {

  if (!image_use %in% names(det_obj$per_fov)) {
    stop("No necrotic-core-like result found for FOV: ", image_use)
  }

  df <- det_obj$per_fov[[image_use]]$cell_df
  poly_df <- det_obj$per_fov[[image_use]]$polygons_df
  clus_sum <- det_obj$per_fov[[image_use]]$cluster_summary %>%
    dplyr::filter(area_pass_cutoff)

  severity_label <- get_fov_severity_label(df)
  n_core <- ifelse(nrow(clus_sum) > 0, nrow(clus_sum), 0)

  p1 <- ggplot(df, aes(x = x, y = y)) +
    geom_point(
      data = df %>% dplyr::filter(!is.na(.data[[after_col]])),
      aes(color = .data[[after_col]]),
      size = 0.18,
      alpha = 0.90
    ) +
    scale_color_manual(
      values = cell_colors[moma_levels],
      breaks = moma_levels,
      drop = FALSE
    ) +
    coord_fixed() +
    theme_void(base_size = 12) +
    labs(
      color = "Myeloid annotation",
      title = "Myeloid annotation with core-like boundaries"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )

  if (nrow(poly_df) > 0) {
    p1 <- p1 +
      geom_polygon(
        data = poly_df,
        aes(x = X, y = Y, group = group_path),
        inherit.aes = FALSE,
        fill = "firebrick2",
        alpha = 0.18,
        color = NA
      ) +
      geom_path(
        data = poly_df,
        aes(x = X, y = Y, group = group_path),
        inherit.aes = FALSE,
        color = "firebrick4",
        linewidth = 0.6
      )
  }

  if (nrow(clus_sum) > 0) {
    p1 <- p1 +
      geom_text(
        data = clus_sum,
        aes(x = centroid_x, y = centroid_y, label = core_id),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "bold",
        color = "firebrick4"
      )
  }

  p2 <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "grey85", size = 0.12, alpha = 0.80) +
    geom_point(
      data = df %>% dplyr::filter(NC_like_region),
      color = "firebrick3",
      size = 0.22,
      alpha = 0.95
    ) +
    coord_fixed() +
    theme_void(base_size = 12) +
    labs(
      title = "Cells assigned to core-like regions"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  if (nrow(poly_df) > 0) {
    p2 <- p2 +
      geom_polygon(
        data = poly_df,
        aes(x = X, y = Y, group = group_path),
        inherit.aes = FALSE,
        fill = "firebrick2",
        alpha = 0.20,
        color = NA
      ) +
      geom_path(
        data = poly_df,
        aes(x = X, y = Y, group = group_path),
        inherit.aes = FALSE,
        color = "firebrick4",
        linewidth = 0.7
      )
  }

  if (nrow(clus_sum) > 0) {
    p2 <- p2 +
      geom_text(
        data = clus_sum,
        aes(x = centroid_x, y = centroid_y, label = core_id),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "bold",
        color = "firebrick4"
      )
  }

  p_all <- (p1 | p2) +
    patchwork::plot_annotation(
      title = "Computational identification of necrotic-core-like regions",
      subtitle = paste0(
        "Section: ", image_use,
        " | Severity: ", severity_label,
        " | Core-like regions above area cutoff: ", n_core
      ),
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )

  print(p_all)

  if (save_plot) {
    if (is.null(filename_base)) {
      filename_base <- paste0(OUT_PREFIX, "_", sanitize_filename(image_use), "_necrotic_core_like_publication")
    }
    save_plot_both(p_all, filename_base, width = width, height = height)
  }

  invisible(p_all)
}

# ------------------------------------------------------------
# FOV 柱状图：Y 轴为每个 FOV 的最大坏死核心样区域面积
# 颜色由 severity 定义；高于 cutoff 为亮色，低于 cutoff 为偏灰色
# ------------------------------------------------------------
plot_fov_max_core_area_barplot <- function(
    fov_summary,
    area_cutoff = 20000,
    filename_base = paste0(OUT_PREFIX, "_necrotic_core_like_FOV_max_area_barplot_publication"),
    width = 11,
    height = 6.5
) {

  if (is.null(fov_summary) || nrow(fov_summary) == 0) {
    warning("fov_summary is empty; skip FOV max core area barplot.")
    return(NULL)
  }

  plot_df <- fov_summary %>%
    dplyr::mutate(
      severity = ifelse(is.na(severity) | severity == "", "not annotated", as.character(severity)),
      max_core_area = ifelse(is.na(max_core_area), 0, as.numeric(max_core_area)),
      NC_area_positive = max_core_area >= area_cutoff,
      area_status = ifelse(NC_area_positive, "above", "below"),
      fill_key = paste0(severity, "__", area_status)
    ) %>%
    dplyr::arrange(max_core_area, FOV) %>%
    dplyr::mutate(
      FOV_ordered = factor(FOV, levels = unique(FOV))
    )

  fill_values <- make_nc_area_fill_palette(unique(plot_df$severity))
  fill_values <- fill_values[intersect(names(fill_values), unique(plot_df$fill_key))]

  fill_labels <- names(fill_values)
  fill_labels <- gsub("__above", " above cutoff", fill_labels)
  fill_labels <- gsub("__below", " below cutoff", fill_labels)
  names(fill_labels) <- names(fill_values)

  cutoff_label_x <- plot_df$FOV_ordered[nrow(plot_df)]

  p <- ggplot(plot_df, aes(x = FOV_ordered, y = max_core_area, fill = fill_key)) +
    geom_col(width = 0.72, color = "grey30", linewidth = 0.15) +
    geom_hline(
      yintercept = area_cutoff,
      linetype = "dashed",
      linewidth = 0.6,
      color = "black"
    ) +
    annotate(
      "text",
      x = cutoff_label_x,
      y = area_cutoff,
      label = paste0("Area cutoff = ", scales::comma(round(area_cutoff, 1))),
      hjust = 1.05,
      vjust = -0.45,
      size = 3.5,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = fill_values,
      breaks = names(fill_values),
      labels = fill_labels,
      name = "Severity / cutoff status",
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.10))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "FOVs ranked by maximum core-like area",
      y = "Maximum necrotic-core-like area",
      title = "FOV-level burden of necrotic-core-like regions"
    ) +
    theme_classic(base_family = "Arial", base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(t = 20, r = 30, b = 10, l = 10, unit = "pt")
    )

  print(p)
  save_plot_both(p, filename_base, width = width, height = height)

  return(p)
}




# ============================================================
# 6.2 APOBEC3A 空间热图和疾病进展柱状图
#     新增输出，不改动已有图片文件名
# ============================================================

get_cell_fov_map_from_images <- function(obj) {

  if (!"images" %in% slotNames(obj) || length(obj@images) == 0) {
    return(NULL)
  }

  fov_map_list <- lapply(names(obj@images), function(image_use) {
    coord <- tryCatch(
      as.data.frame(GetTissueCoordinates(obj, image = image_use)),
      error = function(e) NULL
    )

    if (is.null(coord) || !"cell" %in% colnames(coord)) {
      return(NULL)
    }

    tibble::tibble(
      cell = as.character(coord$cell),
      FOV_for_plot = image_use
    )
  })

  fov_map <- dplyr::bind_rows(fov_map_list)
  if (nrow(fov_map) == 0) return(NULL)

  fov_map %>%
    dplyr::distinct(cell, .keep_all = TRUE)
}

plot_one_fov_APOBEC3A_heatmap <- function(
    obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    a3a_count_col = "APOBEC3A_count_transfer",
    save_plot = TRUE,
    filename_base = NULL,
    width = 7,
    height = 6,
    point_size_bg = 0.12,
    point_size_fg = 0.35,
    alpha_bg = 0.70,
    alpha_fg = 0.95,
    use_log10 = TRUE,
    show_only_positive = TRUE
) {

  df <- get_one_fov_df_for_nc(
    obj = obj,
    image_use = image_use,
    after_col = after_col,
    a3a_count_col = a3a_count_col
  )

  if (!a3a_count_col %in% colnames(df)) {
    stop("Column not found: ", a3a_count_col)
  }

  df$APOBEC3A_count_plot <- as.numeric(df[[a3a_count_col]])
  df$APOBEC3A_count_plot[is.na(df$APOBEC3A_count_plot)] <- 0

  if (use_log10) {
    df$APOBEC3A_heat <- log10(df$APOBEC3A_count_plot + 1)
    legend_title <- "log10(APOBEC3A count + 1)"
  } else {
    df$APOBEC3A_heat <- df$APOBEC3A_count_plot
    legend_title <- "APOBEC3A count"
  }

  severity_label <- get_fov_severity_label(df)

  if (show_only_positive) {
    df_fg <- df %>% dplyr::filter(APOBEC3A_count_plot > 0)
  } else {
    df_fg <- df
  }

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(
      color = "grey85",
      size = point_size_bg,
      alpha = alpha_bg
    ) +
    geom_point(
      data = df_fg,
      aes(color = APOBEC3A_heat),
      size = point_size_fg,
      alpha = alpha_fg
    ) +
    scale_color_viridis_c(
      option = "magma",
      name = legend_title
    ) +
    coord_fixed() +
    theme_void(base_size = 12) +
    labs(
      title = "Spatial expression of APOBEC3A",
      subtitle = paste0("Section: ", image_use, " | Severity: ", severity_label)
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )

  print(p)

  if (save_plot) {
    if (is.null(filename_base)) {
      filename_base <- paste0(
        OUT_PREFIX, "_",
        sanitize_filename(image_use),
        "_APOBEC3A_spatial_heatmap"
      )
    }
    save_plot_both(p, filename_base, width = width, height = height)
  }

  invisible(p)
}

plot_one_fov_APOBEC3A_heatmap_with_core <- function(
    obj,
    det_obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    a3a_count_col = "APOBEC3A_count_transfer",
    save_plot = TRUE,
    filename_base = NULL,
    width = 7,
    height = 6,
    use_log10 = TRUE,
    show_only_positive = TRUE,
    core_outline_color = "grey35",
    core_outline_linetype = "dashed",
    core_outline_linewidth = 0.7
) {

  df <- get_one_fov_df_for_nc(
    obj = obj,
    image_use = image_use,
    after_col = after_col,
    a3a_count_col = a3a_count_col
  )

  if (!a3a_count_col %in% colnames(df)) {
    stop("Column not found: ", a3a_count_col)
  }

  df$APOBEC3A_count_plot <- as.numeric(df[[a3a_count_col]])
  df$APOBEC3A_count_plot[is.na(df$APOBEC3A_count_plot)] <- 0

  if (use_log10) {
    df$APOBEC3A_heat <- log10(df$APOBEC3A_count_plot + 1)
    legend_title <- "log10(APOBEC3A count + 1)"
  } else {
    df$APOBEC3A_heat <- df$APOBEC3A_count_plot
    legend_title <- "APOBEC3A count"
  }

  if (show_only_positive) {
    df_fg <- df %>% dplyr::filter(APOBEC3A_count_plot > 0)
  } else {
    df_fg <- df
  }

  poly_df <- data.frame()

  if (!is.null(det_obj) && image_use %in% names(det_obj$per_fov)) {
    poly_df <- det_obj$per_fov[[image_use]]$polygons_df
  }

  severity_label <- get_fov_severity_label(df)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(
      color = "grey85",
      size = 0.12,
      alpha = 0.70
    ) +
    geom_point(
      data = df_fg,
      aes(color = APOBEC3A_heat),
      size = 0.35,
      alpha = 0.95
    ) +
    scale_color_viridis_c(
      option = "magma",
      name = legend_title
    ) +
    coord_fixed() +
    theme_void(base_size = 12) +
    labs(
      title = "Spatial expression of APOBEC3A",
      subtitle = paste0("Section: ", image_use, " | Severity: ", severity_label)
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )

  if (!is.null(poly_df) && nrow(poly_df) > 0) {
    p <- p +
      geom_path(
        data = poly_df,
        aes(x = X, y = Y, group = group_path),
        inherit.aes = FALSE,
        color = core_outline_color,
        linetype = core_outline_linetype,
        linewidth = core_outline_linewidth
      )
  }

  print(p)

  if (save_plot) {
    if (is.null(filename_base)) {
      filename_base <- paste0(
        OUT_PREFIX, "_",
        sanitize_filename(image_use),
        "_APOBEC3A_spatial_heatmap_with_core"
      )
    }
    save_plot_both(p, filename_base, width = width, height = height)
  }

  invisible(p)
}

plot_A3A_overall_expression_and_positive_rate_by_disease <- function(
    obj,
    after_col = "Xenium_MoMa_scPred_Celltype",
    a3a_count_col = "APOBEC3A_count_transfer",
    score_col = "Xenium_MoMa_scPred_score",
    use_high_confidence_only = FALSE,
    confidence_cutoff = 0.5,
    filename_base = paste0(OUT_PREFIX, "_APOBEC3A_overall_myeloid_by_disease"),
    width = 7,
    height = 5.5
) {

  meta <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  severity_col <- get_severity_col(meta)

  if (is.na(severity_col)) {
    stop(
      "Cannot find severity / disease column in metadata. Available columns:\n",
      paste(colnames(meta), collapse = "\n")
    )
  }

  if (!after_col %in% colnames(meta)) {
    stop("Cannot find myeloid annotation column: ", after_col)
  }

  if (!a3a_count_col %in% colnames(meta)) {
    stop("Cannot find APOBEC3A count column: ", a3a_count_col)
  }

  fov_map <- get_cell_fov_map_from_images(obj)

  if (!is.null(fov_map)) {
    meta <- meta %>%
      dplyr::left_join(fov_map, by = "cell")
  } else {
    fov_col <- intersect(
      c("FOV", "fov", "image", "Image", "sample", "Sample", "orig.ident"),
      colnames(meta)
    )

    if (length(fov_col) > 0) {
      meta$FOV_for_plot <- as.character(meta[[fov_col[1]]])
    } else {
      meta$FOV_for_plot <- "all_cells"
      warning("No FOV information found; using one fallback group.")
    }
  }

  meta$Severity <- as.character(meta[[severity_col]])
  meta$CellType <- as.character(meta[[after_col]])
  meta$APOBEC3A_count <- as.numeric(meta[[a3a_count_col]])
  meta$APOBEC3A_count[is.na(meta$APOBEC3A_count)] <- 0
  meta$APOBEC3A_detected <- meta$APOBEC3A_count > 0

  if (score_col %in% colnames(meta)) {
    meta$TransferScore <- as.numeric(meta[[score_col]])
  } else {
    meta$TransferScore <- NA_real_
  }

  plot_meta <- meta %>%
    dplyr::filter(
      !is.na(Severity),
      Severity != "",
      !is.na(FOV_for_plot),
      FOV_for_plot != "",
      !is.na(CellType),
      CellType %in% moma_levels
    )

  if (use_high_confidence_only) {
    plot_meta <- plot_meta %>%
      dplyr::filter(
        !is.na(TransferScore),
        TransferScore >= confidence_cutoff
      )
  }

  if (nrow(plot_meta) == 0) {
    stop("No myeloid cells available for APOBEC3A disease-stage plot.")
  }

  observed_severity <- unique(plot_meta$Severity)

  severity_order <- c(
    severity_levels[severity_levels %in% observed_severity],
    setdiff(observed_severity, severity_levels)
  )

  plot_meta$Severity <- factor(plot_meta$Severity, levels = severity_order)

  fov_summary <- plot_meta %>%
    dplyr::group_by(FOV_for_plot, Severity) %>%
    dplyr::summarise(
      n_myeloid_cells = dplyr::n(),
      n_APOBEC3A_positive = sum(APOBEC3A_detected, na.rm = TRUE),
      APOBEC3A_positive_fraction = mean(APOBEC3A_detected, na.rm = TRUE),
      mean_APOBEC3A_count = mean(APOBEC3A_count, na.rm = TRUE),
      mean_log10_APOBEC3A_count = mean(log10(APOBEC3A_count + 1), na.rm = TRUE),
      .groups = "drop"
    )

  expression_summary <- fov_summary %>%
    dplyr::group_by(Severity) %>%
    dplyr::summarise(
      n_FOV = dplyr::n(),
      mean_value = mean(mean_log10_APOBEC3A_count, na.rm = TRUE),
      sd_value = sd(mean_log10_APOBEC3A_count, na.rm = TRUE),
      sem_value = sd_value / sqrt(n_FOV),
      total_myeloid_cells = sum(n_myeloid_cells, na.rm = TRUE),
      .groups = "drop"
    )

  detection_summary <- fov_summary %>%
    dplyr::group_by(Severity) %>%
    dplyr::summarise(
      n_FOV = dplyr::n(),
      mean_value = mean(APOBEC3A_positive_fraction, na.rm = TRUE),
      sd_value = sd(APOBEC3A_positive_fraction, na.rm = TRUE),
      sem_value = sd_value / sqrt(n_FOV),
      total_myeloid_cells = sum(n_myeloid_cells, na.rm = TRUE),
      .groups = "drop"
    )

  p_expr <- ggplot(
    expression_summary,
    aes(x = Severity, y = mean_value, fill = Severity)
  ) +
    geom_col(
      width = 0.65,
      color = "black",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(
        ymin = mean_value - sem_value,
        ymax = mean_value + sem_value
      ),
      width = 0.20,
      linewidth = 0.35
    ) +
    geom_point(
      data = fov_summary,
      aes(x = Severity, y = mean_log10_APOBEC3A_count),
      position = position_jitter(width = 0.12, height = 0),
      shape = 21,
      size = 2,
      color = "black",
      fill = "white",
      stroke = 0.3,
      alpha = 0.90,
      inherit.aes = FALSE
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = "Disease stage",
      y = "Mean log10(APOBEC3A count + 1)",
      title = "APOBEC3A expression in myeloid cells across disease stages"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "none"
    )

  p_detect <- ggplot(
    detection_summary,
    aes(x = Severity, y = mean_value, fill = Severity)
  ) +
    geom_col(
      width = 0.65,
      color = "black",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_value - sem_value, 0),
        ymax = mean_value + sem_value
      ),
      width = 0.20,
      linewidth = 0.35
    ) +
    geom_point(
      data = fov_summary,
      aes(x = Severity, y = APOBEC3A_positive_fraction),
      position = position_jitter(width = 0.12, height = 0),
      shape = 21,
      size = 2,
      color = "black",
      fill = "white",
      stroke = 0.3,
      alpha = 0.90,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_classic(base_size = 12) +
    labs(
      x = "Disease stage",
      y = "APOBEC3A-positive fraction",
      title = "APOBEC3A-positive myeloid cells across disease stages"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "none"
    )

  print(p_expr)
  print(p_detect)

  stat_res <- tibble::tibble(
    metric = c(
      "mean_log10_APOBEC3A_count",
      "APOBEC3A_positive_fraction"
    ),
    n_stage = dplyr::n_distinct(fov_summary$Severity),
    p_kruskal = c(
      ifelse(
        dplyr::n_distinct(fov_summary$Severity) >= 2,
        kruskal.test(
          fov_summary$mean_log10_APOBEC3A_count ~ fov_summary$Severity
        )$p.value,
        NA_real_
      ),
      ifelse(
        dplyr::n_distinct(fov_summary$Severity) >= 2,
        kruskal.test(
          fov_summary$APOBEC3A_positive_fraction ~ fov_summary$Severity
        )$p.value,
        NA_real_
      )
    )
  )

  write.csv(
    fov_summary,
    file.path(TABLE_DIR, paste0(filename_base, "_FOV_level_values.csv")),
    row.names = FALSE
  )

  write.csv(
    expression_summary,
    file.path(TABLE_DIR, paste0(filename_base, "_expression_bar_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    detection_summary,
    file.path(TABLE_DIR, paste0(filename_base, "_positive_fraction_bar_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    stat_res,
    file.path(TABLE_DIR, paste0(filename_base, "_Kruskal_test.csv")),
    row.names = FALSE
  )

  save_plot_both(
    p_expr,
    paste0(filename_base, "_expression_barplot"),
    width = width,
    height = height
  )

  save_plot_both(
    p_detect,
    paste0(filename_base, "_positive_fraction_barplot"),
    width = width,
    height = height
  )

  message("========== APOBEC3A disease-stage FOV-level values ==========")
  print(fov_summary)

  message("========== APOBEC3A expression summary ==========")
  print(expression_summary)

  message("========== APOBEC3A positive fraction summary ==========")
  print(detection_summary)

  message("========== APOBEC3A Kruskal-Wallis tests ==========")
  print(stat_res)

  return(list(
    expression_plot = p_expr,
    positive_fraction_plot = p_detect,
    fov_summary = fov_summary,
    expression_summary = expression_summary,
    detection_summary = detection_summary,
    stat_res = stat_res
  ))
}


# ------------------------------------------------------------
# 参数打印：便于确认当前使用哪一套阈值
# ------------------------------------------------------------
print_nc_parameters <- function(
    local_params = NC_LOCAL_PARAMS,
    cluster_params = NC_CLUSTER_PARAMS,
    area_cutoff = NC_CORE_AREA_CUTOFF,
    image_names = NC_IMAGE_NAMES,
    random_seed = NC_RANDOM_SEED
) {
  message("========== Current necrotic-core-like region parameters ==========")
  message("NC_IMAGE_NAMES = ", ifelse(is.null(image_names), "ALL FOVs", paste(image_names, collapse = ", ")))
  message("NC_CORE_AREA_CUTOFF = ", area_cutoff)
  message("NC_RANDOM_SEED = ", random_seed)
  message("---------- Local enrichment seed parameters ----------")
  print(local_params)
  message("---------- DBSCAN / polygon / cluster filter parameters ----------")
  print(cluster_params)
  invisible(TRUE)
}

# ------------------------------------------------------------
# 诊断 1：local enrichment seed 是否存在
# 用于定位 min_total / min_foam / min_foam_fraction / min_z / max_p 是否过严
# ------------------------------------------------------------
diagnose_necrotic_core_like_thresholds_one_fov <- function(
    obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    radius = 80,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.60,
    min_z = 2,
    max_p = 0.05,
    seed = 123
) {
  check_nc_required_packages()

  df <- get_one_fov_df_for_nc(
    obj = obj,
    image_use = image_use,
    after_col = after_col,
    score_col = score_col
  )

  myeloid_df <- df %>% dplyr::filter(is_myeloid_transfer)

  if (nrow(myeloid_df) < 5) {
    return(tibble::tibble(
      FOV = image_use,
      n_total_cells = nrow(df),
      n_myeloid = nrow(myeloid_df),
      n_foam = sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
      global_foam_fraction = NA_real_,
      local_total_q50 = NA_real_, local_total_q75 = NA_real_, local_total_q90 = NA_real_, local_total_q95 = NA_real_, local_total_max = NA_real_,
      local_foam_q50 = NA_real_, local_foam_q75 = NA_real_, local_foam_q90 = NA_real_, local_foam_q95 = NA_real_, local_foam_max = NA_real_,
      local_foam_fraction_q50 = NA_real_, local_foam_fraction_q75 = NA_real_, local_foam_fraction_q90 = NA_real_, local_foam_fraction_q95 = NA_real_, local_foam_fraction_max = NA_real_,
      local_z_q50 = NA_real_, local_z_q75 = NA_real_, local_z_q90 = NA_real_, local_z_q95 = NA_real_, local_z_max = NA_real_,
      local_p_min = NA_real_, local_p_q05 = NA_real_, local_p_q10 = NA_real_,
      pass_total = 0L, pass_foam = 0L, pass_fraction = 0L, pass_z = 0L, pass_p = 0L,
      pass_total_foam_fraction = 0L, pass_all_seed = 0L
    ))
  }

  coords <- as.matrix(myeloid_df[, c("x", "y")])
  nn <- dbscan::frNN(coords, eps = radius)

  obs_total <- lengths(nn$id)
  obs_foam <- sapply(nn$id, function(idx) sum(myeloid_df$is_foam_transfer[idx], na.rm = TRUE))
  obs_frac <- ifelse(obs_total > 0, obs_foam / obs_total, NA_real_)

  set.seed(seed)
  perm_mat <- replicate(n_perm, {
    shuffled <- sample(myeloid_df$is_foam_transfer)
    sapply(nn$id, function(idx) sum(shuffled[idx], na.rm = TRUE))
  })
  if (is.null(dim(perm_mat))) perm_mat <- matrix(perm_mat, ncol = 1)

  perm_mean <- rowMeans(perm_mat, na.rm = TRUE)
  perm_sd <- apply(perm_mat, 1, sd, na.rm = TRUE)
  perm_sd[is.na(perm_sd) | perm_sd == 0] <- 1e-8

  obs_z <- (obs_foam - perm_mean) / perm_sd
  obs_p <- (rowSums(perm_mat >= obs_foam) + 1) / (n_perm + 1)

  pass_total <- obs_total >= min_total
  pass_foam <- obs_foam >= min_foam
  pass_fraction <- obs_frac >= min_foam_fraction
  pass_z <- obs_z >= min_z
  pass_p <- obs_p <= max_p
  pass_total_foam_fraction <- pass_total & pass_foam & pass_fraction
  pass_all_seed <- pass_total & pass_foam & pass_fraction & pass_z & pass_p

  tibble::tibble(
    FOV = image_use,
    n_total_cells = nrow(df),
    n_myeloid = nrow(myeloid_df),
    n_foam = sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
    global_foam_fraction = mean(myeloid_df$is_foam_transfer, na.rm = TRUE),

    local_total_q50 = as.numeric(stats::quantile(obs_total, 0.50, na.rm = TRUE)),
    local_total_q75 = as.numeric(stats::quantile(obs_total, 0.75, na.rm = TRUE)),
    local_total_q90 = as.numeric(stats::quantile(obs_total, 0.90, na.rm = TRUE)),
    local_total_q95 = as.numeric(stats::quantile(obs_total, 0.95, na.rm = TRUE)),
    local_total_max = max(obs_total, na.rm = TRUE),

    local_foam_q50 = as.numeric(stats::quantile(obs_foam, 0.50, na.rm = TRUE)),
    local_foam_q75 = as.numeric(stats::quantile(obs_foam, 0.75, na.rm = TRUE)),
    local_foam_q90 = as.numeric(stats::quantile(obs_foam, 0.90, na.rm = TRUE)),
    local_foam_q95 = as.numeric(stats::quantile(obs_foam, 0.95, na.rm = TRUE)),
    local_foam_max = max(obs_foam, na.rm = TRUE),

    local_foam_fraction_q50 = as.numeric(stats::quantile(obs_frac, 0.50, na.rm = TRUE)),
    local_foam_fraction_q75 = as.numeric(stats::quantile(obs_frac, 0.75, na.rm = TRUE)),
    local_foam_fraction_q90 = as.numeric(stats::quantile(obs_frac, 0.90, na.rm = TRUE)),
    local_foam_fraction_q95 = as.numeric(stats::quantile(obs_frac, 0.95, na.rm = TRUE)),
    local_foam_fraction_max = max(obs_frac, na.rm = TRUE),

    local_z_q50 = as.numeric(stats::quantile(obs_z, 0.50, na.rm = TRUE)),
    local_z_q75 = as.numeric(stats::quantile(obs_z, 0.75, na.rm = TRUE)),
    local_z_q90 = as.numeric(stats::quantile(obs_z, 0.90, na.rm = TRUE)),
    local_z_q95 = as.numeric(stats::quantile(obs_z, 0.95, na.rm = TRUE)),
    local_z_max = max(obs_z, na.rm = TRUE),

    local_p_min = min(obs_p, na.rm = TRUE),
    local_p_q05 = as.numeric(stats::quantile(obs_p, 0.05, na.rm = TRUE)),
    local_p_q10 = as.numeric(stats::quantile(obs_p, 0.10, na.rm = TRUE)),

    pass_total = sum(pass_total, na.rm = TRUE),
    pass_foam = sum(pass_foam, na.rm = TRUE),
    pass_fraction = sum(pass_fraction, na.rm = TRUE),
    pass_z = sum(pass_z, na.rm = TRUE),
    pass_p = sum(pass_p, na.rm = TRUE),
    pass_total_foam_fraction = sum(pass_total_foam_fraction, na.rm = TRUE),
    pass_all_seed = sum(pass_all_seed, na.rm = TRUE)
  )
}

diagnose_necrotic_core_like_thresholds_all_fovs <- function(
    obj,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    image_names = names(obj@images),
    radius = 80,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.60,
    min_z = 2,
    max_p = 0.05,
    seed = 123,
    output_prefix = paste0(OUT_PREFIX, "_necrotic_core_like_threshold_diagnostics")
) {
  message("========== Diagnose local enrichment thresholds before final core detection ==========")
  message("radius = ", radius, "; min_total = ", min_total, "; min_foam = ", min_foam,
          "; min_foam_fraction = ", min_foam_fraction, "; min_z = ", min_z, "; max_p = ", max_p)

  diag_tbl <- dplyr::bind_rows(lapply(image_names, function(image_use) {
    message("Diagnosing local threshold for FOV: ", image_use)
    diagnose_necrotic_core_like_thresholds_one_fov(
      obj = obj,
      image_use = image_use,
      after_col = after_col,
      score_col = score_col,
      radius = radius,
      n_perm = n_perm,
      min_total = min_total,
      min_foam = min_foam,
      min_foam_fraction = min_foam_fraction,
      min_z = min_z,
      max_p = max_p,
      seed = seed
    )
  }))

  diag_tbl <- diag_tbl %>%
    dplyr::arrange(dplyr::desc(pass_all_seed), dplyr::desc(local_foam_max), dplyr::desc(local_foam_fraction_max))

  message("========== Key diagnostic columns: local enrichment ==========")
  print(
    diag_tbl %>%
      dplyr::select(
        FOV, n_myeloid, n_foam, global_foam_fraction,
        local_total_max, local_foam_max, local_foam_fraction_max,
        local_z_max, local_p_min,
        pass_total, pass_foam, pass_fraction, pass_z, pass_p,
        pass_total_foam_fraction, pass_all_seed
      ) %>%
      as.data.frame()
  )

  diag_out <- file.path(TABLE_DIR, paste0(output_prefix, ".csv"))
  write.csv(diag_tbl, diag_out, row.names = FALSE)
  message("Saved local threshold diagnostic table: ", diag_out)

  message("---------- How to interpret local diagnostics ----------")
  message("pass_total == 0: increase radius or decrease min_total")
  message("pass_foam == 0: decrease min_foam")
  message("pass_fraction == 0: decrease min_foam_fraction")
  message("pass_total_foam_fraction > 0 but pass_all_seed == 0: relax min_z or max_p")

  return(diag_tbl)
}

# ------------------------------------------------------------
# 诊断 2：DBSCAN / polygon / cluster filter 是哪一步筛没
# 用于定位 dbscan_eps / dbscan_minPts / min_cluster_* / concavity 是否不合适
# ------------------------------------------------------------
diagnose_nc_dbscan_polygon_one_fov <- function(
    obj,
    image_use,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    radius = 80,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.60,
    min_z = 2,
    max_p = 0.05,
    dbscan_eps = 80,
    dbscan_minPts = 10,
    min_cluster_cells = 25,
    min_cluster_foam = 15,
    min_cluster_foam_fraction = 0.60,
    concavity = 2,
    seed = 123
) {
  check_nc_required_packages()

  df <- get_one_fov_df_for_nc(
    obj = obj,
    image_use = image_use,
    after_col = after_col,
    score_col = score_col
  )

  myeloid_df <- df %>% dplyr::filter(is_myeloid_transfer)

  if (nrow(myeloid_df) < 5) {
    return(tibble::tibble(
      FOV = image_use,
      dbscan_eps = dbscan_eps,
      dbscan_minPts = dbscan_minPts,
      n_myeloid = nrow(myeloid_df),
      n_foam = sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
      n_candidate_seed = 0L,
      n_dbscan_clusters = 0L,
      raw_cluster = NA_integer_,
      n_seed_points = NA_integer_,
      polygon_area = NA_real_,
      n_cells_in_polygon = NA_integer_,
      n_myeloid_in_polygon = NA_integer_,
      n_foam_in_polygon = NA_integer_,
      foam_fraction_in_polygon = NA_real_,
      pass_min_cluster_cells = FALSE,
      pass_min_cluster_foam = FALSE,
      pass_min_cluster_foam_fraction = FALSE,
      pass_all_cluster_filter = FALSE,
      reason = "too_few_myeloid"
    ))
  }

  coords <- as.matrix(myeloid_df[, c("x", "y")])
  nn <- dbscan::frNN(coords, eps = radius)

  obs_total <- lengths(nn$id)
  obs_foam <- sapply(nn$id, function(idx) sum(myeloid_df$is_foam_transfer[idx], na.rm = TRUE))
  obs_frac <- ifelse(obs_total > 0, obs_foam / obs_total, NA_real_)

  set.seed(seed)
  perm_mat <- replicate(n_perm, {
    shuffled <- sample(myeloid_df$is_foam_transfer)
    sapply(nn$id, function(idx) sum(shuffled[idx], na.rm = TRUE))
  })
  if (is.null(dim(perm_mat))) perm_mat <- matrix(perm_mat, ncol = 1)

  perm_mean <- rowMeans(perm_mat, na.rm = TRUE)
  perm_sd <- apply(perm_mat, 1, sd, na.rm = TRUE)
  perm_sd[is.na(perm_sd) | perm_sd == 0] <- 1e-8

  obs_z <- (obs_foam - perm_mean) / perm_sd
  obs_p <- (rowSums(perm_mat >= obs_foam) + 1) / (n_perm + 1)

  seed_flag <- (
    obs_total >= min_total &
      obs_foam >= min_foam &
      obs_frac >= min_foam_fraction &
      obs_z >= min_z &
      obs_p <= max_p
  )

  myeloid_df$NC_local_total <- as.integer(obs_total)
  myeloid_df$NC_local_foam_count <- as.integer(obs_foam)
  myeloid_df$NC_local_foam_fraction <- as.numeric(obs_frac)
  myeloid_df$NC_local_foam_z <- as.numeric(obs_z)
  myeloid_df$NC_local_foam_p <- as.numeric(obs_p)
  myeloid_df$NC_candidate_seed <- seed_flag

  candidate_df <- myeloid_df %>% dplyr::filter(NC_candidate_seed)

  message("========== DBSCAN / polygon diagnostic: ", image_use, " ==========")
  message("n_myeloid = ", nrow(myeloid_df), "; n_foam = ", sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
          "; n_candidate_seed = ", nrow(candidate_df))

  if (nrow(candidate_df) < dbscan_minPts) {
    return(tibble::tibble(
      FOV = image_use,
      dbscan_eps = dbscan_eps,
      dbscan_minPts = dbscan_minPts,
      n_myeloid = nrow(myeloid_df),
      n_foam = sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
      n_candidate_seed = nrow(candidate_df),
      n_dbscan_clusters = 0L,
      raw_cluster = NA_integer_,
      n_seed_points = NA_integer_,
      polygon_area = NA_real_,
      n_cells_in_polygon = NA_integer_,
      n_myeloid_in_polygon = NA_integer_,
      n_foam_in_polygon = NA_integer_,
      foam_fraction_in_polygon = NA_real_,
      pass_min_cluster_cells = FALSE,
      pass_min_cluster_foam = FALSE,
      pass_min_cluster_foam_fraction = FALSE,
      pass_all_cluster_filter = FALSE,
      reason = "candidate_seed_less_than_dbscan_minPts"
    ))
  }

  db <- dbscan::dbscan(
    as.matrix(candidate_df[, c("x", "y")]),
    eps = dbscan_eps,
    minPts = dbscan_minPts
  )

  candidate_df$raw_cluster <- db$cluster
  n_noise <- sum(candidate_df$raw_cluster == 0)
  cluster_ids <- sort(setdiff(unique(candidate_df$raw_cluster), 0))
  n_clusters <- length(cluster_ids)

  message("dbscan_eps = ", dbscan_eps, "; dbscan_minPts = ", dbscan_minPts,
          "; n_noise = ", n_noise, "; n_dbscan_clusters = ", n_clusters)

  if (n_clusters == 0) {
    return(tibble::tibble(
      FOV = image_use,
      dbscan_eps = dbscan_eps,
      dbscan_minPts = dbscan_minPts,
      n_myeloid = nrow(myeloid_df),
      n_foam = sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
      n_candidate_seed = nrow(candidate_df),
      n_dbscan_clusters = 0L,
      raw_cluster = NA_integer_,
      n_seed_points = NA_integer_,
      polygon_area = NA_real_,
      n_cells_in_polygon = NA_integer_,
      n_myeloid_in_polygon = NA_integer_,
      n_foam_in_polygon = NA_integer_,
      foam_fraction_in_polygon = NA_real_,
      pass_min_cluster_cells = FALSE,
      pass_min_cluster_foam = FALSE,
      pass_min_cluster_foam_fraction = FALSE,
      pass_all_cluster_filter = FALSE,
      reason = "no_dbscan_cluster"
    ))
  }

  all_pts_sf <- sf::st_as_sf(df, coords = c("x", "y"), crs = NA)
  out_list <- list()

  for (cl in cluster_ids) {
    cl_pts <- candidate_df %>% dplyr::filter(raw_cluster == cl)
    reason <- "tested"
    polygon_area <- NA_real_
    n_cells_in_polygon <- NA_integer_
    n_myeloid_in_polygon <- NA_integer_
    n_foam_in_polygon <- NA_integer_
    foam_fraction_in_polygon <- NA_real_

    if (nrow(cl_pts) < 3) {
      reason <- "cluster_less_than_3_points"
    } else {
      core_id <- paste0(image_use, "_rawcluster", cl)
      poly_sf <- tryCatch(
        points_to_core_polygon(cl_pts[, c("x", "y")], core_id = core_id, concavity = concavity),
        error = function(e) {
          message("polygon failed: ", image_use, " cluster ", cl, " | ", e$message)
          NULL
        }
      )

      if (is.null(poly_sf) || nrow(poly_sf) == 0) {
        reason <- "polygon_failed"
      } else {
        poly_sf <- sf::st_make_valid(poly_sf)
        polygon_area <- tryCatch(as.numeric(sf::st_area(poly_sf)), error = function(e) NA_real_)
        inside_list <- sf::st_intersects(all_pts_sf, poly_sf, sparse = TRUE)
        inside_flag <- lengths(inside_list) > 0

        sub_all <- df[inside_flag, , drop = FALSE]
        sub_myeloid <- sub_all[sub_all$is_myeloid_transfer, , drop = FALSE]
        sub_foam <- sub_all[sub_all$is_foam_transfer, , drop = FALSE]

        n_cells_in_polygon <- nrow(sub_all)
        n_myeloid_in_polygon <- nrow(sub_myeloid)
        n_foam_in_polygon <- nrow(sub_foam)
        foam_fraction_in_polygon <- ifelse(
          n_myeloid_in_polygon > 0,
          n_foam_in_polygon / n_myeloid_in_polygon,
          NA_real_
        )
      }
    }

    pass_min_cluster_cells <- !is.na(n_myeloid_in_polygon) && n_myeloid_in_polygon >= min_cluster_cells
    pass_min_cluster_foam <- !is.na(n_foam_in_polygon) && n_foam_in_polygon >= min_cluster_foam
    pass_min_cluster_foam_fraction <- !is.na(foam_fraction_in_polygon) && foam_fraction_in_polygon >= min_cluster_foam_fraction
    pass_all_cluster_filter <- pass_min_cluster_cells && pass_min_cluster_foam && pass_min_cluster_foam_fraction

    if (reason == "tested" && !pass_all_cluster_filter) {
      reason <- paste0(
        "failed_cluster_filter: cells=", n_myeloid_in_polygon,
        ", foam=", n_foam_in_polygon,
        ", frac=", round(foam_fraction_in_polygon, 3)
      )
    }
    if (reason == "tested" && pass_all_cluster_filter) reason <- "pass"

    out_list[[as.character(cl)]] <- tibble::tibble(
      FOV = image_use,
      dbscan_eps = dbscan_eps,
      dbscan_minPts = dbscan_minPts,
      n_myeloid = nrow(myeloid_df),
      n_foam = sum(myeloid_df$is_foam_transfer, na.rm = TRUE),
      n_candidate_seed = nrow(candidate_df),
      n_dbscan_clusters = n_clusters,
      raw_cluster = cl,
      n_seed_points = nrow(cl_pts),
      polygon_area = polygon_area,
      n_cells_in_polygon = n_cells_in_polygon,
      n_myeloid_in_polygon = n_myeloid_in_polygon,
      n_foam_in_polygon = n_foam_in_polygon,
      foam_fraction_in_polygon = foam_fraction_in_polygon,
      pass_min_cluster_cells = pass_min_cluster_cells,
      pass_min_cluster_foam = pass_min_cluster_foam,
      pass_min_cluster_foam_fraction = pass_min_cluster_foam_fraction,
      pass_all_cluster_filter = pass_all_cluster_filter,
      reason = reason
    )
  }

  res <- dplyr::bind_rows(out_list) %>%
    dplyr::arrange(
      dplyr::desc(pass_all_cluster_filter),
      dplyr::desc(n_seed_points),
      dplyr::desc(n_foam_in_polygon),
      dplyr::desc(polygon_area)
    )

  print(as.data.frame(res))
  return(res)
}

diagnose_nc_dbscan_polygon_all_fovs <- function(
    obj,
    after_col = "Xenium_MoMa_scPred_Celltype",
    score_col = "Xenium_MoMa_scPred_score",
    image_names = names(obj@images),
    radius = 80,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.60,
    min_z = 2,
    max_p = 0.05,
    dbscan_eps = 80,
    dbscan_minPts = 10,
    min_cluster_cells = 25,
    min_cluster_foam = 15,
    min_cluster_foam_fraction = 0.60,
    concavity = 2,
    seed = 123,
    output_prefix = paste0(OUT_PREFIX, "_necrotic_core_like_DBSCAN_polygon_diagnostics")
) {
  message("========== Diagnose DBSCAN / polygon / cluster filter before final core detection ==========")
  message("dbscan_eps = ", dbscan_eps, "; dbscan_minPts = ", dbscan_minPts,
          "; min_cluster_cells = ", min_cluster_cells,
          "; min_cluster_foam = ", min_cluster_foam,
          "; min_cluster_foam_fraction = ", min_cluster_foam_fraction,
          "; concavity = ", concavity)

  res <- dplyr::bind_rows(lapply(image_names, function(image_use) {
    diagnose_nc_dbscan_polygon_one_fov(
      obj = obj,
      image_use = image_use,
      after_col = after_col,
      score_col = score_col,
      radius = radius,
      n_perm = n_perm,
      min_total = min_total,
      min_foam = min_foam,
      min_foam_fraction = min_foam_fraction,
      min_z = min_z,
      max_p = max_p,
      dbscan_eps = dbscan_eps,
      dbscan_minPts = dbscan_minPts,
      min_cluster_cells = min_cluster_cells,
      min_cluster_foam = min_cluster_foam,
      min_cluster_foam_fraction = min_cluster_foam_fraction,
      concavity = concavity,
      seed = seed
    )
  }))

  out_csv <- file.path(TABLE_DIR, paste0(output_prefix, ".csv"))
  write.csv(res, out_csv, row.names = FALSE)
  message("Saved DBSCAN / polygon diagnostic table: ", out_csv)

  message("========== Passed clusters ==========")
  print(
    res %>%
      dplyr::filter(pass_all_cluster_filter) %>%
      dplyr::arrange(dplyr::desc(polygon_area)) %>%
      as.data.frame()
  )

  message("========== Failed clusters, top by seed points ==========")
  print(
    res %>%
      dplyr::filter(!pass_all_cluster_filter) %>%
      dplyr::arrange(dplyr::desc(n_seed_points), dplyr::desc(n_foam_in_polygon)) %>%
      head(30) %>%
      as.data.frame()
  )

  message("---------- How to interpret DBSCAN / polygon diagnostics ----------")
  message("n_candidate_seed > 0 but n_dbscan_clusters == 0: increase dbscan_eps or decrease dbscan_minPts")
  message("polygon_area is NA or reason == polygon_failed: increase concavity or install/check concaveman")
  message("failed_cluster_filter with enough foam but low fraction: decrease min_cluster_foam_fraction")
  message("failed_cluster_filter with few cells/foam: decrease min_cluster_cells or min_cluster_foam")

  return(res)
}

# ------------------------------------------------------------
# 面积分布 summary：用于决定 NC_CORE_AREA_CUTOFF
# ------------------------------------------------------------
summarise_nc_core_area_distribution <- function(cluster_summary) {
  if (is.null(cluster_summary) || nrow(cluster_summary) == 0) {
    message("No candidate core-like regions were detected before area cutoff. Area summary cannot be calculated.")
    return(tibble::tibble(
      n_candidate_core_regions = 0L,
      mean_core_area = NA_real_,
      median_core_area = NA_real_,
      area_q25 = NA_real_,
      area_q75 = NA_real_,
      area_q90 = NA_real_,
      max_core_area = NA_real_
    ))
  }

  area_summary <- cluster_summary %>%
    dplyr::summarise(
      n_candidate_core_regions = dplyr::n(),
      n_area_cutoff_pass = sum(area_pass_cutoff, na.rm = TRUE),
      mean_core_area = mean(polygon_area, na.rm = TRUE),
      median_core_area = median(polygon_area, na.rm = TRUE),
      area_q25 = as.numeric(stats::quantile(polygon_area, 0.25, na.rm = TRUE)),
      area_q75 = as.numeric(stats::quantile(polygon_area, 0.75, na.rm = TRUE)),
      area_q90 = as.numeric(stats::quantile(polygon_area, 0.90, na.rm = TRUE)),
      max_core_area = max(polygon_area, na.rm = TRUE)
    )

  message("========== Necrotic-core-like candidate area distribution ==========")
  print(area_summary)

  message("========== Necrotic-core-like candidate area distribution by FOV ==========")
  print(
    cluster_summary %>%
      dplyr::group_by(FOV, severity) %>%
      dplyr::summarise(
        n_candidate_core_regions = dplyr::n(),
        n_area_cutoff_pass = sum(area_pass_cutoff, na.rm = TRUE),
        mean_core_area = mean(polygon_area, na.rm = TRUE),
        median_core_area = median(polygon_area, na.rm = TRUE),
        max_core_area = max(polygon_area, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(max_core_area)) %>%
      as.data.frame()
  )

  return(area_summary)
}

# ============================================================
# 7. 正式运行：读取 sc_ref 和 Xenium
# ============================================================

sc_ref <- prepare_sc_ref_for_transfer(
  sc_ref_path = file.path(SC_ROOT, "Result", "sub_integrated_data_Final.rds"),
  ref_label_source_col = "Celltype_raw"
)

xen_path <- file.path(SPATIAL_ROOT, "GSE315246_xenium.obj.integrated.rds")

if (!file.exists(xen_path)) {
  stop("找不到 Xenium 对象: ", xen_path)
}

message("Loading Xenium object: ", xen_path)
xen <- readRDS(xen_path)

if (!inherits(xen, "Seurat")) {
  stop("xen 不是 Seurat object。当前 class: ", paste(class(xen), collapse = ", "))
}

xen <- repair_fov_slots_safely(
  xen,
  object_name = "xen"
)

xen <- set_default_assay_safely(
  xen,
  preferred = c("Xenium", "RNA", "Spatial", "SCT")
)

# 如果原对象只有 disease 而没有 severity，自动复制一列 severity
if (!"severity" %in% colnames(xen@meta.data) && "disease" %in% colnames(xen@meta.data)) {
  xen@meta.data$severity <- as.character(xen@meta.data$disease)
}

# ============================================================
# 8. 选择 Xenium 中需要映射的细胞
# ============================================================

all_xen_cells <- Cells(xen)

if ("predicted.id" %in% colnames(xen@meta.data)) {

  predicted_id_vec <- as.character(xen@meta.data$predicted.id)
  names(predicted_id_vec) <- rownames(xen@meta.data)

  xen_myeloid_cells <- names(predicted_id_vec)[predicted_id_vec == "Myeloid"]
  xen_myeloid_cells <- intersect(xen_myeloid_cells, all_xen_cells)

  message("Xenium predicted.id found.")
  message("Total Xenium cells: ", length(all_xen_cells))
  message("Xenium Myeloid cells for transfer: ", length(xen_myeloid_cells))

  if (length(xen_myeloid_cells) < 20) {
    warning(
      "Xenium 中 predicted.id == 'Myeloid' 的细胞少于 20 个；",
      "将退回到对所有细胞做 transfer，但非髓系细胞会被强行分到 Mo/Ma 亚型。"
    )
    xen_myeloid_cells <- all_xen_cells
  }

} else {

  xen_myeloid_cells <- all_xen_cells

  warning(
    "xen@meta.data 中没有 predicted.id；",
    "将对所有 Xenium 细胞做 Mo/Ma label transfer。",
    "注意：非髓系细胞也会被强行分配为 Mo/Ma 亚型。"
  )
}

# ============================================================
# 9. 记录注释前使用哪一列
# ============================================================

before_col <- choose_first_meta_col(
  xen,
  c(
    "myeloid_substate",
    "xen_myeloid_substate",
    "Xenium_myeloid_substate",
    "dominant_myeloid_program",
    "predicted.id",
    "celltype",
    "Celltype",
    "seurat_clusters"
  )
)

if (is.na(before_col)) {
  stop("找不到可用于注释前展示的 metadata 列。请检查 colnames(xen@meta.data)。")
}

message("Before annotation column: ", before_col)

# ============================================================
# 10. 执行 label transfer
# ============================================================

transfer_res <- transfer_moma_labels_to_spatial_keep_fov(
  ref_obj = sc_ref,
  query_obj = xen,
  ref_label_col = "Celltype_transfer",
  query_assay = DefaultAssay(xen),
  out_prefix = "Xenium_MoMa",
  dims_use = 1:30,
  restrict_query_cells = xen_myeloid_cells
)

xen <- transfer_res$query_obj
xen_transfer_pred <- transfer_res$pred

after_col <- "Xenium_MoMa_scPred_Celltype"
score_col <- "Xenium_MoMa_scPred_score"

# 统一转移后标签顺序
xen@meta.data[[after_col]] <- factor(
  as.character(xen@meta.data[[after_col]]),
  levels = moma_levels
)

# ============================================================
# 11. 统计 APOBEC3A 检出
# ============================================================

xen_counts <- get_assay_data_safe(
  xen,
  assay = DefaultAssay(xen),
  layer = "counts"
)

if (TARGET_GENE %in% rownames(xen_counts)) {
  xen$APOBEC3A_count_transfer <- as.numeric(xen_counts[TARGET_GENE, Cells(xen)])
  xen$APOBEC3A_detected_transfer <- xen$APOBEC3A_count_transfer > 0
} else {
  xen$APOBEC3A_count_transfer <- NA_real_
  xen$APOBEC3A_detected_transfer <- NA
  warning("Xenium counts 中没有 ", TARGET_GENE)
}

# ============================================================
# 12. 保存转移结果摘要，不保存完整 xen 对象
# ============================================================

xen_meta_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_label_transfer_metadata.csv")
)

if (SAVE_FULL_METADATA) {
  write.csv(
    xen@meta.data,
    xen_meta_out,
    row.names = TRUE
  )
  message("Saved full Xenium transfer metadata: ", xen_meta_out)
} else {
  message("Skipped full Xenium transfer metadata; set AST_SAVE_FULL_METADATA=1 to write it.")
}

transfer_meta_cols <- unique(c(
  "severity",
  "disease",
  before_col,
  after_col,
  score_col,
  grep("^Xenium_MoMa_prediction.score.", colnames(xen@meta.data), value = TRUE),
  "APOBEC3A_count_transfer",
  "APOBEC3A_detected_transfer"
))

transfer_meta_cols <- intersect(transfer_meta_cols, colnames(xen@meta.data))

xen_transfer_metadata <- xen@meta.data[, transfer_meta_cols, drop = FALSE] %>%
  tibble::rownames_to_column("cell")

xen_transfer_meta_csv <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_transfer_metadata_only.csv")
)

xen_transfer_meta_rds <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_transfer_metadata_only.rds")
)

if (SAVE_CELL_LEVEL_TABLES) {
  write.csv(
    xen_transfer_metadata,
    xen_transfer_meta_csv,
    row.names = FALSE
  )
  message("Saved compact transfer metadata CSV: ", xen_transfer_meta_csv)
} else {
  message("Skipped compact transfer metadata CSV; set AST_SAVE_CELL_LEVEL_TABLES=1 to write it.")
}

if (SAVE_TRANSFER_RDS) {
  saveRDS(
    xen_transfer_metadata,
    xen_transfer_meta_rds
  )

  saveRDS(
    xen_transfer_pred,
    file.path(TABLE_DIR, paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_transfer_prediction_raw.rds"))
  )

  message("Saved transfer metadata RDS: ", xen_transfer_meta_rds)
} else {
  message("Skipped transfer RDS outputs; set AST_SAVE_TRANSFER_RDS=1 to write them.")
}

# ============================================================
# 13. 输出整体统计表
# ============================================================

xen_transfer_summary <- xen@meta.data %>%
  mutate(
    Xenium_MoMa_scPred_Celltype = as.character(.data[[after_col]]),
    Xenium_MoMa_scPred_score = as.numeric(.data[[score_col]]),
    high_confidence = Xenium_MoMa_scPred_score >= 0.5
  ) %>%
  filter(!is.na(Xenium_MoMa_scPred_Celltype)) %>%
  group_by(Xenium_MoMa_scPred_Celltype, high_confidence) %>%
  summarise(
    n_cells = n(),
    median_prediction_score = median(Xenium_MoMa_scPred_score, na.rm = TRUE),
    mean_prediction_score = mean(Xenium_MoMa_scPred_score, na.rm = TRUE),
    n_APOBEC3A_positive = sum(APOBEC3A_detected_transfer, na.rm = TRUE),
    APOBEC3A_detection_rate = mean(APOBEC3A_detected_transfer, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(high_confidence), desc(n_cells))

xen_summary_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_label_transfer_summary.csv")
)

write.csv(
  xen_transfer_summary,
  xen_summary_out,
  row.names = FALSE
)

message("Saved Xenium transfer summary: ", xen_summary_out)
print(xen_transfer_summary)

message("Transferred myeloid cell type counts:")
print(table(xen@meta.data[[after_col]], useNA = "ifany"))

message("Transferred myeloid cell type counts by confidence:")
print(table(
  Pred = xen@meta.data[[after_col]],
  HighConf = as.numeric(xen@meta.data[[score_col]]) >= 0.5,
  useNA = "ifany"
))

# ============================================================
# 14. UMAP 对照图
# ============================================================

p_compare_all <- plot_before_after_dim(
  obj = xen,
  before_col = before_col,
  after_col = after_col,
  title_prefix = "Xenium cell annotation overview",
  cells = NULL,
  filename_base = paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_before_after_all_cells_UMAP_publication"),
  width = 14,
  height = 6
)

transferred_cells <- rownames(xen@meta.data)[!is.na(xen@meta.data[[after_col]])]

p_compare_transferred <- plot_before_after_dim(
  obj = xen,
  before_col = before_col,
  after_col = after_col,
  title_prefix = "Xenium myeloid cell annotation overview",
  cells = transferred_cells,
  filename_base = paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_before_after_transferred_cells_UMAP_publication"),
  width = 14,
  height = 6
)

p_confidence <- plot_transfer_score_dim(
  obj = xen,
  score_col = score_col,
  title = "Myeloid annotation transfer confidence",
  filename_base = paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_MoMa_transfer_confidence_UMAP_publication"),
  width = 7,
  height = 6
)

# ============================================================
# 15. 所有 FOV 空间 before/after 图
# ============================================================

all_fov_summary <- plot_all_xenium_fovs_with_severity_subtitle(
  obj = xen,
  before_col = before_col,
  after_col = after_col,
  score_col = score_col,
  only_transferred_after = TRUE,
  max_points = 200000
)

message("========== All FOV summary ==========")
print(all_fov_summary)

message("========== FOVs ranked by transferred myeloid cells ==========")
print(
  all_fov_summary %>%
    arrange(desc(n_transferred))
)

message("========== FOVs ranked by APOBEC3A+ transferred myeloid cells ==========")
print(
  all_fov_summary %>%
    arrange(desc(n_APOBEC3A_positive_transferred))
)

# ============================================================
# 16. 不同进展时期 Monocyte / Macrophage / LAM/Foam Cell 比例堆叠柱状图
# ============================================================

composition_res <- plot_moma_composition_by_severity(
  obj = xen,
  after_col = after_col,
  score_col = score_col,
  use_high_confidence_only = FALSE,
  confidence_cutoff = 0.5,
  filename_base = paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity_publication"),
  width = 9,
  height = 6.5
)

composition_res_highconf <- plot_moma_composition_by_severity(
  obj = xen,
  after_col = after_col,
  score_col = score_col,
  use_high_confidence_only = TRUE,
  confidence_cutoff = 0.5,
  filename_base = paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity_highconf_publication"),
  width = 9,
  height = 6.5
)


# ============================================================
# 16.1 坏死核心样区域识别、诊断、细胞标记和 FOV 面积柱状图
#     已有图片文件名不改；以下均为新增输出
# ============================================================

# 统一选择需要分析的 FOV
if (is.null(NC_IMAGE_NAMES)) {
  NC_IMAGE_NAMES_USE <- names(xen@images)
} else {
  NC_IMAGE_NAMES_USE <- intersect(NC_IMAGE_NAMES, names(xen@images))
  if (length(NC_IMAGE_NAMES_USE) == 0) {
    stop("NC_IMAGE_NAMES 与 names(xen@images) 没有交集。请检查 FOV 名称。")
  }
}

print_nc_parameters(
  local_params = NC_LOCAL_PARAMS,
  cluster_params = NC_CLUSTER_PARAMS,
  area_cutoff = NC_CORE_AREA_CUTOFF,
  image_names = NC_IMAGE_NAMES_USE,
  random_seed = NC_RANDOM_SEED
)

# ------------------------------------------------------------
# 16.1.1 诊断 Step 1：local seed 是否存在
# ------------------------------------------------------------
if (isTRUE(NC_RUN_THRESHOLD_DIAGNOSTICS)) {
  nc_diag <- do.call(
    diagnose_necrotic_core_like_thresholds_all_fovs,
    c(
      list(
        obj = xen,
        after_col = after_col,
        score_col = score_col,
        image_names = NC_IMAGE_NAMES_USE,
        seed = NC_RANDOM_SEED,
        output_prefix = paste0(OUT_PREFIX, "_necrotic_core_like_threshold_diagnostics")
      ),
      NC_LOCAL_PARAMS
    )
  )
} else {
  nc_diag <- NULL
  message("Skipped local threshold diagnostics because NC_RUN_THRESHOLD_DIAGNOSTICS = FALSE")
}

# ------------------------------------------------------------
# 16.1.2 诊断 Step 2：DBSCAN / polygon / cluster filter 哪一步筛没
# ------------------------------------------------------------
if (isTRUE(NC_RUN_DBSCAN_POLYGON_DIAGNOSTICS)) {
  nc_dbscan_diag <- do.call(
    diagnose_nc_dbscan_polygon_all_fovs,
    c(
      list(
        obj = xen,
        after_col = after_col,
        score_col = score_col,
        image_names = NC_IMAGE_NAMES_USE,
        seed = NC_RANDOM_SEED,
        output_prefix = paste0(OUT_PREFIX, "_necrotic_core_like_DBSCAN_polygon_diagnostics")
      ),
      NC_LOCAL_PARAMS,
      NC_CLUSTER_PARAMS
    )
  )
} else {
  nc_dbscan_diag <- NULL
  message("Skipped DBSCAN / polygon diagnostics because NC_RUN_DBSCAN_POLYGON_DIAGNOSTICS = FALSE")
}

# ------------------------------------------------------------
# 16.1.3 正式识别：写回 NC_like_region / NC_region_id 等 metadata
# ------------------------------------------------------------
if (isTRUE(NC_RUN_FINAL_DETECTION)) {
  nc_res <- do.call(
    detect_necrotic_core_like_all_fovs,
    c(
      list(
        obj = xen,
        after_col = after_col,
        score_col = score_col,
        image_names = NC_IMAGE_NAMES_USE,
        min_core_area_cutoff = NC_CORE_AREA_CUTOFF,
        polygon_buffer = NC_POLYGON_BUFFER,
        seed = NC_RANDOM_SEED
      ),
      NC_LOCAL_PARAMS,
      NC_CLUSTER_PARAMS
    )
  )

  xen <- nc_res$obj

} else {
  stop("NC_RUN_FINAL_DETECTION = FALSE。当前脚本后续步骤需要 nc_res；如只想诊断，请把后续绘图和 A3A 统计注释掉。")
}

# ------------------------------------------------------------
# 16.1.4 保存诊断和识别结果
# ------------------------------------------------------------

# 保存 FOV-level summary
nc_fov_summary_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_necrotic_core_like_FOV_summary.csv")
)
write.csv(nc_res$fov_summary, nc_fov_summary_out, row.names = FALSE)

# 保存 core-region-level summary，包括低于/高于 area cutoff 的候选区域
nc_cluster_summary_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_necrotic_core_like_cluster_summary.csv")
)
write.csv(nc_res$cluster_summary, nc_cluster_summary_out, row.names = FALSE)

# 保存 area distribution summary，方便决定 NC_CORE_AREA_CUTOFF
nc_core_area_summary <- summarise_nc_core_area_distribution(nc_res$cluster_summary)
nc_core_area_summary_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_necrotic_core_like_area_distribution_summary.csv")
)
write.csv(nc_core_area_summary, nc_core_area_summary_out, row.names = FALSE)

# 保存带坏死核心样区域标记的 metadata
nc_meta_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_xenium_metadata_with_necrotic_core_like.csv")
)
if (SAVE_CELL_LEVEL_TABLES || SAVE_FULL_METADATA) {
  write.csv(xen@meta.data, nc_meta_out, row.names = TRUE)
  message("Saved Xenium metadata with necrotic-core-like columns: ", nc_meta_out)
} else {
  message("Skipped Xenium metadata with necrotic-core-like columns; set AST_SAVE_CELL_LEVEL_TABLES=1 to write it.")
}

message("Saved necrotic-core-like FOV summary: ", nc_fov_summary_out)
message("Saved necrotic-core-like cluster summary: ", nc_cluster_summary_out)
message("Saved necrotic-core-like area distribution summary: ", nc_core_area_summary_out)

message("========== Necrotic-core-like FOV summary ==========")
print(nc_res$fov_summary)

message("========== Necrotic-core-like region summary ==========")
print(nc_res$cluster_summary)

# ------------------------------------------------------------
# 16.1.5 所有 FOV 的坏死核心样区域空间图
# ------------------------------------------------------------
if (isTRUE(NC_PLOT_ALL_FOV_NC)) {
  for (image_use in names(nc_res$per_fov)) {
    tryCatch(
      {
        plot_one_fov_necrotic_core_like(
          det_obj = nc_res,
          image_use = image_use,
          after_col = after_col,
          save_plot = TRUE,
          filename_base = paste0(
            OUT_PREFIX, "_",
            sanitize_filename(image_use),
            "_necrotic_core_like_publication"
          ),
          width = 13,
          height = 6
        )
      },
      error = function(e) {
        warning("Necrotic-core-like FOV plot failed for ", image_use, ": ", e$message)
      }
    )
  }
} else {
  message("Skipped all FOV NC plots because NC_PLOT_ALL_FOV_NC = FALSE")
}

# 新增：FOV 最大坏死核心样区域面积柱状图
p_nc_fov_area <- plot_fov_max_core_area_barplot(
  fov_summary = nc_res$fov_summary,
  area_cutoff = NC_CORE_AREA_CUTOFF,
  filename_base = paste0(OUT_PREFIX, "_necrotic_core_like_FOV_max_area_barplot_publication"),
  width = 11,
  height = 6.5
)

# 新增：A3A 在坏死核心样区域内外的简要统计表
nc_a3a_summary <- xen@meta.data %>%
  dplyr::filter(!is.na(.data[[after_col]])) %>%
  dplyr::mutate(
    transferred_celltype = as.character(.data[[after_col]]),
    NC_like_region = as.logical(NC_like_region),
    NC_FOV_has_core = as.logical(NC_FOV_has_core)
  ) %>%
  dplyr::group_by(NC_like_region, transferred_celltype) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_APOBEC3A_positive = sum(APOBEC3A_detected_transfer, na.rm = TRUE),
    APOBEC3A_detection_rate = mean(APOBEC3A_detected_transfer, na.rm = TRUE),
    mean_APOBEC3A_count = mean(APOBEC3A_count_transfer, na.rm = TRUE),
    median_APOBEC3A_count = median(APOBEC3A_count_transfer, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(NC_like_region), transferred_celltype)

nc_a3a_summary_out <- file.path(
  TABLE_DIR,
  paste0(OUT_PREFIX, "_APOBEC3A_by_necrotic_core_like_region_summary.csv")
)
write.csv(nc_a3a_summary, nc_a3a_summary_out, row.names = FALSE)
message("Saved APOBEC3A by necrotic-core-like region summary: ", nc_a3a_summary_out)
print(nc_a3a_summary)

# 新增图片路径汇总，用于第 18 节 manifest 合并
nc_fov_plot_manifest <- tibble(
  plot_name = paste0("Necrotic-core-like FOV annotation: ", nc_res$fov_summary$FOV),
  pdf = file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX, "_",
      sanitize_filename(nc_res$fov_summary$FOV),
      "_necrotic_core_like_publication.pdf"
    )
  ),
  png = file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX, "_",
      sanitize_filename(nc_res$fov_summary$FOV),
      "_necrotic_core_like_publication.png"
    )
  )
)

nc_barplot_manifest <- tibble(
  plot_name = "FOV maximum necrotic-core-like area barplot",
  pdf = file.path(
    FIG_DIR,
    paste0(OUT_PREFIX, "_necrotic_core_like_FOV_max_area_barplot_publication.pdf")
  ),
  png = file.path(
    FIG_DIR,
    paste0(OUT_PREFIX, "_necrotic_core_like_FOV_max_area_barplot_publication.png")
  )
)

nc_plot_manifest <- bind_rows(nc_barplot_manifest, nc_fov_plot_manifest)

# ============================================================
# 16.2 APOBEC3A 空间热图和疾病进展柱状图
#     新增输出，不改动已有图片文件名
# ============================================================

if (is.null(A3A_HEATMAP_IMAGE_NAMES)) {
  A3A_HEATMAP_IMAGE_NAMES_USE <- names(xen@images)
} else {
  A3A_HEATMAP_IMAGE_NAMES_USE <- intersect(A3A_HEATMAP_IMAGE_NAMES, names(xen@images))
  if (length(A3A_HEATMAP_IMAGE_NAMES_USE) == 0) {
    warning("A3A_HEATMAP_IMAGE_NAMES 与 names(xen@images) 没有交集，跳过 A3A 空间热图。")
  }
}

a3a_heatmap_manifest <- tibble(
  plot_name = character(),
  pdf = character(),
  png = character()
)

a3a_heatmap_core_manifest <- tibble(
  plot_name = character(),
  pdf = character(),
  png = character()
)

if (exists("A3A_HEATMAP_IMAGE_NAMES_USE") && length(A3A_HEATMAP_IMAGE_NAMES_USE) > 0) {

  if (isTRUE(A3A_RUN_SPATIAL_HEATMAPS)) {
    for (image_use in A3A_HEATMAP_IMAGE_NAMES_USE) {
      tryCatch(
        {
          plot_one_fov_APOBEC3A_heatmap(
            obj = xen,
            image_use = image_use,
            after_col = after_col,
            a3a_count_col = "APOBEC3A_count_transfer",
            save_plot = TRUE,
            filename_base = paste0(
              OUT_PREFIX, "_",
              sanitize_filename(image_use),
              "_APOBEC3A_spatial_heatmap"
            ),
            width = 7,
            height = 6,
            use_log10 = TRUE,
            show_only_positive = TRUE
          )
        },
        error = function(e) {
          warning("APOBEC3A spatial heatmap failed for ", image_use, ": ", e$message)
        }
      )
    }

    a3a_heatmap_manifest <- tibble(
      plot_name = paste0("APOBEC3A spatial heatmap: ", A3A_HEATMAP_IMAGE_NAMES_USE),
      pdf = file.path(
        FIG_DIR,
        paste0(
          OUT_PREFIX, "_",
          sanitize_filename(A3A_HEATMAP_IMAGE_NAMES_USE),
          "_APOBEC3A_spatial_heatmap.pdf"
        )
      ),
      png = file.path(
        FIG_DIR,
        paste0(
          OUT_PREFIX, "_",
          sanitize_filename(A3A_HEATMAP_IMAGE_NAMES_USE),
          "_APOBEC3A_spatial_heatmap.png"
        )
      )
    )
  } else {
    message("Skipped APOBEC3A spatial heatmaps because A3A_RUN_SPATIAL_HEATMAPS = FALSE")
  }

  if (isTRUE(A3A_RUN_SPATIAL_HEATMAPS_WITH_CORE)) {
    for (image_use in A3A_HEATMAP_IMAGE_NAMES_USE) {
      tryCatch(
        {
          plot_one_fov_APOBEC3A_heatmap_with_core(
            obj = xen,
            det_obj = nc_res,
            image_use = image_use,
            after_col = after_col,
            a3a_count_col = "APOBEC3A_count_transfer",
            save_plot = TRUE,
            filename_base = paste0(
              OUT_PREFIX, "_",
              sanitize_filename(image_use),
              "_APOBEC3A_spatial_heatmap_with_core"
            ),
            width = 7,
            height = 6,
            use_log10 = TRUE,
            show_only_positive = TRUE,
            core_outline_color = "grey35",
            core_outline_linetype = "dashed",
            core_outline_linewidth = 0.7
          )
        },
        error = function(e) {
          warning("APOBEC3A spatial heatmap with core failed for ", image_use, ": ", e$message)
        }
      )
    }

    a3a_heatmap_core_manifest <- tibble(
      plot_name = paste0("APOBEC3A spatial heatmap with core boundary: ", A3A_HEATMAP_IMAGE_NAMES_USE),
      pdf = file.path(
        FIG_DIR,
        paste0(
          OUT_PREFIX, "_",
          sanitize_filename(A3A_HEATMAP_IMAGE_NAMES_USE),
          "_APOBEC3A_spatial_heatmap_with_core.pdf"
        )
      ),
      png = file.path(
        FIG_DIR,
        paste0(
          OUT_PREFIX, "_",
          sanitize_filename(A3A_HEATMAP_IMAGE_NAMES_USE),
          "_APOBEC3A_spatial_heatmap_with_core.png"
        )
      )
    )
  } else {
    message("Skipped APOBEC3A spatial heatmaps with core because A3A_RUN_SPATIAL_HEATMAPS_WITH_CORE = FALSE")
  }
}

a3a_disease_plot_manifest <- tibble(
  plot_name = character(),
  pdf = character(),
  png = character()
)

if (isTRUE(A3A_RUN_DISEASE_BARPLOTS)) {
  a3a_overall_pair_res <- plot_A3A_overall_expression_and_positive_rate_by_disease(
    obj = xen,
    after_col = after_col,
    a3a_count_col = "APOBEC3A_count_transfer",
    score_col = score_col,
    use_high_confidence_only = FALSE,
    confidence_cutoff = 0.5,
    filename_base = paste0(
      OUT_PREFIX,
      "_APOBEC3A_overall_myeloid_expression_and_positive_fraction_by_disease_stage_publication"
    ),
    width = 7,
    height = 5.5
  )

  a3a_disease_plot_manifest <- tibble(
    plot_name = c(
      "APOBEC3A expression in all myeloid cells across disease stages",
      "APOBEC3A-positive myeloid fraction across disease stages"
    ),
    pdf = file.path(
      FIG_DIR,
      paste0(c(
        paste0(OUT_PREFIX, "_APOBEC3A_overall_myeloid_expression_and_positive_fraction_by_disease_stage_publication_expression_barplot"),
        paste0(OUT_PREFIX, "_APOBEC3A_overall_myeloid_expression_and_positive_fraction_by_disease_stage_publication_positive_fraction_barplot")
      ), ".pdf")
    ),
    png = file.path(
      FIG_DIR,
      paste0(c(
        paste0(OUT_PREFIX, "_APOBEC3A_overall_myeloid_expression_and_positive_fraction_by_disease_stage_publication_expression_barplot"),
        paste0(OUT_PREFIX, "_APOBEC3A_overall_myeloid_expression_and_positive_fraction_by_disease_stage_publication_positive_fraction_barplot")
      ), ".png")
    )
  )
} else {
  message("Skipped APOBEC3A disease-stage barplots because A3A_RUN_DISEASE_BARPLOTS = FALSE")
}

a3a_extra_plot_manifest <- bind_rows(
  a3a_heatmap_manifest,
  a3a_heatmap_core_manifest,
  a3a_disease_plot_manifest
)


# ============================================================
# 17. APOBEC3A+ 在新注释里的分布
# ============================================================

xen_a3a_transfer_summary <- xen@meta.data %>%
  filter(!is.na(.data[[after_col]])) %>%
  mutate(
    transferred_celltype = as.character(.data[[after_col]]),
    high_confidence = as.numeric(.data[[score_col]]) >= 0.5
  ) %>%
  group_by(transferred_celltype, high_confidence) %>%
  summarise(
    n_cells = n(),
    n_APOBEC3A_positive = sum(APOBEC3A_detected_transfer, na.rm = TRUE),
    APOBEC3A_detection_rate = mean(APOBEC3A_detected_transfer, na.rm = TRUE),
    median_prediction_score = median(as.numeric(.data[[score_col]]), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(APOBEC3A_detection_rate), desc(n_cells))

write.csv(
  xen_a3a_transfer_summary,
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_transfer_a3a_summary.csv")),
  row.names = FALSE
)

print(xen_a3a_transfer_summary)

p_a3a_by_new_label <- xen_a3a_transfer_summary %>%
  filter(high_confidence) %>%
  mutate(
    transferred_celltype = factor(
      transferred_celltype,
      levels = moma_levels
    )
  ) %>%
  ggplot(aes(x = transferred_celltype, y = APOBEC3A_detection_rate, fill = transferred_celltype)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = cell_colors[moma_levels],
    breaks = moma_levels,
    drop = FALSE
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "APOBEC3A-positive fraction",
    title = "APOBEC3A-positive fraction across myeloid annotations"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

print(p_a3a_by_new_label)

save_plot_both(
  p_a3a_by_new_label,
  paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_APOBEC3A_fraction_by_transferred_MoMa_state_publication"),
  width = 7,
  height = 5
)

# ============================================================
# 18. 保存图片路径汇总
# ============================================================

fov_plot_manifest <- tibble(
  plot_name = paste0("FOV publication annotation: ", all_fov_summary$FOV),
  pdf = file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX, "_",
      sanitize_filename(all_fov_summary$FOV),
      "_FOV_publication_annotation.pdf"
    )
  ),
  png = file.path(
    FIG_DIR,
    paste0(
      OUT_PREFIX, "_",
      sanitize_filename(all_fov_summary$FOV),
      "_FOV_publication_annotation.png"
    )
  )
)

base_plot_manifest <- tibble(
  plot_name = c(
    "All cells UMAP before vs after",
    "Transferred myeloid cells UMAP before vs after",
    "Transfer confidence UMAP",
    "Myeloid composition by severity",
    "Myeloid composition by severity high confidence",
    "APOBEC3A fraction by transferred myeloid state"
  ),
  pdf = file.path(
    FIG_DIR,
    paste0(c(
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_before_after_all_cells_UMAP_publication"),
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_before_after_transferred_cells_UMAP_publication"),
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_MoMa_transfer_confidence_UMAP_publication"),
      paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity_publication"),
      paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity_highconf_publication"),
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_APOBEC3A_fraction_by_transferred_MoMa_state_publication")
    ), ".pdf")
  ),
  png = file.path(
    FIG_DIR,
    paste0(c(
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_before_after_all_cells_UMAP_publication"),
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_before_after_transferred_cells_UMAP_publication"),
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_MoMa_transfer_confidence_UMAP_publication"),
      paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity_publication"),
      paste0(OUT_PREFIX, "_Xenium_MoMa_composition_by_severity_highconf_publication"),
      paste0(OUT_PREFIX, "_Xenium_CelltypeRaw_APOBEC3A_fraction_by_transferred_MoMa_state_publication")
    ), ".png")
  )
)

empty_plot_manifest <- tibble(
  plot_name = character(),
  pdf = character(),
  png = character()
)

plot_manifest <- bind_rows(
  base_plot_manifest,
  fov_plot_manifest,
  if (exists("nc_plot_manifest")) nc_plot_manifest else empty_plot_manifest,
  if (exists("a3a_extra_plot_manifest")) a3a_extra_plot_manifest else empty_plot_manifest
) %>%
  mutate(
    pdf_exists = file.exists(pdf),
    png_exists = file.exists(png)
  )

write.csv(
  plot_manifest,
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_xenium_moma_CelltypeRaw_transfer_plot_manifest_publication.csv")),
  row.names = FALSE
)

message("========== Plot manifest ==========")
print(plot_manifest)

message("========== Xenium Mo/Ma Celltype_raw publication-ready label transfer finished ==========")
# ============================================================
# GeoMx WTA ROI-level APOBEC3A expression
# 合并所有 CD45/CD4 subset，不分 facet，不分 marker
# 只比较 Plaque 和 Adventitia ROI
# 按 Normal / Mild / Moderate / Severe 顺序作图
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# 0. 路径设置
# ------------------------------------------------------------

ROOT <- env_path("AST_ROOT", ROOT)
OUTPUT_ROOT <- env_path("AST_OUTPUT_ROOT", OUTPUT_ROOT)
TABLE_DIR <- file.path(OUTPUT_ROOT, "tables")
FIG_DIR <- file.path(OUTPUT_ROOT, "figures")
OUT_PREFIX <- env_path("AST_OUT_PREFIX", OUT_PREFIX)

GEOMX_OUT_DIR <- file.path(FIG_DIR, "GeoMx_WTA_ROI_level_myeloid_A3A")
dir.create(GEOMX_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. 自动寻找 GeoMx ROI 表
# ------------------------------------------------------------

roi_table_candidates <- c(
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_ROI_level_myeloid_program_scores.csv")),
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_ROI_clean_for_A3A_myeloid_program.csv")),
  file.path(TABLE_DIR, "geomx_wta_apobec_roi_counts.csv"),
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_roi_a3a_by_grade_localisation_subset.csv")),
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_ROI_A3A_by_grade_localisation_subset.csv"))
)

ROI_TABLE <- roi_table_candidates[file.exists(roi_table_candidates)][1]

if (is.na(ROI_TABLE) || !file.exists(ROI_TABLE)) {
  stop(
    "没有找到 GeoMx ROI 表。\n请确认以下文件至少存在一个：\n",
    paste(roi_table_candidates, collapse = "\n")
  )
}

message("Using GeoMx ROI table: ", ROI_TABLE)

geomx_roi <- read.csv(
  ROI_TABLE,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()

message("ROI table columns:")
print(colnames(geomx_roi))

# ------------------------------------------------------------
# 2. 兼容不同列名
# ------------------------------------------------------------

# APOBEC3A_count
if (!"APOBEC3A_count" %in% colnames(geomx_roi)) {
  a3a_candidates <- grep(
    "APOBEC3A.*count|^APOBEC3A$|A3A.*count",
    colnames(geomx_roi),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(a3a_candidates) > 0) {
    geomx_roi$APOBEC3A_count <- as.numeric(geomx_roi[[a3a_candidates[1]]])
    message("Use column as APOBEC3A_count: ", a3a_candidates[1])
  } else {
    stop("ROI 表中找不到 APOBEC3A_count 或 APOBEC3A 相关列。")
  }
}

# grade
if (!"grade" %in% colnames(geomx_roi)) {
  grade_candidates <- grep(
    "grade|severity|lesion",
    colnames(geomx_roi),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(grade_candidates) > 0) {
    geomx_roi$grade <- as.character(geomx_roi[[grade_candidates[1]]])
    message("Use column as grade: ", grade_candidates[1])
  } else {
    stop("ROI 表中找不到 grade / severity / lesion 相关列。")
  }
}

# localisation
if (!"localisation" %in% colnames(geomx_roi)) {
  loc_candidates <- grep(
    "localisation|localization|location|region",
    colnames(geomx_roi),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(loc_candidates) > 0) {
    geomx_roi$localisation <- as.character(geomx_roi[[loc_candidates[1]]])
    message("Use column as localisation: ", loc_candidates[1])
  } else {
    stop("ROI 表中找不到 localisation / localization / location / region 相关列。")
  }
}

# roi ID
if (!"roi" %in% colnames(geomx_roi)) {
  roi_candidates <- grep(
    "^roi$|ROI|AOI|sample|gsm|dcc",
    colnames(geomx_roi),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(roi_candidates) > 0) {
    geomx_roi$roi <- as.character(geomx_roi[[roi_candidates[1]]])
    message("Use column as roi: ", roi_candidates[1])
  } else {
    geomx_roi$roi <- paste0("ROI_", seq_len(nrow(geomx_roi)))
  }
}

# positive
if (!"APOBEC3A_positive" %in% colnames(geomx_roi)) {
  geomx_roi$APOBEC3A_positive <- as.numeric(geomx_roi$APOBEC3A_count) > 0
}

# ------------------------------------------------------------
# 3. 清洗数据
# ------------------------------------------------------------

grade_levels <- c("Normal", "Mild", "Moderate", "Severe")
localisation_keep <- c("Plaque", "Adventitia")

geomx_roi_clean <- geomx_roi %>%
  mutate(
    roi = as.character(roi),
    grade = as.character(grade),
    localisation = as.character(localisation),
    APOBEC3A_count = as.numeric(APOBEC3A_count),
    APOBEC3A_positive = as.logical(APOBEC3A_positive),
    APOBEC3A_log1p = log1p(APOBEC3A_count),

    grade = case_when(
      grepl("^normal$", grade, ignore.case = TRUE) ~ "Normal",
      grepl("^mild$", grade, ignore.case = TRUE) ~ "Mild",
      grepl("^moderate$", grade, ignore.case = TRUE) ~ "Moderate",
      grepl("^severe$", grade, ignore.case = TRUE) ~ "Severe",
      TRUE ~ grade
    ),

    localisation = case_when(
      grepl("^plaque$", localisation, ignore.case = TRUE) ~ "Plaque",
      grepl("^adventitia$", localisation, ignore.case = TRUE) ~ "Adventitia",
      TRUE ~ localisation
    )
  ) %>%
  filter(
    grade %in% grade_levels,
    localisation %in% localisation_keep,
    !is.na(APOBEC3A_count)
  ) %>%
  mutate(
    grade = factor(grade, levels = grade_levels),
    localisation = factor(localisation, levels = localisation_keep)
  )

if (nrow(geomx_roi_clean) == 0) {
  stop(
    "过滤后没有剩余 ROI。\n",
    "请检查 grade 是否包含 Normal/Mild/Moderate/Severe，",
    "localisation 是否包含 Plaque/Adventitia。"
  )
}

# ------------------------------------------------------------
# 4. 输出用于作图的数据
# ------------------------------------------------------------

write.csv(
  geomx_roi_clean,
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_plot_data.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. 统计 summary
# ------------------------------------------------------------

geomx_summary <- geomx_roi_clean %>%
  group_by(grade, localisation) %>%
  summarise(
    n_roi = n(),
    n_APOBEC3A_positive = sum(APOBEC3A_positive, na.rm = TRUE),
    APOBEC3A_detection_rate = mean(APOBEC3A_positive, na.rm = TRUE),
    mean_APOBEC3A_count = mean(APOBEC3A_count, na.rm = TRUE),
    median_APOBEC3A_count = median(APOBEC3A_count, na.rm = TRUE),
    mean_APOBEC3A_log1p = mean(APOBEC3A_log1p, na.rm = TRUE),
    median_APOBEC3A_log1p = median(APOBEC3A_log1p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(localisation, grade)

write.csv(
  geomx_summary,
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_summary.csv")),
  row.names = FALSE
)

print(geomx_summary)

# ------------------------------------------------------------
# 6. 统计检验
# ------------------------------------------------------------

stat_list <- list()

# 6.1 所有 Plaque + Adventitia ROI 合并后，不同 grade 比较
if (dplyr::n_distinct(geomx_roi_clean$grade) >= 2) {
  stat_list[[length(stat_list) + 1]] <- tibble(
    comparison = "All Plaque and Adventitia ROIs",
    test = "Kruskal-Wallis",
    variable = "APOBEC3A_log1p by lesion grade",
    p_value = tryCatch(
      kruskal.test(APOBEC3A_log1p ~ grade, data = geomx_roi_clean)$p.value,
      error = function(e) NA_real_
    )
  )
}

# 6.2 Plaque 内部不同 grade 比较
plaque_df <- geomx_roi_clean %>% filter(localisation == "Plaque")
if (nrow(plaque_df) > 0 && dplyr::n_distinct(plaque_df$grade) >= 2) {
  stat_list[[length(stat_list) + 1]] <- tibble(
    comparison = "Plaque ROIs only",
    test = "Kruskal-Wallis",
    variable = "APOBEC3A_log1p by lesion grade",
    p_value = tryCatch(
      kruskal.test(APOBEC3A_log1p ~ grade, data = plaque_df)$p.value,
      error = function(e) NA_real_
    )
  )
}

# 6.3 Adventitia 内部不同 grade 比较
adv_df <- geomx_roi_clean %>% filter(localisation == "Adventitia")
if (nrow(adv_df) > 0 && dplyr::n_distinct(adv_df$grade) >= 2) {
  stat_list[[length(stat_list) + 1]] <- tibble(
    comparison = "Adventitia ROIs only",
    test = "Kruskal-Wallis",
    variable = "APOBEC3A_log1p by lesion grade",
    p_value = tryCatch(
      kruskal.test(APOBEC3A_log1p ~ grade, data = adv_df)$p.value,
      error = function(e) NA_real_
    )
  )
}

# 6.4 同一 grade 下 Plaque vs Adventitia 比较
grade_pair_tests <- lapply(levels(geomx_roi_clean$grade), function(g) {
  tmp <- geomx_roi_clean %>% filter(grade == g)

  if (nrow(tmp) == 0 || dplyr::n_distinct(tmp$localisation) < 2) {
    return(NULL)
  }

  tibble(
    comparison = paste0(g, ": Plaque vs Adventitia"),
    test = "Wilcoxon rank-sum",
    variable = "APOBEC3A_log1p",
    p_value = tryCatch(
      wilcox.test(APOBEC3A_log1p ~ localisation, data = tmp)$p.value,
      error = function(e) NA_real_
    )
  )
})

stat_df <- bind_rows(stat_list, bind_rows(grade_pair_tests)) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write.csv(
  stat_df,
  file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_tests.csv")),
  row.names = FALSE
)

print(stat_df)

# ------------------------------------------------------------
# 7. 作图：一张图，不分 facet，不分 marker
# ------------------------------------------------------------

p_a3a_plaque_adventitia <- ggplot(
  geomx_roi_clean,
  aes(x = grade, y = APOBEC3A_log1p, fill = localisation)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.78,
    color = "black",
    linewidth = 0.45
  ) +
  geom_point(
    aes(group = localisation),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.75
    ),
    size = 2.1,
    alpha = 0.8,
    color = "black"
  ) +
  labs(
    title = "GeoMx WTA ROI-level APOBEC3A expression",
    subtitle = "Plaque and adventitia ROIs pooled across CD45/CD4-defined subsets",
    x = "Lesion grade",
    y = "log1p(APOBEC3A count)",
    fill = "ROI localisation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, size = 17),
    plot.subtitle = element_text(hjust = 0, size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_a3a_plaque_adventitia)

ggsave(
  file.path(GEOMX_OUT_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_expression_one_panel.pdf")),
  p_a3a_plaque_adventitia,
  width = 8.5,
  height = 5.5,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(GEOMX_OUT_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_expression_one_panel.png")),
  p_a3a_plaque_adventitia,
  width = 8.5,
  height = 5.5,
  dpi = 320,
  bg = "white"
)

# ------------------------------------------------------------
# 8. 可选：检测率柱状图，同样一张图
# ------------------------------------------------------------

p_a3a_detection <- ggplot(
  geomx_summary,
  aes(x = grade, y = APOBEC3A_detection_rate, fill = localisation)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  geom_text(
    aes(label = paste0(n_APOBEC3A_positive, "/", n_roi)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "GeoMx WTA ROI-level APOBEC3A detection rate",
    subtitle = "Plaque and adventitia ROIs pooled across CD45/CD4-defined subsets",
    x = "Lesion grade",
    y = "APOBEC3A-positive ROI fraction",
    fill = "ROI localisation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, size = 17),
    plot.subtitle = element_text(hjust = 0, size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_a3a_detection)

ggsave(
  file.path(GEOMX_OUT_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_detection_rate_one_panel.pdf")),
  p_a3a_detection,
  width = 8.5,
  height = 5.5,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(GEOMX_OUT_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_detection_rate_one_panel.png")),
  p_a3a_detection,
  width = 8.5,
  height = 5.5,
  dpi = 320,
  bg = "white"
)

# ------------------------------------------------------------
# 9. 完成提示
# ------------------------------------------------------------

message("========== Finished ==========")
message("Output figure directory: ", GEOMX_OUT_DIR)

message("Main output figures:")
message(file.path(GEOMX_OUT_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_expression_one_panel.png")))
message(file.path(GEOMX_OUT_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_detection_rate_one_panel.png")))

message("Main output tables:")
message(file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_plot_data.csv")))
message(file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_summary.csv")))
message(file.path(TABLE_DIR, paste0(OUT_PREFIX, "_geomx_Plaque_Adventitia_ROI_A3A_tests.csv")))

# ============================================================
# 17. Optional exploratory FOV7 / FOV11 / FOV10 analysis:
#     1) p_xen_map 风格同框展示 + core-like region 虚线边界
#     2) 分析 ISG-like Myeloid 与坏死核心样区域的空间关系
#
# 使用前提：
# - 已经运行完第二份代码，至少已经生成 nc_res 和 xen
# - 如果 nc_res 不在当前环境，本脚本会尝试读取保存的 nc_res RDS
# - 如果也没有 RDS，但 detect_necrotic_core_like_all_fovs() 函数存在，则只对 fov7/11/10 重新跑 NC detection
# - 默认不运行；设置 AST_RUN_EXPLORATORY_ISG_CORE_RELATION=1 后启用
# ============================================================

if (RUN_EXPLORATORY_ISG_CORE_RELATION) {

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(Matrix)
  library(sf)
})

# -----------------------------
# 0. 基础路径与参数
# -----------------------------

ROOT <- if (exists("ROOT")) ROOT else "/public3/DSC/single_cell/spatial/atherosclerosis_spatial_publication"
SPATIAL_ROOT <- if (exists("SPATIAL_ROOT")) SPATIAL_ROOT else "/public3/DSC/single_cell/spatial"
TABLE_DIR <- if (exists("TABLE_DIR")) TABLE_DIR else file.path(ROOT, "results", "tables")
FIG_DIR <- if (exists("FIG_DIR")) FIG_DIR else file.path(ROOT, "results", "figures")
OUT_PREFIX <- if (exists("OUT_PREFIX")) OUT_PREFIX else "myeloid_story"
TARGET_GENE <- if (exists("TARGET_GENE")) TARGET_GENE else "APOBEC3A"

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# 目标 FOV 顺序：按你要求 fov7, fov11, fov10 同框展示
TARGET_FOV_REQUEST <- c("fov.7", "fov.11", "fov.10")

# 判定“临近坏死核心”的距离阈值。
# 建议先用 120，因为你 NC_LOCAL_PARAMS 里 radius = 120。
# 如果图上觉得太宽/太窄，可试 80 / 100 / 150。
ADJACENT_DISTANCE <- 120

# permutation 次数，用于判断 ISG-like Myeloid 是否比随机 myeloid 更靠近 core
N_PERM <- 1000
PERM_SEED <- 20260611

# 如果没有 spatial_myeloid_program，但有 Xenium 表达矩阵，则用 ISG signature score top 25% 作为 ISG-like fallback
ISG_FALLBACK_QUANTILE <- 0.75

# 输出文件前缀
OUT_TAG <- paste0(OUT_PREFIX, "_FOV7_11_10_ISG_core_relation")

# -----------------------------
# 1. 小工具函数
# -----------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

save_plot_both2 <- function(p, filename_base, width = 10, height = 6.5) {
  pdf_file <- file.path(FIG_DIR, paste0(filename_base, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename_base, ".png"))

  ggsave(
    pdf_file,
    p,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    png_file,
    p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )

  message("Saved: ", pdf_file)
  message("Saved: ", png_file)

  invisible(c(pdf_file, png_file))
}

sanitize_filename2 <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-\\.]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

theme_story2 <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold")
    )
}

resolve_fov_names <- function(requested, available) {
  norm <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[^a-z0-9]", "", x)
    x
  }

  available_norm <- norm(available)

  resolved <- vapply(requested, function(req) {
    req_norm <- norm(req)

    hit <- which(available_norm == req_norm)
    if (length(hit) > 0) return(available[hit[1]])

    # 支持输入 7 / 10 / 11
    req_num <- gsub("[^0-9]", "", req)
    if (nzchar(req_num)) {
      candidate_norms <- c(
        paste0("fov", req_num),
        paste0("fov0", req_num)
      )
      hit2 <- which(available_norm %in% candidate_norms)
      if (length(hit2) > 0) return(available[hit2[1]])
    }

    NA_character_
  }, character(1))

  if (any(is.na(resolved))) {
    stop(
      "以下 FOV 没有在 xen@images 或 nc_res$per_fov 中找到：",
      paste(requested[is.na(resolved)], collapse = ", "),
      "\n当前可用 FOV：\n",
      paste(available, collapse = "\n")
    )
  }

  unique(resolved)
}

choose_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

to_logical_robust <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x) || is.integer(x)) return(x > 0)
  x2 <- tolower(as.character(x))
  x2 %in% c("true", "t", "1", "yes", "y")
}

get_assay_data_any <- function(object, assay = NULL, layer = "data") {
  if (is.null(assay)) assay <- DefaultAssay(object)

  out <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e1) {
      tryCatch(
        GetAssayData(object, assay = assay, slot = layer),
        error = function(e2) NULL
      )
    }
  )

  out
}

score_gene_set_simple <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }

  x <- as.matrix(mat[genes, , drop = FALSE])

  if (nrow(x) == 1) {
    score <- as.numeric(scale(x[1, ]))
  } else {
    z <- t(scale(t(x)))
    z[is.nan(z) | is.infinite(z)] <- NA_real_
    score <- colMeans(z, na.rm = TRUE)
  }

  score[is.nan(score) | is.infinite(score)] <- NA_real_
  names(score) <- colnames(mat)
  score
}

# 如果 metadata 里没有 ISG-like 程序，就根据 Xenium data/counts 补一个 fallback score
prepare_isg_score_if_needed <- function(obj) {
  program_cols <- intersect(
    c("spatial_myeloid_program", "manual_myeloid_program", "myeloid_substate"),
    colnames(obj@meta.data)
  )

  has_isg_program <- FALSE

  if (length(program_cols) > 0) {
    txt <- apply(
      obj@meta.data[, program_cols, drop = FALSE],
      1,
      function(z) paste(as.character(z), collapse = ";")
    )
    has_isg_program <- any(grepl("ISG", txt, ignore.case = TRUE), na.rm = TRUE)
  }

  if (has_isg_program) {
    message("Found existing ISG-like program column. No fallback ISG score needed.")
    return(obj)
  }

  if ("ISG_like_score_fallback" %in% colnames(obj@meta.data)) {
    message("Found existing ISG_like_score_fallback in metadata.")
    return(obj)
  }

  assay_use <- if ("Xenium" %in% names(obj@assays)) "Xenium" else DefaultAssay(obj)

  mat <- get_assay_data_any(obj, assay = assay_use, layer = "data")
  if (is.null(mat)) {
    mat <- get_assay_data_any(obj, assay = assay_use, layer = "counts")
  }

  if (is.null(mat)) {
    warning("Cannot extract Xenium expression matrix. ISG fallback score skipped.")
    obj@meta.data$ISG_like_score_fallback <- NA_real_
    return(obj)
  }

  isg_genes <- c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3",
    "IFI6", "IFI27", "MX1", "OAS1", "OAS2",
    "STAT1", "RSAD2"
  )

  present_isg <- intersect(isg_genes, rownames(mat))

  if (length(present_isg) < 2) {
    warning(
      "Fewer than 2 ISG signature genes are present in Xenium matrix: ",
      paste(present_isg, collapse = ", "),
      "\nISG fallback score may be unstable."
    )
  }

  obj@meta.data$ISG_like_score_fallback <- score_gene_set_simple(mat, isg_genes)

  message(
    "Computed ISG_like_score_fallback using genes: ",
    paste(present_isg, collapse = ", ")
  )

  obj
}

get_one_fov_df_simple <- function(obj, image_use) {
  coord <- GetTissueCoordinates(obj, image = image_use)
  coord <- as.data.frame(coord)

  if (!"cell" %in% colnames(coord)) {
    coord$cell <- rownames(coord)
  }

  if (!"x" %in% colnames(coord) && "imagecol" %in% colnames(coord)) {
    coord$x <- coord$imagecol
  }

  if (!"y" %in% colnames(coord) && "imagerow" %in% colnames(coord)) {
    coord$y <- coord$imagerow
  }

  if (!all(c("cell", "x", "y") %in% colnames(coord))) {
    stop(
      "FOV ", image_use, " 的坐标表没有 cell/x/y。当前列名：",
      paste(colnames(coord), collapse = ", ")
    )
  }

  coord$cell <- as.character(coord$cell)

  meta <- obj@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")

  meta$cell <- as.character(meta$cell)

  df <- coord %>%
    left_join(meta, by = "cell") %>%
    mutate(
      image = image_use,
      FOV = image_use
    )

  df
}

add_common_flags <- function(df, after_col = "Xenium_MoMa_scPred_Celltype") {
  # severity
  severity_col <- choose_first_col(
    df,
    c(
      "severity", "disease", "Disease", "grade", "Grade",
      "category", "Category", "condition", "Condition",
      "sample_type", "Sample_Type", "lesion_type", "plaque_type"
    )
  )

  if (is.na(severity_col)) {
    df$severity <- "not annotated"
  } else {
    df$severity <- as.character(df[[severity_col]])
    df$severity[is.na(df$severity) | df$severity == ""] <- "not annotated"
  }

  # myeloid definition
  pred_col <- choose_first_col(df, c("predicted.id", "celltype", "Celltype", "seurat_clusters"))

  is_myeloid_pred <- rep(FALSE, nrow(df))
  if (!is.na(pred_col)) {
    is_myeloid_pred <- as.character(df[[pred_col]]) == "Myeloid"
  }

  is_myeloid_transfer <- rep(FALSE, nrow(df))
  if (after_col %in% colnames(df)) {
    is_myeloid_transfer <- as.character(df[[after_col]]) %in% c(
      "Monocyte", "Macrophage", "LAM/Foam Cell", "LAM"
    )
  }

  if ("is_myeloid_transfer" %in% colnames(df)) {
    is_myeloid_transfer <- is_myeloid_transfer | to_logical_robust(df$is_myeloid_transfer)
  }

  if ("mapped_myeloid" %in% colnames(df)) {
    is_myeloid_pred <- is_myeloid_pred | to_logical_robust(df$mapped_myeloid)
  }

  df$is_myeloid_any <- is_myeloid_pred | is_myeloid_transfer

  # A3A detection
  a3a_detect_col <- choose_first_col(
    df,
    c("APOBEC3A_detected", "APOBEC3A_detected_transfer", "A3A_detected")
  )
  a3a_count_col <- choose_first_col(
    df,
    c("APOBEC3A_count", "APOBEC3A_count_transfer", "A3A_count")
  )

  if (!is.na(a3a_detect_col)) {
    df$APOBEC3A_detected_plot <- to_logical_robust(df[[a3a_detect_col]])
  } else if (!is.na(a3a_count_col)) {
    df$APOBEC3A_detected_plot <- as.numeric(df[[a3a_count_col]]) > 0
  } else {
    df$APOBEC3A_detected_plot <- FALSE
  }

  df$is_A3A_positive_myeloid <- df$is_myeloid_any & df$APOBEC3A_detected_plot

  # Foam/LAM myeloid
  program_cols_for_foam <- intersect(
    c("spatial_myeloid_program", "manual_myeloid_program", "myeloid_substate"),
    colnames(df)
  )

  foam_txt <- rep("", nrow(df))
  if (length(program_cols_for_foam) > 0) {
    foam_txt <- apply(
      df[, program_cols_for_foam, drop = FALSE],
      1,
      function(z) paste(as.character(z), collapse = ";")
    )
  }

  foam_from_transfer <- rep(FALSE, nrow(df))
  if (after_col %in% colnames(df)) {
    foam_from_transfer <- as.character(df[[after_col]]) %in% c("LAM", "LAM/Foam Cell")
  }

  df$is_foam_lam_myeloid <- df$is_myeloid_any & (
    grepl("Foam|LAM", foam_txt, ignore.case = TRUE) | foam_from_transfer
  )

  # ISG-like myeloid
  program_cols_for_isg <- intersect(
    c("spatial_myeloid_program", "manual_myeloid_program", "myeloid_substate"),
    colnames(df)
  )

  isg_from_text <- rep(FALSE, nrow(df))

  if (length(program_cols_for_isg) > 0) {
    isg_txt <- apply(
      df[, program_cols_for_isg, drop = FALSE],
      1,
      function(z) paste(as.character(z), collapse = ";")
    )
    isg_from_text <- grepl("ISG", isg_txt, ignore.case = TRUE)
  }

  # 如果原始 program 里确实有 ISG，就用原始 program；
  # 如果完全没有 ISG program，就用 fallback score top 25%。
  if (sum(isg_from_text, na.rm = TRUE) > 0) {
    df$is_ISG_like_myeloid <- df$is_myeloid_any & isg_from_text
  } else if ("ISG_like_score_fallback" %in% colnames(df)) {
    myeloid_scores <- df$ISG_like_score_fallback[df$is_myeloid_any]
    cutoff <- suppressWarnings(
      quantile(myeloid_scores, probs = ISG_FALLBACK_QUANTILE, na.rm = TRUE)
    )
    if (!is.finite(cutoff)) cutoff <- Inf

    df$is_ISG_like_myeloid <- df$is_myeloid_any &
      !is.na(df$ISG_like_score_fallback) &
      df$ISG_like_score_fallback >= cutoff
  } else {
    df$is_ISG_like_myeloid <- FALSE
  }

  # p_xen_map 风格的叠加组：A3A+ 优先级最高
  df$map_group <- dplyr::case_when(
    df$is_A3A_positive_myeloid ~ "A3A+ Myeloid",
    df$is_ISG_like_myeloid ~ "ISG-like Myeloid",
    df$is_foam_lam_myeloid ~ "Foam/LAM Myeloid",
    TRUE ~ NA_character_
  )

  df
}

polygon_sf_to_df2 <- function(poly_sf) {
  if (is.null(poly_sf) || nrow(poly_sf) == 0) {
    return(tibble(
      X = numeric(),
      Y = numeric(),
      core_id = character(),
      group_path = character(),
      image = character()
    ))
  }

  out_list <- lapply(seq_len(nrow(poly_sf)), function(i) {
    cc <- as.data.frame(sf::st_coordinates(poly_sf[i, ]))
    if (nrow(cc) == 0) return(NULL)

    cc$core_id <- as.character(poly_sf$core_id[i])

    grp_cols <- intersect(c("L1", "L2", "L3"), colnames(cc))
    if (length(grp_cols) == 0) {
      cc$group_path <- cc$core_id
    } else {
      grp_val <- apply(cc[, grp_cols, drop = FALSE], 1, paste, collapse = "_")
      cc$group_path <- paste0(cc$core_id, "__", grp_val)
    }

    if ("image" %in% colnames(poly_sf)) {
      cc$image <- as.character(poly_sf$image[i])
    } else {
      cc$image <- NA_character_
    }

    cc
  })

  bind_rows(out_list)
}

# -----------------------------
# 2. 准备 xen 和 nc_res
# -----------------------------

if (!exists("xen")) {
  xen_path <- file.path(SPATIAL_ROOT, "GSE315246_xenium.obj.integrated.rds")
  if (!file.exists(xen_path)) {
    stop("当前环境没有 xen 对象，也找不到默认 xen_path: ", xen_path)
  }
  message("Loading xen from: ", xen_path)
  xen <- readRDS(xen_path)
}

xen <- prepare_isg_score_if_needed(xen)

available_fovs_from_xen <- names(xen@images)

# 如果当前环境已经有 nc_res，则顺手保存一份，方便以后单独运行
nc_res_rds <- file.path(TABLE_DIR, paste0(OUT_PREFIX, "_necrotic_core_like_nc_res.rds"))

if (exists("nc_res")) {
  tryCatch(
    saveRDS(nc_res, nc_res_rds),
    error = function(e) warning("保存 nc_res RDS 失败：", e$message)
  )
}

if (!exists("nc_res")) {
  if (file.exists(nc_res_rds)) {
    message("Loading existing nc_res from: ", nc_res_rds)
    nc_res <- readRDS(nc_res_rds)
  } else if (exists("detect_necrotic_core_like_all_fovs")) {
    message("nc_res not found. Re-running necrotic-core-like detection only for target FOVs.")

    target_fovs_tmp <- resolve_fov_names(TARGET_FOV_REQUEST, available_fovs_from_xen)

    after_col_tmp <- if (exists("after_col")) after_col else "Xenium_MoMa_scPred_Celltype"
    score_col_tmp <- if (exists("score_col")) score_col else "Xenium_MoMa_scPred_score"

    local_params_tmp <- if (exists("NC_LOCAL_PARAMS")) {
      NC_LOCAL_PARAMS
    } else {
      list(
        radius = 120,
        n_perm = 200,
        min_total = 20,
        min_foam = 10,
        min_foam_fraction = 0.45,
        min_z = 2,
        max_p = 0.05
      )
    }

    cluster_params_tmp <- if (exists("NC_CLUSTER_PARAMS")) {
      NC_CLUSTER_PARAMS
    } else {
      list(
        dbscan_eps = 120,
        dbscan_minPts = 5,
        min_cluster_cells = 12,
        min_cluster_foam = 8,
        min_cluster_foam_fraction = 0.45,
        concavity = 4
      )
    }

    nc_res <- do.call(
      detect_necrotic_core_like_all_fovs,
      c(
        list(
          obj = xen,
          after_col = after_col_tmp,
          score_col = score_col_tmp,
          image_names = target_fovs_tmp,
          min_core_area_cutoff = if (exists("NC_CORE_AREA_CUTOFF")) NC_CORE_AREA_CUTOFF else 20000,
          polygon_buffer = if (exists("NC_POLYGON_BUFFER")) NC_POLYGON_BUFFER else 20,
          seed = if (exists("NC_RANDOM_SEED")) NC_RANDOM_SEED else 123
        ),
        local_params_tmp,
        cluster_params_tmp
      )
    )

    xen <- nc_res$obj
    saveRDS(nc_res, nc_res_rds)
  } else {
    stop(
      "没有找到 nc_res，也没有找到保存的 nc_res RDS：", nc_res_rds,
      "\n请先运行第二份代码的 16.1.3 正式识别步骤，或把本代码接在 nc_res 生成之后运行。"
    )
  }
}

if (!is.null(nc_res$obj)) {
  xen <- nc_res$obj
}

available_fovs <- unique(c(names(xen@images), names(nc_res$per_fov)))
TARGET_FOVS <- resolve_fov_names(TARGET_FOV_REQUEST, available_fovs)

message("Target FOVs resolved as: ", paste(TARGET_FOVS, collapse = ", "))

after_col_use <- if (exists("after_col")) after_col else "Xenium_MoMa_scPred_Celltype"

# -----------------------------
# 3. 提取三个 FOV 的细胞表和 core polygon
# -----------------------------

cell_list <- lapply(TARGET_FOVS, function(img) {
  if (!is.null(nc_res$per_fov[[img]]) && !is.null(nc_res$per_fov[[img]]$cell_df)) {
    df <- nc_res$per_fov[[img]]$cell_df
    if (!"image" %in% colnames(df)) df$image <- img
    if (!"FOV" %in% colnames(df)) df$FOV <- img
    df
  } else {
    get_one_fov_df_simple(xen, img)
  }
})

xen_fov_df <- bind_rows(cell_list) %>%
  add_common_flags(after_col = after_col_use)

# 保留指定顺序
xen_fov_df$image <- factor(as.character(xen_fov_df$image), levels = TARGET_FOVS)

# 每个 FOV 一个 severity label
fov_severity <- xen_fov_df %>%
  mutate(image = as.character(image)) %>%
  group_by(image) %>%
  summarise(
    severity = names(sort(table(severity), decreasing = TRUE))[1],
    .groups = "drop"
  )

facet_levels <- fov_severity %>%
  mutate(facet_label = paste0(image, " | Severity: ", severity)) %>%
  arrange(match(image, TARGET_FOVS)) %>%
  pull(facet_label)

xen_fov_df <- xen_fov_df %>%
  mutate(
    image = as.character(image),
    facet_label = paste0(image, " | Severity: ", severity),
    facet_label = factor(facet_label, levels = facet_levels)
  )

# core polygons
core_sf_list <- lapply(TARGET_FOVS, function(img) {
  if (is.null(nc_res$per_fov[[img]])) return(NULL)
  poly <- nc_res$per_fov[[img]]$polygons_sf
  if (is.null(poly) || nrow(poly) == 0) return(NULL)
  poly$image <- img
  poly
})

core_sf_list <- core_sf_list[!vapply(core_sf_list, is.null, logical(1))]

if (length(core_sf_list) > 0) {
  core_sf <- do.call(rbind, core_sf_list)
  core_poly_df <- polygon_sf_to_df2(core_sf)
} else {
  core_sf <- NULL
  core_poly_df <- tibble(
    X = numeric(),
    Y = numeric(),
    core_id = character(),
    group_path = character(),
    image = character()
  )
  warning("这三个 FOV 中没有 accepted core-like polygon。图中不会出现 core 虚线边界。")
}

core_poly_df <- core_poly_df %>%
  left_join(fov_severity, by = "image") %>%
  mutate(
    facet_label = paste0(image, " | Severity: ", severity),
    facet_label = factor(facet_label, levels = facet_levels)
  )

# -----------------------------
# 4. 图 1：p_xen_map 风格 + core-like region 虚线边界
# -----------------------------

xen_map_background <- xen_fov_df
xen_map_overlay <- xen_fov_df %>%
  filter(!is.na(map_group))

p_fov_core_map <- ggplot() +
  geom_point(
    data = xen_map_background,
    aes(x = x, y = y),
    color = "grey86",
    size = 0.08,
    alpha = 0.5
  ) +
  geom_point(
    data = xen_map_overlay,
    aes(x = x, y = y, color = map_group),
    size = 0.22,
    alpha = 0.9
  ) +
  geom_path(
    data = core_poly_df,
    aes(x = X, y = Y, group = group_path),
    color = "black",
    linewidth = 0.45,
    linetype = "22",
    alpha = 0.95
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  scale_y_reverse() +
  scale_color_manual(
    values = c(
      "A3A+ Myeloid" = "#B6424B",
      "ISG-like Myeloid" = "#2ca02c",
      "Foam/LAM Myeloid" = "#9467bd"
    ),
    breaks = c("A3A+ Myeloid", "ISG-like Myeloid", "Foam/LAM Myeloid")
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Myeloid state",
    title = "Xenium FOVs show APOBEC3A-positive and ISG-like myeloid niches around core-like regions",
    subtitle = paste0(
      "Dashed boundaries indicate LAM/Foam Cell-enriched necrotic-core-like regions; adjacent distance cutoff = ",
      ADJACENT_DISTANCE
    )
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9),
    legend.position = "bottom"
  )

print(p_fov_core_map)

save_plot_both2(
  p_fov_core_map,
  paste0(OUT_TAG, "_spatial_map_with_core_boundary"),
  width = 11.5,
  height = 6.5
)

# -----------------------------
# 5. 计算 ISG-like Myeloid 到 core-like region 的空间关系
# -----------------------------

compute_core_relation_one_fov <- function(df_img, core_sf_all, image_use, adjacent_distance = 120) {
  df_img <- df_img %>%
    mutate(
      has_core = FALSE,
      inside_core_polygon = FALSE,
      distance_to_core = NA_real_,
      spatial_relation_to_core = "No accepted core"
    )

  if (is.null(core_sf_all) || nrow(core_sf_all) == 0) {
    return(df_img)
  }

  poly_img <- core_sf_all[as.character(core_sf_all$image) == image_use, , drop = FALSE]

  if (is.null(poly_img) || nrow(poly_img) == 0) {
    return(df_img)
  }

  pts_sf <- sf::st_as_sf(
    df_img,
    coords = c("x", "y"),
    crs = sf::st_crs(poly_img),
    remove = FALSE
  )

  inside_list <- sf::st_within(pts_sf, poly_img, sparse = TRUE)
  inside_flag <- lengths(inside_list) > 0

  dist_mat <- as.matrix(sf::st_distance(pts_sf, poly_img))
  min_dist <- apply(dist_mat, 1, min, na.rm = TRUE)
  min_dist[!is.finite(min_dist)] <- NA_real_

  df_img$has_core <- TRUE
  df_img$inside_core_polygon <- inside_flag
  df_img$distance_to_core <- as.numeric(min_dist)

  df_img$spatial_relation_to_core <- dplyr::case_when(
    inside_flag ~ "Core-colocalized",
    !is.na(min_dist) & min_dist <= adjacent_distance ~ "Peri-core adjacent",
    !is.na(min_dist) & min_dist > adjacent_distance ~ "Distant from core",
    TRUE ~ "No accepted core"
  )

  df_img
}

relation_df <- bind_rows(lapply(TARGET_FOVS, function(img) {
  df_img <- xen_fov_df %>%
    filter(as.character(image) == img)

  compute_core_relation_one_fov(
    df_img = df_img,
    core_sf_all = core_sf,
    image_use = img,
    adjacent_distance = ADJACENT_DISTANCE
  )
}))

relation_levels <- c(
  "Core-colocalized",
  "Peri-core adjacent",
  "Distant from core",
  "No accepted core"
)

relation_df <- relation_df %>%
  mutate(
    spatial_relation_to_core = factor(
      spatial_relation_to_core,
      levels = relation_levels
    ),
    myeloid_group_for_core = ifelse(
      is_ISG_like_myeloid,
      "ISG-like Myeloid",
      "Other Myeloid"
    ),
    proximal_to_core = spatial_relation_to_core %in% c(
      "Core-colocalized",
      "Peri-core adjacent"
    )
  )

# 保存细胞级距离表
cell_relation_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_cell_level_core_distance.csv"))
write.csv(relation_df, cell_relation_out, row.names = FALSE)
message("Saved cell-level ISG/core relation table: ", cell_relation_out)

# -----------------------------
# 6. 统计：ISG-like Myeloid 是共定位、临近还是无关？
# -----------------------------

summarise_relation_one <- function(df_img, image_label = NULL) {
  d <- df_img %>%
    filter(is_myeloid_any)

  if (nrow(d) == 0) {
    return(tibble(
      image = image_label %||% NA_character_,
      severity = NA_character_,
      n_myeloid = 0L,
      n_isg = 0L,
      n_other_myeloid = 0L,
      n_core_regions = 0L,
      isg_core_n = 0L,
      isg_adjacent_n = 0L,
      isg_distant_n = 0L,
      isg_core_fraction = NA_real_,
      isg_proximal_fraction = NA_real_,
      other_proximal_fraction = NA_real_,
      proximal_odds_ratio = NA_real_,
      fisher_p_proximal_enrichment = NA_real_
    ))
  }

  sev <- names(sort(table(d$severity), decreasing = TRUE))[1]
  img <- image_label %||% unique(as.character(d$image))[1]

  n_core_regions <- if (!is.null(core_sf) && nrow(core_sf) > 0) {
    sum(as.character(core_sf$image) == img)
  } else {
    0L
  }

  isg <- d$is_ISG_like_myeloid
  prox <- d$proximal_to_core
  core <- d$spatial_relation_to_core == "Core-colocalized"
  adjacent <- d$spatial_relation_to_core == "Peri-core adjacent"
  distant <- d$spatial_relation_to_core == "Distant from core"

  n_isg <- sum(isg, na.rm = TRUE)
  n_other <- sum(!isg, na.rm = TRUE)

  isg_prox <- sum(isg & prox, na.rm = TRUE)
  isg_nonprox <- sum(isg & !prox, na.rm = TRUE)
  other_prox <- sum(!isg & prox, na.rm = TRUE)
  other_nonprox <- sum(!isg & !prox, na.rm = TRUE)

  fisher_p <- NA_real_
  odds_ratio <- NA_real_

  if (n_isg > 0 && n_other > 0 && (isg_prox + other_prox) > 0 && (isg_nonprox + other_nonprox) > 0) {
    ft <- tryCatch(
      fisher.test(
        matrix(
          c(isg_prox, isg_nonprox, other_prox, other_nonprox),
          nrow = 2,
          byrow = TRUE
        ),
        alternative = "greater"
      ),
      error = function(e) NULL
    )

    if (!is.null(ft)) {
      fisher_p <- ft$p.value
      odds_ratio <- unname(ft$estimate)
    }
  }

  tibble(
    image = img,
    severity = sev,
    n_myeloid = nrow(d),
    n_isg = n_isg,
    n_other_myeloid = n_other,
    n_core_regions = n_core_regions,
    isg_core_n = sum(isg & core, na.rm = TRUE),
    isg_adjacent_n = sum(isg & adjacent, na.rm = TRUE),
    isg_distant_n = sum(isg & distant, na.rm = TRUE),
    isg_core_fraction = ifelse(n_isg > 0, sum(isg & core, na.rm = TRUE) / n_isg, NA_real_),
    isg_proximal_fraction = ifelse(n_isg > 0, isg_prox / n_isg, NA_real_),
    other_proximal_fraction = ifelse(n_other > 0, other_prox / n_other, NA_real_),
    proximal_odds_ratio = odds_ratio,
    fisher_p_proximal_enrichment = fisher_p
  )
}

relation_summary_fov <- bind_rows(lapply(TARGET_FOVS, function(img) {
  summarise_relation_one(
    relation_df %>% filter(as.character(image) == img),
    image_label = img
  )
}))

relation_summary_all <- summarise_relation_one(
  relation_df,
  image_label = "ALL_SELECTED_FOVS"
)

relation_summary <- bind_rows(relation_summary_fov, relation_summary_all)

# permutation：ISG-like Myeloid 的 median distance 是否小于随机 myeloid
run_distance_permutation_one <- function(df_img, image_label, n_perm = 1000, seed = 1) {
  set.seed(seed)

  d <- df_img %>%
    filter(
      is_myeloid_any,
      has_core,
      is.finite(distance_to_core)
    )

  if (nrow(d) == 0 || sum(d$is_ISG_like_myeloid) < 3) {
    return(tibble(
      image = image_label,
      observed_median_distance_ISG = NA_real_,
      median_null_distance = NA_real_,
      permutation_p_nearer_than_random = NA_real_,
      n_perm = n_perm
    ))
  }

  n_isg <- sum(d$is_ISG_like_myeloid)
  obs <- median(d$distance_to_core[d$is_ISG_like_myeloid], na.rm = TRUE)

  null <- replicate(
    n_perm,
    median(sample(d$distance_to_core, size = n_isg, replace = FALSE), na.rm = TRUE)
  )

  p_near <- (sum(null <= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)

  tibble(
    image = image_label,
    observed_median_distance_ISG = obs,
    median_null_distance = median(null, na.rm = TRUE),
    permutation_p_nearer_than_random = p_near,
    n_perm = n_perm
  )
}

perm_fov <- bind_rows(lapply(seq_along(TARGET_FOVS), function(i) {
  img <- TARGET_FOVS[i]
  run_distance_permutation_one(
    relation_df %>% filter(as.character(image) == img),
    image_label = img,
    n_perm = N_PERM,
    seed = PERM_SEED + i
  )
}))

perm_all <- run_distance_permutation_one(
  relation_df,
  image_label = "ALL_SELECTED_FOVS",
  n_perm = N_PERM,
  seed = PERM_SEED + 999
)

perm_summary <- bind_rows(perm_fov, perm_all)

relation_summary <- relation_summary %>%
  left_join(perm_summary, by = "image") %>%
  mutate(
    relationship_call = case_when(
      n_core_regions == 0 ~ "No accepted core detected",
      n_isg < 3 ~ "Too few ISG-like myeloid cells",
      isg_core_fraction >= 0.30 &
        !is.na(proximal_odds_ratio) &
        proximal_odds_ratio > 1.5 &
        !is.na(fisher_p_proximal_enrichment) &
        fisher_p_proximal_enrichment <= 0.05 ~ "Core-colocalized / enriched inside core-like region",
      isg_proximal_fraction >= 0.50 &
        (
          (!is.na(proximal_odds_ratio) & proximal_odds_ratio > 1.5) |
            (!is.na(permutation_p_nearer_than_random) & permutation_p_nearer_than_random <= 0.05)
        ) ~ "Peri-core adjacent / enriched near core-like region",
      isg_proximal_fraction <= 0.20 &
        (is.na(proximal_odds_ratio) | proximal_odds_ratio <= 1.2) &
        (is.na(permutation_p_nearer_than_random) | permutation_p_nearer_than_random > 0.05) ~ "No obvious spatial association",
      TRUE ~ "Ambiguous / weak trend"
    )
  )

relation_summary_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_summary_relation_call.csv"))
write.csv(relation_summary, relation_summary_out, row.names = FALSE)

message("Saved ISG/core relation summary: ", relation_summary_out)
message("========== ISG-like Myeloid vs core-like region relationship ==========")
print(relation_summary)

# -----------------------------
# 7. 图 2：ISG-like Myeloid 与 core 的空间关系图
# -----------------------------

isg_relation_plot_df <- relation_df %>%
  filter(is_ISG_like_myeloid)

p_isg_core_relation_map <- ggplot() +
  geom_point(
    data = xen_fov_df,
    aes(x = x, y = y),
    color = "grey88",
    size = 0.07,
    alpha = 0.45
  ) +
  geom_path(
    data = core_poly_df,
    aes(x = X, y = Y, group = group_path),
    color = "black",
    linewidth = 0.48,
    linetype = "22",
    alpha = 0.95
  ) +
  geom_point(
    data = isg_relation_plot_df,
    aes(x = x, y = y, color = spatial_relation_to_core),
    size = 0.35,
    alpha = 0.95
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  scale_y_reverse() +
  scale_color_manual(
    values = c(
      "Core-colocalized" = "#B6424B",
      "Peri-core adjacent" = "#F58518",
      "Distant from core" = "#2ca02c",
      "No accepted core" = "#8C8C8C"
    ),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "ISG-like myeloid\nrelation to core",
    title = "Spatial relationship between ISG-like myeloid cells and core-like regions",
    subtitle = paste0(
      "Core-colocalized = inside dashed boundary; peri-core adjacent = outside but within ",
      ADJACENT_DISTANCE,
      " coordinate units"
    )
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9),
    legend.position = "bottom"
  )

print(p_isg_core_relation_map)

save_plot_both2(
  p_isg_core_relation_map,
  paste0(OUT_TAG, "_ISG_myeloid_spatial_relation_map"),
  width = 11.5,
  height = 6.5
)

# -----------------------------
# 8. 图 3：ISG-like Myeloid core/adjacent/distant 比例柱状图
# -----------------------------

isg_relation_fraction <- relation_df %>%
  filter(is_ISG_like_myeloid) %>%
  count(image, severity, facet_label, spatial_relation_to_core, name = "n_cells") %>%
  group_by(image, severity, facet_label) %>%
  mutate(
    total_ISG_like_myeloid = sum(n_cells),
    fraction = n_cells / total_ISG_like_myeloid
  ) %>%
  ungroup()

isg_relation_fraction_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_ISG_relation_fraction.csv"))
write.csv(isg_relation_fraction, isg_relation_fraction_out, row.names = FALSE)
message("Saved ISG relation fraction table: ", isg_relation_fraction_out)

p_isg_relation_bar <- ggplot(
  isg_relation_fraction,
  aes(x = facet_label, y = fraction, fill = spatial_relation_to_core)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(
    values = c(
      "Core-colocalized" = "#B6424B",
      "Peri-core adjacent" = "#F58518",
      "Distant from core" = "#2ca02c",
      "No accepted core" = "#8C8C8C"
    ),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Fraction of ISG-like myeloid cells",
    fill = "Relation to core",
    title = "ISG-like myeloid cells are classified by proximity to core-like regions"
  ) +
  theme_story2(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_isg_relation_bar)

save_plot_both2(
  p_isg_relation_bar,
  paste0(OUT_TAG, "_ISG_core_relation_fraction_barplot"),
  width = 8.8,
  height = 5.2
)

# -----------------------------
# 9. 图 4：ISG-like vs Other Myeloid 到 core 的距离分布
# -----------------------------

distance_plot_df <- relation_df %>%
  filter(
    is_myeloid_any,
    has_core,
    is.finite(distance_to_core)
  ) %>%
  mutate(
    myeloid_group_for_core = factor(
      myeloid_group_for_core,
      levels = c("Other Myeloid", "ISG-like Myeloid")
    ),
    log10_distance_to_core = log10(distance_to_core + 1)
  )

distance_test_summary <- distance_plot_df %>%
  group_by(image, severity, facet_label) %>%
  summarise(
    n_ISG = sum(myeloid_group_for_core == "ISG-like Myeloid"),
    n_other = sum(myeloid_group_for_core == "Other Myeloid"),
    median_distance_ISG = median(distance_to_core[myeloid_group_for_core == "ISG-like Myeloid"], na.rm = TRUE),
    median_distance_other = median(distance_to_core[myeloid_group_for_core == "Other Myeloid"], na.rm = TRUE),
    wilcox_p_ISG_nearer = {
      x <- distance_to_core[myeloid_group_for_core == "ISG-like Myeloid"]
      y <- distance_to_core[myeloid_group_for_core == "Other Myeloid"]
      if (length(x) >= 3 && length(y) >= 3) {
        tryCatch(wilcox.test(x, y, alternative = "less")$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

distance_test_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_distance_test_ISG_vs_other_myeloid.csv"))
write.csv(distance_test_summary, distance_test_out, row.names = FALSE)
message("Saved distance test table: ", distance_test_out)

p_distance <- ggplot(
  distance_plot_df,
  aes(x = myeloid_group_for_core, y = log10_distance_to_core, fill = myeloid_group_for_core)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, linewidth = 0.25) +
  geom_jitter(width = 0.15, size = 0.25, alpha = 0.35) +
  facet_wrap(~ facet_label, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Other Myeloid" = "#BDBDBD",
      "ISG-like Myeloid" = "#2ca02c"
    )
  ) +
  labs(
    x = NULL,
    y = "log10(distance to nearest core-like region + 1)",
    fill = "Myeloid group",
    title = "ISG-like myeloid proximity to necrotic-core-like regions",
    subtitle = "Lower distance indicates stronger core colocalization or peri-core localization"
  ) +
  theme_story2(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

print(p_distance)

save_plot_both2(
  p_distance,
  paste0(OUT_TAG, "_distance_to_core_ISG_vs_other_myeloid"),
  width = 10.5,
  height = 5.8
)

# -----------------------------
# 10. 终端输出一个简明结论
# -----------------------------

message("\n========== Final relationship calls ==========")

relation_summary %>%
  select(
    image,
    severity,
    n_isg,
    n_core_regions,
    isg_core_fraction,
    isg_proximal_fraction,
    other_proximal_fraction,
    proximal_odds_ratio,
    fisher_p_proximal_enrichment,
    permutation_p_nearer_than_random,
    relationship_call
  ) %>%
  print(n = Inf)

message("\n主要输出：")
message("1) ", file.path(FIG_DIR, paste0(OUT_TAG, "_spatial_map_with_core_boundary.pdf/png")))
message("2) ", file.path(FIG_DIR, paste0(OUT_TAG, "_ISG_myeloid_spatial_relation_map.pdf/png")))
message("3) ", file.path(FIG_DIR, paste0(OUT_TAG, "_ISG_core_relation_fraction_barplot.pdf/png")))
message("4) ", file.path(FIG_DIR, paste0(OUT_TAG, "_distance_to_core_ISG_vs_other_myeloid.pdf/png")))
message("5) ", relation_summary_out)
message("6) ", cell_relation_out)

} else {
  message("Skipped optional FOV7/FOV10/FOV11 ISG-core exploratory analysis; set AST_RUN_EXPLORATORY_ISG_CORE_RELATION=1 to run it.")
}
