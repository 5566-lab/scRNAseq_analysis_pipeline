#!/usr/bin/env Rscript

# ============================================================
# One-click Xenium core-region workflow
# Core-colocalized / Peri-core adjacent / Other region
# + APOBEC3A-like myeloid score heatmap and statistics
#
# 关键修复：
# 1) 只在“核心多边形生成”阶段使用 NC_POLYGON_BUFFER；
# 2) 后续三分类只计算细胞点到 accepted core polygon 的距离，不再 st_buffer()；
# 3) 不重复 DBSCAN，不重复 build polygon；
# 4) 不再出现 classify_spatial_region_one_fov() 被重复定义后覆盖的情况；
# 5) 所有 sf 点/面均不写 crs = NA，避免 sf 版本兼容问题。
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(sf)
  library(dbscan)
})

try(sf::sf_use_s2(FALSE), silent = TRUE)

# ============================================================
# 0. 路径和参数区：一般只改这里
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

XENIUM_RDS <- file.path(SPATIAL_ROOT, "GSE315246_xenium.obj.integrated.rds")
SC_REF_RDS <- file.path(SC_ROOT, "Result", "sub_integrated_data_Final.rds")

OUTPUT_ROOT <- env_path("AST_OUTPUT_ROOT", file.path(ROOT, "results"))
TABLE_DIR <- file.path(OUTPUT_ROOT, "tables")
FIG_DIR <- file.path(OUTPUT_ROOT, "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PREFIX <- env_path("AST_OUT_PREFIX", "myeloid_story")
OUT_TAG <- paste0(OUT_PREFIX, "_core_region3_ISG_from_scratch")

SAVE_CELL_LEVEL_TABLES <- env_flag("AST_SAVE_CELL_LEVEL_TABLES", FALSE)
SAVE_PER_FOV_PLOTS <- env_flag("AST_SAVE_PER_FOV_PLOTS", FALSE)
RUN_FINAL_ALL_FOV_SUMMARY <- env_flag("AST_RUN_FINAL_ALL_FOV_SUMMARY", FALSE)

# NULL = 自动分析 xen@images 中所有 FOV。
# 如果只想分析部分 FOV，改成例如：c("fov.7", "fov.10", "fov.11")
TARGET_FOV_REQUEST <- NULL

# Xenium 中已有的 Mo/Ma label-transfer 列。
# 如果没有该列，脚本会尝试从 SC_REF_RDS 重新做一次 label transfer。
AFTER_COL <- "Xenium_MoMa_scPred_Celltype"
SCORE_COL <- "Xenium_MoMa_scPred_score"
DO_LABEL_TRANSFER_IF_MISSING <- TRUE
SC_REF_LABEL_COL <- "Celltype_raw"
TRANSFER_DIMS <- 1:30

# APOBEC3A-like myeloid score：优先读取已有缓存 signature；没有则使用内置基因。
SC_A3A_CACHE_TAG <- paste0(OUT_PREFIX, "_sc_ref_A3Apos_within_celltype_signature")
BAD_GENES_FOR_MYELOID_SCORE <- c(
  "MYH11", "ACTA2", "TAGLN", "MYLK",
  "COL1A1", "COL3A1", "FBN1", "AEBP1",
  "PECAM1", "APOLD1", "THBS1"
)

# ------------------------------------------------------------
# core-like region 识别参数
# ------------------------------------------------------------
NC_RANDOM_SEED <- 123

# 局部泡沫/LAM 富集 seed 参数
# 这组参数相对宽松，适合先保证 fov.10/fov.11 可检出；
# 若过松，可提高 min_foam_fraction / min_z，降低 max_p。
NC_LOCAL_PARAMS <- list(
  radius = 100,
  n_perm = 200,
  min_total = 20,
  min_foam = 5,
  min_foam_fraction = 0.30,
  min_z = 1.20,
  max_p = 0.15
)

# seed 聚类及 polygon 后过滤参数
NC_CLUSTER_PARAMS <- list(
  dbscan_eps = 120,
  dbscan_minPts = 4,
  min_cluster_cells = 12,
  min_cluster_foam = 8,
  min_cluster_foam_fraction = 0.10,
  concavity = 4
)

# 关键参数：
# 只在 build_core_polygon 阶段用一次。
# 设 0 表示完全不扩张核心 polygon。
NC_POLYGON_BUFFER <- 20

# accepted core-like region 的最小面积过滤。
NC_CORE_AREA_CUTOFF <- 20000

# Peri-core adjacent 定义：
# 非 core 内细胞，且到 accepted core polygon 的距离 <= ADJACENT_DISTANCE。
# 注意：这里不再 st_buffer()，只计算距离。
ADJACENT_DISTANCE <- 120

# ============================================================
# 1. 标签、颜色和 signature
# ============================================================

TARGET_GENE <- "APOBEC3A"
LAM_DISPLAY_LABEL <- "LAM/Foam Cell"

region_levels3 <- c(
  "Core-colocalized",
  "Peri-core adjacent",
  "Other region"
)

region_colors3 <- c(
  "Core-colocalized" = "#B6424B",
  "Peri-core adjacent" = "#F58518",
  "Other region" = "#8C8C8C"
)

myeloid_signatures <- list(
  Mono_classical = c("FCN1", "S100A8", "S100A9", "VCAN", "SELL", "CCR2", "LYZ", "LST1", "CD14"),
  Macrophage_core = c("CD68", "CSF1R", "AIF1", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "MSR1"),
  Foam_LAM = c("APOE", "APOC1", "TREM2", "LPL", "SPP1", "GPNMB", "PLA2G7", "LGALS3", "FABP5", "LIPA"),
  ISG_like = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27", "MX1", "OAS1", "OAS2", "STAT1", "RSAD2", "MARCO", "LPL", "C7", "TNF"),
  Inflammatory = c("IL1B", "S100A8", "S100A9", "CXCL8", "CCL2", "CCL3", "CCL4", "NFKBIA", "NAMPT"),
  Resident_TRM = c("LYVE1", "SELENOP", "MRC1", "FOLR2", "F13A1", "C1QA", "C1QB", "C1QC")
)

# ============================================================
# 2. 通用函数
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

cat2 <- function(...) {
  cat(..., "\n")
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-\\.]", "_", x)
  x <- gsub("_+", "_", x)
  x
}

theme_story <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}

save_plot_both <- function(p, filename_base, width = 8, height = 6, dpi = 320) {
  pdf_file <- file.path(FIG_DIR, paste0(filename_base, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename_base, ".png"))

  ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )

  message("Saved PDF: ", pdf_file)
  message("Saved PNG: ", png_file)

  invisible(c(pdf_file, png_file))
}

choose_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_severity_col <- function(meta_df) {
  choose_first_col(
    meta_df,
    c(
      "Severity", "severity",
      "Status", "status",
      "disease", "Disease",
      "grade", "Grade",
      "category", "Category",
      "condition", "Condition",
      "sample_type", "Sample_Type",
      "lesion_type", "plaque_type"
    )
  )
}

get_fov_severity_label <- function(df) {
  severity_col <- get_severity_col(df)

  if (is.na(severity_col)) {
    return("not annotated")
  }

  sev_vec <- as.character(df[[severity_col]])
  sev_vec <- sev_vec[!is.na(sev_vec) & sev_vec != ""]

  if (length(sev_vec) == 0) {
    return("not annotated")
  }

  sev_tab <- sort(table(sev_vec), decreasing = TRUE)

  if (length(sev_tab) == 1) {
    return(names(sev_tab)[1])
  }

  paste0(names(sev_tab), " n=", as.integer(sev_tab), collapse = "; ")
}

set_default_assay_safely <- function(obj, preferred = c("Xenium", "RNA", "Spatial", "SCT", "integrated")) {
  if (!inherits(obj, "Seurat")) {
    stop("输入对象不是 Seurat object。当前 class: ", paste(class(obj), collapse = ", "))
  }

  assay_names <- names(obj@assays)
  if (length(assay_names) == 0) {
    stop("Seurat object 中没有 assay。")
  }

  chosen <- intersect(preferred, assay_names)
  if (length(chosen) > 0) {
    DefaultAssay(obj) <- chosen[1]
  } else {
    DefaultAssay(obj) <- assay_names[1]
    warning("没有找到 preferred assay，使用第一个 assay: ", assay_names[1])
  }

  message("DefaultAssay set to: ", DefaultAssay(obj))
  obj
}

join_layers_safely <- function(obj, assay = NULL, object_name = "object") {
  if (is.null(assay)) assay <- DefaultAssay(obj)

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
      paste(class(assay_obj), collapse = "/")
    )
  }

  obj
}

strip_images_copy <- function(obj) {
  obj2 <- obj
  if ("images" %in% slotNames(obj2)) {
    obj2@images <- list()
  }
  obj2
}

get_assay_data_or_null <- function(object, assay = NULL, layer = "data") {
  if (is.null(assay)) assay <- DefaultAssay(object)

  tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e1) {
      tryCatch(
        GetAssayData(object, assay = assay, slot = layer),
        error = function(e2) NULL
      )
    }
  )
}

get_best_assay_matrix <- function(object, assay = NULL, prefer_layers = c("data", "counts")) {
  if (is.null(assay)) assay <- DefaultAssay(object)

  for (ly in prefer_layers) {
    mat <- get_assay_data_or_null(object, assay = assay, layer = ly)
    if (!is.null(mat)) {
      return(list(mat = mat, layer = ly, assay = assay))
    }
  }

  stop("无法从 assay=", assay, " 提取 data 或 counts layer/slot。")
}

natural_fov_order <- function(x) {
  num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
  num[is.na(num)] <- Inf
  x[order(num, x)]
}

resolve_fov_names <- function(requested, available) {
  available <- as.character(available)

  if (is.null(requested) || length(requested) == 0) {
    return(natural_fov_order(available))
  }

  norm <- function(x) {
    x <- tolower(as.character(x))
    gsub("[^a-z0-9]", "", x)
  }

  available_norm <- norm(available)

  resolved <- vapply(requested, function(req) {
    req_norm <- norm(req)

    hit <- which(available_norm == req_norm)
    if (length(hit) > 0) {
      return(available[hit[1]])
    }

    req_num <- gsub("[^0-9]", "", req)
    if (nzchar(req_num)) {
      candidate_norms <- c(
        paste0("fov", req_num),
        paste0("fov0", req_num)
      )
      hit2 <- which(available_norm %in% candidate_norms)
      if (length(hit2) > 0) {
        return(available[hit2[1]])
      }
    }

    NA_character_
  }, character(1))

  if (any(is.na(resolved))) {
    stop(
      "以下 FOV 没有找到：",
      paste(requested[is.na(resolved)], collapse = ", "),
      "\n当前可用 FOV:\n",
      paste(available, collapse = "\n")
    )
  }

  unique(resolved)
}

rename_lam_labels <- function(x) {
  x <- as.character(x)
  x[x == "LAM"] <- LAM_DISPLAY_LABEL
  x[x == "Foam cells"] <- LAM_DISPLAY_LABEL
  x[x == "Foam cell"] <- LAM_DISPLAY_LABEL
  x[x == "Foam cells1"] <- LAM_DISPLAY_LABEL
  x[x == "Foam cells2"] <- LAM_DISPLAY_LABEL
  x
}

broad_myeloid_label <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    is.na(x) | x == "" ~ NA_character_,
    grepl("LAM|Foam", x, ignore.case = TRUE) ~ LAM_DISPLAY_LABEL,
    grepl("Mac", x, ignore.case = TRUE) ~ "Macrophage",
    grepl("Mono", x, ignore.case = TRUE) ~ "Monocyte",
    TRUE ~ x
  )
}

winsorize_vec <- function(x, probs = c(0.01, 0.99)) {
  q <- suppressWarnings(quantile(x, probs = probs, na.rm = TRUE))
  if (any(!is.finite(q))) return(x)
  pmin(pmax(x, q[1]), q[2])
}

# ============================================================
# 3. 读取 Xenium object
# ============================================================

if (!file.exists(XENIUM_RDS)) {
  stop("找不到 Xenium RDS: ", XENIUM_RDS)
}

message("Loading Xenium object: ", XENIUM_RDS)
xen <- readRDS(XENIUM_RDS)

if (!inherits(xen, "Seurat")) {
  stop("Xenium object 不是 Seurat object。当前 class: ", paste(class(xen), collapse = ", "))
}

xen <- set_default_assay_safely(xen, preferred = c("Xenium", "RNA", "Spatial", "SCT", "integrated"))
xen <- join_layers_safely(xen, assay = DefaultAssay(xen), object_name = "xen")

if (length(xen@images) == 0) {
  stop("xen@images 为空，没有 FOV 坐标，无法画空间图。")
}

TARGET_FOVS <- resolve_fov_names(TARGET_FOV_REQUEST, names(xen@images))
message("Target FOVs: ", paste(TARGET_FOVS, collapse = ", "))

# ============================================================
# 4. 如没有 Mo/Ma transfer 列，自动补 label transfer
# ============================================================

prepare_sc_ref_for_transfer <- function(sc_ref_path = SC_REF_RDS, ref_label_col = SC_REF_LABEL_COL) {
  if (!file.exists(sc_ref_path)) {
    stop("找不到单细胞参考对象: ", sc_ref_path)
  }

  message("Loading single-cell reference: ", sc_ref_path)
  sc_ref <- readRDS(sc_ref_path)

  if (!inherits(sc_ref, "Seurat")) {
    stop("sc_ref 不是 Seurat object。当前 class: ", paste(class(sc_ref), collapse = ", "))
  }

  if (!ref_label_col %in% colnames(sc_ref@meta.data)) {
    fallback <- choose_first_col(
      sc_ref@meta.data,
      c("Celltype_raw", "Celltype_raw1", "celltype", "CellType", "predicted.id", "seurat_clusters")
    )

    if (is.na(fallback)) {
      stop(
        "sc_ref@meta.data 中找不到 ", ref_label_col,
        "，也找不到可替代细胞类型列。当前列:\n",
        paste(colnames(sc_ref@meta.data), collapse = "\n")
      )
    }

    warning("找不到 ", ref_label_col, "，改用 ", fallback, " 做 transfer label。")
    ref_label_col <- fallback
  }

  sc_ref <- set_default_assay_safely(sc_ref, preferred = c("RNA", "SCT", "integrated"))
  sc_ref <- join_layers_safely(sc_ref, assay = DefaultAssay(sc_ref), object_name = "sc_ref")

  raw_label <- rename_lam_labels(sc_ref@meta.data[[ref_label_col]])
  broad_label <- broad_myeloid_label(raw_label)

  keep <- !is.na(broad_label) & broad_label %in% c("Monocyte", "Macrophage", LAM_DISPLAY_LABEL)

  if (sum(keep) < 50) {
    stop("sc_ref 中可用于 Mo/Ma transfer 的细胞太少: ", sum(keep))
  }

  sc_ref <- subset(sc_ref, cells = rownames(sc_ref@meta.data)[keep])
  sc_ref$Celltype_transfer <- broad_label[keep]

  sc_ref <- tryCatch(
    NormalizeData(sc_ref, verbose = FALSE),
    error = function(e) {
      message("NormalizeData skipped for sc_ref: ", e$message)
      sc_ref
    }
  )

  sc_ref <- tryCatch(
    FindVariableFeatures(sc_ref, verbose = FALSE),
    error = function(e) {
      message("FindVariableFeatures skipped for sc_ref: ", e$message)
      sc_ref
    }
  )

  message("Reference transfer label counts:")
  print(table(sc_ref$Celltype_transfer, useNA = "ifany"))

  sc_ref
}

run_label_transfer_if_needed <- function(xen) {
  if (AFTER_COL %in% colnames(xen@meta.data)) {
    message("Found existing transfer column: ", AFTER_COL)
    return(xen)
  }

  if (!isTRUE(DO_LABEL_TRANSFER_IF_MISSING)) {
    stop("xen@meta.data 中没有 ", AFTER_COL, "，且 DO_LABEL_TRANSFER_IF_MISSING = FALSE。")
  }

  message("No ", AFTER_COL, " found. Start label transfer from sc_ref.")

  sc_ref <- prepare_sc_ref_for_transfer()

  query <- strip_images_copy(xen)
  query <- set_default_assay_safely(query, preferred = c("Xenium", "RNA", "Spatial", "SCT", "integrated"))
  query <- join_layers_safely(query, assay = DefaultAssay(query), object_name = "xen_query_no_images")

  query <- tryCatch(
    NormalizeData(query, verbose = FALSE),
    error = function(e) {
      message("NormalizeData skipped for query: ", e$message)
      query
    }
  )

  ref_features <- rownames(sc_ref)
  query_features <- rownames(query)
  common_features <- intersect(ref_features, query_features)

  if (length(common_features) < 30) {
    stop(
      "sc_ref 与 Xenium query 共同基因太少，无法稳定 transfer: ",
      length(common_features),
      "\n共同基因：", paste(common_features, collapse = ", ")
    )
  }

  var_features <- intersect(VariableFeatures(sc_ref), common_features)
  if (length(var_features) >= 30) {
    features_use <- var_features
  } else {
    features_use <- common_features
  }

  message("Common transfer features: ", length(common_features))
  message("Features used for FindTransferAnchors: ", length(features_use))

  anchors <- FindTransferAnchors(
    reference = sc_ref,
    query = query,
    reference.assay = DefaultAssay(sc_ref),
    query.assay = DefaultAssay(query),
    normalization.method = "LogNormalize",
    features = features_use,
    dims = TRANSFER_DIMS
  )

  pred <- TransferData(
    anchorset = anchors,
    refdata = sc_ref$Celltype_transfer,
    dims = TRANSFER_DIMS
  )

  common_cells <- intersect(rownames(xen@meta.data), rownames(pred))

  xen@meta.data[[AFTER_COL]] <- NA_character_
  xen@meta.data[[SCORE_COL]] <- NA_real_

  xen@meta.data[common_cells, AFTER_COL] <- as.character(pred[common_cells, "predicted.id"])
  if ("prediction.score.max" %in% colnames(pred)) {
    xen@meta.data[common_cells, SCORE_COL] <- as.numeric(pred[common_cells, "prediction.score.max"])
  }

  message("Transfer result counts:")
  print(table(xen@meta.data[[AFTER_COL]], useNA = "ifany"))

  xen
}

xen <- run_label_transfer_if_needed(xen)

# ============================================================
# 5. 计算 ISG-like score
# ============================================================

load_isg_signature <- function() {
  candidate_cache <- c(
    file.path(TABLE_DIR, paste0(SC_A3A_CACHE_TAG, "_updated_ISG_like_signature_A3Apos_within_celltype.csv")),
    list.files(TABLE_DIR, pattern = "updated_ISG_like_signature_A3Apos_within_celltype\\.csv$", full.names = TRUE)
  )
  candidate_cache <- unique(candidate_cache[file.exists(candidate_cache)])

  if (length(candidate_cache) > 0) {
    cache_file <- candidate_cache[1]
    message("Read cached ISG-like signature: ", cache_file)

    sig_df <- read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE)
    if ("gene" %in% colnames(sig_df)) {
      sig <- unique(as.character(sig_df$gene))
      sig <- sig[!is.na(sig) & sig != ""]
      sig <- setdiff(sig, BAD_GENES_FOR_MYELOID_SCORE)
      if (length(sig) >= 2) {
        return(sig)
      }
    }
  }

  message("No usable cached ISG-like signature found. Use built-in myeloid_signatures$ISG_like.")
  sig <- unique(myeloid_signatures$ISG_like)
  sig <- setdiff(sig, BAD_GENES_FOR_MYELOID_SCORE)
  sig
}

score_gene_set <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))

  if (length(genes) < 2) {
    stop("用于 score 的 genes 少于 2 个。当前匹配基因: ", paste(genes, collapse = ", "))
  }

  x <- as.matrix(mat[genes, , drop = FALSE])

  if (nrow(x) == 1) {
    score <- as.numeric(scale(x[1, ]))
  } else {
    gene_sd <- apply(x, 1, sd, na.rm = TRUE)
    keep <- is.finite(gene_sd) & gene_sd > 0

    if (sum(keep) < 2) {
      stop("用于 score 的基因方差太低，无法计算 z-score。")
    }

    x <- x[keep, , drop = FALSE]
    z <- t(scale(t(x)))
    z[is.nan(z) | is.infinite(z)] <- NA_real_
    score <- colMeans(z, na.rm = TRUE)
  }

  score[is.nan(score) | is.infinite(score)] <- NA_real_
  names(score) <- colnames(mat)
  score
}

calculate_isg_score <- function(obj, signature_genes) {
  assay_use <- if ("Xenium" %in% names(obj@assays)) "Xenium" else DefaultAssay(obj)

  mat_info <- get_best_assay_matrix(obj, assay = assay_use, prefer_layers = c("data", "counts"))
  mat <- mat_info$mat

  if (mat_info$layer == "counts") {
    message("Use counts layer for score; applying log1p(counts).")
    mat <- log1p(mat)
  } else {
    message("Use data layer for score.")
  }

  present_genes <- intersect(signature_genes, rownames(mat))
  missing_genes <- setdiff(signature_genes, rownames(mat))

  message("Total signature genes: ", length(signature_genes))
  message("Detected signature genes in Xenium: ", length(present_genes))
  message("Detected genes: ", paste(present_genes, collapse = ", "))

  if (length(missing_genes) > 0) {
    message("Missing genes: ", paste(missing_genes, collapse = ", "))
  }

  score <- score_gene_set(mat, present_genes)

  obj@meta.data$ISG_like_score <- NA_real_
  obj@meta.data$A3A_like_myeloid_score <- NA_real_

  common_cells <- intersect(rownames(obj@meta.data), names(score))
  obj@meta.data[common_cells, "ISG_like_score"] <- score[common_cells]
  obj@meta.data[common_cells, "A3A_like_myeloid_score"] <- score[common_cells]

  message("ISG_like_score summary:")
  print(summary(obj@meta.data$ISG_like_score))

  obj
}

signature_genes <- load_isg_signature()
message("Final ISG-like signature used:")
print(signature_genes)

xen <- calculate_isg_score(xen, signature_genes)

# ============================================================
# 6. 提取 FOV cell table，并标记 myeloid / foam
# ============================================================

get_one_fov_df <- function(obj, image_use) {
  if (!image_use %in% names(obj@images)) {
    stop("Image/FOV not found in xen@images: ", image_use)
  }

  img_obj <- obj@images[[image_use]]

  coords <- tryCatch(
    Seurat::GetTissueCoordinates(img_obj),
    error = function(e1) {
      tryCatch(
        Seurat::GetTissueCoordinates(obj, image = image_use),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(coords)) {
    stop("无法提取 FOV 坐标: ", image_use)
  }

  coords <- as.data.frame(coords)

  if (!"cell" %in% colnames(coords)) {
    if ("cells" %in% colnames(coords)) {
      coords$cell <- coords$cells
    } else if (!is.null(rownames(coords))) {
      coords$cell <- rownames(coords)
    } else {
      stop("坐标表中找不到 cell/cells，也没有 rownames: ", image_use)
    }
  }

  x_col <- choose_first_col(coords, c("x", "X", "x_centroid", "centroid_x", "global_x", "imagecol", "pxl_col_in_fullres"))
  y_col <- choose_first_col(coords, c("y", "Y", "y_centroid", "centroid_y", "global_y", "imagerow", "pxl_row_in_fullres"))

  if (is.na(x_col) || is.na(y_col)) {
    stop(
      "无法识别坐标 x/y 列: ", image_use,
      "\n当前坐标列: ", paste(colnames(coords), collapse = ", ")
    )
  }

  coords <- coords %>%
    transmute(
      cell = as.character(cell),
      x = as.numeric(.data[[x_col]]),
      y = as.numeric(.data[[y_col]])
    )

  common_cells <- intersect(coords$cell, rownames(obj@meta.data))

  if (length(common_cells) == 0) {
    stop("FOV 坐标 cell 与 xen@meta.data 完全无法匹配: ", image_use)
  }

  coords <- coords %>%
    filter(cell %in% common_cells)

  meta <- obj@meta.data[coords$cell, , drop = FALSE] %>%
    as.data.frame(check.names = FALSE)

  meta <- meta[, setdiff(colnames(meta), c("cell", "x", "y", "image")), drop = FALSE]

  out <- bind_cols(coords, meta)
  out$image <- image_use
  out
}

add_myeloid_and_foam_flags <- function(df, ann_col = AFTER_COL) {
  if (!ann_col %in% colnames(df)) {
    candidates <- c(
      "Xenium_MoMa_scPred_Celltype",
      "transferred_celltype",
      "Celltype_transfer",
      "Celltype_raw1",
      "Celltype_raw",
      "predicted.id",
      "celltype",
      "CellType"
    )
    hit <- intersect(candidates, colnames(df))

    if (length(hit) == 0) {
      stop(
        "找不到髓系注释列。当前可用列:\n",
        paste(colnames(df), collapse = "\n")
      )
    }

    ann_col <- hit[1]
    warning("找不到 AFTER_COL，改用注释列: ", ann_col)
  }

  ann <- as.character(df[[ann_col]])
  ann <- rename_lam_labels(ann)
  broad <- broad_myeloid_label(ann)

  df$CellType_transfer_raw <- ann
  df$CellType_transfer_broad <- broad

  df$is_myeloid_transfer <- !is.na(ann) &
    grepl("Mono|Mac|LAM|Foam|myeloid", ann, ignore.case = TRUE)

  df$is_foam_transfer <- !is.na(ann) &
    grepl("LAM|Foam", ann, ignore.case = TRUE)

  df
}

message("Extracting selected FOV coordinates and metadata...")
cell_df <- bind_rows(lapply(TARGET_FOVS, function(img) {
  message("Extracting FOV: ", img)
  get_one_fov_df(xen, img)
}))

cell_df <- add_myeloid_and_foam_flags(cell_df, ann_col = AFTER_COL)

message("Cell counts by FOV:")
print(table(cell_df$image, useNA = "ifany"))

message("Transferred myeloid counts by FOV:")
print(table(cell_df$image, cell_df$is_myeloid_transfer, useNA = "ifany"))

message("Foam/LAM counts by FOV:")
print(table(cell_df$image, cell_df$is_foam_transfer, useNA = "ifany"))

# ============================================================
# 7. core-like polygon 构建工具
# ============================================================

polygon_sf_to_df <- function(poly_sf) {
  if (is.null(poly_sf) || nrow(poly_sf) == 0) {
    return(tibble(
      X = numeric(),
      Y = numeric(),
      group_path = character(),
      core_id = character(),
      image = character()
    ))
  }

  out_list <- lapply(seq_len(nrow(poly_sf)), function(i) {
    cc <- as.data.frame(sf::st_coordinates(poly_sf[i, ]))

    if (nrow(cc) == 0) return(NULL)

    if (!"core_id" %in% colnames(poly_sf)) {
      core_id_i <- paste0("core", i)
    } else {
      core_id_i <- as.character(poly_sf$core_id[i])
    }

    cc$core_id <- core_id_i

    grp_cols <- intersect(c("L1", "L2", "L3"), colnames(cc))
    if (length(grp_cols) == 0) {
      cc$group_path <- core_id_i
    } else {
      grp_val <- apply(cc[, grp_cols, drop = FALSE], 1, paste, collapse = "_")
      cc$group_path <- paste0(core_id_i, "__", grp_val)
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

points_to_core_polygon <- function(df_pts, core_id, concavity = 4) {
  df_pts <- df_pts %>%
    transmute(
      x = as.numeric(x),
      y = as.numeric(y)
    ) %>%
    filter(is.finite(x), is.finite(y)) %>%
    distinct()

  if (nrow(df_pts) < 3) return(NULL)

  sf_pts <- sf::st_as_sf(
    df_pts,
    coords = c("x", "y"),
    remove = FALSE
  )

  geom <- tryCatch(
    {
      if (requireNamespace("concaveman", quietly = TRUE) && nrow(df_pts) >= 4) {
        tmp <- concaveman::concaveman(sf_pts, concavity = concavity, length_threshold = 0)
        sf::st_geometry(tmp)
      } else {
        sf::st_sfc(sf::st_convex_hull(sf::st_union(sf::st_geometry(sf_pts))))
      }
    },
    error = function(e) {
      message("concave hull failed for ", core_id, "; fallback to convex hull. Error: ", e$message)
      sf::st_sfc(sf::st_convex_hull(sf::st_union(sf::st_geometry(sf_pts))))
    }
  )

  poly <- sf::st_sf(core_id = core_id, geometry = geom)
  poly <- suppressWarnings(sf::st_make_valid(poly))

  poly <- tryCatch(
    suppressWarnings(sf::st_collection_extract(poly, "POLYGON")),
    error = function(e) poly
  )

  if (is.null(poly) || nrow(poly) == 0) return(NULL)

  geom_type <- as.character(sf::st_geometry_type(poly, by_geometry = FALSE))
  if (!grepl("POLYGON", geom_type)) return(NULL)

  poly$core_id <- core_id
  poly <- poly[, "core_id", drop = FALSE]

  poly
}

nc_empty_cluster_summary <- function() {
  tibble(
    FOV = character(),
    severity = character(),
    core_id = character(),
    area_pass_cutoff = logical(),
    n_cells_in_polygon = integer(),
    n_myeloid_in_polygon = integer(),
    n_foam_in_polygon = integer(),
    foam_fraction_in_polygon = numeric(),
    polygon_area_raw = numeric(),
    polygon_area = numeric(),
    polygon_buffer_used = numeric(),
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
    mean_core_area = NA_real_,
    total_core_area = 0,
    n_total_cells = n_total_cells,
    n_myeloid_cells = n_myeloid_cells,
    n_core_like_cells = 0L,
    n_myeloid_core_like_cells = 0L,
    n_foam_core_like_cells = 0L
  )
}

# ============================================================
# 8. core-like region 识别
#    注意：polygon_buffer 在这里真正传入并使用一次
# ============================================================

validate_polygon_buffer <- function(polygon_buffer) {
  if (is.null(polygon_buffer)) return(0)
  if (!is.numeric(polygon_buffer) || length(polygon_buffer) != 1) {
    stop("polygon_buffer 必须是长度为 1 的 numeric。当前: ", paste(polygon_buffer, collapse = ", "))
  }
  if (is.na(polygon_buffer) || !is.finite(polygon_buffer)) {
    stop("polygon_buffer 不能是 NA/NaN/Inf。")
  }
  if (polygon_buffer < 0) {
    stop("polygon_buffer 不能小于 0。")
  }
  as.numeric(polygon_buffer)
}

detect_necrotic_core_like_one_fov <- function(
    df,
    image_use,
    radius = 120,
    n_perm = 200,
    min_total = 20,
    min_foam = 10,
    min_foam_fraction = 0.45,
    min_z = 2,
    max_p = 0.05,
    dbscan_eps = 120,
    dbscan_minPts = 5,
    min_cluster_cells = 12,
    min_cluster_foam = 8,
    min_cluster_foam_fraction = 0.45,
    min_core_area_cutoff = 20000,
    concavity = 4,
    polygon_buffer = 20,
    seed = 123
) {
  polygon_buffer <- validate_polygon_buffer(polygon_buffer)

  df <- df %>%
    mutate(
      x = as.numeric(x),
      y = as.numeric(y)
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

  myeloid_df <- df %>%
    filter(
      is_myeloid_transfer,
      is.finite(x),
      is.finite(y)
    )

  if (nrow(myeloid_df) < min_total) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = tibble(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = tibble(),
      cluster_summary = nc_empty_cluster_summary(),
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
  obs_foam <- sapply(nn$id, function(idx) {
    sum(myeloid_df$is_foam_transfer[idx], na.rm = TRUE)
  })

  obs_frac <- ifelse(obs_total > 0, obs_foam / obs_total, NA_real_)

  set.seed(seed)

  perm_mat <- replicate(n_perm, {
    shuffled <- sample(myeloid_df$is_foam_transfer)
    sapply(nn$id, function(idx) {
      sum(shuffled[idx], na.rm = TRUE)
    })
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

  candidate_df <- myeloid_df %>%
    filter(NC_candidate_seed)

  if (nrow(candidate_df) < dbscan_minPts) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = tibble(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = tibble(),
      cluster_summary = nc_empty_cluster_summary(),
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
  candidate_df <- candidate_df %>%
    filter(raw_cluster > 0)

  if (nrow(candidate_df) == 0) {
    return(list(
      cell_df = df,
      polygons_sf = NULL,
      polygons_df = tibble(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = tibble(),
      cluster_summary = nc_empty_cluster_summary(),
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

  valid_all <- is.finite(df$x) & is.finite(df$y)
  valid_idx <- which(valid_all)

  all_pts_sf <- sf::st_as_sf(
    df[valid_all, , drop = FALSE],
    coords = c("x", "y"),
    remove = FALSE
  )

  candidate_poly_list <- list()
  accepted_poly_list <- list()
  summary_list <- list()
  core_counter <- 0L

  for (cl in raw_clusters) {
    cl_pts <- candidate_df %>%
      filter(raw_cluster == cl)

    if (nrow(cl_pts) < 3) next

    core_counter <- core_counter + 1L
    core_id <- paste0(image_use, "_core", core_counter)

    poly_raw <- tryCatch(
      points_to_core_polygon(
        cl_pts[, c("x", "y")],
        core_id = core_id,
        concavity = concavity
      ),
      error = function(e) {
        message("points_to_core_polygon failed for ", core_id, ": ", e$message)
        NULL
      }
    )

    if (is.null(poly_raw) || nrow(poly_raw) == 0) next

    poly_raw <- suppressWarnings(sf::st_make_valid(poly_raw))

    raw_area_val <- tryCatch(
      as.numeric(sf::st_area(poly_raw)),
      error = function(e) NA_real_
    )

    # 关键修复：
    # polygon_buffer 参数只在这里使用一次；
    # 后续 cell relation / peri-core 定义不再 st_buffer。
    poly_sf <- poly_raw
    if (polygon_buffer > 0) {
      poly_sf <- suppressWarnings(sf::st_buffer(poly_sf, dist = polygon_buffer))
      poly_sf <- suppressWarnings(sf::st_make_valid(poly_sf))
    }

    poly_sf$core_id <- core_id
    poly_sf$image <- image_use
    poly_sf$polygon_buffer_used <- polygon_buffer

    inside_list <- sf::st_within(all_pts_sf, poly_sf, sparse = TRUE)
    inside_flag_valid <- lengths(inside_list) > 0

    sub_all <- df[valid_idx[inside_flag_valid], , drop = FALSE]
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

    area_val <- tryCatch(
      as.numeric(sf::st_area(poly_sf)),
      error = function(e) NA_real_
    )

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
      polygon_area_raw = raw_area_val,
      polygon_area = area_val,
      polygon_buffer_used = polygon_buffer,
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
      polygons_df = tibble(),
      candidate_polygons_sf = NULL,
      candidate_polygons_df = tibble(),
      cluster_summary = nc_empty_cluster_summary(),
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

  candidate_polygons_sf <- NULL
  candidate_polygons_df <- tibble()

  if (length(candidate_poly_list) > 0) {
    candidate_polygons_sf <- do.call(rbind, candidate_poly_list)
    candidate_polygons_df <- polygon_sf_to_df(candidate_polygons_sf)
  }

  accepted_polygons_sf <- NULL
  accepted_polygons_df <- tibble()

  if (length(accepted_poly_list) > 0) {
    accepted_polygons_sf <- do.call(rbind, accepted_poly_list)
    accepted_polygons_df <- polygon_sf_to_df(accepted_polygons_sf)

    inside_all <- sf::st_within(all_pts_sf, accepted_polygons_sf, sparse = TRUE)
    region_id_vec <- rep(NA_character_, nrow(df))

    for (j in seq_along(valid_idx)) {
      hit <- inside_all[[j]]
      if (length(hit) > 0) {
        region_id_vec[valid_idx[j]] <- accepted_polygons_sf$core_id[hit[1]]
      }
    }

    df$NC_region_id <- region_id_vec
    df$NC_like_region <- !is.na(df$NC_region_id)
  }

  max_core_area <- max(cluster_summary$polygon_area, na.rm = TRUE)
  if (!is.finite(max_core_area)) max_core_area <- 0

  accepted_area <- cluster_summary$polygon_area[cluster_summary$area_pass_cutoff]
  mean_core_area <- mean(accepted_area, na.rm = TRUE)
  if (!is.finite(mean_core_area)) mean_core_area <- NA_real_

  total_core_area <- sum(accepted_area, na.rm = TRUE)

  fov_summary <- tibble(
    FOV = image_use,
    severity = severity_label,
    NC_like_positive = max_core_area >= min_core_area_cutoff,
    NC_area_cutoff = min_core_area_cutoff,
    n_candidate_core_regions = nrow(cluster_summary),
    n_core_regions = sum(cluster_summary$area_pass_cutoff, na.rm = TRUE),
    max_core_area = max_core_area,
    mean_core_area = mean_core_area,
    total_core_area = total_core_area,
    n_total_cells = nrow(df),
    n_myeloid_cells = nrow(myeloid_df),
    n_core_like_cells = sum(df$NC_like_region, na.rm = TRUE),
    n_myeloid_core_like_cells = sum(df$NC_like_region & df$is_myeloid_transfer, na.rm = TRUE),
    n_foam_core_like_cells = sum(df$NC_like_region & df$is_foam_transfer, na.rm = TRUE)
  )

  list(
    cell_df = df,
    polygons_sf = accepted_polygons_sf,
    polygons_df = accepted_polygons_df,
    candidate_polygons_sf = candidate_polygons_sf,
    candidate_polygons_df = candidate_polygons_df,
    cluster_summary = cluster_summary,
    fov_summary = fov_summary
  )
}

detect_necrotic_core_like_selected_fovs <- function(cell_df, image_names) {
  per_fov <- list()

  message("NC_POLYGON_BUFFER used exactly once during polygon construction: ", NC_POLYGON_BUFFER)

  for (i in seq_along(image_names)) {
    image_use <- image_names[i]
    message("Detecting core-like regions in: ", image_use)

    df_one <- cell_df %>%
      filter(as.character(image) == as.character(image_use))

    per_fov[[image_use]] <- do.call(
      detect_necrotic_core_like_one_fov,
      c(
        list(
          df = df_one,
          image_use = image_use,
          min_core_area_cutoff = NC_CORE_AREA_CUTOFF,
          polygon_buffer = NC_POLYGON_BUFFER,
          seed = NC_RANDOM_SEED + i
        ),
        NC_LOCAL_PARAMS,
        NC_CLUSTER_PARAMS
      )
    )
  }

  fov_summary <- bind_rows(lapply(per_fov, function(x) x$fov_summary))
  cluster_summary <- bind_rows(lapply(per_fov, function(x) x$cluster_summary))

  list(
    per_fov = per_fov,
    fov_summary = fov_summary,
    cluster_summary = cluster_summary
  )
}

nc_res <- detect_necrotic_core_like_selected_fovs(cell_df, TARGET_FOVS)

nc_fov_summary_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_NC_like_fov_summary.csv"))
nc_cluster_summary_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_NC_like_cluster_summary.csv"))

write.csv(nc_res$fov_summary, nc_fov_summary_out, row.names = FALSE)
write.csv(nc_res$cluster_summary, nc_cluster_summary_out, row.names = FALSE)

message("Saved FOV summary: ", nc_fov_summary_out)
message("Saved cluster summary: ", nc_cluster_summary_out)

message("========== Core-like FOV summary ==========")
print(nc_res$fov_summary, n = Inf)

message("========== Core-like cluster summary ==========")
print(nc_res$cluster_summary, n = Inf)

accepted_area_all <- nc_res$cluster_summary$polygon_area[nc_res$cluster_summary$area_pass_cutoff]
message("Accepted core-like regions average size: ", mean(accepted_area_all, na.rm = TRUE))

# ============================================================
# 9. 构建 accepted core polygon 和 facet 信息
# ============================================================

cell_df_with_nc <- bind_rows(lapply(TARGET_FOVS, function(img) {
  nc_res$per_fov[[img]]$cell_df
}))

cell_df_with_nc$image <- as.character(cell_df_with_nc$image)

fov_severity <- bind_rows(lapply(TARGET_FOVS, function(img) {
  df_img <- cell_df_with_nc %>% filter(as.character(image) == as.character(img))
  tibble(image = as.character(img), severity = get_fov_severity_label(df_img))
}))

facet_levels <- fov_severity %>%
  mutate(facet_label = paste0(image, " | Severity: ", severity)) %>%
  arrange(match(image, TARGET_FOVS)) %>%
  pull(facet_label)

cell_df_with_nc <- cell_df_with_nc %>%
  left_join(fov_severity, by = "image", suffix = c("", ".fov")) %>%
  mutate(
    facet_label = paste0(image, " | Severity: ", severity.fov),
    facet_label = factor(facet_label, levels = facet_levels),
    ISG_score_raw = as.numeric(ISG_like_score),
    ISG_score_plot = ifelse(
      is_myeloid_transfer & is.finite(ISG_score_raw),
      ISG_score_raw,
      NA_real_
    )
  )

cell_df_with_nc$ISG_score_plot <- winsorize_vec(cell_df_with_nc$ISG_score_plot, probs = c(0.01, 0.99))

core_sf_list <- lapply(TARGET_FOVS, function(img) {
  poly <- nc_res$per_fov[[img]]$polygons_sf
  if (is.null(poly) || nrow(poly) == 0) return(NULL)
  poly$image <- as.character(img)
  if (!"core_id" %in% colnames(poly)) {
    poly$core_id <- paste0(img, "_core", seq_len(nrow(poly)))
  }
  poly
})

core_sf_list <- core_sf_list[!vapply(core_sf_list, is.null, logical(1))]

if (length(core_sf_list) > 0) {
  core_sf <- do.call(rbind, core_sf_list)
  core_poly_df <- polygon_sf_to_df(core_sf)
} else {
  core_sf <- NULL
  core_poly_df <- tibble(
    X = numeric(),
    Y = numeric(),
    group_path = character(),
    core_id = character(),
    image = character(),
    facet_label = factor(levels = facet_levels)
  )
  warning("目标 FOV 中没有 accepted core-like polygon。空间图不会显示虚线边界。")
}

if (nrow(core_poly_df) > 0) {
  core_poly_df <- core_poly_df %>%
    mutate(image = as.character(image)) %>%
    left_join(fov_severity, by = "image") %>%
    mutate(
      facet_label = paste0(image, " | Severity: ", severity),
      facet_label = factor(facet_label, levels = facet_levels)
    )
}

# ============================================================
# 10. 三分类：Core-colocalized / Peri-core adjacent / Other
#     注意：这里不再使用 polygon_buffer，不再 st_buffer()
# ============================================================

compute_core_relation_one_fov_no_buffer <- function(
    df_img,
    core_sf_all,
    image_use,
    adjacent_distance = ADJACENT_DISTANCE
) {
  df_img <- df_img %>%
    mutate(
      has_core = FALSE,
      inside_core_polygon = FALSE,
      distance_to_core = NA_real_,
      spatial_relation_to_core = "Other region",
      spatial_region_3class = factor("Other region", levels = region_levels3),
      proximal_to_core = FALSE
    )

  if (is.null(core_sf_all) || nrow(core_sf_all) == 0) {
    return(df_img)
  }

  poly_img <- core_sf_all[
    as.character(core_sf_all$image) == as.character(image_use),
    ,
    drop = FALSE
  ]

  if (is.null(poly_img) || nrow(poly_img) == 0) {
    return(df_img)
  }

  df_img <- df_img %>%
    mutate(
      x = as.numeric(x),
      y = as.numeric(y)
    )

  valid_xy <- is.finite(df_img$x) & is.finite(df_img$y)

  if (sum(valid_xy) == 0) {
    warning("No valid x/y coordinates in ", image_use)
    return(df_img)
  }

  pts_sf <- sf::st_as_sf(
    df_img[valid_xy, , drop = FALSE],
    coords = c("x", "y"),
    remove = FALSE
  )

  poly_union <- suppressWarnings(sf::st_union(poly_img))

  inside_flag <- lengths(sf::st_within(pts_sf, poly_union, sparse = TRUE)) > 0

  dist_to_core <- as.numeric(sf::st_distance(pts_sf, poly_union))
  dist_to_core[!is.finite(dist_to_core)] <- NA_real_

  idx_valid <- which(valid_xy)

  df_img$has_core <- TRUE
  df_img$inside_core_polygon[idx_valid] <- inside_flag
  df_img$distance_to_core[idx_valid] <- dist_to_core

  df_img$spatial_relation_to_core[idx_valid[inside_flag]] <- "Core-colocalized"

  df_img$spatial_relation_to_core[
    idx_valid[!inside_flag & !is.na(dist_to_core) & dist_to_core <= adjacent_distance]
  ] <- "Peri-core adjacent"

  df_img$spatial_relation_to_core[
    idx_valid[!inside_flag & !is.na(dist_to_core) & dist_to_core > adjacent_distance]
  ] <- "Other region"

  df_img$spatial_region_3class <- factor(
    df_img$spatial_relation_to_core,
    levels = region_levels3
  )

  df_img$proximal_to_core <- df_img$spatial_region_3class %in% c(
    "Core-colocalized",
    "Peri-core adjacent"
  )

  df_img
}

relation_df <- bind_rows(lapply(TARGET_FOVS, function(img) {
  df_img <- cell_df_with_nc %>%
    filter(as.character(image) == as.character(img))

  compute_core_relation_one_fov_no_buffer(
    df_img = df_img,
    core_sf_all = core_sf,
    image_use = img,
    adjacent_distance = ADJACENT_DISTANCE
  )
}))

relation_df <- relation_df %>%
  mutate(
    spatial_region_3class = factor(
      as.character(spatial_region_3class),
      levels = region_levels3
    )
  )

message("========== Spatial region counts, all cells ==========")
print(table(relation_df$image, relation_df$spatial_region_3class, useNA = "ifany"))

message("========== Spatial region counts, myeloid cells only ==========")
print(table(
  relation_df$image[relation_df$is_myeloid_transfer],
  relation_df$spatial_region_3class[relation_df$is_myeloid_transfer],
  useNA = "ifany"
))

cell_relation_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_cell_level_core_distance_region3_ISG_score.csv"))
if (SAVE_CELL_LEVEL_TABLES) {
  write.csv(relation_df, cell_relation_out, row.names = FALSE)
  message("Saved cell-level relation table: ", cell_relation_out)
} else {
  message("Skipped cell-level relation table; set AST_SAVE_CELL_LEVEL_TABLES=1 to write it.")
}

spatial_region_summary_all <- relation_df %>%
  group_by(image, spatial_region_3class) %>%
  summarise(
    n_cells = n(),
    n_myeloid = sum(is_myeloid_transfer, na.rm = TRUE),
    n_non_myeloid = sum(!is_myeloid_transfer, na.rm = TRUE),
    mean_distance_to_core = mean(distance_to_core, na.rm = TRUE),
    median_distance_to_core = median(distance_to_core, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(image, spatial_region_3class)

spatial_region_summary_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_spatial_region3_summary_ALL_cells.csv"))
write.csv(spatial_region_summary_all, spatial_region_summary_out, row.names = FALSE)
message("Saved all-cell spatial region summary: ", spatial_region_summary_out)

# ============================================================
# 11. 空间图
# ============================================================

n_fov <- length(TARGET_FOVS)
spatial_width <- ifelse(n_fov <= 3, 11.5, min(24, 4.2 * min(n_fov, 6)))
spatial_height <- ifelse(n_fov <= 3, 6.5, max(6.5, 3.4 * ceiling(n_fov / 4)))

p_spatial_region_3class_all_cells <- ggplot() +
  geom_point(
    data = relation_df,
    aes(x = x, y = y, color = spatial_region_3class),
    size = 0.12,
    alpha = 0.70
  ) +
  geom_path(
    data = core_poly_df,
    aes(x = X, y = Y, group = group_path),
    color = "black",
    linewidth = 0.48,
    linetype = "22",
    alpha = 0.95
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  scale_color_manual(
    values = region_colors3,
    drop = FALSE,
    name = "Spatial region"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Spatial distribution of core-colocalized, peri-core adjacent and other cells",
    subtitle = paste0(
      "Peri-core adjacent = non-core cells within ",
      ADJACENT_DISTANCE,
      " spatial units from accepted core-like regions; dashed boundaries indicate accepted core-like regions"
    )
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 14, r = 22, b = 14, l = 22)
  )

print(p_spatial_region_3class_all_cells)

save_plot_both(
  p_spatial_region_3class_all_cells,
  paste0(OUT_TAG, "_spatial_distribution_region3_ALL_cells"),
  width = spatial_width,
  height = spatial_height
)

p_spatial_region_3class_myeloid <- ggplot() +
  geom_point(
    data = relation_df %>% filter(!is_myeloid_transfer),
    aes(x = x, y = y),
    color = "grey90",
    size = 0.06,
    alpha = 0.30
  ) +
  geom_point(
    data = relation_df %>% filter(is_myeloid_transfer),
    aes(x = x, y = y, color = spatial_region_3class),
    size = 0.24,
    alpha = 0.92
  ) +
  geom_path(
    data = core_poly_df,
    aes(x = X, y = Y, group = group_path),
    color = "black",
    linewidth = 0.48,
    linetype = "22",
    alpha = 0.95
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  scale_color_manual(
    values = region_colors3,
    drop = FALSE,
    name = "Spatial region"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Spatial distribution of transferred myeloid cells across core-related regions",
    subtitle = "Grey background: non-myeloid or unassigned cells; dashed boundaries: accepted core-like regions"
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 14, r = 22, b = 14, l = 22)
  )

print(p_spatial_region_3class_myeloid)

save_plot_both(
  p_spatial_region_3class_myeloid,
  paste0(OUT_TAG, "_spatial_distribution_region3_MYELOID_cells"),
  width = spatial_width,
  height = spatial_height
)

p_myeloid_isg_heatmap <- ggplot() +
  geom_point(
    data = relation_df %>% filter(!is_myeloid_transfer | is.na(ISG_score_plot)),
    aes(x = x, y = y),
    color = "grey86",
    size = 0.07,
    alpha = 0.42
  ) +
  geom_point(
    data = relation_df %>% filter(is_myeloid_transfer, !is.na(ISG_score_plot)),
    aes(x = x, y = y, color = ISG_score_plot),
    size = 0.22,
    alpha = 0.95
  ) +
  geom_path(
    data = core_poly_df,
    aes(x = X, y = Y, group = group_path),
    color = "black",
    linewidth = 0.50,
    linetype = "22",
    alpha = 0.95
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  scale_color_gradientn(
    colors = c("#2166AC", "#F7F7F7", "#B2182B"),
    name = "A3A-like myeloid score"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Spatial heatmap of APOBEC3A-like myeloid score in transferred myeloid cells",
    subtitle = "Non-myeloid or unassigned cells are shown in grey; dashed boundaries indicate accepted core-like regions"
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 14, r = 22, b = 14, l = 22)
  )

print(p_myeloid_isg_heatmap)

save_plot_both(
  p_myeloid_isg_heatmap,
  paste0(OUT_TAG, "_myeloid_ISG_score_spatial_heatmap"),
  width = spatial_width,
  height = spatial_height
)



# ============================================================
# 11B. 恢复图 2：髓系细胞 APOBEC3A-like / APOBEC3A-like myeloid score + core-like 边界
#      说明：这里把 winsorized ISG_score_plot 再 rescale 到 0-1，
#      以匹配灰-白-红色标和 limits = c(0, 1)。
# ============================================================

relation_df$ISG_score_plot01 <- NA_real_
valid_isg_score <- relation_df$is_myeloid_transfer & is.finite(relation_df$ISG_score_plot)

if (sum(valid_isg_score, na.rm = TRUE) > 0) {
  score_range <- range(relation_df$ISG_score_plot[valid_isg_score], na.rm = TRUE)

  if (is.finite(score_range[1]) && is.finite(score_range[2]) && diff(score_range) > 0) {
    relation_df$ISG_score_plot01[valid_isg_score] <- scales::rescale(
      relation_df$ISG_score_plot[valid_isg_score],
      to = c(0, 1),
      from = score_range
    )
  } else {
    relation_df$ISG_score_plot01[valid_isg_score] <- 0.5
  }
}

p_myeloid_isg_region_map <- ggplot() +
  geom_point(
    data = relation_df %>% filter(!is_myeloid_transfer | is.na(ISG_score_plot01)),
    aes(x = x, y = y),
    color = "grey88",
    size = 0.07,
    alpha = 0.42
  ) +
  geom_point(
    data = relation_df %>% filter(is_myeloid_transfer, !is.na(ISG_score_plot01)),
    aes(x = x, y = y, color = ISG_score_plot01),
    size = 0.22,
    alpha = 0.95
  ) +
  geom_path(
    data = core_poly_df,
    aes(x = X, y = Y, group = group_path),
    color = "black",
    linewidth = 0.5,
    linetype = "22",
    alpha = 0.95
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  scale_color_gradientn(
    colors = c("gray", "#F7F7F7", "#B2182B"),
    values = scales::rescale(c(0, 0.25, 1), from = c(0, 1)),
    limits = c(0, 1),
    oob = scales::squish,
    name = "A3A-like myeloid score"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "APOBEC3A-like myeloid score around core-like regions",
    subtitle = "Grey background: non-myeloid or unassigned cells; dashed boundaries: accepted core-like regions"
  ) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(
      face = "bold",
      size = 11,
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13,
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0.5,
      margin = margin(b = 6)
    ),
    legend.position = "right",
    plot.margin = margin(t = 14, r = 22, b = 14, l = 22)
  )

print(p_myeloid_isg_region_map)

save_plot_both(
  p_myeloid_isg_region_map,
  paste0(OUT_TAG, "_myeloid_A3Alike_score_region_map"),
  width = spatial_width,
  height = spatial_height
)

# ============================================================
# 11C. 每个 FOV 单独保存空间分布图
#      同时保留上面的 facet 合并图，方便总览。
# ============================================================

save_single_fov_spatial_plots <- function(img) {
  df_one <- relation_df %>% filter(as.character(image) == as.character(img))
  poly_one <- core_poly_df %>% filter(as.character(image) == as.character(img))

  if (nrow(df_one) == 0) {
    warning("Skip single-FOV plot because no cells found: ", img)
    return(invisible(NULL))
  }

  img_safe <- sanitize_filename(as.character(img))
  title_suffix <- unique(as.character(df_one$facet_label))
  title_suffix <- title_suffix[!is.na(title_suffix)]
  if (length(title_suffix) == 0) title_suffix <- as.character(img)
  title_suffix <- title_suffix[1]

  p_one_all <- ggplot() +
    geom_point(
      data = df_one,
      aes(x = x, y = y, color = spatial_region_3class),
      size = 0.16,
      alpha = 0.72
    ) +
    geom_path(
      data = poly_one,
      aes(x = X, y = Y, group = group_path),
      color = "black",
      linewidth = 0.55,
      linetype = "22",
      alpha = 0.95
    ) +
    scale_y_reverse() +
    coord_cartesian(clip = "off") +
    scale_color_manual(
      values = region_colors3,
      drop = FALSE,
      name = "Spatial region"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Spatial region classification: ", title_suffix),
      subtitle = paste0("Peri-core adjacent = non-core cells within ", ADJACENT_DISTANCE, " spatial units")
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(t = 14, r = 20, b = 14, l = 20)
    )

  save_plot_both(
    p_one_all,
    paste0(OUT_TAG, "_", img_safe, "_spatial_distribution_region3_ALL_cells"),
    width = 6.2,
    height = 5.4
  )

  p_one_myeloid_region <- ggplot() +
    geom_point(
      data = df_one %>% filter(!is_myeloid_transfer),
      aes(x = x, y = y),
      color = "grey90",
      size = 0.07,
      alpha = 0.30
    ) +
    geom_point(
      data = df_one %>% filter(is_myeloid_transfer),
      aes(x = x, y = y, color = spatial_region_3class),
      size = 0.30,
      alpha = 0.92
    ) +
    geom_path(
      data = poly_one,
      aes(x = X, y = Y, group = group_path),
      color = "black",
      linewidth = 0.55,
      linetype = "22",
      alpha = 0.95
    ) +
    scale_y_reverse() +
    coord_cartesian(clip = "off") +
    scale_color_manual(
      values = region_colors3,
      drop = FALSE,
      name = "Spatial region"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Transferred myeloid cells across core-related regions: ", title_suffix),
      subtitle = "Grey background: non-myeloid or unassigned cells; dashed boundaries: accepted core-like regions"
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(t = 14, r = 20, b = 14, l = 20)
    )

  save_plot_both(
    p_one_myeloid_region,
    paste0(OUT_TAG, "_", img_safe, "_spatial_distribution_region3_MYELOID_cells"),
    width = 6.2,
    height = 5.4
  )

  p_one_myeloid_score <- ggplot() +
    geom_point(
      data = df_one %>% filter(!is_myeloid_transfer | is.na(ISG_score_plot01)),
      aes(x = x, y = y),
      color = "grey88",
      size = 0.07,
      alpha = 0.42
    ) +
    geom_point(
      data = df_one %>% filter(is_myeloid_transfer, !is.na(ISG_score_plot01)),
      aes(x = x, y = y, color = ISG_score_plot01),
      size = 0.30,
      alpha = 0.95
    ) +
    geom_path(
      data = poly_one,
      aes(x = X, y = Y, group = group_path),
      color = "black",
      linewidth = 0.55,
      linetype = "22",
      alpha = 0.95
    ) +
    scale_y_reverse() +
    coord_cartesian(clip = "off") +
    scale_color_gradientn(
      colors = c("gray", "#F7F7F7", "#B2182B"),
      values = scales::rescale(c(0, 0.25, 1), from = c(0, 1)),
      limits = c(0, 1),
      oob = scales::squish,
      name = "A3A-like myeloid score"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("APOBEC3A-like myeloid score around core-like regions: ", title_suffix),
      subtitle = "Grey background: non-myeloid or unassigned cells; dashed boundaries: accepted core-like regions"
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right",
      plot.margin = margin(t = 14, r = 20, b = 14, l = 20)
    )

  save_plot_both(
    p_one_myeloid_score,
    paste0(OUT_TAG, "_", img_safe, "_myeloid_A3Alike_score_region_map"),
    width = 6.2,
    height = 5.4
  )

  invisible(NULL)
}

if (SAVE_PER_FOV_PLOTS) {
  invisible(lapply(TARGET_FOVS, save_single_fov_spatial_plots))
} else {
  message("Skipped per-FOV spatial plot export; set AST_SAVE_PER_FOV_PLOTS=1 to write it.")
}


# ============================================================
# 12. FOV core area 诊断图
# ============================================================

p_core_area <- nc_res$fov_summary %>%
  mutate(
    FOV = factor(FOV, levels = FOV[order(max_core_area)]),
    NC_like_status = ifelse(NC_like_positive, "Above cutoff", "Below cutoff")
  ) %>%
  ggplot(aes(x = FOV, y = max_core_area, fill = severity)) +
  geom_col(width = 0.75, alpha = 0.90) +
  geom_hline(
    yintercept = NC_CORE_AREA_CUTOFF,
    linetype = "22",
    linewidth = 0.5,
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "Maximum core-like region area",
    fill = "Severity",
    title = "FOV-level maximum core-like region area",
    subtitle = paste0("Dashed line: accepted core area cutoff = ", NC_CORE_AREA_CUTOFF)
  ) +
  theme_story(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_core_area)

save_plot_both(
  p_core_area,
  paste0(OUT_TAG, "_FOV_max_core_area_barplot"),
  width = max(7, min(14, 0.55 * length(TARGET_FOVS) + 4)),
  height = 5.2
)

# ============================================================
# 13. 髓系 APOBEC3A-like myeloid score 统计
# ============================================================

isg_region_test_df <- relation_df %>%
  filter(
    is_myeloid_transfer,
    !is.na(ISG_like_score),
    is.finite(ISG_like_score),
    !is.na(spatial_region_3class)
  ) %>%
  mutate(
    spatial_region_3class = factor(
      as.character(spatial_region_3class),
      levels = region_levels3
    ),
    log10_distance_to_core = log10(distance_to_core + 1)
  )

if (nrow(isg_region_test_df) == 0) {
  stop("No transferred myeloid cells with valid ISG_like_score and spatial_region_3class.")
}

isg_region_cell_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_region3_cell_table.csv"))
write.csv(isg_region_test_df, isg_region_cell_out, row.names = FALSE)
message("Saved APOBEC3A-like myeloid score region cell table: ", isg_region_cell_out)

message("========== Myeloid cells used for ISG region test ==========")
print(table(isg_region_test_df$image, isg_region_test_df$spatial_region_3class, useNA = "ifany"))

summarise_isg_by_region <- function(df, label = "ALL_SELECTED_FOVS") {
  df %>%
    group_by(spatial_region_3class) %>%
    summarise(
      image = label,
      n_cells = n(),
      mean_ISG_score = mean(ISG_like_score, na.rm = TRUE),
      median_ISG_score = median(ISG_like_score, na.rm = TRUE),
      sd_ISG_score = sd(ISG_like_score, na.rm = TRUE),
      q25_ISG_score = quantile(ISG_like_score, 0.25, na.rm = TRUE),
      q75_ISG_score = quantile(ISG_like_score, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(image, everything())
}

run_kruskal_one <- function(df, label = "ALL_SELECTED_FOVS") {
  df <- df %>%
    filter(!is.na(spatial_region_3class), !is.na(ISG_like_score)) %>%
    droplevels()

  n_group <- n_distinct(df$spatial_region_3class)

  if (nrow(df) < 3 || n_group < 2) {
    return(tibble(
      image = label,
      n_cells = nrow(df),
      n_groups = n_group,
      kruskal_chisq = NA_real_,
      kruskal_df = NA_real_,
      kruskal_p = NA_real_
    ))
  }

  kt <- tryCatch(
    kruskal.test(ISG_like_score ~ spatial_region_3class, data = df),
    error = function(e) NULL
  )

  if (is.null(kt)) {
    return(tibble(
      image = label,
      n_cells = nrow(df),
      n_groups = n_group,
      kruskal_chisq = NA_real_,
      kruskal_df = NA_real_,
      kruskal_p = NA_real_
    ))
  }

  tibble(
    image = label,
    n_cells = nrow(df),
    n_groups = n_group,
    kruskal_chisq = unname(kt$statistic),
    kruskal_df = unname(kt$parameter),
    kruskal_p = kt$p.value
  )
}

run_pairwise_wilcox_one <- function(df, label = "ALL_SELECTED_FOVS") {
  df <- df %>%
    filter(!is.na(spatial_region_3class), !is.na(ISG_like_score)) %>%
    droplevels()

  groups <- levels(droplevels(df$spatial_region_3class))
  groups <- groups[groups %in% unique(as.character(df$spatial_region_3class))]

  if (length(groups) < 2) {
    return(tibble(
      image = character(),
      group1 = character(),
      group2 = character(),
      n1 = integer(),
      n2 = integer(),
      median1 = numeric(),
      median2 = numeric(),
      median_diff_group1_minus_group2 = numeric(),
      wilcox_p = numeric(),
      wilcox_p_adj_BH = numeric()
    ))
  }

  comps <- combn(groups, 2, simplify = FALSE)

  out <- lapply(comps, function(cp) {
    g1 <- cp[1]
    g2 <- cp[2]

    x <- df$ISG_like_score[df$spatial_region_3class == g1]
    y <- df$ISG_like_score[df$spatial_region_3class == g2]

    p <- if (length(x) >= 3 && length(y) >= 3) {
      tryCatch(
        wilcox.test(x, y, alternative = "two.sided")$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }

    tibble(
      image = label,
      group1 = g1,
      group2 = g2,
      n1 = length(x),
      n2 = length(y),
      median1 = median(x, na.rm = TRUE),
      median2 = median(y, na.rm = TRUE),
      median_diff_group1_minus_group2 = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
      wilcox_p = p
    )
  }) %>%
    bind_rows()

  out %>%
    mutate(wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH"))
}

isg_summary_fov <- bind_rows(lapply(TARGET_FOVS, function(img) {
  summarise_isg_by_region(
    isg_region_test_df %>% filter(as.character(image) == img),
    label = img
  )
}))

isg_summary_all <- summarise_isg_by_region(isg_region_test_df, label = "ALL_SELECTED_FOVS")
isg_summary <- bind_rows(isg_summary_fov, isg_summary_all)

isg_kruskal_fov <- bind_rows(lapply(TARGET_FOVS, function(img) {
  run_kruskal_one(
    isg_region_test_df %>% filter(as.character(image) == img),
    label = img
  )
}))

isg_kruskal_all <- run_kruskal_one(isg_region_test_df, label = "ALL_SELECTED_FOVS")
isg_kruskal <- bind_rows(isg_kruskal_fov, isg_kruskal_all)

isg_pairwise_fov <- bind_rows(lapply(TARGET_FOVS, function(img) {
  run_pairwise_wilcox_one(
    isg_region_test_df %>% filter(as.character(image) == img),
    label = img
  )
}))

isg_pairwise_all <- run_pairwise_wilcox_one(isg_region_test_df, label = "ALL_SELECTED_FOVS")
isg_pairwise <- bind_rows(isg_pairwise_fov, isg_pairwise_all)

isg_summary_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_region3_summary.csv"))
isg_kruskal_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_region3_kruskal.csv"))
isg_pairwise_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_region3_pairwise_wilcox_BH.csv"))

write.csv(isg_summary, isg_summary_out, row.names = FALSE)
write.csv(isg_kruskal, isg_kruskal_out, row.names = FALSE)
write.csv(isg_pairwise, isg_pairwise_out, row.names = FALSE)

message("Saved APOBEC3A-like myeloid score region summary: ", isg_summary_out)
message("Saved APOBEC3A-like myeloid score Kruskal-Wallis test: ", isg_kruskal_out)
message("Saved APOBEC3A-like myeloid score pairwise Wilcoxon test: ", isg_pairwise_out)

message("========== APOBEC3A-like myeloid score by spatial region ==========")
print(isg_summary, n = Inf)

message("========== Kruskal-Wallis test ==========")
print(isg_kruskal, n = Inf)

message("========== Pairwise Wilcoxon test, BH-adjusted ==========")
print(isg_pairwise, n = Inf)

# ============================================================
# 14. APOBEC3A-like myeloid score 统计图
# ============================================================

p_isg_region_violin_fov <- ggplot(
  isg_region_test_df,
  aes(x = spatial_region_3class, y = ISG_like_score, fill = spatial_region_3class)
) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.90, linewidth = 0.25) +
  geom_jitter(width = 0.12, size = 0.20, alpha = 0.20) +
  facet_wrap(~ facet_label, scales = "free_y") +
  scale_fill_manual(values = region_colors3, drop = FALSE) +
  labs(
    x = NULL,
    y = "APOBEC3A-like myeloid score",
    fill = "Spatial region",
    title = "Differential APOBEC3A-like myeloid score across core-colocalized, peri-core and other myeloid regions",
    subtitle = "Statistics: Kruskal-Wallis across three regions; pairwise Wilcoxon rank-sum with BH correction"
  ) +
  theme_story(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_isg_region_violin_fov)

save_plot_both(
  p_isg_region_violin_fov,
  paste0(OUT_TAG, "_myeloid_ISG_score_region3_violin_boxplot_by_FOV"),
  width = spatial_width,
  height = spatial_height
)

p_isg_region_violin_all <- ggplot(
  isg_region_test_df,
  aes(x = spatial_region_3class, y = ISG_like_score, fill = spatial_region_3class)
) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.90, linewidth = 0.25) +
  geom_jitter(width = 0.14, size = 0.18, alpha = 0.18) +
  scale_fill_manual(values = region_colors3, drop = FALSE) +
  labs(
    x = NULL,
    y = "APOBEC3A-like myeloid score",
    fill = "Spatial region",
    title = "APOBEC3A-like myeloid score is compared across spatial regions",
    subtitle = "All selected FOVs combined"
  ) +
  theme_story(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_isg_region_violin_all)

save_plot_both(
  p_isg_region_violin_all,
  paste0(OUT_TAG, "_myeloid_ISG_score_region3_violin_boxplot_ALL_SELECTED_FOVS"),
  width = 6.6,
  height = 5.4
)

# ============================================================
# 15. 可选：如果 xen 里有 UMAP，则画 APOBEC3A-like myeloid score UMAP 热图
# ============================================================

plot_myeloid_isg_umap_if_available <- function(obj) {
  red_use <- intersect(
    c("umap", "UMAP", "integrated.umap", "harmony.umap", "ref.umap"),
    names(obj@reductions)
  )

  if (length(red_use) == 0) {
    message("No UMAP reduction found in xen object. Skip APOBEC3A-like myeloid score UMAP heatmap.")
    return(invisible(NULL))
  }

  red_use <- red_use[1]
  emb <- as.data.frame(Embeddings(obj, reduction = red_use))

  if (ncol(emb) < 2) {
    message("UMAP reduction has < 2 dimensions. Skip APOBEC3A-like myeloid score UMAP heatmap.")
    return(invisible(NULL))
  }

  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
  emb$cell <- rownames(emb)

  meta <- obj@meta.data %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("cell")

  plot_df <- emb %>%
    left_join(meta, by = "cell") %>%
    add_myeloid_and_foam_flags(ann_col = AFTER_COL) %>%
    mutate(
      ISG_score_raw = as.numeric(ISG_like_score),
      ISG_score_plot = ifelse(
        is_myeloid_transfer & is.finite(ISG_score_raw),
        ISG_score_raw,
        NA_real_
      )
    )

  plot_df$ISG_score_plot <- winsorize_vec(plot_df$ISG_score_plot, probs = c(0.01, 0.99))

  p_umap <- ggplot() +
    geom_point(
      data = plot_df %>% filter(!is_myeloid_transfer | is.na(ISG_score_plot)),
      aes(x = UMAP_1, y = UMAP_2),
      color = "grey86",
      size = 0.08,
      alpha = 0.35
    ) +
    geom_point(
      data = plot_df %>% filter(is_myeloid_transfer, !is.na(ISG_score_plot)),
      aes(x = UMAP_1, y = UMAP_2, color = ISG_score_plot),
      size = 0.18,
      alpha = 0.90
    ) +
    scale_color_gradientn(
      colors = c("#2166AC", "#F7F7F7", "#B2182B"),
      name = "A3A-like myeloid score"
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = "UMAP heatmap of APOBEC3A-like myeloid score in transferred myeloid cells",
      subtitle = paste0("Reduction: ", red_use, "; non-myeloid or unassigned cells are shown in grey")
    ) +
    theme_classic(base_size = 9) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  print(p_umap)

  save_plot_both(
    p_umap,
    paste0(OUT_TAG, "_myeloid_ISG_score_UMAP_heatmap"),
    width = 6.5,
    height = 5.8
  )

  invisible(p_umap)
}

plot_myeloid_isg_umap_if_available(xen)

# ============================================================
# 16. 最终输出提示
# ============================================================

final_print <- isg_kruskal %>%
  select(
    image,
    n_cells,
    n_groups,
    kruskal_chisq,
    kruskal_df,
    kruskal_p
  )

final_print_out <- file.path(TABLE_DIR, paste0(OUT_TAG, "_final_ISG_score_region3_test_summary.csv"))
write.csv(final_print, final_print_out, row.names = FALSE)

message("\n========== Final APOBEC3A-like myeloid score region test summary ==========")
print(final_print, n = Inf)

message("\n========== Main outputs ==========")
message("1) ", file.path(FIG_DIR, paste0(OUT_TAG, "_spatial_distribution_region3_ALL_cells.pdf/png")))
message("2) ", file.path(FIG_DIR, paste0(OUT_TAG, "_spatial_distribution_region3_MYELOID_cells.pdf/png")))
message("3) ", file.path(FIG_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_spatial_heatmap.pdf/png")))
message("4) ", file.path(FIG_DIR, paste0(OUT_TAG, "_FOV_max_core_area_barplot.pdf/png")))
message("5) ", file.path(FIG_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_region3_violin_boxplot_by_FOV.pdf/png")))
message("6) ", file.path(FIG_DIR, paste0(OUT_TAG, "_myeloid_ISG_score_region3_violin_boxplot_ALL_SELECTED_FOVS.pdf/png")))
message("7) ", isg_summary_out)
message("8) ", isg_kruskal_out)
message("9) ", isg_pairwise_out)
message("10) ", cell_relation_out)
message("11) ", nc_fov_summary_out)
message("12) ", nc_cluster_summary_out)
message("13) ", final_print_out)


# ============================================================
# 15. FINAL ADD-ON:
#     Use ALL FOVs for final spatial-region statistics
#     + combined p_isg_region_bar_all
#     + one p_spatial_region_myeloid per FOV
#
# 关键点：
# 1) 不再使用前面只包含 3 个 FOV 的 relation_df 做最终统计；
# 2) 这里单独构建 relation_df_final；
# 3) TARGET_FOVS_FINAL 默认 = names(xen@images)，即所有 FOV；
# 4) 不覆盖前面的 TARGET_FOVS / relation_df / nc_res；
# 5) core polygon 仍然调用你前面已有 detect_necrotic_core_like_selected_fovs()，
#    不额外 st_buffer，不重复定义 classify 函数。
# ============================================================

if (RUN_FINAL_ALL_FOV_SUMMARY) {

# ------------------------------------------------------------
# 15.0 最终统计使用全部 FOV
# ------------------------------------------------------------

if (!exists("xen")) {
  stop("找不到 xen 对象。要统计所有 FOV，必须在前面已经读取 xen。")
}

if (length(xen@images) == 0) {
  stop("xen@images 为空，无法获取所有 FOV。")
}

if (exists("natural_fov_order")) {
  TARGET_FOVS_FINAL <- natural_fov_order(names(xen@images))
} else {
  TARGET_FOVS_FINAL <- names(xen@images)
}

message("========== FINAL SECTION: use ALL FOVs ==========")
message("Number of FOVs used in final analysis: ", length(TARGET_FOVS_FINAL))
message("FOVs used: ", paste(TARGET_FOVS_FINAL, collapse = ", "))

OUT_TAG_FINAL <- paste0(OUT_TAG, "_ALL_FOV_FINAL")

# ------------------------------------------------------------
# 15.1 为所有 FOV 重新提取 cell table
# ------------------------------------------------------------

message("Extracting all FOV coordinates and metadata for final analysis...")

cell_df_final <- bind_rows(lapply(TARGET_FOVS_FINAL, function(img) {
  message("Extracting final FOV: ", img)
  get_one_fov_df(xen, img)
}))

cell_df_final <- add_myeloid_and_foam_flags(cell_df_final, ann_col = AFTER_COL)

message("Final cell counts by FOV:")
print(table(cell_df_final$image, useNA = "ifany"))

message("Final transferred myeloid counts by FOV:")
print(table(cell_df_final$image, cell_df_final$is_myeloid_transfer, useNA = "ifany"))

message("Final Foam/LAM counts by FOV:")
print(table(cell_df_final$image, cell_df_final$is_foam_transfer, useNA = "ifany"))

# ------------------------------------------------------------
# 15.2 对所有 FOV 重新识别 core-like region
#     注意：这里调用的是你前面已经定义好的函数
# ------------------------------------------------------------

message("Detecting core-like regions for ALL FOVs in final analysis...")

nc_res_final <- detect_necrotic_core_like_selected_fovs(
  cell_df = cell_df_final,
  image_names = TARGET_FOVS_FINAL
)

nc_fov_summary_final_out <- file.path(
  TABLE_DIR,
  paste0(OUT_TAG_FINAL, "_NC_like_fov_summary.csv")
)

nc_cluster_summary_final_out <- file.path(
  TABLE_DIR,
  paste0(OUT_TAG_FINAL, "_NC_like_cluster_summary.csv")
)

write.csv(nc_res_final$fov_summary, nc_fov_summary_final_out, row.names = FALSE)
write.csv(nc_res_final$cluster_summary, nc_cluster_summary_final_out, row.names = FALSE)

message("Saved final all-FOV FOV summary: ", nc_fov_summary_final_out)
message("Saved final all-FOV cluster summary: ", nc_cluster_summary_final_out)

# ------------------------------------------------------------
# 15.3 构建所有 FOV 的 cell_df_with_nc_final 和 facet 信息
# ------------------------------------------------------------

cell_df_with_nc_final <- bind_rows(lapply(TARGET_FOVS_FINAL, function(img) {
  nc_res_final$per_fov[[img]]$cell_df
}))

cell_df_with_nc_final$image <- as.character(cell_df_with_nc_final$image)

fov_severity_final <- bind_rows(lapply(TARGET_FOVS_FINAL, function(img) {
  df_img <- cell_df_with_nc_final %>%
    filter(as.character(image) == as.character(img))

  tibble(
    image = as.character(img),
    severity_fov = get_fov_severity_label(df_img)
  )
}))

facet_levels_final <- fov_severity_final %>%
  mutate(
    facet_label = paste0(image, " | Severity: ", severity_fov)
  ) %>%
  arrange(match(image, TARGET_FOVS_FINAL)) %>%
  pull(facet_label)

cell_df_with_nc_final <- cell_df_with_nc_final %>%
  select(-any_of(c("severity_fov", "facet_label"))) %>%
  left_join(fov_severity_final, by = "image") %>%
  mutate(
    facet_label = paste0(image, " | Severity: ", severity_fov),
    facet_label = factor(facet_label, levels = facet_levels_final)
  )

# 准备 ISG_score_plot
if (!"ISG_score_plot" %in% colnames(cell_df_with_nc_final)) {
  if ("ISG_like_score" %in% colnames(cell_df_with_nc_final)) {
    cell_df_with_nc_final$ISG_score_plot <- cell_df_with_nc_final$ISG_like_score
  } else if ("A3A_like_myeloid_score" %in% colnames(cell_df_with_nc_final)) {
    cell_df_with_nc_final$ISG_score_plot <- cell_df_with_nc_final$A3A_like_myeloid_score
  } else {
    stop(
      "cell_df_with_nc_final 中没有 ISG_like_score / A3A_like_myeloid_score，无法生成 ISG_score_plot。\n",
      "当前列名：\n",
      paste(colnames(cell_df_with_nc_final), collapse = ", ")
    )
  }
}

cell_df_with_nc_final$ISG_score_plot <- as.numeric(cell_df_with_nc_final$ISG_score_plot)
cell_df_with_nc_final$ISG_score_plot[!cell_df_with_nc_final$is_myeloid_transfer] <- NA_real_
cell_df_with_nc_final$ISG_score_plot <- winsorize_vec(
  cell_df_with_nc_final$ISG_score_plot,
  probs = c(0.01, 0.99)
)

# ------------------------------------------------------------
# 15.4 构建所有 FOV 的 accepted core polygon 边界
# ------------------------------------------------------------

core_sf_list_final <- lapply(TARGET_FOVS_FINAL, function(img) {
  poly <- nc_res_final$per_fov[[img]]$polygons_sf

  if (is.null(poly) || nrow(poly) == 0) {
    return(NULL)
  }

  poly$image <- as.character(img)

  if (!"core_id" %in% colnames(poly)) {
    poly$core_id <- paste0(img, "_core", seq_len(nrow(poly)))
  }

  poly
})

core_sf_list_final <- core_sf_list_final[
  !vapply(core_sf_list_final, is.null, logical(1))
]

if (length(core_sf_list_final) > 0) {
  core_sf_final <- do.call(rbind, core_sf_list_final)
  core_poly_df_final <- polygon_sf_to_df(core_sf_final)
} else {
  core_sf_final <- NULL
  core_poly_df_final <- tibble(
    X = numeric(),
    Y = numeric(),
    group_path = character(),
    core_id = character(),
    image = character(),
    facet_label = factor(character(), levels = facet_levels_final)
  )

  warning("所有 FOV 中没有 accepted core-like polygon。空间图不会显示虚线边界。")
}

if (nrow(core_poly_df_final) > 0) {
  core_poly_df_final <- core_poly_df_final %>%
    mutate(image = as.character(image)) %>%
    select(-any_of(c("severity", "severity.fov", "severity_fov", "facet_label"))) %>%
    left_join(fov_severity_final, by = "image") %>%
    mutate(
      facet_label = paste0(image, " | Severity: ", severity_fov),
      facet_label = factor(facet_label, levels = facet_levels_final)
    )
}

# ------------------------------------------------------------
# 15.5 重新构建所有 FOV 的 relation_df_final
# ------------------------------------------------------------

relation_df_final <- bind_rows(lapply(TARGET_FOVS_FINAL, function(img) {
  df_img <- cell_df_with_nc_final %>%
    filter(as.character(image) == as.character(img))

  compute_core_relation_one_fov_no_buffer(
    df_img = df_img,
    core_sf_all = core_sf_final,
    image_use = img,
    adjacent_distance = ADJACENT_DISTANCE
  )
}))

relation_df_final <- relation_df_final %>%
  mutate(
    spatial_region_3class = factor(
      as.character(spatial_region_3class),
      levels = c(
        "Core-colocalized",
        "Peri-core adjacent",
        "Other region"
      )
    )
  )

message("========== FINAL all-FOV spatial region counts, all cells ==========")
print(table(
  relation_df_final$image,
  relation_df_final$spatial_region_3class,
  useNA = "ifany"
))

message("========== FINAL all-FOV spatial region counts, myeloid cells only ==========")
print(table(
  relation_df_final$image[relation_df_final$is_myeloid_transfer],
  relation_df_final$spatial_region_3class[relation_df_final$is_myeloid_transfer],
  useNA = "ifany"
))

relation_df_final_out <- file.path(
  TABLE_DIR,
  paste0(OUT_TAG_FINAL, "_cell_level_core_distance_region3_ISG_score.csv")
)

write.csv(relation_df_final, relation_df_final_out, row.names = FALSE)
message("Saved final all-FOV cell-level relation table: ", relation_df_final_out)

# ------------------------------------------------------------
# 15.6 准备所有 FOV 合并后的 APOBEC3A-like myeloid score 统计数据
# ------------------------------------------------------------

isg_region_test_df_final <- relation_df_final %>%
  filter(
    is_myeloid_transfer,
    !is.na(ISG_score_plot),
    is.finite(ISG_score_plot),
    !is.na(spatial_region_3class)
  ) %>%
  mutate(
    spatial_region_3class = factor(
      as.character(spatial_region_3class),
      levels = c(
        "Core-colocalized",
        "Peri-core adjacent",
        "Other region"
      )
    )
  )

if (nrow(isg_region_test_df_final) == 0) {
  stop("isg_region_test_df_final 为空：没有可用于统计的髓系细胞 APOBEC3A-like myeloid score。")
}

isg_region_test_final_out <- file.path(
  TABLE_DIR,
  paste0(OUT_TAG_FINAL, "_myeloid_ISG_score_spatial_region_cell_level_ALL_FOVS.csv")
)

write.csv(isg_region_test_df_final, isg_region_test_final_out, row.names = FALSE)
message("Saved final all-FOV ISG test table: ", isg_region_test_final_out)

isg_region_bar_summary_final <- isg_region_test_df_final %>%
  group_by(spatial_region_3class) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(ISG_score_plot, na.rm = TRUE),
    median_score = median(ISG_score_plot, na.rm = TRUE),
    sd_score = sd(ISG_score_plot, na.rm = TRUE),
    se_score = sd_score / sqrt(n_cells),
    q25_score = as.numeric(quantile(ISG_score_plot, 0.25, na.rm = TRUE)),
    q75_score = as.numeric(quantile(ISG_score_plot, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    spatial_region_3class = factor(
      as.character(spatial_region_3class),
      levels = c(
        "Core-colocalized",
        "Peri-core adjacent",
        "Other region"
      )
    )
  ) %>%
  arrange(spatial_region_3class)

isg_region_bar_summary_final_out <- file.path(
  TABLE_DIR,
  paste0(OUT_TAG_FINAL, "_myeloid_ISG_score_spatial_region_bar_summary_ALL_FOVS.csv")
)

write.csv(
  isg_region_bar_summary_final,
  isg_region_bar_summary_final_out,
  row.names = FALSE
)

message("Saved final all-FOV bar summary: ", isg_region_bar_summary_final_out)
print(isg_region_bar_summary_final, n = Inf)

# ------------------------------------------------------------
# 15.7 统计检验：所有 FOV 合并
# ------------------------------------------------------------

n_groups_available_final <- length(unique(as.character(
  isg_region_test_df_final$spatial_region_3class
)))

if (n_groups_available_final >= 2) {
  isg_region_kruskal_final <- kruskal.test(
    ISG_score_plot ~ spatial_region_3class,
    data = isg_region_test_df_final
  )

  isg_region_kruskal_final_df <- tibble(
    comparison_scope = "ALL_FOVS",
    fov_count = length(TARGET_FOVS_FINAL),
    fovs = paste(TARGET_FOVS_FINAL, collapse = "_"),
    n_cells = nrow(isg_region_test_df_final),
    n_groups = n_groups_available_final,
    kruskal_chisq = unname(isg_region_kruskal_final$statistic),
    kruskal_df = unname(isg_region_kruskal_final$parameter),
    kruskal_p = isg_region_kruskal_final$p.value
  )

  isg_region_kruskal_final_out <- file.path(
    TABLE_DIR,
    paste0(OUT_TAG_FINAL, "_myeloid_ISG_score_spatial_region_kruskal_ALL_FOVS.csv")
  )

  write.csv(
    isg_region_kruskal_final_df,
    isg_region_kruskal_final_out,
    row.names = FALSE
  )

  message("Saved final all-FOV Kruskal test: ", isg_region_kruskal_final_out)
  print(isg_region_kruskal_final_df, n = Inf)

  pair_list_final <- combn(
    x = c("Core-colocalized", "Peri-core adjacent", "Other region"),
    m = 2,
    simplify = FALSE
  )

  pair_list_final <- pair_list_final[vapply(pair_list_final, function(x) {
    all(x %in% unique(as.character(
      isg_region_test_df_final$spatial_region_3class
    )))
  }, logical(1))]

  isg_region_pairwise_final <- bind_rows(lapply(pair_list_final, function(pair_i) {
    g1 <- pair_i[1]
    g2 <- pair_i[2]

    x1 <- isg_region_test_df_final %>%
      filter(as.character(spatial_region_3class) == g1) %>%
      pull(ISG_score_plot)

    x2 <- isg_region_test_df_final %>%
      filter(as.character(spatial_region_3class) == g2) %>%
      pull(ISG_score_plot)

    wt <- wilcox.test(x1, x2)

    tibble(
      group1 = g1,
      group2 = g2,
      n1 = length(x1),
      n2 = length(x2),
      median1 = median(x1, na.rm = TRUE),
      median2 = median(x2, na.rm = TRUE),
      median_diff_group1_minus_group2 = median1 - median2,
      wilcox_p = wt$p.value
    )
  })) %>%
    mutate(
      wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH")
    )

  isg_region_pairwise_final_out <- file.path(
    TABLE_DIR,
    paste0(OUT_TAG_FINAL, "_myeloid_ISG_score_spatial_region_pairwise_wilcox_BH_ALL_FOVS.csv")
  )

  write.csv(
    isg_region_pairwise_final,
    isg_region_pairwise_final_out,
    row.names = FALSE
  )

  message("Saved final all-FOV pairwise Wilcoxon test: ", isg_region_pairwise_final_out)
  print(isg_region_pairwise_final, n = Inf)
} else {
  warning("Final all-FOV available spatial regions fewer than 2. Skip statistical tests.")
}

# ------------------------------------------------------------
# 15.8 画所有 FOV 合并后的 p_isg_region_bar_all
# ------------------------------------------------------------

y_max <- max(
  isg_region_bar_summary_final$mean_score +
    isg_region_bar_summary_final$se_score,
  na.rm = TRUE
)

y_min <- min(
  isg_region_bar_summary_final$mean_score -
    isg_region_bar_summary_final$se_score,
  na.rm = TRUE
)

y_range <- y_max - y_min

if (!is.finite(y_range) || y_range == 0) {
  y_range <- 1
}

y_lower <- min(0, y_min - 0.08 * y_range)
y_upper <- y_max + 0.22 * y_range

p_isg_region_bar_all <- ggplot(
  isg_region_bar_summary_final,
  aes(
    x = spatial_region_3class,
    y = mean_score,
    fill = spatial_region_3class
  )
) +
  geom_col(
    width = 0.68,
    color = "black",
    linewidth = 0.25,
    alpha = 0.92
  ) +
  geom_errorbar(
    aes(
      ymin = mean_score - se_score,
      ymax = mean_score + se_score
    ),
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = paste0("n=", n_cells),
      y = mean_score + se_score + 0.04 * y_range
    ),
    size = 3.2,
    vjust = 0
  ) +
  scale_fill_manual(
    values = c(
      "Core-colocalized" = "#B6424B",
      "Peri-core adjacent" = "#F58518",
      "Other region" = "#2ca02c"
    ),
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  coord_cartesian(
    ylim = c(y_lower, y_upper),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Mean APOBEC3A-like myeloid score",
    fill = "Spatial region",
    title = "APOBEC3A-like myeloid score across spatial regions"
  ) +
  theme_story(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(t = 16, r = 20, b = 16, l = 20)
  )

print(p_isg_region_bar_all)

save_plot_both(
  p_isg_region_bar_all,
  paste0(OUT_TAG_FINAL, "_myeloid_ISG_score_spatial_region_barplot_ALL_FOVS"),
  width = 5,
  height = 4
)

# ------------------------------------------------------------
# 15.9 对每个 FOV 单独画一张 p_spatial_region_myeloid
# ------------------------------------------------------------

per_fov_region_summary_final <- relation_df_final %>%
  filter(is_myeloid_transfer) %>%
  group_by(image, spatial_region_3class) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(ISG_score_plot, na.rm = TRUE),
    median_score = median(ISG_score_plot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(image, spatial_region_3class)

per_fov_region_summary_final_out <- file.path(
  TABLE_DIR,
  paste0(OUT_TAG_FINAL, "_myeloid_spatial_region_summary_PER_FOV.csv")
)

write.csv(
  per_fov_region_summary_final,
  per_fov_region_summary_final_out,
  row.names = FALSE
)

message("Saved final per-FOV myeloid region summary: ", per_fov_region_summary_final_out)

for (img in TARGET_FOVS_FINAL) {
  df_plot <- relation_df_final %>%
    filter(as.character(image) == as.character(img))

  core_poly_df_plot <- core_poly_df_final %>%
    filter(as.character(image) == as.character(img))

  sev_label <- unique(as.character(df_plot$severity_fov))
  sev_label <- sev_label[!is.na(sev_label) & sev_label != ""]
  if (length(sev_label) == 0) {
    sev_label <- "not annotated"
  } else {
    sev_label <- sev_label[1]
  }

  p_spatial_region_myeloid <- ggplot() +
    geom_point(
      data = df_plot %>% filter(!is_myeloid_transfer),
      aes(x = x, y = y),
      color = "grey90",
      size = 0.06,
      alpha = 0.30
    ) +
    geom_point(
      data = df_plot %>% filter(is_myeloid_transfer),
      aes(x = x, y = y, color = spatial_region_3class),
      size = 0.24,
      alpha = 0.92
    ) +
    geom_path(
      data = core_poly_df_plot,
      aes(x = X, y = Y, group = group_path),
      color = "black",
      linewidth = 0.48,
      linetype = "22",
      alpha = 0.95
    ) +
    scale_y_reverse() +
    coord_cartesian(clip = "off") +
    scale_color_manual(
      values = c(
        "Core-colocalized" = "#B6424B",
        "Peri-core adjacent" = "#F58518",
        "Other region" = "#8C8C8C"
      ),
      drop = FALSE,
      name = "Spatial region"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(img, ": spatial distribution of transferred myeloid cells"),
      subtitle = paste0(
        "Severity: ", sev_label,
        "; grey background: non-myeloid or unassigned cells"
      )
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(t = 14, r = 22, b = 14, l = 22)
    )

  print(p_spatial_region_myeloid)

  save_plot_both(
    p_spatial_region_myeloid,
    paste0(
      OUT_TAG_FINAL,
      "_",
      sanitize_filename(img),
      "_spatial_region_MYELOID"
    ),
    width = 5.2,
    height = 4.8
  )
}

message("========== FINAL all-FOV added outputs ==========")
message("FOV count used: ", length(TARGET_FOVS_FINAL))
message("Combined barplot: ", file.path(
  FIG_DIR,
  paste0(
    OUT_TAG_FINAL,
    "_myeloid_ISG_score_spatial_region_barplot_ALL_FOVS.pdf"
  )
))
message("Per-FOV myeloid spatial plots saved with suffix: _spatial_region_MYELOID")
message("========== Done ==========")

} else {
  message("Skipped final all-FOV add-on; set AST_RUN_FINAL_ALL_FOV_SUMMARY=1 to run it.")
}
