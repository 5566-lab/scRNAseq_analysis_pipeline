#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

SPATIAL_ROOT <- Sys.getenv(
  "AST_SPATIAL_ROOT",
  unset = "/public3/DSC/single_cell/spatial"
)
ROOT <- Sys.getenv(
  "AST_ROOT",
  unset = normalizePath("workflow/04_spatial", mustWork = TRUE)
)
OUTPUT_ROOT <- Sys.getenv("AST_OUTPUT_ROOT", unset = file.path(ROOT, "results"))
TABLE_DIR <- file.path(OUTPUT_ROOT, "tables")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

target_gene <- "APOBEC3A"

detect_from_matrix <- function(mat, sample_cells = NULL) {
  if (!target_gene %in% rownames(mat)) {
    return(list(
      assay_has_APOBEC3A = FALSE,
      n_units = ncol(mat),
      n_APOBEC3A_positive = NA_integer_,
      APOBEC3A_total_counts = NA_real_,
      APOBEC3A_detection_rate = NA_real_
    ))
  }
  if (!is.null(sample_cells)) {
    sample_cells <- intersect(sample_cells, colnames(mat))
    mat <- mat[, sample_cells, drop = FALSE]
  }
  x <- mat[target_gene, , drop = TRUE]
  list(
    assay_has_APOBEC3A = TRUE,
    n_units = length(x),
    n_APOBEC3A_positive = sum(x > 0),
    APOBEC3A_total_counts = sum(x),
    APOBEC3A_detection_rate = mean(x > 0)
  )
}

as_row <- function(dataset, platform, sample_id, disease_or_source, analysis_unit,
                   status, details, reason = "") {
  c(
    dataset = dataset,
    platform = platform,
    sample_id = sample_id,
    disease_or_source = disease_or_source,
    analysis_unit = analysis_unit,
    status = status,
    details,
    reason = reason
  )
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

rows <- list()

message("Checking GSE315246 Xenium")
xen <- readRDS(file.path(SPATIAL_ROOT, "GSE315246_xenium.obj.integrated.rds"))
xen_counts <- GetAssayData(xen, assay = "Xenium", layer = "counts")
for (image_name in Images(xen)) {
  cells <- Cells(xen[[image_name]])
  det <- detect_from_matrix(xen_counts, cells)
  disease <- names(sort(table(xen$disease[cells]), decreasing = TRUE))[1]
  myeloid_cells <- intersect(cells, rownames(xen@meta.data)[xen$predicted.id == "Myeloid"])
  myeloid_det <- detect_from_matrix(xen_counts, myeloid_cells)
  detail <- c(
    assay_has_APOBEC3A = det$assay_has_APOBEC3A,
    n_units = det$n_units,
    n_APOBEC3A_positive = det$n_APOBEC3A_positive,
    APOBEC3A_total_counts = det$APOBEC3A_total_counts,
    APOBEC3A_detection_rate = det$APOBEC3A_detection_rate,
    n_myeloid_units = myeloid_det$n_units,
    n_myeloid_APOBEC3A_positive = myeloid_det$n_APOBEC3A_positive,
    myeloid_APOBEC3A_detection_rate = myeloid_det$APOBEC3A_detection_rate
  )
  rows[[length(rows) + 1]] <- as_row(
    "GSE315246", "Xenium", image_name, disease, "cell/FOV",
    ifelse(det$n_APOBEC3A_positive > 0, "APOBEC3A_detected", "APOBEC3A_not_detected"),
    detail
  )
}

message("Checking GSE314851 Visium FFPE")
vis <- readRDS(file.path(SPATIAL_ROOT, "GSE314851_Visium_FFPE_integrated.rds"))
vis_counts <- GetAssayData(vis, assay = "Spatial", layer = "counts")
for (sample_id in sort(unique(vis$sample))) {
  cells <- rownames(vis@meta.data)[vis$sample == sample_id]
  det <- detect_from_matrix(vis_counts, cells)
  severity <- names(sort(table(vis$category[cells]), decreasing = TRUE))[1]
  detail <- c(
    assay_has_APOBEC3A = det$assay_has_APOBEC3A,
    n_units = det$n_units,
    n_APOBEC3A_positive = det$n_APOBEC3A_positive,
    APOBEC3A_total_counts = det$APOBEC3A_total_counts,
    APOBEC3A_detection_rate = det$APOBEC3A_detection_rate,
    n_myeloid_units = NA_integer_,
    n_myeloid_APOBEC3A_positive = NA_integer_,
    myeloid_APOBEC3A_detection_rate = NA_real_
  )
  rows[[length(rows) + 1]] <- as_row(
    "GSE314851", "10x Visium FFPE", sample_id, severity, "spot/sample",
    ifelse(det$n_APOBEC3A_positive > 0, "APOBEC3A_detected", "APOBEC3A_not_detected"),
    detail
  )
}

message("Checking GSE243179 erosion-plaque Visium")
gse243 <- readRDS(file.path(SPATIAL_ROOT, "GSE243179", "merge.rds"))
spatial_layers <- grep("^counts", Layers(gse243[["Spatial"]]), value = TRUE)
sample_ids <- Images(gse243)
for (i in seq_along(spatial_layers)) {
  sample_id <- sample_ids[i]
  mat <- LayerData(gse243, assay = "Spatial", layer = spatial_layers[i])
  det <- detect_from_matrix(mat)
  detail <- c(
    assay_has_APOBEC3A = det$assay_has_APOBEC3A,
    n_units = det$n_units,
    n_APOBEC3A_positive = det$n_APOBEC3A_positive,
    APOBEC3A_total_counts = det$APOBEC3A_total_counts,
    APOBEC3A_detection_rate = det$APOBEC3A_detection_rate,
    n_myeloid_units = NA_integer_,
    n_myeloid_APOBEC3A_positive = NA_integer_,
    myeloid_APOBEC3A_detection_rate = NA_real_
  )
  rows[[length(rows) + 1]] <- as_row(
    "GSE243179", "10x Visium FFPE", sample_id, "erosion plaque", "spot/sample",
    ifelse(det$n_APOBEC3A_positive > 0, "APOBEC3A_detected", "APOBEC3A_not_detected"),
    detail,
    "APOBEC3A is nearly absent in this erosion-plaque Visium object"
  )
}

read_cosmx_expr_apobec3a <- function(tarfile_path, member_name) {
  tmp <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(tmp), add = TRUE)
  untar(tarfile_path, files = member_name, exdir = dirname(tmp))
  extracted <- file.path(dirname(tmp), member_name)
  on.exit(unlink(extracted), add = TRUE)
  header <- readLines(gzfile(extracted), n = 1)
  fields <- strsplit(header, ",", fixed = TRUE)[[1]]
  keep_cols <- fields %in% c("fov", "cell_ID", target_gene)
  has_target <- target_gene %in% fields
  col_classes <- ifelse(keep_cols, NA, "NULL")
  dat <- read.csv(gzfile(extracted), colClasses = col_classes, check.names = FALSE)
  if (!has_target) {
    return(list(assay_has_APOBEC3A = FALSE, n_units = nrow(dat), n_fovs = length(unique(dat$fov)),
                n_APOBEC3A_positive = NA_integer_, APOBEC3A_total_counts = NA_real_,
                APOBEC3A_detection_rate = NA_real_))
  }
  x <- dat[[target_gene]]
  list(
    assay_has_APOBEC3A = TRUE,
    n_units = nrow(dat),
    n_fovs = length(unique(dat$fov)),
    n_APOBEC3A_positive = sum(x > 0),
    APOBEC3A_total_counts = sum(x),
    APOBEC3A_detection_rate = mean(x > 0)
  )
}

message("Checking GSE277441 raw CosMx")
cosmx_tar <- file.path(SPATIAL_ROOT, "GSE277441_RAW.tar")
cosmx_members <- untar(cosmx_tar, list = TRUE)
cosmx_expr <- grep("exprMat_file[.]csv[.]gz$", cosmx_members, value = TRUE)
cosmx_meta <- read.csv(file.path(TABLE_DIR, "raw_gse_sample_metadata.csv"), check.names = FALSE)
for (member in cosmx_expr) {
  gsm <- sub("^(GSM[0-9]+).*", "\\1", basename(member))
  det <- read_cosmx_expr_apobec3a(cosmx_tar, member)
  severity <- cosmx_meta$lesion_severity[match(gsm, cosmx_meta$gsm)]
  detail <- c(
    assay_has_APOBEC3A = det$assay_has_APOBEC3A,
    n_units = det$n_units,
    n_APOBEC3A_positive = det$n_APOBEC3A_positive,
    APOBEC3A_total_counts = det$APOBEC3A_total_counts,
    APOBEC3A_detection_rate = det$APOBEC3A_detection_rate,
    n_myeloid_units = NA_integer_,
    n_myeloid_APOBEC3A_positive = NA_integer_,
    myeloid_APOBEC3A_detection_rate = NA_real_
  )
  rows[[length(rows) + 1]] <- as_row(
    "GSE277441", "NanoString CosMx SMI", gsm, severity, "cell/raw sample",
    ifelse(!det$assay_has_APOBEC3A, "APOBEC3A_not_in_panel",
      ifelse(det$n_APOBEC3A_positive > 0, "APOBEC3A_detected", "APOBEC3A_not_detected")
    ),
    detail,
    paste0("n_fovs=", det$n_fovs)
  )
}

parse_mtx_target <- function(tarfile_path, features_member, matrix_member) {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  untar(tarfile_path, files = c(features_member, matrix_member), exdir = tmpdir)
  features_path <- file.path(tmpdir, features_member)
  matrix_path <- file.path(tmpdir, matrix_member)
  features <- read.delim(gzfile(features_path), header = FALSE, stringsAsFactors = FALSE)
  gene_idx <- which(features[[2]] == target_gene)
  dims <- readLines(gzfile(matrix_path), n = 3)
  while (startsWith(dims[length(dims)], "%")) {
    dims <- c(dims, readLines(gzfile(matrix_path), n = 1))
  }
  # Reopen after simple header read; parse line by line for the target row.
  con <- gzfile(matrix_path, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0 || !startsWith(line, "%")) break
  }
  dim_vals <- as.integer(strsplit(line, " ", fixed = TRUE)[[1]])
  n_spots <- dim_vals[2]
  if (length(gene_idx) == 0) {
    return(list(assay_has_APOBEC3A = FALSE, n_units = n_spots, n_APOBEC3A_positive = NA_integer_,
                APOBEC3A_total_counts = NA_real_, APOBEC3A_detection_rate = NA_real_))
  }
  positive_spots <- integer()
  total_counts <- 0
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    vals <- strsplit(line, " ", fixed = TRUE)[[1]]
    if (as.integer(vals[1]) == gene_idx) {
      count <- as.numeric(vals[3])
      if (count > 0) positive_spots <- c(positive_spots, as.integer(vals[2]))
      total_counts <- total_counts + count
    }
  }
  list(
    assay_has_APOBEC3A = TRUE,
    n_units = n_spots,
    n_APOBEC3A_positive = length(unique(positive_spots)),
    APOBEC3A_total_counts = total_counts,
    APOBEC3A_detection_rate = length(unique(positive_spots)) / n_spots
  )
}

message("Checking GSE283269 raw DIT Visium")
dit_tar <- file.path(SPATIAL_ROOT, "GSE283269_RAW.tar")
dit_members <- untar(dit_tar, list = TRUE)
dit_meta <- cosmx_meta
dit_features <- grep("features[.]tsv[.]gz$", dit_members, value = TRUE)
dit_features <- dit_features[!grepl("raw_features", dit_features)]
for (features_member in dit_features) {
  gsm <- sub("^(GSM[0-9]+).*", "\\1", basename(features_member))
  prefix <- sub("_features[.]tsv[.]gz$", "", features_member)
  matrix_member <- paste0(prefix, "_matrix.mtx.gz")
  if (!matrix_member %in% dit_members) next
  det <- parse_mtx_target(dit_tar, features_member, matrix_member)
  source <- dit_meta$tissue_source[match(gsm, dit_meta$gsm)]
  title <- dit_meta$title[match(gsm, dit_meta$gsm)]
  qc_note <- ifelse(grepl("Sample2d1", title, ignore.case = TRUE),
                    "source study removed Sample2d1 after QC", "")
  detail <- c(
    assay_has_APOBEC3A = det$assay_has_APOBEC3A,
    n_units = det$n_units,
    n_APOBEC3A_positive = det$n_APOBEC3A_positive,
    APOBEC3A_total_counts = det$APOBEC3A_total_counts,
    APOBEC3A_detection_rate = det$APOBEC3A_detection_rate,
    n_myeloid_units = NA_integer_,
    n_myeloid_APOBEC3A_positive = NA_integer_,
    myeloid_APOBEC3A_detection_rate = NA_real_
  )
  rows[[length(rows) + 1]] <- as_row(
    "GSE283269", "10x Visium FFPE", gsm, source, "spot/raw sample",
    ifelse(det$n_APOBEC3A_positive > 0, "APOBEC3A_detected", "APOBEC3A_not_detected"),
    detail,
    qc_note
  )
}

message("Checking GSE277170 GeoMx WTA DCC resolved through Bruker PKC")
geomx_apobec_path <- file.path(TABLE_DIR, "geomx_wta_apobec_roi_counts.csv")
if (!file.exists(geomx_apobec_path)) {
  system2("python3", file.path(ROOT, "scripts", "08_geomx_wta_apobec3a.py"))
}
geomx_apobec <- read.csv(geomx_apobec_path, check.names = FALSE)
for (i in seq_len(nrow(geomx_apobec))) {
  is_ntc <- as_bool(geomx_apobec$is_no_template_control[i])
  is_negative_control <- as_bool(geomx_apobec$is_negative_control_roi[i])
  apobec_count <- as.numeric(geomx_apobec$APOBEC3A_count[i])
  status <- if (is_ntc) {
    "control_no_template"
  } else if (apobec_count > 0) {
    "APOBEC3A_detected"
  } else {
    "APOBEC3A_not_detected"
  }
  reason <- paste0(
    "Bruker Hs_R_NGS_WTA_v1.0.pkc resolved DCC RTS_ID counts",
    "; localisation=", geomx_apobec$localisation[i],
    "; subset=", geomx_apobec$subset[i],
    "; APOBEC3A_B_count=", geomx_apobec$APOBEC3A_B_count[i],
    "; APOBEC3B_count=", geomx_apobec$APOBEC3B_count[i]
  )
  if (is_ntc) {
    reason <- paste(reason, "No Template Control", sep = "; ")
  } else if (is_negative_control) {
    reason <- paste(reason, "ROI annotated Negative control", sep = "; ")
  }
  detail <- c(
    assay_has_APOBEC3A = TRUE,
    n_units = 1L,
    n_APOBEC3A_positive = ifelse(!is_ntc && apobec_count > 0, 1L, 0L),
    APOBEC3A_total_counts = apobec_count,
    APOBEC3A_detection_rate = ifelse(!is_ntc && apobec_count > 0, 1, 0),
    n_myeloid_units = NA_integer_,
    n_myeloid_APOBEC3A_positive = NA_integer_,
    myeloid_APOBEC3A_detection_rate = NA_real_
  )
  rows[[length(rows) + 1]] <- as_row(
    "GSE277170", "NanoString GeoMx DSP WTA", geomx_apobec$gsm[i],
    ifelse(is.na(geomx_apobec$grade[i]) || geomx_apobec$grade[i] == "NA", "not annotated", geomx_apobec$grade[i]),
    "ROI/DCC",
    status,
    detail,
    reason
  )
}

out <- bind_rows(lapply(rows, as.data.frame.list, stringsAsFactors = FALSE)) %>%
  mutate(
    n_units = as.integer(n_units),
    n_APOBEC3A_positive = as.integer(n_APOBEC3A_positive),
    APOBEC3A_total_counts = as.numeric(APOBEC3A_total_counts),
    APOBEC3A_detection_rate = as.numeric(APOBEC3A_detection_rate),
    n_myeloid_units = as.integer(n_myeloid_units),
    n_myeloid_APOBEC3A_positive = as.integer(n_myeloid_APOBEC3A_positive),
    myeloid_APOBEC3A_detection_rate = as.numeric(myeloid_APOBEC3A_detection_rate)
  )

write.csv(out, file.path(TABLE_DIR, "apobec3a_all_sample_detection_inventory.csv"), row.names = FALSE)

summary <- out %>%
  group_by(dataset, platform, status) %>%
  summarise(
    n_samples = n(),
    total_units = sum(n_units, na.rm = TRUE),
    APOBEC3A_positive_units = sum(n_APOBEC3A_positive, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary, file.path(TABLE_DIR, "apobec3a_all_dataset_detection_summary.csv"), row.names = FALSE)

message("Wrote APOBEC3A all-sample inventory")
print(summary)
