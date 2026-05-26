#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(GEOquery)
})

source("R/utils/config.R")
source("R/utils/seurat_io.R")

cfg <- load_config()
set.seed(cfg$project$seed)

load_gse260657 <- function(cfg) {
  data_dir <- cfg$inputs$gse260657_dir
  files <- list.files(data_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  objects <- list()

  for (file_path in files) {
    patient_num <- as.numeric(sub(".*human_([0-9]+)\\.txt\\.gz", "\\1", basename(file_path)))
    project_name <- paste0("GSE260657_Human_", patient_num)
    plaque <- if (!is.na(patient_num) && patient_num <= 7) {
      "Carotid Stable Plaque"
    } else {
      "Carotid Unstable Plaque"
    }

    message_step("Loading ", project_name)
    counts <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    obj <- create_seurat_from_counts(counts, project_name, cfg)
    obj$Patient_ID <- paste0("Human_", patient_num)
    obj$AC_PA <- plaque
    obj$Source_GSE <- "GSE260657"
    objects[[project_name]] <- obj
  }

  if (length(objects) == 0) {
    return(NULL)
  }
  merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), project = "GSE260657_Combined")
}

load_gse247238 <- function(cfg) {
  patient_db <- list(
    "1" = list(type = "Carotid Stable Plaque", age = 70, sex = "Male"),
    "2" = list(type = "Carotid Stable Plaque", age = 72, sex = "Female"),
    "3" = list(type = "Carotid Unstable Plaque", age = 74, sex = "Male"),
    "4" = list(type = "Carotid Unstable Plaque", age = 81, sex = "Male"),
    "5" = list(type = "Carotid Stable Plaque", age = 71, sex = "Male"),
    "6" = list(type = "Carotid Unstable Plaque", age = 80, sex = "Male"),
    "7" = list(type = "Carotid Stable Plaque", age = 78, sex = "Male"),
    "8" = list(type = "Carotid Stable Plaque", age = 69, sex = "Male"),
    "9" = list(type = "Carotid Unstable Plaque", age = 82, sex = "Male"),
    "10" = list(type = "Carotid Stable Plaque", age = 71, sex = "Male")
  )

  gse_info <- tryCatch(
    GEOquery::getGEO("GSE247238", GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) NULL
  )
  gsm_to_title <- character()
  if (!is.null(gse_info)) {
    gsm_to_title <- Biobase::pData(gse_info[[1]])$title
    names(gsm_to_title) <- Biobase::pData(gse_info[[1]])$geo_accession
  }

  data_dir <- cfg$inputs$gse247238_dir
  matrix_files <- list.files(data_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = FALSE)
  sample_ids <- str_remove(matrix_files, "_matrix\\.mtx\\.gz$")
  objects <- list()

  for (sample in sample_ids) {
    gsm_id <- strsplit(sample, "_")[[1]][1]
    message_step("Loading ", sample)
    counts <- ReadMtx(
      mtx = file.path(data_dir, paste0(sample, "_matrix.mtx.gz")),
      cells = file.path(data_dir, paste0(sample, "_barcodes.tsv.gz")),
      features = file.path(data_dir, paste0(sample, "_features.tsv.gz")),
      feature.column = 1
    )
    obj <- create_seurat_from_counts(counts, sample, cfg)
    car_id <- "Unknown"
    plaque <- "Unknown"
    age <- NA
    sex <- "Unknown"

    if (gsm_id %in% names(gsm_to_title)) {
      car_match <- str_extract(gsm_to_title[[gsm_id]], "CAR[0-9]+")
      if (!is.na(car_match)) {
        car_id <- str_remove(car_match, "CAR")
        if (car_id %in% names(patient_db)) {
          plaque <- patient_db[[car_id]]$type
          age <- patient_db[[car_id]]$age
          sex <- patient_db[[car_id]]$sex
        }
      }
    }

    obj$GSM_ID <- gsm_id
    obj$Source_GSE <- "GSE247238"
    obj$Patient_ID <- paste0("Patient_", car_id)
    obj$AC_PA <- plaque
    obj$Age <- age
    obj$Sex <- sex
    objects[[sample]] <- obj
  }

  if (length(objects) == 0) {
    return(NULL)
  }
  merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), project = "GSE247238_Combined")
}

load_existing_rds <- function(path, source_gse) {
  if (!file.exists(path)) {
    return(NULL)
  }
  obj <- readRDS(path)
  obj$Source_GSE <- source_gse
  obj
}

objects <- list()
objects$GSE260657 <- load_gse260657(cfg)
objects$GSE247238 <- load_gse247238(cfg)
objects$GSE210152 <- load_existing_rds(cfg$inputs$gse210152_rds, "GSE210152")
objects <- objects[!vapply(objects, is.null, logical(1))]

out <- project_path(cfg, cfg$outputs$prepared_objects_rds)
safe_save_rds(objects, out)
message_step("Saved prepared object list: ", out)

