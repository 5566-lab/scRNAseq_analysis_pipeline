#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

option_list <- list(
  make_option(c("-c", "--config"), default = "configs/config.yaml"),
  make_option(c("--check-inputs"), action = "store_true", default = FALSE, dest = "check_inputs")
)
opt <- parse_args(OptionParser(option_list = option_list))
config_path <- normalizePath(opt$config, mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(config_path), ".."), mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)
failures <- character()

check <- function(ok, message) {
  if (isTRUE(ok)) {
    cat("PASS  ", message, "\n", sep = "")
  } else {
    cat("FAIL  ", message, "\n", sep = "")
    failures <<- c(failures, message)
  }
}

repo_file <- function(...) file.path(repo_root, ...)

r_files <- list.files(repo_root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
r_parse_ok <- vapply(r_files, function(path) {
  tryCatch({
    parse(path)
    TRUE
  }, error = function(e) {
    message("R parse error in ", path, ": ", conditionMessage(e))
    FALSE
  })
}, logical(1))
check(all(r_parse_ok), paste("R syntax:", length(r_files), "files"))

python_files <- list.files(repo_root, pattern = "\\.py$", recursive = TRUE, full.names = TRUE)
python_check <- paste(
  "import ast, pathlib, sys;",
  "p=pathlib.Path(sys.argv[1]);",
  "ast.parse(p.read_text(encoding='utf-8'), filename=str(p))"
)
python_ok <- vapply(python_files, function(path) {
  identical(system2("python3", c("-c", shQuote(python_check), shQuote(path))), 0L)
}, logical(1))
check(all(python_ok), paste("Python syntax:", length(python_files), "files"))

membership <- read.csv(
  repo_file("data", "gene_sets", "MMI_GSE5099_GSE11864_Top200_membership.csv"),
  stringsAsFactors = FALSE
)
mature <- membership$gene_symbol[membership$direction == "mature_positive"]
immature <- membership$gene_symbol[membership$direction == "immature_negative"]
check(length(mature) == 200L, "MMI mature-positive set contains 200 genes")
check(length(immature) == 200L, "MMI monocyte-high set contains 200 genes")
check(!anyDuplicated(mature) && !anyDuplicated(immature), "MMI sets contain no within-set duplicates")
check(!length(intersect(mature, immature)), "MMI directions do not overlap")
check(identical(unique(membership$top_n_final_per_direction), 200L), "MMI membership is frozen at Top200")

mmi_table <- read.csv(repo_file("data", "gene_sets", "mmi_top200_signatures.csv"), stringsAsFactors = FALSE)
frozen_mature <- mmi_table$gene[mmi_table$gene_set == "MMI macrophage-maturation signature"]
frozen_immature <- mmi_table$gene[mmi_table$gene_set == "MMI monocyte-associated signature"]
check(setequal(mature, frozen_mature), "S6 MMI mature signature matches ranked membership")
check(setequal(immature, frozen_immature), "S6 MMI monocyte signature matches ranked membership")

mpi_table <- read.csv(repo_file("data", "gene_sets", "mpi_signatures.csv"), stringsAsFactors = FALSE)
mpi_counts <- table(mpi_table$gene_set)
check(unname(mpi_counts["MPI M1-like signature"]) == 145L, "MPI M1-like signature contains 145 genes")
check(unname(mpi_counts["MPI M2-like signature"]) == 165L, "MPI M2-like signature contains 165 genes")

manifest <- read.delim(
  repo_file("metadata", "rna_editing_sample_manifest.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
check(sum(manifest$Condition == "WT") == 6L, "RNA-editing manifest contains 6 WT samples")
check(sum(manifest$Condition == "KO") == 6L, "RNA-editing manifest contains 6 KO samples")
check(all(table(manifest$Batch, manifest$Condition) == 3L), "Each clone contains 3 WT and 3 KO samples")
check(!any(grepl("THP_1", manifest$Sample_Name)), "Primary manifest excludes THP-1")
check(all(manifest$JACUSA_Condition[manifest$Condition == "WT"] == "cond1"), "JACUSA cond1 is WT")
check(all(manifest$JACUSA_Condition[manifest$Condition == "KO"] == "cond2"), "JACUSA cond2 is KO")

trajectory_text <- paste(vapply(
  list.files(repo_file("workflow", "03_trajectory"), pattern = "\\.R$", full.names = TRUE),
  function(path) paste0(readLines(path, warn = FALSE), collapse = "\n"),
  character(1)
), collapse = "\n")
check(grepl("data_pseudotime\\.rds", trajectory_text), "Trajectory figures use the frozen publication checkpoint")
check(!grepl("readRDS\\([^\\n]*cds_MM_foam\\.rds", trajectory_text), "Trajectory figures do not read the drifted CDS checkpoint")

bulk_text <- paste(readLines(repo_file("workflow", "05_bulk_rnaseq", "analyze_bulk_rnaseq.R"), warn = FALSE), collapse = "\n")
check(grepl("~ clone \\+ condition", bulk_text), "Combined bulk model includes clone")
check(!grepl("ComBat_seq", bulk_text), "DESeq2 does not consume ComBat-seq-adjusted counts")

editing_text <- paste(readLines(repo_file("workflow", "06_rna_editing", "02_analyze_editing.R"), warn = FALSE), collapse = "\n")
check(grepl("final_fdr", editing_text), "RNA-editing tables export BH-FDR")
check(grepl("mean_ko - mean_wt", editing_text), "RNA-editing delta is KO minus WT")

required_keys <- c("scrna", "scoring", "spatial", "bulk_rnaseq", "rna_editing", "supplementary_tables")
check(all(required_keys %in% names(cfg)), "Configuration contains all publication modules")

if (opt$check_inputs) {
  input_paths <- c(
    cfg$scrna$annotated_rds,
    cfg$scrna$scored_rds,
    cfg$scrna$pseudotime_rds,
    cfg$spatial$xenium_rds,
    cfg$spatial$visium_rds,
    cfg$bulk_rnaseq$count_table,
    cfg$rna_editing$jacusa_jar,
    cfg$rna_editing$reference_fasta,
    cfg$rna_editing$annotation_gtf,
    cfg$rna_editing$annotation_gff3,
    cfg$rna_editing$kegg_rds,
    cfg$supplementary_tables$public7_ngs_result_root
  )
  path_ok <- file.exists(input_paths)
  if (any(!path_ok)) cat("Missing configured inputs:\n", paste(input_paths[!path_ok], collapse = "\n"), "\n")
  check(all(path_ok), "Configured primary inputs exist")
}

if (length(failures)) {
  stop(length(failures), " validation check(s) failed:\n- ", paste(failures, collapse = "\n- "), call. = FALSE)
}
cat("\nRepository validation passed.\n")
