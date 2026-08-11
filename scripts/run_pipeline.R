#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

options(error = function() {
  traceback(2)
  quit(status = 1, save = "no")
})

option_list <- list(
  make_option(c("-c", "--config"), default = "configs/config.yaml"),
  make_option(
    c("-s", "--stage"), default = "all",
    help = paste(
      "Comma-separated stages:",
      "scrna,scoring,trajectory,pathway,spatial,bulk,rna-editing,supplementary,all"
    )
  ),
  make_option(c("--dry-run"), action = "store_true", default = FALSE, dest = "dry_run")
)
opt <- parse_args(OptionParser(option_list = option_list))
config_path <- normalizePath(opt$config, mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(config_path), ".."), mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)

resolve_path <- function(path) {
  if (grepl("^/", path)) path else file.path(repo_root, path)
}

run_command <- function(command, args = character(), env = character()) {
  display <- paste(c(env, shQuote(command), vapply(args, shQuote, character(1))), collapse = " ")
  message("[pipeline] ", display)
  if (opt$dry_run) return(invisible(0L))
  status <- system2(command, args = args, env = env)
  if (!identical(status, 0L)) stop("Command failed with status ", status, ": ", display)
  invisible(status)
}

run_r <- function(script, args = character(), env = character()) {
  run_command("Rscript", c(file.path(repo_root, script), args), env)
}

run_python <- function(script, args = character(), env = character()) {
  run_command("python3", c(file.path(repo_root, script), args), env)
}

common_env <- c(
  paste0("PIPELINE_REPO_ROOT=", repo_root),
  paste0("SCRNA_WORKSPACE_ROOT=", cfg$project$workspace_root),
  paste0("SCRNA_SCORED_RDS=", cfg$scrna$scored_rds),
  paste0("SCORING_OUTPUT_DIR=", resolve_path(cfg$scoring$output_dir)),
  paste0("SCORING_GEO_CACHE=", cfg$scrna$geo_cache),
  paste0("SCRNA_MONOCLE_DIR=", file.path(repo_root, "results", "trajectory")),
  paste0("SCRNA_PSEUDOTIME_RDS=", cfg$scrna$pseudotime_rds)
)

stage_scrna <- function() {
  scripts <- if (file.exists(cfg$inputs$existing_integrated_rds)) {
    message("[pipeline] Using configured annotated checkpoint; raw object preparation is skipped.")
    c(
      "workflow/01_scrna/02_integrate_cluster.R",
      "workflow/01_scrna/03_hdWGCNA_mo_ma.R"
    )
  } else {
    c(
      "workflow/01_scrna/01_prepare_seurat_objects.R",
      "workflow/01_scrna/02_integrate_cluster.R",
      "workflow/01_scrna/03_hdWGCNA_mo_ma.R"
    )
  }
  for (script in scripts) run_r(script, c("--config", config_path), common_env)
}

stage_scoring <- function() {
  generated_scored <- resolve_path(cfg$outputs$hdwgcnna_rds)
  scoring_input <- if (file.exists(generated_scored)) generated_scored else cfg$scrna$scored_rds
  scoring_env <- c(
    common_env[!grepl("^SCRNA_SCORED_RDS=", common_env)],
    paste0("SCRNA_SCORED_RDS=", scoring_input)
  )
  scripts <- c(
    "workflow/02_scoring/01_build_external_rankings_and_sensitivity.R",
    "workflow/02_scoring/02_freeze_top200_consensus.R",
    "workflow/02_scoring/03_score_and_plot_mpi_mmi.R"
  )
  for (script in scripts) run_r(script, env = scoring_env)
}

stage_trajectory <- function() {
  dir.create(file.path(repo_root, "results", "trajectory"), recursive = TRUE, showWarnings = FALSE)
  scripts <- c(
    "workflow/03_trajectory/01_plot_global_pseudotime.R",
    "workflow/03_trajectory/02_plot_selected_fate_branches.R",
    "workflow/03_trajectory/03_plot_region_pseudotime.R",
    "workflow/03_trajectory/04_plot_expression_umap.R"
  )
  for (script in scripts) run_r(script, env = common_env)
}

stage_pathway <- function() {
  run_r(
    "workflow/01_scrna/04_pathway_analysis.R",
    c("--config", config_path),
    common_env
  )
}

stage_spatial <- function() {
  spatial_env <- c(
    common_env,
    paste0("AST_ROOT=", file.path(repo_root, "workflow", "04_spatial")),
    paste0("AST_OUTPUT_ROOT=", file.path(repo_root, "results", "spatial")),
    paste0("AST_SPATIAL_ROOT=", cfg$spatial$spatial_root),
    paste0("AST_SC_ROOT=", cfg$project$workspace_root)
  )
  scripts <- c(
    "workflow/04_spatial/00_inventory.R",
    "workflow/04_spatial/01_visium_publication_figures.R",
    "workflow/04_spatial/02_xenium_publication_figures.R",
    "workflow/04_spatial/03_cross_platform_summary.R",
    "workflow/04_spatial/04_raw_archive_manifest.py",
    "workflow/04_spatial/08_geomx_wta_apobec3a.py",
    "workflow/04_spatial/05_raw_archive_overview.R",
    "workflow/04_spatial/06_apobec3a_coexpression.R",
    "workflow/04_spatial/07_apobec3a_all_sample_inventory.R",
    "workflow/04_spatial/09_myeloid_apobec3a_integrated_story.R",
    "workflow/04_spatial/10_apobec3a_like_myeloid_core_region.R"
  )
  for (script in scripts) {
    if (grepl("\\.py$", script)) run_python(script, env = spatial_env) else run_r(script, env = spatial_env)
  }
}

stage_bulk <- function() {
  run_r(
    "workflow/05_bulk_rnaseq/analyze_bulk_rnaseq.R",
    c("--config", config_path, "--contrast", "all"),
    common_env
  )
}

stage_rna_editing <- function() {
  result_root <- file.path(repo_root, "results", "rna_editing")
  expression_root <- file.path(repo_root, "results", "bulk_rnaseq", "counts")
  contrasts <- c("combined", "clone13", "clone37")
  for (contrast in contrasts) {
    jacusa_dir <- file.path(result_root, paste0("jacusa2_", contrast))
    analysis_dir <- file.path(result_root, contrast)
    expected <- if (contrast == "combined") 6L else 3L
    editing_env <- c(
      common_env,
      paste0("RNA_EDITING_CONTRAST=", contrast),
      paste0("RNA_EDITING_ALIGNMENT_DIR=", cfg$rna_editing$alignment_dir),
      paste0("RNA_EDITING_JACUSA_JAR=", cfg$rna_editing$jacusa_jar),
      paste0("RNA_EDITING_REFERENCE_FASTA=", cfg$rna_editing$reference_fasta),
      paste0("RNA_EDITING_ANNOTATION_GTF=", cfg$rna_editing$annotation_gtf),
      paste0("RNA_EDITING_ANNOTATION_GFF3=", cfg$rna_editing$annotation_gff3),
      paste0("RNA_EDITING_KEGG_RDS=", cfg$rna_editing$kegg_rds),
      paste0("RNA_EDITING_JACUSA_OUTPUT_DIR=", jacusa_dir),
      paste0("RNA_EDITING_EXPRESSION_OUTPUT_DIR=", expression_root),
      paste0("RNA_EDITING_THREADS=", cfg$rna_editing$threads),
      paste0("RNA_EDITING_MIN_COVERAGE=", cfg$rna_editing$min_coverage),
      paste0("RNA_EDITING_ANALYSIS_OUTPUT_DIR=", analysis_dir),
      paste0("RNA_EDITING_SAMPLE_MANIFEST=", file.path(jacusa_dir, "jacusa_sample_manifest.tsv")),
      paste0("RNA_EDITING_EXPECTED_REPS=", expected),
      paste0("RNA_EDITING_WT_MIN_REPS=", if (contrast == "combined") cfg$rna_editing$wt_min_covered_replicates else 2L),
      paste0("RNA_EDITING_KO_MIN_REPS=", if (contrast == "combined") cfg$rna_editing$ko_min_covered_replicates else 2L),
      paste0("RNA_EDITING_PVALUE_CUTOFF=", cfg$rna_editing$pvalue_cutoff),
      paste0("RNA_EDITING_DELTA_CUTOFF=", cfg$rna_editing$delta_editing_cutoff),
      "RNA_EDITING_SCORE_VALUES=1",
      paste0("RNA_EDITING_RNASEQ_RESULT_ROOT=", resolve_path(cfg$bulk_rnaseq$output_dir))
    )
    run_python("workflow/06_rna_editing/01_call_jacusa2_and_featurecounts.py", env = editing_env)
    run_r("workflow/06_rna_editing/02_analyze_editing.R", env = editing_env)
  }
  run_r(
    "workflow/06_rna_editing/03_compare_contrasts.R",
    env = c(common_env, paste0("RNA_EDITING_RESULT_ROOT=", result_root))
  )
}

stage_supplementary <- function() {
  env <- c(
    common_env,
    paste0("SUPPLEMENTARY_SOURCE_ROOT=", cfg$supplementary_tables$source_workspace),
    paste0("PUBLIC7_NGS_RESULT_ROOT=", cfg$supplementary_tables$public7_ngs_result_root),
    paste0("SUPPLEMENTARY_OUTPUT_DIR=", resolve_path(cfg$supplementary_tables$output_dir))
  )
  run_python("workflow/07_supplementary_tables/create_tables.py", env = env)
}

stage_functions <- list(
  scrna = stage_scrna,
  scoring = stage_scoring,
  trajectory = stage_trajectory,
  pathway = stage_pathway,
  spatial = stage_spatial,
  bulk = stage_bulk,
  `rna-editing` = stage_rna_editing,
  supplementary = stage_supplementary
)
requested <- trimws(strsplit(opt$stage, ",", fixed = TRUE)[[1]])
if (identical(requested, "all")) requested <- names(stage_functions)
unknown <- setdiff(requested, names(stage_functions))
if (length(unknown)) stop("Unknown stage(s): ", paste(unknown, collapse = ", "))

for (stage in requested) {
  message("\n[pipeline] ===== ", stage, " =====")
  stage_functions[[stage]]()
}
message("[pipeline] Completed requested stages: ", paste(requested, collapse = ", "))
