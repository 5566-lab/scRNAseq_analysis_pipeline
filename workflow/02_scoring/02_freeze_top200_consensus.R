#!/usr/bin/env Rscript

# Build an exact 200 + 200 external MMI gene set from GSE5099 and GSE11864.
# Selection is source-balanced, deterministic, and independent of carotid
# expression values and subtype labels.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

repo_root <- Sys.getenv(
  "PIPELINE_REPO_ROOT",
  unset = normalizePath(".", mustWork = TRUE)
)
ranking_dir <- Sys.getenv(
  "SCORING_OUTPUT_DIR",
  unset = file.path(repo_root, "results", "scoring")
)
gene_set_dir <- file.path(repo_root, "data", "gene_sets")
dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 200L
gse5099_file <- file.path(
  ranking_dir, "GSE5099_complete_external_effect_ranking.csv"
)
gse11864_file <- file.path(
  ranking_dir, "GSE11864_complete_external_effect_ranking.csv"
)
membership_file <- file.path(
  gene_set_dir, "MMI_GSE5099_GSE11864_Top200_membership.csv"
)
conflict_file <- file.path(
  gene_set_dir, "MMI_GSE5099_GSE11864_Top200_direction_conflicts.csv"
)
summary_file <- file.path(
  gene_set_dir, "MMI_GSE5099_GSE11864_Top200_selection_summary.csv"
)

required_columns <- c(
  "gene_symbol", "logFC", "moderated_t",
  "rank_macrophage", "rank_monocyte"
)
read_ranking <- function(path, source_name) {
  ranking <- read.csv(path, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, colnames(ranking))
  if (length(missing_columns)) {
    stop(
      source_name, " ranking is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  ranking %>%
    select(all_of(required_columns)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(
      source_gene_count = n(),
      macrophage_percentile = rank_macrophage / source_gene_count,
      monocyte_percentile = rank_monocyte / source_gene_count
    )
}

gse5099 <- read_ranking(gse5099_file, "GSE5099")
gse11864 <- read_ranking(gse11864_file, "GSE11864")

common <- inner_join(
  gse5099, gse11864,
  by = "gene_symbol",
  suffix = c("_GSE5099", "_GSE11864")
)

discordant <- common %>%
  filter(sign(moderated_t_GSE5099) != sign(moderated_t_GSE11864)) %>%
  transmute(
    gene_symbol,
    GSE5099_logFC = logFC_GSE5099,
    GSE5099_moderated_t = moderated_t_GSE5099,
    GSE11864_logFC = logFC_GSE11864,
    GSE11864_moderated_t = moderated_t_GSE11864,
    action = "excluded_before_consensus_ranking",
    reason = "opposite_differential_direction_between_sources"
  ) %>%
  arrange(gene_symbol)

build_direction <- function(common_table, direction_label, top_n) {
  if (direction_label == "mature_positive") {
    candidates <- common_table %>%
      filter(moderated_t_GSE5099 > 0, moderated_t_GSE11864 > 0) %>%
      mutate(
        GSE5099_rank_in_direction = rank_macrophage_GSE5099,
        GSE11864_rank_in_direction = rank_macrophage_GSE11864,
        GSE5099_percentile_in_direction = macrophage_percentile_GSE5099,
        GSE11864_percentile_in_direction = macrophage_percentile_GSE11864
      )
  } else {
    candidates <- common_table %>%
      filter(moderated_t_GSE5099 < 0, moderated_t_GSE11864 < 0) %>%
      mutate(
        GSE5099_rank_in_direction = rank_monocyte_GSE5099,
        GSE11864_rank_in_direction = rank_monocyte_GSE11864,
        GSE5099_percentile_in_direction = monocyte_percentile_GSE5099,
        GSE11864_percentile_in_direction = monocyte_percentile_GSE11864
      )
  }

  candidates %>%
    mutate(
      direction = direction_label,
      consensus_mean_percentile = (
        GSE5099_percentile_in_direction +
          GSE11864_percentile_in_direction
      ) / 2,
      consensus_worst_source_percentile = pmax(
        GSE5099_percentile_in_direction,
        GSE11864_percentile_in_direction
      )
    ) %>%
    arrange(
      consensus_mean_percentile,
      consensus_worst_source_percentile,
      gene_symbol
    ) %>%
    slice_head(n = top_n) %>%
    mutate(
      consensus_rank = row_number(),
      top_n_final_per_direction = top_n,
      selected_in_GSE5099_individual_top200 =
        GSE5099_rank_in_direction <= top_n,
      selected_in_GSE11864_individual_top200 =
        GSE11864_rank_in_direction <= top_n
    ) %>%
    transmute(
      top_n_final_per_direction,
      direction,
      consensus_rank,
      gene_symbol,
      consensus_mean_percentile,
      consensus_worst_source_percentile,
      GSE5099_rank_in_direction,
      GSE11864_rank_in_direction,
      GSE5099_percentile_in_direction,
      GSE11864_percentile_in_direction,
      selected_in_GSE5099_individual_top200,
      selected_in_GSE11864_individual_top200,
      GSE5099_logFC = logFC_GSE5099,
      GSE11864_logFC = logFC_GSE11864,
      GSE5099_moderated_t = moderated_t_GSE5099,
      GSE11864_moderated_t = moderated_t_GSE11864
    )
}

membership <- bind_rows(
  build_direction(common, "mature_positive", top_n),
  build_direction(common, "immature_negative", top_n)
)

if (
  sum(membership$direction == "mature_positive") != top_n ||
    sum(membership$direction == "immature_negative") != top_n ||
    anyDuplicated(membership[c("direction", "gene_symbol")]) ||
    length(intersect(
      membership$gene_symbol[membership$direction == "mature_positive"],
      membership$gene_symbol[membership$direction == "immature_negative"]
    ))
) {
  stop("Final consensus gene-set validation failed.")
}

selection_summary <- tibble(
  parameter = c(
    "external_sources",
    "genes_ranked_in_both_sources",
    "direction_discordant_genes_excluded",
    "direction_concordant_macrophage_candidates",
    "direction_concordant_monocyte_candidates",
    "final_mature_positive_genes",
    "final_immature_negative_genes",
    "consensus_ranking",
    "selection_uses_carotid_expression_or_subtype_results"
  ),
  value = c(
    "GSE5099;GSE11864",
    nrow(common),
    nrow(discordant),
    sum(common$moderated_t_GSE5099 > 0 & common$moderated_t_GSE11864 > 0),
    sum(common$moderated_t_GSE5099 < 0 & common$moderated_t_GSE11864 < 0),
    top_n,
    top_n,
    "mean_of_source_specific_directional_percentile_ranks",
    "FALSE"
  )
)

write.csv(membership, membership_file, row.names = FALSE)
write.csv(discordant, conflict_file, row.names = FALSE)
write.csv(selection_summary, summary_file, row.names = FALSE)

message(
  "Finished: exact ", top_n, " mature-positive and ", top_n,
  " monocyte-high genes written to ", membership_file
)
