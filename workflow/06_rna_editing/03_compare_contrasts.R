#!/usr/bin/env Rscript

# Compare direction and magnitude of editing changes across clone13, clone37,
# and the batch-adjusted combined analysis. All contrasts use delta = KO - WT.

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(purrr)
  library(tidyr)
})

root <- Sys.getenv("RNA_EDITING_RESULT_ROOT", "results/rna_editing")
out_dir <- file.path(root, "contrast_consistency")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_path <- function(contrast) {
  file.path(
    root, contrast, "Score_1", "Final_Results_Score1_PureLimma",
    "Integrated_Results_Score1_PureLimma.xlsx"
  )
}

read_contrast <- function(contrast) {
  path <- result_path(contrast)
  if (!file.exists(path)) stop("Missing editing result: ", path)
  sheets <- getSheetNames(path)
  map_dfr(sheets[grepl("^(C_to_U|A_to_I)_limma$", sheets)], function(sheet) {
    read.xlsx(path, sheet = sheet) %>%
      mutate(
        contrast = contrast,
        editing_type = sub("_limma$", "", sheet),
        site_id = paste(seqnames, start, end, editing_type, sep = ":")
      ) %>%
      select(
        contrast, editing_type, site_id, seqnames, start, end,
        gene_name, delta_ratio, final_p_value, final_fdr,
        significance, significance_fdr
      )
  })
}

all_sites <- map_dfr(c("clone13", "clone37", "combined"), read_contrast)
write.csv(all_sites, file.path(out_dir, "all_contrast_sites.csv"), row.names = FALSE)

pairwise_summary <- function(left, right) {
  x <- all_sites %>%
    filter(contrast == left) %>%
    select(site_id, editing_type, delta_left = delta_ratio)
  y <- all_sites %>%
    filter(contrast == right) %>%
    select(site_id, delta_right = delta_ratio)
  joined <- inner_join(x, y, by = "site_id") %>%
    mutate(
      comparison = paste(left, right, sep = "_vs_"),
      same_direction = sign(delta_left) == sign(delta_right)
    )
  write.csv(
    joined,
    file.path(out_dir, paste0(left, "_vs_", right, "_shared_sites.csv")),
    row.names = FALSE
  )
  joined %>%
    group_by(comparison, editing_type) %>%
    summarise(
      n_shared = n(),
      pearson_r = cor(delta_left, delta_right, method = "pearson", use = "complete.obs"),
      spearman_rho = cor(delta_left, delta_right, method = "spearman", use = "complete.obs"),
      same_direction_fraction = mean(same_direction, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_table <- bind_rows(
  pairwise_summary("clone13", "clone37"),
  pairwise_summary("clone13", "combined"),
  pairwise_summary("clone37", "combined")
)
write.csv(summary_table, file.path(out_dir, "pairwise_consistency_summary.csv"), row.names = FALSE)
message("Saved RNA-editing consistency analysis to: ", out_dir)
