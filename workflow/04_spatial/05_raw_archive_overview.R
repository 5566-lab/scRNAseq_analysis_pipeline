#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(yaml)
})

ROOT <- Sys.getenv(
  "AST_ROOT",
  unset = normalizePath("workflow/04_spatial", mustWork = TRUE)
)
source(file.path(ROOT, "R", "utils.R"))
cfg <- yaml::read_yaml(file.path(ROOT, "config", "config.yaml"))
fig_dir <- cfg$outputs$figures
tab_dir <- cfg$outputs$tables

summary_tab <- read.csv(file.path(tab_dir, "raw_archive_sample_summary.csv"), check.names = FALSE)
sample_tab <- read.csv(file.path(tab_dir, "raw_gse_sample_metadata.csv"), check.names = FALSE)
qc_tab <- read.csv(file.path(tab_dir, "raw_archive_qc_summary.csv"), check.names = FALSE)

p_summary <- summary_tab %>%
  pivot_longer(c(n_geo_samples, n_raw_files), names_to = "metric", values_to = "n") %>%
  ggplot(aes(x = gse, y = n, fill = platform_hint)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = "Count", fill = "Platform", title = "Additional raw spatial archives available for extension analyses") +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
save_figure(p_summary, file.path(fig_dir, "FigS1A_Raw_archive_inventory.pdf"), width = 7.2, height = 4.2)

severity_tab <- sample_tab %>%
  mutate(
    across(any_of(c("lesion_severity", "grade", "tissue_source", "source_name")), ~ na_if(.x, "")),
    across(any_of(c("lesion_severity", "grade", "tissue_source", "source_name")), ~ na_if(.x, "NA")),
    severity = case_when(
      gse == "GSE277441" ~ lesion_severity,
      gse == "GSE277170" ~ grade,
      gse == "GSE283269" ~ tissue_source,
      TRUE ~ coalesce(lesion_severity, grade, tissue_source, source_name)
    ),
    severity = ifelse(is.na(severity) | severity == "", "Not annotated", severity)
  ) %>%
  count(gse, severity, name = "n_samples")
write.csv(severity_tab, file.path(tab_dir, "raw_archive_sample_annotation_summary.csv"), row.names = FALSE)

p_annot <- ggplot(severity_tab, aes(x = gse, y = n_samples, fill = severity)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  labs(x = NULL, y = "GEO samples / ROIs", fill = "Annotation", title = "Raw archive disease-stage and source annotations") +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
save_figure(p_annot, file.path(fig_dir, "FigS1B_Raw_archive_sample_annotations.pdf"), width = 8.2, height = 4.4)

cosmx_qc <- qc_tab %>% filter(platform == "CosMx") %>% mutate(n_cells = as.numeric(n_cells))
if (nrow(cosmx_qc) > 0) {
  p_cosmx <- ggplot(cosmx_qc, aes(x = gsm, y = n_cells, fill = gsm)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    labs(x = NULL, y = "Cells in metadata file", fill = "GSM", title = "CosMx raw metadata cell counts") +
    theme_pub(9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  save_figure(p_cosmx, file.path(fig_dir, "FigS1C_CosMx_raw_cell_counts.pdf"), width = 5.6, height = 3.8)
}

visium_qc <- qc_tab %>%
  filter(platform == "Visium") %>%
  mutate(
    spots_under_tissue = as.numeric(spots_under_tissue),
    median_genes_per_spot = as.numeric(median_genes_per_spot)
  )
if (nrow(visium_qc) > 0) {
  p_visium <- visium_qc %>%
    select(gsm, spots_under_tissue, median_genes_per_spot) %>%
    pivot_longer(-gsm, names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = gsm, y = value, fill = metric)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = NULL, y = "Metric value", fill = "Metric", title = "DIT Visium raw metrics from Space Ranger summaries") +
    theme_pub(9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  save_figure(p_visium, file.path(fig_dir, "FigS1D_DIT_Visium_raw_metrics.pdf"), width = 8.4, height = 4.2)
}

geomx_qc <- qc_tab %>%
  filter(platform == "GeoMx") %>%
  mutate(
    Aligned = as.numeric(Aligned),
    total_probe_counts = as.numeric(total_probe_counts),
    grade = sample_tab$grade[match(gsm, sample_tab$gsm)],
    localisation = sample_tab$localisation[match(gsm, sample_tab$gsm)]
  ) %>%
  filter(!is.na(Aligned), !is.na(total_probe_counts), !is.na(grade), grade != "NA")
if (nrow(geomx_qc) > 0) {
  geomx_summary <- geomx_qc %>%
    group_by(grade, localisation) %>%
    summarise(n_roi = n(), median_aligned_reads = median(Aligned), median_probe_counts = median(total_probe_counts), .groups = "drop")
  write.csv(geomx_summary, file.path(tab_dir, "geomx_roi_qc_by_grade_localisation.csv"), row.names = FALSE)
  p_geomx <- ggplot(geomx_summary, aes(x = localisation, y = n_roi, fill = grade)) +
    geom_col(position = "dodge", width = 0.72, color = "white", linewidth = 0.2) +
    labs(x = NULL, y = "ROI count", fill = "Grade", title = "GeoMx ROI coverage by plaque severity and localization") +
    theme_pub(9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
  save_figure(p_geomx, file.path(fig_dir, "FigS1E_GeoMx_ROI_coverage.pdf"), width = 6.8, height = 4.2)
}

message("Raw archive overview figures complete")
