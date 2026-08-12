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

visium_stats <- read.csv(file.path(tab_dir, "visium_signature_severity_kruskal.csv"), check.names = FALSE) %>%
  mutate(dataset = "Visium_FFPE")
xenium_stats <- read.csv(file.path(tab_dir, "xenium_signature_severity_kruskal.csv"), check.names = FALSE) %>%
  mutate(dataset = "Xenium")

combined <- bind_rows(visium_stats, xenium_stats) %>%
  mutate(
    severe_minus_mild = severe_mean - mild_mean,
    direction = case_when(
      severe_minus_mild > 0 ~ "Higher in severe",
      severe_minus_mild < 0 ~ "Lower in severe",
      TRUE ~ "No change"
    )
  )
write.csv(combined, file.path(tab_dir, "cross_platform_signature_severity_summary.csv"), row.names = FALSE)

p_cross <- combined %>%
  mutate(signature = factor(signature, levels = unique(signature))) %>%
  ggplot(aes(x = signature, y = severe_minus_mild, fill = dataset)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, color = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_manual(values = c(Visium_FFPE = "#4C78A8", Xenium = "#B6424B")) +
  labs(x = NULL, y = "Mean score difference: Severe - Mild", fill = "Dataset", title = "Cross-platform severity trend of spatial plaque programs") +
  theme_pub()
save_figure(p_cross, file.path(fig_dir, "Fig7A_CrossPlatform_signature_severity_trends.pdf"), width = 7.2, height = 4.8)

visium_comp <- read.csv(file.path(tab_dir, "visium_celltype_composition_by_severity.csv"), check.names = FALSE) %>%
  mutate(dataset = "Visium_FFPE")
xenium_comp <- read.csv(file.path(tab_dir, "xenium_celltype_composition_by_severity.csv"), check.names = FALSE) %>%
  mutate(dataset = "Xenium")

shared_types <- c("Myeloid", "SMCPericyte", "ModSMC", "Endothelium", "TCells", "BCells", "Fibroblast1", "Fibroblast2")
comp_combined <- bind_rows(visium_comp, xenium_comp) %>%
  filter(cell_type %in% shared_types) %>%
  mutate(group = factor(group, levels = severity_levels))
write.csv(comp_combined, file.path(tab_dir, "cross_platform_celltype_composition_summary.csv"), row.names = FALSE)

p_comp <- ggplot(comp_combined, aes(x = group, y = frac, color = dataset, group = dataset)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ cell_type, ncol = 4, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c(Visium_FFPE = "#4C78A8", Xenium = "#B6424B")) +
  labs(x = NULL, y = "Fraction", color = "Dataset", title = "Cross-platform cell-composition shifts across plaque severity") +
  theme_pub(9)
save_figure(p_comp, file.path(fig_dir, "Fig7B_CrossPlatform_celltype_composition_trends.pdf"), width = 9.8, height = 6.8)

summary_text <- c(
  "# Spatial Atherosclerosis Analysis Summary",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "## Main evidence chain",
  "",
  "- Xenium provides cell-level validation of plaque-severity-associated myeloid, foam-cell, inflammatory, and APOBEC3A-related programs.",
  "- Visium provides whole-transcriptome spatial spot-level validation across 16 samples and enables module-level spatial maps.",
  "- Cross-platform tables compare severe-versus-mild score shifts and cell-composition trends.",
  "",
  "## Key output tables",
  "",
  "- `cross_platform_signature_severity_summary.csv`",
  "- `cross_platform_celltype_composition_summary.csv`",
  "- `xenium_myeloid_neighbor_summary.csv`",
  "- `visium_marker_expression_by_severity.csv`",
  "- `xenium_marker_expression_by_celltype_severity.csv`",
  "",
  "## Interpretation guardrails",
  "",
  "- The current figures support spatial association and disease-severity trends, not direct causality.",
  "- Xenium panel signatures are limited by the targeted gene panel; Visium signatures are broader and should be used for pathway-level claims.",
  "- Representative spatial maps use object coordinates and should be paired with histology overlays in the final manuscript if journal figure style requires it."
)
writeLines(summary_text, file.path(tab_dir, "analysis_summary.md"))
message("Cross-platform summary complete")
