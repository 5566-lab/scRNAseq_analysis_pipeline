#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(patchwork)
})

base_dir <- Sys.getenv(
  "SCRNA_MONOCLE_DIR",
  unset = "/public3/DSC/single_cell/Result/figer_new/monocle3"
)
out_dir <- file.path(base_dir, "Core_vs_Adjacent_from_main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(Sys.getenv(
  "SCRNA_PSEUDOTIME_RDS",
  unset = file.path(base_dir, "data_pseudotime.rds")
))
meta <- as.data.frame(SummarizedExperiment::colData(cds))
meta$pseudotime <- monocle3::pseudotime(cds)
meta$Region <- factor(
  meta$Sample_Type,
  levels = c("Proximal Adjacent", "Atherosclerotic Core")
)

plot_df <- meta %>%
  filter(is.finite(pseudotime), !is.na(Region), !is.na(Celltype_raw1))

if (nlevels(droplevels(plot_df$Region)) != 2L) {
  stop("Both Proximal Adjacent and Atherosclerotic Core cells are required.")
}

# Derive the branch point from the original interactive selections preserved in
# data_pseudotime metadata. Later checkpoints contain a different subset_2
# selection and must not be used to reconstruct the published figure.
if (!all(c("subset_1", "subset_2") %in% colnames(meta))) {
  stop("The saved cds_MM_foam object does not contain subset_1/subset_2 labels.")
}
fate1_cells <- rownames(meta)[as.character(meta$subset_1) == "Fate 1"]
fate2_cells <- rownames(meta)[as.character(meta$subset_2) == "Fate 2"]
shared_cells <- intersect(fate1_cells, fate2_cells)
if (length(shared_cells) == 0L) {
  stop("The original Fate 1 and Fate 2 selections do not share a trajectory trunk.")
}
branch_t <- max(meta[shared_cells, "pseudotime"], na.rm = TRUE)

pt_min <- min(plot_df$pseudotime)
pt_max <- max(plot_df$pseudotime)
pt_grid <- seq(pt_min, pt_max, length.out = 1000)

dominant_type <- function(region_name) {
  region_df <- filter(plot_df, Region == region_name)
  cell_types <- unique(region_df$Celltype_raw1)
  density_matrix <- sapply(cell_types, function(cell_type) {
    values <- region_df$pseudotime[region_df$Celltype_raw1 == cell_type]
    if (length(values) < 5L) return(rep(0, length(pt_grid)))
    fitted <- density(values, from = pt_min, to = pt_max, n = length(pt_grid))
    fitted$y * length(values) / nrow(region_df)
  })
  data.frame(
    pseudotime = pt_grid,
    Dominant_Type = cell_types[apply(density_matrix, 1, which.max)],
    Region = region_name
  )
}

dominant_df <- bind_rows(lapply(levels(plot_df$Region), dominant_type)) %>%
  mutate(Region = factor(Region, levels = levels(plot_df$Region)))

# Use run-length groups so repeated dominance windows receive separate labels.
dominant_df <- dominant_df %>%
  group_by(Region) %>%
  mutate(run = cumsum(Dominant_Type != lag(Dominant_Type, default = first(Dominant_Type)))) %>%
  ungroup()

track_labels <- dominant_df %>%
  count(Region, run, Dominant_Type, name = "n_grid") %>%
  filter(n_grid >= 45L) %>%
  left_join(
    dominant_df %>% group_by(Region, run, Dominant_Type) %>%
      summarise(mid_x = median(pseudotime), .groups = "drop"),
    by = c("Region", "run", "Dominant_Type")
  )

region_colors <- c(
  "Proximal Adjacent" = "#2F6B9A",
  "Atherosclerotic Core" = "#D9772B"
)
celltype_colors <- c(
  "Classical Mono" = "#D97B72", "Inflammatory Mono" = "#C44E52",
  "ISG+ Mono" = "#4EA8C7", "Non-classical Mono" = "#B565A7",
  "Transitional Mac" = "#8172B3", "Foam cells1" = "#79A943",
  "Foam cells2" = "#55AE8A", "LAM" = "#8E82B8",
  "CX3CR1+ TRM" = "#BE9B2F", "LYVE1+ TRM" = "#4C956C"
)

density_values <- lapply(split(plot_df$pseudotime, plot_df$Region), density)
max_density <- max(vapply(density_values, function(x) max(x$y), numeric(1)))

p_density <- ggplot(plot_df, aes(pseudotime, fill = Region, color = Region)) +
  geom_density(alpha = 0.30, linewidth = 1.25, adjust = 1) +
  geom_vline(xintercept = branch_t, linetype = "dashed", color = "grey45", linewidth = 0.9) +
  annotate("text", x = branch_t - 0.55, y = max_density * 0.84,
           label = "Bifurcation", angle = 90, color = "grey30",
           size = 4.5, fontface = "italic") +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  scale_x_continuous(limits = c(pt_min, pt_max), expand = c(0, 0)) +
  labs(title = "Pseudotime Distribution by Plaque Region", y = "Cell Density") +
  theme_classic(base_size = 15) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.line.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.position = "top", legend.title = element_blank()
  )

p_tracks <- ggplot(dominant_df, aes(pseudotime, Region, fill = Dominant_Type)) +
  geom_tile(height = 0.8) +
  geom_vline(xintercept = branch_t, linetype = "dashed", color = "black", linewidth = 0.9) +
  geom_text(data = track_labels,
            aes(mid_x, Region, label = Dominant_Type), inherit.aes = FALSE,
            size = 3.5, fontface = "bold", check_overlap = TRUE) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_x_continuous(limits = c(pt_min, pt_max), expand = c(0, 0)) +
  labs(x = "Pseudotime", fill = "Dominant Cell Type") +
  theme_classic(base_size = 15) +
  theme(
    text = element_text(family = "Arial"),
    axis.title.y = element_blank(), axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom", legend.title = element_text(face = "bold")
  )

final_plot <- p_density / p_tracks + plot_layout(heights = c(4, 1.35))

ggsave(file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Density.pdf"),
       final_plot, width = 12, height = 7.5, device = cairo_pdf)
ggsave(file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Density.png"),
       final_plot, width = 12, height = 7.5, dpi = 400)

celltype_summary <- plot_df %>%
  group_by(Region, Celltype_raw1) %>%
  summarise(
    n_cells = n(),
    median_pseudotime = median(pseudotime),
    q25_pseudotime = quantile(pseudotime, 0.25),
    q75_pseudotime = quantile(pseudotime, 0.75),
    .groups = "drop"
  )

celltype_order <- plot_df %>%
  group_by(Celltype_raw1) %>%
  summarise(median_pseudotime = median(pseudotime), .groups = "drop") %>%
  arrange(median_pseudotime) %>%
  pull(Celltype_raw1)

# Put early states at the top and later states at the bottom.
ridge_df <- plot_df %>%
  mutate(
    Celltype_raw1 = factor(
      Celltype_raw1,
      levels = rev(celltype_order)
    )
  )

# Shared colors for the count-based ridgeline variants.
ridge_region_colors <- c(
  "Proximal Adjacent" = "#4C72B0",
  "Atherosclerotic Core" = "#D55E00"
)

# 按伪时序 bin 计数（非核密度归一化），突出细胞量变化。
ridge_binwidth <- (pt_max - pt_min) / 45
region_n <- plot_df %>% count(Region, name = "region_n")
ridge_df_norm <- left_join(ridge_df, region_n, by = "Region") %>%
  mutate(weight_by_region = 1 / region_n)

p_ridges <- ggplot(
  ridge_df,
  aes(
    x = pseudotime,
    y = Celltype_raw1,
    fill = Region,
    color = Region,
    group = interaction(Celltype_raw1, Region)
  )
) +
  geom_density_ridges(
    stat = "binline",
    binwidth = ridge_binwidth,
    position = "identity",
    scale = 1.35,
    rel_min_height = 0.005,
    alpha = 0.52,
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = branch_t,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = ridge_region_colors) +
  scale_color_manual(values = ridge_region_colors) +
  scale_x_continuous(
    limits = c(pt_min, pt_max),
    expand = expansion(mult = c(0.005, 0.015))
  ) +
  labs(
    title = "Pseudotime distribution of MPS cell states",
    subtitle = "Bin-count ridgelines (absolute cell counts)",
    x = "Pseudotime",
    y = NULL,
    fill = "Region",
    color = "Region"
  ) +
  theme_ridges(font_size = 13, grid = TRUE, center_axis_labels = TRUE) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 18,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 8, 0),
    legend.background = element_rect(
      fill = "white", color = "grey80",
      linewidth = 0.4, size = 0.4
    ),
    legend.box.background = element_rect(
      fill = "white", color = "grey85",
      linewidth = 0.4, size = 0.4
    ),
    legend.box.margin = margin(2, 2, 2, 2)
  )

# 按区域内细胞数归一化后的 bin 计数（方便比较两个区域分布差异，不受总细胞数影响）。
p_ridges_norm <- ggplot(
  ridge_df_norm,
  aes(
    x = pseudotime,
    y = Celltype_raw1,
    fill = Region,
    color = Region,
    group = interaction(Celltype_raw1, Region),
    weight = weight_by_region
  )
) +
  geom_density_ridges(
    stat = "binline",
    binwidth = ridge_binwidth,
    position = "identity",
    scale = 1.35,
    rel_min_height = 0.005,
    alpha = 0.52,
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = branch_t,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = ridge_region_colors) +
  scale_color_manual(values = ridge_region_colors) +
  scale_x_continuous(
    limits = c(pt_min, pt_max),
    expand = expansion(mult = c(0.005, 0.015))
  ) +
  labs(
    title = "Pseudotime distribution of MPS cell states",
    subtitle = "Bin-count ridgelines (within-region normalized)",
    x = "Pseudotime",
    y = NULL,
    fill = "Region",
    color = "Region"
  ) +
  theme_ridges(font_size = 13, grid = TRUE, center_axis_labels = TRUE) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = 18,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 8, 0),
    legend.background = element_rect(
      fill = "white", color = "grey80",
      linewidth = 0.4, size = 0.4
    ),
    legend.box.background = element_rect(
      fill = "white", color = "grey85",
      linewidth = 0.4, size = 0.4
    ),
    legend.box.margin = margin(2, 2, 2, 2)
  )

ggsave(
  file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Ridgelines_overlay_Counts.pdf"),
  p_ridges, width = 14, height = 8.5, device = cairo_pdf
)
ggsave(
  file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Ridgelines_overlay_Counts.png"),
  p_ridges, width = 14, height = 8.5, dpi = 400, bg = "white"
)

ggsave(
  file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Ridgelines_overlay_Normalized.pdf"),
  p_ridges_norm, width = 14, height = 8.5, device = cairo_pdf
)
ggsave(
  file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Ridgelines_overlay_Normalized.png"),
  p_ridges_norm, width = 14, height = 8.5, dpi = 400, bg = "white"
)

write.csv(
  celltype_summary,
  file.path(out_dir, "Core_vs_Adjacent_celltype_pseudotime_summary.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# 丰度缩放山脊图（Abundance-scaled ridgelines）
# 每个细胞类型内部，把核心/旁组织的密度曲线高度按“该类型在这片区域的占比”缩放，
# 从而同时保留 伪时间轴 + 山脊形状，并直观反映丰度差异（如 FC1 核心 91.7% vs 旁 8.3%）。
# 用占比而非原始细胞数，天然规避“两个样品细胞总数差异大”的问题。
# ---------------------------------------------------------------------------
abun_grid <- seq(pt_min, pt_max, length.out = 500)
region_levels <- c("Atherosclerotic Core", "Proximal Adjacent")

abun_list <- list()
# 每个区域归一化到相同的总面积：
# height 正比于“该细胞类型在该区域内部所占比例”，而非“细胞类型内部的核心/旁占比”。
# 这样核心与旁总面积相当，比较的是相对富集/耗竭，避免核心因样本量大而全面占优。
region_totals <- plot_df %>% count(Region) %>% tibble::deframe()
for (ct in celltype_order) {
  for (rg in region_levels) {
    vals <- plot_df$pseudotime[plot_df$Celltype_raw1 == ct & plot_df$Region == rg]
    n <- length(vals)
    if (n < 5) next
    prop <- n / region_totals[[rg]]  # 该细胞类型在该区域内的占比 (0-1)
    d <- density(vals, from = pt_min, to = pt_max, n = length(abun_grid))
    abun_list[[length(abun_list) + 1L]] <- data.frame(
      x = abun_grid,
      y = ct,
      height = d$y * prop,  # 面积 = 区域内占比
      Region = rg
    )
  }
}
abun_df <- do.call(rbind, abun_list)
# 每个细胞类型各自归一化（放开共用 Y 轴）：每个类型内最高的山脊高度≈1，
# 只保留“类型内部核心 vs 旁谁高”的对比，不再横向比较不同类型的丰度。
abun_df <- abun_df %>%
  group_by(y) %>%
  mutate(height = height / max(height, na.rm = TRUE) * 0.9) %>%
  ungroup()
abun_df$y <- factor(abun_df$y, levels = rev(celltype_order))
abun_df$Region <- factor(abun_df$Region, levels = region_levels)

# Publication-style warm/cool pairing used only by the abundance-scaled panel.
abun_region_fills <- c(
  "Proximal Adjacent" = "#4DBBD5",
  "Atherosclerotic Core" = "#E64B35"
)
abun_region_lines <- c(
  "Proximal Adjacent" = "#287F96",
  "Atherosclerotic Core" = "#B6382A"
)

p_abun <- ggplot() +
  geom_ridgeline(
    data = filter(abun_df, Region == "Atherosclerotic Core"),
    aes(x, y, height = height, fill = Region, color = Region),
    alpha = 0.42, linewidth = 0.55, min_height = 0
  ) +
  geom_ridgeline(
    data = filter(abun_df, Region == "Proximal Adjacent"),
    aes(x, y, height = height, fill = Region, color = Region),
    alpha = 0.42, linewidth = 0.55, min_height = 0
  ) +
  geom_vline(xintercept = branch_t, linetype = "22", color = "grey38", linewidth = 0.7) +
  annotate(
    "text",
    x = branch_t - 0.45,
    y = "LAM",
    label = "Bifurcation",
    angle = 90,
    color = "grey32",
    size = 4,
    family = "Arial",
    fontface = "italic"
  ) +
  scale_fill_manual(values = abun_region_fills) +
  scale_color_manual(values = abun_region_lines) +
  scale_x_continuous(limits = c(pt_min, pt_max), expand = expansion(mult = c(0.005, 0.015))) +
  labs(
    title = "Pseudotime Distribution Across MPS Cell Types",
    subtitle = "Ridge area scaled by cell-type abundance within each region",
    x = "Pseudotime",
    y = NULL
  ) +
  theme_ridges(font_size = 12, grid = TRUE, center_axis_labels = TRUE) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17, margin = margin(b = 3)),
    plot.subtitle = element_text(hjust = 0.5, size = 10.5, color = "grey32", margin = margin(b = 5)),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 5)),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.text.y = element_text(size = 10.5, face = "bold", color = "black", margin = margin(r = 3)),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 10.5, face = "bold"),
    legend.key.width = unit(0.52, "cm"),
    legend.key.height = unit(0.36, "cm"),
    legend.spacing.x = unit(0.25, "cm"),
    legend.margin = margin(0, 0, 5, 0),
    legend.background = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(8, 12, 8, 10)
  )

ggsave(
  file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Ridgelines_abundance.pdf"),
  p_abun, width = 13, height = 8.8, device = cairo_pdf
)
ggsave(
  file.path(out_dir, "F2.6_Core_vs_Adjacent_Pseudotime_Ridgelines_abundance.png"),
  p_abun, width = 13, height = 8.8, dpi = 400, bg = "white",
  device = png, type = "cairo"
)

message(sprintf(
  "Saved regional density and cell-type ridgeline figures; branch point = %.3f",
  branch_t
))
