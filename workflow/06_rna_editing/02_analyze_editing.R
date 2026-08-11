# ===================================================================
# 1. 基础库加载与全局设置
# ===================================================================
options(clusterProfiler.check.update = FALSE)
options(warn = 1) # 强制警告实时打印，不延迟

suppressMessages({
  library(plyranges)
  library(magrittr)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(openxlsx)
  library(JACUSA2helper)
  library(limma)
  library(stringr)
  library(DOSE)
  library(ggrepel)
})

# ===================================================================
# 2. 设置路径和全局变量
# ===================================================================
env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

# 与上游 JACUSA2 完全一致：cond1 = WT，cond2 = KO；已排除 THP_1。
jacusa_output_dir <- env_path(
  "RNA_EDITING_JACUSA_OUTPUT_DIR",
  "results/rna_editing/jacusa2_combined"
)
input_file <- env_path(
  "RNA_EDITING_JACUSA_OUT",
  file.path(jacusa_output_dir, "Aggregated_WTcond1_KOcond2_no_THP1.out")
)
sample_manifest_file <- env_path(
  "RNA_EDITING_SAMPLE_MANIFEST",
  file.path(jacusa_output_dir, "jacusa_sample_manifest.tsv")
)
gtf_file <- env_path(
  "RNA_EDITING_ANNOTATION_GFF3",
  "/public7/DSC_Public7/DSC/DSC/20251004_NGS/RawData/Homo_sapiens.GRCh38.104.gff3"
)

main_output_dir <- env_path(
  "RNA_EDITING_ANALYSIS_OUTPUT_DIR",
  "results/rna_editing/combined"
)
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

PVAL_CUTOFF <- as.numeric(Sys.getenv("RNA_EDITING_PVALUE_CUTOFF", "0.1"))
DELTA_CUTOFF <- as.numeric(Sys.getenv("RNA_EDITING_DELTA_CUTOFF", "0.05"))
WT_MIN_COVERED_REPS <- as.integer(Sys.getenv("RNA_EDITING_WT_MIN_REPS", "4"))
KO_MIN_COVERED_REPS <- as.integer(Sys.getenv("RNA_EDITING_KO_MIN_REPS", "3"))
EXPECTED_REPLICATES <- as.integer(Sys.getenv("RNA_EDITING_EXPECTED_REPS", "6"))

kegg_offline_file <- env_path(
  "RNA_EDITING_KEGG_RDS",
  "/public7/DSC_Public7/DSC/DSC/20251004_NGS/result/M013_result_Score1_Hybrid/kegg_hsa_offline.rds"
)
if (file.exists(kegg_offline_file)) {
  cat("\n✅ 成功找到离线 KEGG 文件，全管线将启用超高速离线富集模式！\n")
  kegg_local <- readRDS(kegg_offline_file)
} else {
  stop("❌ 错误：离线 KEGG 文件不存在，请检查路径！")
}

# ===================================================================
# 3. 定义辅助函数 (纯净版，剔除旧算法)
# ===================================================================

calculate_ratios <- function(rep_tibble) {
  rep_tibble$A <- as.numeric(rep_tibble$A); rep_tibble$C <- as.numeric(rep_tibble$C)
  rep_tibble$G <- as.numeric(rep_tibble$G); rep_tibble$T <- as.numeric(rep_tibble$T)
  total <- rowSums(rep_tibble[c("A","C","G","T")], na.rm = TRUE)
  safe_div <- function(num, den) ifelse(den == 0, 0, num / den)
  dplyr::tibble(A_ratio = safe_div(rep_tibble$A, total), C_ratio = safe_div(rep_tibble$C, total),
                G_ratio = safe_div(rep_tibble$G, total), T_ratio = safe_div(rep_tibble$T, total))
}

unnest_ratios <- function(condition_ratios_tibble, condition_name) {
  rep_names <- names(condition_ratios_tibble)
  condition_ratios_tibble$query_idx <- 1:nrow(condition_ratios_tibble)
  long_df <- tidyr::pivot_longer(condition_ratios_tibble, cols = all_of(rep_names), names_to = "replicate", values_to = "ratios_tibble")
  long_df <- tidyr::unnest(long_df, cols = c(ratios_tibble))
  long_df$condition <- condition_name
  return(long_df)
}

# ---------- 全局 PCA 分析模块 ----------
plot_global_pca <- function(filtered_data, sample_info, output_dir, score_label) {
  cat("\n========== 📊 开始绘制全局 RNA 编辑 PCA 图 ==========\n")
  refs <- as.character(mcols(filtered_data)$ref)

  get_editing_ratio <- function(sample_bases) {
    mat <- as.matrix(sample_bases[, c("A", "C", "G", "T")])
    total <- rowSums(mat, na.rm = TRUE)
    total[total == 0] <- 1
    ref_idx <- cbind(1:nrow(mat), match(refs, c("A", "C", "G", "T")))
    ref_counts <- mat[ref_idx]
    return((total - ref_counts) / total)
  }

  wt_ratios <- do.call(cbind, lapply(filtered_data$bases$cond1, get_editing_ratio))  # cond1 = WT
  ko_ratios <- do.call(cbind, lapply(filtered_data$bases$cond2, get_editing_ratio))  # cond2 = KO
  colnames(wt_ratios) <- paste0("cond1_", names(filtered_data$bases$cond1))
  colnames(ko_ratios) <- paste0("cond2_", names(filtered_data$bases$cond2))
  ratio_matrix <- cbind(wt_ratios, ko_ratios)

  ratio_matrix <- na.omit(ratio_matrix)
  ratio_matrix <- ratio_matrix[apply(ratio_matrix, 1, var) > 0, ]

  if(nrow(ratio_matrix) < 10) return(NULL)

  overall_plot_dir <- file.path(output_dir, "Overall_Plots")
  dir.create(overall_plot_dir, recursive = TRUE, showWarnings = FALSE)

  generate_pca_plot <- function(mat, title, subtitle, filename) {
    pca_res <- prcomp(t(mat), scale. = TRUE)
    sum_pca <- summary(pca_res)
    pca_df <- as.data.frame(pca_res$x)
    pca_df$Sample_ID <- rownames(pca_df)
    pca_df <- dplyr::left_join(pca_df, sample_info, by = "Sample_ID")

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Batch)) +
      geom_point(size = 5, alpha = 0.85) +
      geom_text_repel(aes(label = Real_Sample_Name), size = 3.5, show.legend = FALSE, max.overlaps = 20) +
      scale_color_manual(values = c("WT" = "#377EB8", "KO" = "#E41A1C")) +
      labs(title = title, subtitle = subtitle, x = paste0("PC1 (", round(sum_pca$importance[2,"PC1"]*100,1), "%)"), y = paste0("PC2 (", round(sum_pca$importance[2,"PC2"]*100,1), "%)")) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(overall_plot_dir, filename), plot = p, width = 9, height = 7)
    return(p)
  }

  # 图1: Raw PCA
  generate_pca_plot(ratio_matrix, "PCA - Raw Editing Profiles", "Before Batch Correction", paste0("PCA_Raw_Score", score_label, ".pdf"))

  # 图2: 去批次 PCA
  if("Batch" %in% colnames(sample_info) && length(unique(sample_info$Batch)) > 1) {
    matched_batches <- sample_info$Batch[match(colnames(ratio_matrix), sample_info$Sample_ID)]
    clean_matrix <- limma::removeBatchEffect(ratio_matrix, batch = matched_batches)
    generate_pca_plot(clean_matrix, "PCA - Batch Corrected", "After limma::removeBatchEffect", paste0("PCA_BatchCorrected_Score", score_label, ".pdf"))
  }
}

run_enrichment_GO <- function(gene_list, gene_prefix, plot_dir, analysis_type) {
  if (length(gene_list) == 0) return(character(0))
  hg <- tryCatch({ suppressMessages(bitr(gene_list, fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db", drop = TRUE)) }, error = function(e) NULL)
  if (is.null(hg) || nrow(hg) == 0) return(character(0))

  go_results_df <- data.frame()
  for (ont in c("BP", "CC", "MF")) {
    try({
      go_ont <- suppressMessages(enrichGO(hg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = ont, pAdjustMethod = 'BH', pvalueCutoff = PVAL_CUTOFF, qvalueCutoff = 1.0))
      if (!is.null(go_ont) && nrow(go_ont) > 0) {
        go_simp <- clusterProfiler::simplify(go_ont, cutoff = 0.7, by = "p.adjust", select_fun = min)
        if (nrow(go_simp) > 0) {
          df <- as.data.frame(go_simp); df$ONTOLOGY <- ont
          go_results_df <- rbind(go_results_df, df)
          if (ont %in% c("BP", "CC")) {
            ggsave(file.path(plot_dir, paste0(analysis_type, "_", gene_prefix, "_go_", ont, "_limma_dotplot.pdf")), plot = dotplot(go_simp, showCategory=15, title=paste("GO", ont, "-", gene_prefix)), width = 9, height = 7)
            ggsave(file.path(plot_dir, paste0(analysis_type, "_", gene_prefix, "_go_", ont, "_limma_emap.pdf")), plot = emapplot(pairwise_termsim(go_simp), showCategory=15), width = 9, height = 7)
          }
        }
      }
    }, silent = TRUE)
  }
  if (nrow(go_results_df) > 0) {
    write.csv(go_results_df, file=file.path(plot_dir, paste0(analysis_type, "_", gene_prefix, "_go_limma.csv")), row.names = FALSE)
    return(head(go_results_df$Description[order(go_results_df$p.adjust)], 20))
  }
  return(character(0))
}

run_enrichment_KEGG <- function(gene_list, gene_prefix, plot_dir, analysis_type) {
  if (length(gene_list) == 0) return(character(0))
  hg <- tryCatch({ suppressMessages(bitr(gene_list, fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db", drop = TRUE)) }, error = function(e) NULL)
  if (is.null(hg) || nrow(hg) == 0) return(character(0))

  try({
    kegg <- suppressMessages(enricher(hg$ENTREZID, TERM2GENE=kegg_local$KEGGPATHID2EXTID, TERM2NAME=kegg_local$KEGGPATHID2NAME, pvalueCutoff=PVAL_CUTOFF, pAdjustMethod='BH', qvalueCutoff=1.0))
    if (!is.null(kegg) && nrow(kegg) > 0) {
      kegg_df <- as.data.frame(kegg)
      if (nrow(kegg_df) > 0) {
        write.csv(kegg_df, file=file.path(plot_dir, paste0(analysis_type, "_", gene_prefix, "_kegg_limma.csv")), row.names = FALSE)
        ggsave(file.path(plot_dir, paste0(analysis_type, "_", gene_prefix, "_kegg_limma_dotplot.pdf")), plot = dotplot(kegg, showCategory=20, title=paste("KEGG -", gene_prefix)), width = 10, height = 8)
        return(head(kegg_df$Description[order(kegg_df$p.adjust)], 20))
      }
    }
  }, silent = TRUE)
  return(character(0))
}


# ===================================================================
# 3.5 后处理与补图函数：从 M013 迁移并适配 PureLimma / 多 Score
# ===================================================================

plot_overall_base_substitution <- function(pooled_df, overall_plot_dir, score_label) {
  tryCatch({
    if (is.null(pooled_df) || nrow(pooled_df) == 0) return(NULL)
    dir.create(overall_plot_dir, recursive = TRUE, showWarnings = FALSE)

    count_table <- table(KO = pooled_df$ref2ko_type, WT = pooled_df$ref2wt_type)
    heatmap_df <- as.data.frame(count_table)

    heatmap_plot <- ggplot(heatmap_df, aes(x = WT, y = KO, fill = Freq)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_text(aes(label = ifelse(Freq == 0, "", Freq)), color = "black", size = 2.5) +
      scale_fill_gradientn(colours = c("white", "#00A6CA", "#F29E2E", "#D7191C")) +
      labs(
        title = paste0("Base Substitution Comparison (Score >= ", score_label, ")"),
        x = "WT base substitution",
        y = "KO base substitution",
        fill = "Count"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    ggsave(
      filename = file.path(overall_plot_dir, paste0("overall_heatmap_ggplot_Score", score_label, ".pdf")),
      plot = heatmap_plot,
      width = 14,
      height = 11
    )

    bar_df <- pooled_df %>%
      dplyr::select(query_idx, ref2wt_type, ref2ko_type) %>%
      tidyr::pivot_longer(
        cols = c("ref2wt_type", "ref2ko_type"),
        names_to = "Condition",
        values_to = "Substitution"
      ) %>%
      dplyr::filter(Substitution != "no change") %>%
      dplyr::count(Condition, Substitution) %>%
      dplyr::mutate(Condition = ifelse(Condition == "ref2wt_type", "WT", "KO"))

    if (nrow(bar_df) > 0) {
      barplot_gg <- ggplot(bar_df, aes(x = Substitution, y = n, fill = Condition)) +
        geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
        scale_fill_manual(values = c("WT" = "#377EB8", "KO" = "#E41A1C")) +
        theme_minimal() +
        labs(
          title = paste0("Overall Base Substitution Counts (Score >= ", score_label, ")"),
          x = "Substitution Type",
          y = "Count"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          panel.grid.minor = element_blank()
        )

      ggsave(
        filename = file.path(overall_plot_dir, "overall_base_substitution_barplot.pdf"),
        plot = barplot_gg,
        width = 10,
        height = 6
      )
    }

    cat("  ✅ Overall base substitution heatmap/barplot 已补齐\n")
  }, error = function(e) {
    cat(paste0("  ⚠️ Overall base substitution plots 生成失败: ", e$message, "\n"))
  })
}

redraw_kegg_emap_from_csv <- function(target_dir, show_category = 15) {
  tryCatch({
    cat("\n▶ 后处理：扫描 KEGG CSV 并重绘 Emap 网络图...\n")

    csv_files <- list.files(
      target_dir,
      pattern = "(kegg_limma\\.csv$|_KEGG\\.csv$)",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(csv_files) == 0) {
      cat("  ⚠️ 未找到 KEGG CSV 文件，跳过 Redrawn_Emap。\n")
      return(NULL)
    }

    for (f in csv_files) {
      df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0 || !"geneID" %in% names(df) || !"Description" %in% names(df)) next

      df <- df %>% dplyr::filter(!is.na(geneID), geneID != "")
      if (nrow(df) < 2) next

      if ("p.adjust" %in% names(df)) {
        df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
        df_plot <- df[order(df$p.adjust), , drop = FALSE]
      } else if ("pvalue" %in% names(df)) {
        df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
        df$p.adjust <- df$pvalue
        df_plot <- df[order(df$pvalue), , drop = FALSE]
      } else {
        df$p.adjust <- 1
        df_plot <- df
      }

      df_plot <- head(df_plot, show_category)
      if (nrow(df_plot) < 2) next

      if ("ID" %in% names(df_plot)) {
        rownames(df_plot) <- make.unique(as.character(df_plot$ID))
      } else {
        df_plot$ID <- paste0("Term_", seq_len(nrow(df_plot)))
        rownames(df_plot) <- df_plot$ID
      }

      gene_sets <- strsplit(as.character(df_plot$geneID), "/")
      names(gene_sets) <- rownames(df_plot)
      all_genes <- unique(unlist(gene_sets))
      all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

      dummy_enrich_obj <- methods::new(
        "enrichResult",
        result = df_plot,
        pvalueCutoff = 1.0,
        pAdjustMethod = "none",
        qvalueCutoff = 1.0,
        gene = all_genes,
        universe = character(0),
        geneSets = gene_sets,
        organism = "UNKNOWN",
        keytype = "UNKNOWN",
        ontology = "KEGG",
        readable = FALSE
      )

      emap_obj <- tryCatch(pairwise_termsim(dummy_enrich_obj), error = function(e) NULL)
      if (is.null(emap_obj)) next

      base_name <- tools::file_path_sans_ext(basename(f))
      out_file <- file.path(dirname(f), paste0("Redrawn_Emap_", base_name, ".pdf"))

      p_emap <- emapplot(emap_obj, showCategory = min(show_category, nrow(df_plot)), color = "p.adjust", layout = "kk") +
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
        labs(title = paste("KEGG Emap -", base_name))

      ggsave(out_file, plot = p_emap, width = 10, height = 8)
      cat("  ✅ 已生成: ", out_file, "\n", sep = "")
    }
  }, error = function(e) {
    cat(paste0("  ⚠️ KEGG Emap 后处理失败: ", e$message, "\n"))
  })
}

plot_editing_mrna_relationship <- function(all_sites, log2fc_df, cross_plot_dir) {
  tryCatch({
    cat("\n▶ 后处理：绘制 RNA editing 与 mRNA 表达四象限关系图...\n")

    if (is.null(all_sites) || nrow(all_sites) == 0 || is.null(log2fc_df) || nrow(log2fc_df) == 0) {
      cat("  ⚠️ all_sites 或 log2fc_df 为空，跳过四象限关系图。\n")
      return(NULL)
    }

    required_cols <- c("gene_name", "Editing_Type", "significance", "delta_ratio")
    if (!all(required_cols %in% names(all_sites)) || !all(c("Gene_Symbol", "Transcriptome_Log2FC") %in% names(log2fc_df))) {
      cat("  ⚠️ 缺少必要列，跳过四象限关系图。\n")
      return(NULL)
    }

    dir.create(cross_plot_dir, showWarnings = FALSE, recursive = TRUE)

    correlation_data <- all_sites %>%
      dplyr::filter(significance %in% c("Up", "Down")) %>%
      dplyr::inner_join(log2fc_df, by = c("gene_name" = "Gene_Symbol")) %>%
      dplyr::mutate(
        delta_ratio = suppressWarnings(as.numeric(delta_ratio)),
        Transcriptome_Log2FC = suppressWarnings(as.numeric(Transcriptome_Log2FC))
      ) %>%
      dplyr::filter(!is.na(delta_ratio), !is.na(Transcriptome_Log2FC)) %>%
      dplyr::mutate(
        Relationship = dplyr::case_when(
          (delta_ratio > 0 & Transcriptome_Log2FC > 0) | (delta_ratio < 0 & Transcriptome_Log2FC < 0) ~ "Positive (Synergistic)",
          (delta_ratio > 0 & Transcriptome_Log2FC < 0) | (delta_ratio < 0 & Transcriptome_Log2FC > 0) ~ "Negative (Antagonistic)",
          TRUE ~ "Neutral"
        )
      )

    if (nrow(correlation_data) < 5) {
      cat("  ⚠️ 可用于 editing-mRNA 关联的点过少，跳过。\n")
      return(NULL)
    }

    plot_editing_expression_corr <- function(data, edit_type, save_path) {
      sub_data <- data %>% dplyr::filter(Editing_Type == edit_type)
      if (nrow(sub_data) < 5) {
        cat(paste0("  ⚠️ ", edit_type, " 有效点过少，跳过绘图。\n"))
        return(NULL)
      }

      cor_res <- suppressWarnings(cor.test(sub_data$delta_ratio, sub_data$Transcriptome_Log2FC, method = "spearman"))
      p_val_label <- if (cor_res$p.value < 0.001) "p < 0.001" else paste0("p = ", format(cor_res$p.value, digits = 3))

      label_df <- sub_data[order(abs(sub_data$delta_ratio * sub_data$Transcriptome_Log2FC), decreasing = TRUE), , drop = FALSE]
      label_df <- head(label_df, min(10, nrow(label_df)))

      p <- ggplot(sub_data, aes(x = delta_ratio, y = Transcriptome_Log2FC)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
        geom_point(aes(color = Relationship), alpha = 0.7, size = 2.5) +
        geom_smooth(method = "lm", color = "black", linetype = "solid", se = TRUE, alpha = 0.1, linewidth = 0.8) +
        scale_color_manual(values = c(
          "Positive (Synergistic)" = "#D88782",
          "Negative (Antagonistic)" = "#8EB9CB",
          "Neutral" = "gray70"
        )) +
        geom_text_repel(
          data = label_df,
          aes(label = gene_name),
          size = 3.5,
          fontface = "italic",
          box.padding = 0.5,
          segment.color = "grey50",
          max.overlaps = 20
        ) +
        theme_bw() +
        labs(
          title = paste(edit_type, "Delta vs mRNA Log2FC"),
          subtitle = paste0("Spearman R = ", round(cor_res$estimate, 3), ", ", p_val_label),
          x = "Editing Efficiency Delta (KO - WT)",
          y = "mRNA Log2FoldChange (KO / WT)"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
          legend.position = "bottom",
          panel.grid.minor = element_blank()
        )

      ggsave(save_path, plot = p, width = 8, height = 7.5)
      cat("  ✅ 已生成: ", save_path, "\n", sep = "")
      return(p)
    }

    plot_relationship_ratio <- function(data, save_path) {
      ratio_df <- data %>%
        dplyr::group_by(Editing_Type, Relationship) %>%
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(Editing_Type) %>%
        dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
        dplyr::ungroup()

      p <- ggplot(ratio_df, aes(x = Editing_Type, y = Percentage, fill = Relationship)) +
        geom_bar(stat = "identity", position = "stack", width = 0.5, color = "white", linewidth = 0.3) +
        geom_text(
          aes(label = paste0(round(Percentage, 1), "%")),
          position = position_stack(vjust = 0.5),
          color = "white",
          fontface = "bold",
          size = 4.5
        ) +
        scale_fill_manual(values = c(
          "Positive (Synergistic)" = "#D88782",
          "Negative (Antagonistic)" = "#8EB9CB",
          "Neutral" = "gray70"
        )) +
        theme_minimal() +
        labs(
          title = "Directional Impact of Editing on Gene Expression",
          x = "Editing Type",
          y = "Percentage of Significant Target Genes (%)"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text = element_text(size = 12, face = "bold"),
          legend.position = "right"
        )

      ggsave(save_path, plot = p, width = 7.5, height = 6)
      cat("  ✅ 已生成: ", save_path, "\n", sep = "")
    }

    plot_editing_expression_corr(correlation_data, "C_to_U", file.path(cross_plot_dir, "C2U_vs_mRNA_Quadrant_Plot.pdf"))
    plot_editing_expression_corr(correlation_data, "A_to_I", file.path(cross_plot_dir, "A2I_vs_mRNA_Quadrant_Plot.pdf"))
    plot_relationship_ratio(correlation_data, file.path(cross_plot_dir, "Editing_mRNA_Relationship_Barplot.pdf"))
  }, error = function(e) {
    cat(paste0("  ⚠️ editing-mRNA 后处理失败: ", e$message, "\n"))
  })
}

generate_final_dimension_report <- function(all_sites, de_status_df, log2fc_df, cross_plot_dir, parent_output_dir) {
  tryCatch({
    cat("\n▶ 后处理：生成 5 大维度真实靶点汇总表...\n")

    if (is.null(all_sites) || nrow(all_sites) == 0 || is.null(de_status_df) || nrow(de_status_df) == 0) {
      cat("  ⚠️ all_sites 或 de_status_df 为空，跳过 Final_True_Dimensions_Targets_Summary。\n")
      return(NULL)
    }

    dim_csv_files <- list.files(cross_plot_dir, pattern = "^[1-5]_.*\\.csv$", full.names = TRUE)
    if (length(dim_csv_files) == 0) {
      cat("  ⚠️ 未找到 5 大维度 CSV，跳过 Final_True_Dimensions_Targets_Summary。\n")
      return(NULL)
    }

    dim_pathway_list <- list()
    for (f in dim_csv_files) {
      df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0 || !"Description" %in% names(df)) next

      file_base <- basename(f)
      dim_prefix <- stringr::str_extract(file_base, "^[1-5]_[A-Za-z0-9_]+(?=_(GO_BP|KEGG)\\.csv$)")
      if (is.na(dim_prefix) || dim_prefix == "") next
      db_type <- ifelse(grepl("GO_BP", file_base), "GO_BP", "KEGG")
      gene_col <- ifelse("geneID" %in% names(df), "geneID", ifelse("core_enrichment" %in% names(df), "core_enrichment", NA))
      if (is.na(gene_col)) next

      for (i in seq_len(nrow(df))) {
        raw_genes <- as.character(df[[gene_col]][i])
        if (is.na(raw_genes) || raw_genes == "") next
        raw_ids <- stringr::str_split(raw_genes, "/")[[1]]
        raw_ids <- raw_ids[raw_ids != "" & !is.na(raw_ids)]

        symbols <- raw_ids
        if (length(raw_ids) > 0 && all(grepl("^[0-9]+$", raw_ids))) {
          symbols <- tryCatch({
            id_map <- suppressMessages(bitr(raw_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop = TRUE))
            unique(id_map$SYMBOL)
          }, error = function(e) raw_ids)
        }
        symbols <- symbols[!is.na(symbols) & symbols != ""]
        if (length(symbols) == 0) next

        pathway_label <- if ("ID" %in% names(df)) {
          paste0(df$Description[i], " (", df$ID[i], ")")
        } else {
          as.character(df$Description[i])
        }

        for (sym in symbols) {
          dim_pathway_list[[length(dim_pathway_list) + 1]] <- data.frame(
            Gene_Symbol = sym,
            Analysis_Dimension = dim_prefix,
            Database = db_type,
            Pathway = pathway_label,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(dim_pathway_list) > 0) {
      true_dim_pathway_df <- dplyr::bind_rows(dim_pathway_list) %>%
        dplyr::group_by(Gene_Symbol, Analysis_Dimension) %>%
        dplyr::summarise(Enriched_Pathways = paste(unique(Pathway), collapse = " | "), .groups = "drop")
    } else {
      true_dim_pathway_df <- data.frame(Gene_Symbol = character(), Analysis_Dimension = character(), Enriched_Pathways = character())
    }

    sig_edit_genes <- all_sites %>%
      dplyr::filter(significance %in% c("Up", "Down")) %>%
      dplyr::pull(gene_name) %>%
      unique()
    sig_mrna_genes <- de_status_df %>%
      dplyr::filter(mRNA_DE_Status %in% c("significantly_upregulated", "significantly_downregulated")) %>%
      dplyr::pull(Gene_Symbol) %>%
      unique()
    all_cross_genes <- intersect(sig_edit_genes, sig_mrna_genes)

    if (length(all_cross_genes) == 0) {
      cat("  ⚠️ 没有 editing 和 mRNA 双重显著基因，跳过最终维度表。\n")
      return(NULL)
    }

    all_sites_annotated <- all_sites %>%
      dplyr::left_join(de_status_df, by = c("gene_name" = "Gene_Symbol")) %>%
      dplyr::mutate(mRNA_DE_Status = ifelse(is.na(mRNA_DE_Status), "Not_Significant", mRNA_DE_Status))

    final_dimension_table <- all_sites_annotated %>%
      dplyr::filter(gene_name %in% all_cross_genes) %>%
      dplyr::left_join(true_dim_pathway_df, by = c("gene_name" = "Gene_Symbol"), relationship = "many-to-many") %>%
      dplyr::left_join(log2fc_df, by = c("gene_name" = "Gene_Symbol")) %>%
      dplyr::mutate(
        Analysis_Dimension = ifelse(is.na(Analysis_Dimension), "Unclassified_Dimension", Analysis_Dimension),
        Enriched_Pathways = ifelse(is.na(Enriched_Pathways), "No Significant Pathway", Enriched_Pathways),
        Editing_Trend = dplyr::case_when(
          significance == "Up" ~ "Increased in KO",
          significance == "Down" ~ "Decreased in KO",
          TRUE ~ "Not Significant"
        ),
        Transcriptome_Trend = dplyr::case_when(
          mRNA_DE_Status == "significantly_upregulated" ~ "mRNA Upregulated",
          mRNA_DE_Status == "significantly_downregulated" ~ "mRNA Downregulated",
          TRUE ~ "Not Significant"
        ),
        Genomic_Position = paste0(seqnames, ":", start)
      ) %>%
      dplyr::select(
        Gene_Symbol = gene_name,
        Analysis_Dimension,
        Enriched_Pathways,
        Editing_Type,
        Genomic_Position,
        WT_Efficiency = mean_wt,
        KO_Efficiency = mean_ko,
        Delta_Ratio = delta_ratio,
        P_Value = final_p_value,
        Editing_Trend,
        Transcriptome_Trend,
        Transcriptome_Log2FC
      ) %>%
      dplyr::arrange(Analysis_Dimension == "Unclassified_Dimension", Analysis_Dimension, Gene_Symbol, Editing_Type)

    final_report_path <- file.path(parent_output_dir, "Final_True_Dimensions_Targets_Summary.xlsx")
    openxlsx::write.xlsx(final_dimension_table, file = final_report_path, asTable = TRUE, overwrite = TRUE)
    cat("  ✅ 已生成: ", final_report_path, "\n", sep = "")
  }, error = function(e) {
    cat(paste0("  ⚠️ 5 大维度最终报表生成失败: ", e$message, "\n"))
  })
}

run_postprocessing_modules <- function(base_out, parent_output_dir, score_label, all_sites, de_status_df, log2fc_df, cross_plot_dir) {
  cat(paste0("\n========== 后处理模块启动：Score ", score_label, " ==========" , "\n"))
  redraw_kegg_emap_from_csv(base_out)
  plot_editing_mrna_relationship(all_sites, log2fc_df, cross_plot_dir)
  generate_final_dimension_report(all_sites, de_status_df, log2fc_df, cross_plot_dir, parent_output_dir)
  cat(paste0("========== 后处理模块完成：Score ", score_label, " ==========" , "\n"))
}

# ===================================================================
# 4. 纯净 Limma 核心分析函数
# ===================================================================
perform_limma_editing_analysis <- function(filtered_data, gtf_data, analysis_type = "C_to_U", plot_dir = "plots", pooled_types_df = NULL, sample_info = NULL) {

  cat(paste0("  ➤ 正在执行多因素 Limma 计算: ", analysis_type, "...\n"))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  gtf_genes <- gtf_data[gtf_data$type %in% c("gene", "pseudogene", "lnc_RNA", "ncRNA_gene", "miRNA", "snRNA", "snoRNA", "scRNA", "rRNA", "tRNA")]
  if (!("Name" %in% names(mcols(gtf_genes)))) gtf_genes$Name <- gtf_genes$gene_name
  gtf_genes_gr <- GRanges(seqnames(gtf_genes), IRanges(start(gtf_genes), end(gtf_genes)), strand=strand(gtf_genes), gene_name=gtf_genes$Name, gene_type=gtf_genes$type)

  cds_gr <- if(length(gtf_data[gtf_data$type == "CDS"])>0) GRanges(seqnames(gtf_data[gtf_data$type == "CDS"]), IRanges(start(gtf_data[gtf_data$type == "CDS"]), end(gtf_data[gtf_data$type == "CDS"])), region_type = "CDS") else GRanges()
  utr_gr <- c(if(length(gtf_data[gtf_data$type == "five_prime_UTR"])>0) GRanges(gtf_data[gtf_data$type == "five_prime_UTR"], region_type="5UTR") else GRanges(), if(length(gtf_data[gtf_data$type == "three_prime_UTR"])>0) GRanges(gtf_data[gtf_data$type == "three_prime_UTR"], region_type="3UTR") else GRanges())
  if(length(utr_gr)>0) utr_gr$region_type <- "UTR"
  region_gr <- c(cds_gr, utr_gr)

  if (is.null(pooled_types_df)) {
    # Python/JACUSA2 输入顺序：cond1 = WT，cond2 = KO
    rna_wt <- Reduce("+", filtered_data$bases$cond1)
    rna_ko <- Reduce("+", filtered_data$bases$cond2)
    pooled_types_df <- data.frame(
      query_idx = 1:length(filtered_data),
      ref2ko_type = base_sub(rna_ko, mcols(filtered_data)$ref),
      ref2wt_type = base_sub(rna_wt, mcols(filtered_data)$ref)
    )
  }

  cond1_ratios <- mutate(filtered_data$bases$cond1, across(everything(), calculate_ratios))
  cond2_ratios <- mutate(filtered_data$bases$cond2, across(everything(), calculate_ratios))
  all_ratios_long <- bind_rows(unnest_ratios(cond1_ratios, "cond1"), unnest_ratios(cond2_ratios, "cond2"))
  ref_df <- tibble(query_idx = 1:length(filtered_data), ref = mcols(filtered_data)$ref)
  all_ratios_long <- left_join(all_ratios_long, ref_df, by = "query_idx")

  if (analysis_type == "C_to_U") analysis_data <- all_ratios_long %>% filter(ref %in% c("C", "G")) %>% mutate(editing_ratio = ifelse(ref == "C", T_ratio, A_ratio)) else analysis_data <- all_ratios_long %>% filter(ref %in% c("A", "T")) %>% mutate(editing_ratio = ifelse(ref == "A", G_ratio, C_ratio))
  if (nrow(analysis_data) == 0) return(NULL)

  # === 核心多因素 Limma 建模 ===
  wide_ratios <- analysis_data %>% mutate(sample_id = paste(condition, replicate, sep = "_")) %>% dplyr::select(query_idx, sample_id, editing_ratio) %>% tidyr::pivot_wider(names_from = sample_id, values_from = editing_ratio)
  ratio_mat <- as.matrix(wide_ratios[,-1])
  rownames(ratio_mat) <- wide_ratios$query_idx

  keep <- apply(ratio_mat, 1, var, na.rm = TRUE) > 0; keep[is.na(keep)] <- FALSE
  if (sum(keep) == 0) return(NULL)
  ratio_mat_filtered <- ratio_mat[keep, , drop = FALSE]

  sample_ids <- colnames(ratio_mat_filtered)
  if (!is.null(sample_info) && "Condition" %in% colnames(sample_info)) {
    matched_idx <- match(sample_ids, sample_info$Sample_ID)
    if (any(is.na(matched_idx))) {
      stop("❌ ratio_mat 中存在无法匹配到 sample_metadata 的样本列：", paste(sample_ids[is.na(matched_idx)], collapse = ", "))
    }
    groups <- factor(sample_info$Condition[matched_idx], levels = c("WT", "KO"))
  } else {
    # 兜底逻辑仍遵守新版 Python/JACUSA2 输入顺序：cond1 = WT，cond2 = KO
    groups <- factor(ifelse(grepl("^cond1_", sample_ids), "WT", "KO"), levels = c("WT", "KO"))
  }

  if (!is.null(sample_info) && "Batch" %in% colnames(sample_info) && length(unique(sample_info$Batch)) > 1) {
    matched_batches <- sample_info$Batch[match(sample_ids, sample_info$Sample_ID)]
    if (any(is.na(matched_batches))) stop("❌ Batch 信息匹配失败，请检查 sample_metadata。")
    design <- model.matrix(~ factor(matched_batches) + groups)
  } else {
    design <- model.matrix(~ groups)
  }

  fit <- lmFit(ratio_mat_filtered, design)
  fit <- tryCatch({ eBayes(fit, trend = (nrow(ratio_mat_filtered) > 10)) }, error = function(e) { eBayes(fit, trend = FALSE) })

  target_coef <- which(colnames(design) == "groupsKO")
  if (length(target_coef) != 1) target_coef <- ncol(design)
  limma_res <- topTable(fit, coef = target_coef, number = Inf, sort.by = "none")

  stats_df <- tibble(
    query_idx = as.integer(rownames(limma_res)),
    mean_wt = rowMeans(ratio_mat_filtered[, groups == "WT", drop=FALSE], na.rm = TRUE),
    mean_ko = rowMeans(ratio_mat_filtered[, groups == "KO", drop=FALSE], na.rm = TRUE),
    final_p_value = limma_res$P.Value,
    final_fdr = limma_res$adj.P.Val
  ) %>% mutate(delta_ratio = mean_ko - mean_wt, final_fold_change = (mean_ko + 0.001) / (mean_wt + 0.001))
  # ===========================

  coords_df <- data.frame(query_idx = 1:length(filtered_data), seqnames = as.character(seqnames(filtered_data)), start = start(filtered_data), end = end(filtered_data))
  final_df <- left_join(stats_df, coords_df, by = "query_idx")
  query_gr <- GRanges(final_df$seqnames, IRanges(final_df$start, final_df$end))
  suppressWarnings(seqlevelsStyle(query_gr) <- seqlevelsStyle(gtf_genes_gr)[1])

  if (length(region_gr) > 0) {
    suppressWarnings(seqlevelsStyle(region_gr) <- seqlevelsStyle(query_gr))
    ov_reg <- findOverlaps(query_gr, region_gr)
    if (length(ov_reg) > 0) {
      ov_df <- data.frame(query_idx = final_df$query_idx[queryHits(ov_reg)], region_type = region_gr$region_type[subjectHits(ov_reg)]) %>% group_by(query_idx) %>% summarise(region_type = ifelse("CDS" %in% region_type, "CDS", ifelse("UTR" %in% region_type, "UTR", first(region_type)))) %>% ungroup()
      final_df <- final_df %>% left_join(ov_df, by = "query_idx")
    } else final_df$region_type <- NA
  } else final_df$region_type <- NA
  final_df$region_type <- ifelse(is.na(final_df$region_type), "other", final_df$region_type)

  ov_gene <- findOverlaps(query_gr, gtf_genes_gr)
  overlap_df <- data.frame(query_idx = final_df$query_idx[queryHits(ov_gene)], gene_name = gtf_genes_gr$gene_name[subjectHits(ov_gene)], gene_type = gtf_genes_gr$gene_type[subjectHits(ov_gene)], gene_strand = as.character(strand(gtf_genes_gr)[subjectHits(ov_gene)]))
  overlap_df <- overlap_df[!grepl("^ENSG", overlap_df$gene_name), , drop = FALSE]

  merged_df <- final_df %>% left_join(overlap_df, by = "query_idx") %>% left_join(pooled_types_df[!duplicated(pooled_types_df$query_idx), ], by = "query_idx") %>% left_join(ref_df, by="query_idx")

  merged_df_clean <- merged_df %>% dplyr::filter(!is.na(gene_name) | ((ref2wt_type %in% c("C->T", "G->A", "A->G", "T->C")) & (ref2ko_type == "no change")))

  merged_df_plot <- merged_df_clean %>% filter(!is.na(final_fold_change) & !is.na(final_p_value)) %>%
    mutate(
      final_p_value = pmax(final_p_value, .Machine$double.xmin),
      neg_log_p_value = -log10(final_p_value),
      log_fold_change = log2((mean_ko + 0.001) / (mean_wt + 0.001)),
      significance = case_when(final_p_value < PVAL_CUTOFF & delta_ratio >= DELTA_CUTOFF ~ "Up", final_p_value < PVAL_CUTOFF & delta_ratio <= -DELTA_CUTOFF ~ "Down", TRUE ~ "Not Significant"),
      significance_fdr = case_when(final_fdr < PVAL_CUTOFF & delta_ratio >= DELTA_CUTOFF ~ "Up", final_fdr < PVAL_CUTOFF & delta_ratio <= -DELTA_CUTOFF ~ "Down", TRUE ~ "Not Significant")
    )

  merged_df_final <- merged_df_plot %>%
    mutate(Antisense_C_Flag = case_when((gene_strand == '-' & ref2wt_type == 'C->T') | (gene_strand == '+' & ref2wt_type == 'G->A') ~ "The antisense strand is C", TRUE ~ "")) %>%
    dplyr::select(any_of(c("seqnames", "start", "end", "mean_wt", "mean_ko", "delta_ratio", "final_fold_change", "final_p_value", "final_fdr", "log_fold_change", "ref", "gene_name", "gene_type", "gene_strand", "region_type", "ref2wt_type", "ref2ko_type", "Antisense_C_Flag", "significance", "significance_fdr")))

  merged_df_plot <- merged_df_plot %>% mutate(region_type = case_when(region_type %in% c("CDS", "UTR") ~ region_type, is.na(gene_name) ~ "Intergenic", !is.na(gene_name) & gene_type == "protein_coding" ~ "Intronic", TRUE ~ "ncRNA/Other_Exon"))
  region_stats_df <- as.data.frame(table(merged_df_plot$region_type, useNA = "no")); colnames(region_stats_df) <- c("Region_Type", "Count")
  region_stats_df <- subset(region_stats_df, !is.na(Region_Type) & as.character(Region_Type) != "NA" & Count > 0)
  if (sum(region_stats_df$Count) > 0) { region_stats_df$Percentage <- round(region_stats_df$Count / sum(region_stats_df$Count) * 100, 2) } else { region_stats_df$Percentage <- 0 }

  if (nrow(region_stats_df) > 0) {
    color_map <- c("CDS"="#D88782", "UTR"="#8EB9CB", "ncRNA/Other_Exon"="#80B4A2", "Intronic"="#8390A7", "Intergenic"="#E4B19F")
    p_bar <- ggplot(region_stats_df, aes(x = Region_Type, y = Percentage, fill = Region_Type)) + geom_bar(stat = "identity", width = 0.7, color="black") + geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, fontface = "bold") + scale_fill_manual(values = color_map) + theme_minimal() + theme(legend.position = "none")
    ggsave(file.path(plot_dir, paste0(analysis_type, "_region_distribution.pdf")), plot = p_bar, width = 8, height = 6)
  }

  if (nrow(merged_df_plot) > 0) {
    volcano_plot <- ggplot(merged_df_plot, aes(x = log_fold_change, y = neg_log_p_value, color = significance)) + geom_jitter(alpha=0.5, size=1.2) + scale_color_manual(values = c("Up"="#E41A1C", "Down"="#377EB8", "Not Significant"="#D3D3D3")) + theme_bw() + geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed")
    ggsave(file.path(plot_dir, paste0(analysis_type, "_Volcano_Plot.pdf")), plot = volcano_plot, width = 8, height = 7)

    p_scatter <- ggplot(merged_df_plot, aes(x = mean_wt, y = mean_ko, color = significance)) + geom_point(data = filter(merged_df_plot, significance == "Not Significant"), color = "gray80", alpha = 0.4) + geom_point(data = filter(merged_df_plot, significance != "Not Significant"), aes(color = significance), size = 2.5, alpha = 0.8) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") + scale_color_manual(values = c("Up"="#F29E2E", "Down"="#0072B2")) + theme_minimal() + theme(legend.position = "none")
    ggsave(file.path(plot_dir, paste0(analysis_type, "_Scatter_Plot.pdf")), plot = p_scatter, width = 10, height = 7)
  }

  return(list(final_data = merged_df_final, plot_df = merged_df_plot))
}

# ===================================================================
# 5. 主循环控制: Limma 纯血版 Pipeline
# ===================================================================
run_limma_pipeline <- function(wt_vs_ko, gtf_data, score_cutoff, parent_output_dir, base_out, score_label, sample_metadata) {

  overall_plot_dir <- file.path(base_out, "Overall_Plots")
  dir.create(overall_plot_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("  -> 执行 Majority Rule 容缺过滤：cond1/WT 中至少 %d 个 >=10X，cond2/KO 中至少 %d 个 >=10X...\n", WT_MIN_COVERED_REPS, KO_MIN_COVERED_REPS))
  filtered <- wt_vs_ko %>%
    dplyr::filter(score >= score_cutoff) %>%
    dplyr::filter(rowSums(as.matrix(cov$cond1) >= 10) >= WT_MIN_COVERED_REPS & rowSums(as.matrix(cov$cond2) >= 10) >= KO_MIN_COVERED_REPS) %>%
    dplyr::filter(robust(bases))

  cat(sprintf("  -> 过滤完成！存活有效位点数: %d\n", length(filtered)))
  if (length(filtered) == 0) return(NULL)

  plot_global_pca(filtered, sample_metadata, base_out, score_label)

  # Python/JACUSA2 输入顺序：cond1 = WT，cond2 = KO
  rna_wt <- Reduce("+", filtered$bases$cond1)
  rna_ko <- Reduce("+", filtered$bases$cond2)
  pooled_df <- data.frame(
    query_idx = 1:length(filtered),
    ref2ko_type = base_sub(rna_ko, mcols(filtered)$ref),
    ref2wt_type = base_sub(rna_wt, mcols(filtered)$ref)
  )

  # 补回 M013 的 overall base substitution 总览图
  plot_overall_base_substitution(pooled_df, overall_plot_dir, score_label)

  # 执行纯 Limma 分析
  c2u_limma <- perform_limma_editing_analysis(filtered, gtf_data, "C_to_U", file.path(base_out, "C_to_U_Plots"), pooled_df, sample_info = sample_metadata)
  a2i_limma <- perform_limma_editing_analysis(filtered, gtf_data, "A_to_I", file.path(base_out, "A_to_I_Plots"), pooled_df, sample_info = sample_metadata)

  cat("\n    [通路富集] 全部使用去批次纯净版 Limma 结果提取靶标...\n")
  limma_enrichment <- list(C_to_U = list(UP = list(), DOWN = list()), A_to_I = list(UP = list(), DOWN = list()))

  limma_enrichment$C_to_U$UP$GO <- run_enrichment_GO(unique(c2u_limma$plot_df$gene_name[c2u_limma$plot_df$significance == "Up"]), "UP", file.path(base_out, "C_to_U_Plots"), "C_to_U")
  limma_enrichment$C_to_U$UP$KEGG <- run_enrichment_KEGG(unique(c2u_limma$plot_df$gene_name[c2u_limma$plot_df$significance == "Up"]), "UP", file.path(base_out, "C_to_U_Plots"), "C_to_U")
  limma_enrichment$C_to_U$DOWN$GO <- run_enrichment_GO(unique(c2u_limma$plot_df$gene_name[c2u_limma$plot_df$significance == "Down"]), "DOWN", file.path(base_out, "C_to_U_Plots"), "C_to_U")
  limma_enrichment$C_to_U$DOWN$KEGG <- run_enrichment_KEGG(unique(c2u_limma$plot_df$gene_name[c2u_limma$plot_df$significance == "Down"]), "DOWN", file.path(base_out, "C_to_U_Plots"), "C_to_U")

  limma_enrichment$A_to_I$UP$GO <- run_enrichment_GO(unique(a2i_limma$plot_df$gene_name[a2i_limma$plot_df$significance == "Up"]), "UP", file.path(base_out, "A_to_I_Plots"), "A_to_I")
  limma_enrichment$A_to_I$UP$KEGG <- run_enrichment_KEGG(unique(a2i_limma$plot_df$gene_name[a2i_limma$plot_df$significance == "Up"]), "UP", file.path(base_out, "A_to_I_Plots"), "A_to_I")
  limma_enrichment$A_to_I$DOWN$GO <- run_enrichment_GO(unique(a2i_limma$plot_df$gene_name[a2i_limma$plot_df$significance == "Down"]), "DOWN", file.path(base_out, "A_to_I_Plots"), "A_to_I")
  limma_enrichment$A_to_I$DOWN$KEGG <- run_enrichment_KEGG(unique(a2i_limma$plot_df$gene_name[a2i_limma$plot_df$significance == "Down"]), "DOWN", file.path(base_out, "A_to_I_Plots"), "A_to_I")

  excel_out <- file.path(base_out, paste0("Integrated_Results_Score", score_label, "_PureLimma.xlsx"))
  if (!file.exists(excel_out)) {
    wb <- openxlsx::createWorkbook()
    if (!is.null(c2u_limma)) { addWorksheet(wb, "C_to_U_limma"); writeData(wb, "C_to_U_limma", c2u_limma$final_data) }
    if (!is.null(a2i_limma)) { addWorksheet(wb, "A_to_I_limma"); writeData(wb, "A_to_I_limma", a2i_limma$final_data) }
    saveWorkbook(wb, excel_out, overwrite = TRUE)
  }

  return(limma_enrichment)
}

# ===================================================================
# 6. 全局核心数据预加载 & 样本元数据建立
# ===================================================================
gtf_cache_path <- file.path(main_output_dir, "cached_gtf_data.rds")
cat("\n========== 步骤 1: 预加载 GTF 基因注释文件 ==========\n")
if (exists("gtf_data")) { cat("  -> 检测到 GTF 已存在\n")
} else if (file.exists(gtf_cache_path)) { gtf_data <- readRDS(gtf_cache_path)
} else { gtf_data <- import(gtf_file); saveRDS(gtf_data, gtf_cache_path) }

jacusa_cache_path <- file.path(main_output_dir, "cached_jacusa_data_no_THP1_WTcond1.rds")
cat("\n========== 步骤 2: 预加载 JACUSA2 原始分析数据 ==========\n")
if (!file.exists(input_file)) stop("❌ 找不到 JACUSA2 输入文件：", input_file)
if (!file.exists(sample_manifest_file)) stop("❌ 找不到 Python 生成的样本 manifest：", sample_manifest_file)

cache_is_fresh <- file.exists(jacusa_cache_path) && file.info(jacusa_cache_path)$mtime >= file.info(input_file)$mtime
if (exists("wt_vs_ko")) {
  cat("  -> 检测到 JACUSA2 数据已存在于当前 R 会话\n")
} else if (cache_is_fresh) {
  cat("  -> 读取当前 no_THP1/WTcond1 专用缓存\n")
  wt_vs_ko <- readRDS(jacusa_cache_path)
} else {
  cat("  -> 重新读取 JACUSA2 out 文件，并更新 no_THP1/WTcond1 专用缓存\n")
  wt_vs_ko <- read_result(input_file, nThread = 16)
  saveRDS(wt_vs_ko, jacusa_cache_path)
}

cat("\n========== 步骤 3: 预加载 RNA-seq 转录组 DE_Status ==========\n")
rnaseq_base_path <- env_path(
  "RNA_EDITING_RNASEQ_RESULT_ROOT",
  "results/bulk_rnaseq"
)
rnaseq_dirs <- list.dirs(rnaseq_base_path, recursive = FALSE)
ultimate_dirs <- rnaseq_dirs[
  grepl("all_Analysis_Ultimate|^combined$", basename(rnaseq_dirs))
]

de_status_df <- data.frame(Gene_Symbol = character(), mRNA_DE_Status = character(), stringsAsFactors = FALSE)
log2fc_df <- data.frame(Gene_Symbol = character(), Transcriptome_Log2FC = numeric(), stringsAsFactors = FALSE)

if (length(ultimate_dirs) > 0) {
  latest_rnaseq_dir <- ultimate_dirs[order(file.info(ultimate_dirs)$mtime, decreasing = TRUE)][1]
  up_file <- file.path(latest_rnaseq_dir, "upregulated_gene_list.txt")
  down_file <- file.path(latest_rnaseq_dir, "downregulated_gene_list.txt")
  up_gene_list <- if(file.exists(up_file)) trimws(readLines(up_file)) else character(0)
  down_gene_list <- if(file.exists(down_file)) trimws(readLines(down_file)) else character(0)
  deseq_res_file <- file.path(latest_rnaseq_dir, "DESeq2_full_results.csv")
  if (file.exists(deseq_res_file)) {
    deseq_res <- read.csv(deseq_res_file, stringsAsFactors = FALSE)
    if ("SYMBOL" %in% colnames(deseq_res)) {
      if (!length(up_gene_list) && all(c("padj", "log2FoldChange") %in% colnames(deseq_res))) {
        up_gene_list <- deseq_res$SYMBOL[
          !is.na(deseq_res$padj) & deseq_res$padj < 0.05 &
            deseq_res$log2FoldChange >= 1
        ]
        down_gene_list <- deseq_res$SYMBOL[
          !is.na(deseq_res$padj) & deseq_res$padj < 0.05 &
            deseq_res$log2FoldChange <= -1
        ]
      }
      if ("log2FoldChange" %in% colnames(deseq_res)) log2fc_df <- deseq_res %>% dplyr::select(Gene_Symbol = SYMBOL, Transcriptome_Log2FC = log2FoldChange) %>% distinct(Gene_Symbol, .keep_all = TRUE)
      de_status_df <- deseq_res %>% dplyr::select(Gene_Symbol = SYMBOL) %>% distinct() %>% mutate(mRNA_DE_Status = case_when(Gene_Symbol %in% up_gene_list ~ "significantly_upregulated", Gene_Symbol %in% down_gene_list ~ "significantly_downregulated", TRUE ~ "Not_Significant"))
    }
  }
}

# ---------------- 从 Python manifest 构建样本元数据：R 不再手写/猜测顺序 ----------------
build_metadata <- function(jacusa_data, manifest_file) {
  if (!file.exists(manifest_file)) stop("❌ manifest 文件不存在：", manifest_file)
  manifest <- read.delim(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("JACUSA_Condition", "JACUSA_Replicate", "Condition", "Sample_Name", "Batch", "BAM", "Order_In_JACUSA_Command")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) stop("❌ manifest 缺少列：", paste(missing_cols, collapse = ", "))

  if (any(grepl("THP_1", manifest$Sample_Name) | grepl("THP_1", manifest$BAM))) {
    stop("❌ manifest 中仍包含 THP_1 样本；请先重新运行新版 Python pipeline。")
  }

  manifest <- manifest %>%
    dplyr::mutate(
      Order_In_JACUSA_Command = as.integer(Order_In_JACUSA_Command),
      JACUSA_Condition = as.character(JACUSA_Condition),
      Condition = as.character(Condition)
    ) %>%
    dplyr::arrange(factor(JACUSA_Condition, levels = c("cond1", "cond2")), Order_In_JACUSA_Command)

  cond1_reps <- names(jacusa_data$bases$cond1)
  cond2_reps <- names(jacusa_data$bases$cond2)
  m1 <- manifest %>% dplyr::filter(JACUSA_Condition == "cond1") %>% dplyr::arrange(Order_In_JACUSA_Command)
  m2 <- manifest %>% dplyr::filter(JACUSA_Condition == "cond2") %>% dplyr::arrange(Order_In_JACUSA_Command)

  cat(sprintf("\n🔹 [Cond1 - WT组] 共 %d 个样本 | [Cond2 - KO组] 共 %d 个样本\n", length(cond1_reps), length(cond2_reps)))

  if (nrow(m1) != length(cond1_reps) || nrow(m2) != length(cond2_reps)) {
    stop(sprintf(
      "❌ manifest 与 JACUSA2 out 中的 rep 数不一致：manifest cond1=%d, R cond1=%d; manifest cond2=%d, R cond2=%d。请删除旧缓存并确认 input_file/manifest 同源。",
      nrow(m1), length(cond1_reps), nrow(m2), length(cond2_reps)
    ))
  }
  if (any(m1$Condition != "WT") || any(m2$Condition != "KO")) {
    stop("❌ manifest 条件方向错误：要求 cond1 全部为 WT，cond2 全部为 KO。")
  }
  if (length(cond1_reps) != EXPECTED_REPLICATES || length(cond2_reps) != EXPECTED_REPLICATES) {
    warning(sprintf(
      "当前不是预期的 %d WT vs %d KO：cond1/WT=%d, cond2/KO=%d。",
      EXPECTED_REPLICATES, EXPECTED_REPLICATES,
      length(cond1_reps), length(cond2_reps)
    ))
  }

  m1$JACUSA_Replicate_From_R <- cond1_reps
  m2$JACUSA_Replicate_From_R <- cond2_reps

  meta <- dplyr::bind_rows(m1, m2) %>%
    dplyr::mutate(
      Sample_ID = paste0(JACUSA_Condition, "_", JACUSA_Replicate_From_R),
      Real_Sample_Name = Sample_Name
    ) %>%
    dplyr::select(
      Sample_ID,
      Condition,
      Real_Sample_Name,
      Batch,
      BAM,
      JACUSA_Condition,
      JACUSA_Replicate = JACUSA_Replicate_From_R,
      Order_In_JACUSA_Command
    )

  print(meta)
  write.csv(meta, file.path(main_output_dir, "sample_metadata_from_manifest.csv"), row.names = FALSE)
  return(meta)
}
sample_metadata <- build_metadata(wt_vs_ko, sample_manifest_file)


# ===================================================================
# 7. 核心循环: 默认只运行论文采用的 Score 1
# ===================================================================
score_values <- as.numeric(strsplit(
  Sys.getenv("RNA_EDITING_SCORE_VALUES", "1"), ",", fixed = TRUE
)[[1]])
if (any(!is.finite(score_values))) stop("RNA_EDITING_SCORE_VALUES 必须是逗号分隔的数值。")
for (current_score in score_values) {
  score_label <- as.character(current_score)
  cat(paste0("\n🚀 开始处理: FIXED_SCORE = ", score_label, " 🚀\n"))

  parent_output_dir <- file.path(main_output_dir, paste0("Score_", score_label))
  base_out <- file.path(parent_output_dir, paste0("Final_Results_Score", score_label, "_PureLimma"))
  cross_plot_dir <- file.path(base_out, "Cross_Analysis_Plots")

  # 👑 修复点：补回自动创建文件夹的命令！没有它们 R 无法保存文件！
  dir.create(parent_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(cross_plot_dir, recursive = TRUE, showWarnings = FALSE)
  run_limma_pipeline(wt_vs_ko, gtf_data, current_score, parent_output_dir, base_out, score_label, sample_metadata)

  # === 整合并提取脂质靶标 ===
  all_csvs <- list.files(base_out, pattern = "(go|kegg)_limma\\.csv$", recursive = TRUE, full.names = TRUE)
  results_list <- list()
  if (length(all_csvs) > 0) {
    for (f in all_csvs) {
      df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error=function(e) NULL)
      if (is.null(df) || nrow(df) == 0 || !"Description" %in% names(df)) next
      db_match <- ifelse(grepl("go_", f), "GO", "KEGG")
      gene_col <- ifelse("geneID" %in% names(df), "geneID", ifelse("core_enrichment" %in% names(df), "core_enrichment", NA))
      if (is.na(gene_col)) next

      unique_ids <- unique(unlist(str_split(paste(df[[gene_col]], collapse = "/"), "/")))
      unique_ids <- unique_ids[unique_ids != ""]
      id_dict <- NULL
      if (any(grepl("^[0-9]+$", head(unique_ids, 10)))) {
        try({ id_map <- suppressMessages(bitr(unique_ids, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)); id_dict <- setNames(id_map$SYMBOL, id_map$ENTREZID) }, silent=TRUE)
      }

      for (i in 1:nrow(df)) {
        entrez_ids <- str_split(df[i, ][[gene_col]], "/")[[1]]
        symbols_str <- if(!is.null(id_dict)) paste(ifelse(is.na(id_dict[entrez_ids]), entrez_ids, id_dict[entrez_ids]), collapse=", ") else paste(entrez_ids, collapse=", ")
        results_list[[length(results_list) + 1]] <- data.frame(
          Algorithm = "limma", Score_Cutoff = score_label, Analysis_Type = str_extract(f, "(C_to_U|A_to_I)"),
          Direction = str_extract(f, "(UP|DOWN)"), Database = db_match, Pathway_Description = df[i, "Description"],
          P_adjust = df[i, "p.adjust"], Gene_Symbols = symbols_str, stringsAsFactors = FALSE
        )
      }
    }
  }

  site_db_file <- file.path(base_out, paste0("Integrated_Results_Score", score_label, "_PureLimma.xlsx"))
  all_sites <- data.frame()
  if (file.exists(site_db_file)) {
    c2u <- tryCatch(read.xlsx(site_db_file, sheet = "C_to_U_limma"), error=function(e) data.frame())
    a2i <- tryCatch(read.xlsx(site_db_file, sheet = "A_to_I_limma"), error=function(e) data.frame())
    if(nrow(c2u)>0) c2u$Editing_Type <- "C_to_U"
    if(nrow(a2i)>0) a2i$Editing_Type <- "A_to_I"
    all_sites <- bind_rows(c2u, a2i) %>% distinct(seqnames, start, end, ref, gene_name, Editing_Type, .keep_all = TRUE)
  }

  if (length(results_list) > 0) {
    final_results <- bind_rows(results_list) %>% arrange(Database, Analysis_Type, Direction, P_adjust)
    write.xlsx(final_results, file = file.path(parent_output_dir, paste0("ALL_Pathways_Gene_Summary_Limma_Score", score_label, ".xlsx")), asTable = TRUE, overwrite = TRUE)

    target_pathways <- final_results %>% filter(str_detect(tolower(Pathway_Description), "atherosclerosis|lipid"))
    if (nrow(target_pathways) > 0 && nrow(all_sites) > 0) {
      gene_pathway_dict <- target_pathways %>% separate_rows(Gene_Symbols, sep = ",\\s*") %>% rename(gene_name = Gene_Symbols) %>% group_by(gene_name) %>% summarise(Pathways = paste(unique(Pathway_Description), collapse = " ; "))
      final_site_report <- all_sites %>% inner_join(gene_pathway_dict, by = "gene_name") %>% mutate(Log2_FC = round(log_fold_change, 4)) %>% left_join(de_status_df, by = c("gene_name"="Gene_Symbol")) %>% arrange(gene_name, Editing_Type)
      write.xlsx(final_site_report, file = file.path(parent_output_dir, paste0("Target_Atherosclerosis_Lipid_SiteLevel_Report_Limma_Score", score_label, ".xlsx")), asTable = TRUE, overwrite = TRUE)
    }
  }

  # === 5大维度交叉分析 (纯净版) ===
  if(nrow(all_sites) > 0 && nrow(de_status_df) > 0) {
    all_sites_annotated <- all_sites %>% left_join(de_status_df, by = c("gene_name" = "Gene_Symbol"))
    sig_c2u_genes <- all_sites %>% filter(Editing_Type == "C_to_U", significance %in% c("Up", "Down")) %>% pull(gene_name) %>% unique()
    sig_a2i_genes <- all_sites %>% filter(Editing_Type == "A_to_I", significance %in% c("Up", "Down")) %>% pull(gene_name) %>% unique()
    sig_any_edit <- unique(c(sig_c2u_genes, sig_a2i_genes))
    up_mrna <- de_status_df %>% filter(mRNA_DE_Status == "significantly_upregulated") %>% pull(Gene_Symbol) %>% unique()
    down_mrna <- de_status_df %>% filter(mRNA_DE_Status == "significantly_downregulated") %>% pull(Gene_Symbol) %>% unique()
    sig_mrna <- unique(c(up_mrna, down_mrna))

    dims <- list(
      "1_C2U_Edit_and_DE_Sig" = intersect(sig_c2u_genes, sig_mrna),
      "2_A2I_Edit_and_DE_Sig" = intersect(sig_a2i_genes, sig_mrna),
      "3_Any_Edit_and_DE_Sig" = intersect(sig_any_edit, sig_mrna),
      "4_Any_Edit_and_DE_Upregulated" = intersect(sig_any_edit, up_mrna),
      "5_Any_Edit_and_DE_Downregulated" = intersect(sig_any_edit, down_mrna)
    )

    dim_pathways <- list()
    for (d_name in names(dims)) {
      g_list <- dims[[d_name]]
      if(length(g_list) < 5) next
      ids <- suppressMessages(bitr(g_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE))
      if (is.null(ids) || nrow(ids) == 0) next

      go <- tryCatch({ suppressMessages(enrichGO(ids$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05)) }, error=function(e) NULL)
      if(!is.null(go) && nrow(go)>0) { write.csv(as.data.frame(go), file.path(cross_plot_dir, paste0(d_name, "_GO_BP.csv")), row.names=F); ggsave(file.path(cross_plot_dir, paste0(d_name, "_GO_BP.pdf")), plot=dotplot(go, title=d_name), width=10, height=8) }

      kegg <- tryCatch({ suppressMessages(enricher(ids$ENTREZID, TERM2GENE=kegg_local$KEGGPATHID2EXTID, TERM2NAME=kegg_local$KEGGPATHID2NAME, pvalueCutoff=0.05)) }, error=function(e) NULL)
      if(!is.null(kegg) && nrow(kegg)>0) { write.csv(as.data.frame(kegg), file.path(cross_plot_dir, paste0(d_name, "_KEGG.csv")), row.names=F); ggsave(file.path(cross_plot_dir, paste0(d_name, "_KEGG.pdf")), plot=dotplot(kegg, title=d_name), width=10, height=8) }
    }
  }

  # === M013 后处理模块迁移版：KEGG Emap重绘 + editing/mRNA关系图 + 5大维度终表 ===
  run_postprocessing_modules(base_out, parent_output_dir, score_label, all_sites, de_status_df, log2fc_df, cross_plot_dir)
}

cat("\n🎉🎉🎉 纯净版 Limma 多因素分析全流程执行完毕！ 🎉🎉🎉\n")




# ===================================================================
# 终端核对看板：自动读取所有 Score 档位的结果并打印 (UP/DOWN 增强版)
# ===================================================================
cat("\n🚀 正在启动终极核对看板 (区分 UP/DOWN)...\n")

# 定义打印函数
print_top_pathways_v2 <- function(base_dir) {
  # 寻找所有生成的 limma.csv 文件
  csv_files <- list.files(base_dir, pattern = "limma\\.csv$", recursive = TRUE, full.names = TRUE)

  if (length(csv_files) == 0) {
    cat("  ⚠️ 警告：当前目录下未找到 limma.csv 结果文件。\n")
    return()
  }

  for (f in csv_files) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next

    # --- 核心修改：提取路径和文件名信息 ---
    path_info <- basename(dirname(f))
    file_name <- basename(f)

    # 识别富集类型 (GO/KEGG)
    db_type <- ifelse(grepl("go_", f), "GO BP", "KEGG")

    # 识别上下调趋势 (从文件名提取 UP 或 DOWN)
    trend <- ifelse(grepl("_UP_", file_name), "⬆️ UP (编辑增加)",
                    ifelse(grepl("_DOWN_", file_name), "⬇️ DOWN (编辑减少)", "Unknown"))

    # 打印格式化抬头
    cat(sprintf("\n[ %s | %s | %s ]\n", path_info, db_type, trend))
    cat(sprintf("---------------------------------------------------------\n"))

    # 排序并提取 Top 15
    df <- df[order(as.numeric(df$p.adjust)), ]
    top_pathways <- head(df$Description, 15)

    for (i in seq_along(top_pathways)) {
      cat(sprintf("%d. %s\n", i, top_pathways[i]))
    }
  }
}
print_cross_dimension_pathways <- function(base_out) {
  cross_dir <- file.path(base_out, "Cross_Analysis_Plots")
  if (!dir.exists(cross_dir)) return()

  csv_files <- list.files(cross_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(csv_files) > 0) {
    cat("\n✨ [五大维度交叉富集看板]\n")
    for (f in csv_files) {
      df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next

      dim_name <- basename(f)
      cat(sprintf("  ▶ 维度: %s\n", dim_name))

      df <- df[order(as.numeric(df$p.adjust)), ]
      top_n <- head(df$Description, 10) # 维度分析通常通路较少，打前10个

      for (i in seq_along(top_n)) {
        cat(sprintf("    %d. %s\n", i, top_n[i]))
      }
    }
  }
}
# 遍历所有 Score 文件夹进行打印
for (current_score in c(1.0, 0.5, 0)) {
  score_dir <- file.path(main_output_dir, paste0("Score_", current_score), paste0("Final_Results_Score", current_score, "_PureLimma"))

  if (dir.exists(score_dir)) {
    cat(sprintf("\n\n#########################################################\n"))
    cat(sprintf("看板：Score档位 [%s] 的富集结果摘要\n", current_score))
    cat(sprintf("#########################################################\n"))
    print_top_pathways_v2(score_dir)
    # 打印 5 个交叉维度看板
    print_cross_dimension_pathways(score_dir)
  } else {
    cat(sprintf("\n⚠️ 目录不存在，跳过 Score档位 [%s]: %s\n", current_score, score_dir))
  }
}


# Score [1]
# A_to_I_Plots | GO BP | DOWN（编辑减少）
# lysosomal membrane
# → 溶酶体膜
# A_to_I_Plots | KEGG | DOWN（编辑减少）
# Terpenoid backbone biosynthesis
# → 萜类骨架生物合成
# Influenza A
# → 甲型流感
# Measles
# → 麻疹
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# NOD-like receptor signaling pathway
# → NOD 样受体信号通路
# Coronavirus disease - COVID-19
# → 冠状病毒病 COVID-19
# Virion - Ebolavirus, Lyssavirus and Morbillivirus
# → 病毒颗粒：埃博拉病毒、狂犬病毒属和麻疹病毒属
# Bacterial invasion of epithelial cells
# → 细菌侵入上皮细胞
# Leishmaniasis
# → 利什曼病
# Taurine and hypotaurine metabolism
# → 牛磺酸和亚牛磺酸代谢
# Fatty acid biosynthesis
# → 脂肪酸生物合成
# Endocrine and other factor-regulated calcium reabsorption
# → 内分泌及其他因子调控的钙重吸收
# Fatty acid metabolism
# → 脂肪酸代谢
# Huntington disease
# → 亨廷顿病
# Phospholipase D signaling pathway
# → 磷脂酶 D 信号通路
# A_to_I_Plots | GO BP | UP（编辑增加）
# ribonucleoprotein complex biogenesis
# → 核糖核蛋白复合体生成
# cytoplasmic pattern recognition receptor signaling pathway
# → 细胞质模式识别受体信号通路
# regulation of type I interferon production
# → I 型干扰素产生的调控
# type I interferon production
# → I 型干扰素产生
# positive regulation of innate immune response
# → 先天免疫反应的正调控
# positive regulation of response to biotic stimulus
# → 对生物性刺激反应的正调控
# regulation of apoptotic DNA fragmentation
# → 凋亡性 DNA 片段化的调控
# establishment of protein localization to organelle
# → 蛋白质定位至细胞器的建立
# regulation of viral life cycle
# → 病毒生命周期的调控
# protein folding
# → 蛋白质折叠
# regulation of cytoplasmic pattern recognition receptor signaling pathway
# → 细胞质模式识别受体信号通路的调控
# transcription initiation at RNA polymerase II promoter
# → RNA 聚合酶 II 启动子处的转录起始
# mRNA splicing, via spliceosome
# → 通过剪接体进行的 mRNA 剪接
# rRNA metabolic process
# → rRNA 代谢过程
# maturation of SSU-rRNA
# → 小亚基 rRNA 成熟
# A_to_I_Plots | KEGG | UP（编辑增加）
# Tuberculosis
# → 结核病
# Apoptosis
# → 细胞凋亡
# Lysosome biogenesis
# → 溶酶体生物发生
# Protein processing in endoplasmic reticulum
# → 内质网中的蛋白质加工
# RNA polymerase
# → RNA 聚合酶
# N-Glycan biosynthesis
# → N-糖链生物合成
# Fatty acid metabolism
# → 脂肪酸代谢
# Amyotrophic lateral sclerosis
# → 肌萎缩侧索硬化症
# Prion disease
# → 朊病毒病
# Terpenoid backbone biosynthesis
# → 萜类骨架生物合成
# Various types of N-glycan biosynthesis
# → 多种类型的 N-糖链生物合成
# Parkinson disease
# → 帕金森病
# Huntington disease
# → 亨廷顿病
# NOD-like receptor signaling pathway
# → NOD 样受体信号通路
# Shigellosis
# → 志贺菌病
# C_to_U_Plots | GO BP | DOWN（编辑减少）
# phosphate ion transmembrane transport
# → 磷酸根离子跨膜转运
# phosphate transmembrane transporter activity
# → 磷酸盐跨膜转运蛋白活性
# solute:proton symporter activity
# → 溶质:质子同向转运体活性
# regulation of postsynapse organization
# → 突触后结构组织的调控
# phosphate ion transport
# → 磷酸根离子转运
# nuclear speck
# → 核斑
# early endosome
# → 早期内体
# lysosomal membrane
# → 溶酶体膜
# specific granule membrane
# → 特异性颗粒膜
# transition metal ion transmembrane transporter activity
# → 过渡金属离子跨膜转运蛋白活性
# C_to_U_Plots | KEGG | DOWN（编辑减少）
# Glycerophospholipid metabolism
# → 甘油磷脂代谢
# Biosynthesis of nucleotide sugars
# → 核苷酸糖生物合成
# Amino sugar and nucleotide sugar metabolism
# → 氨基糖和核苷酸糖代谢
# MAPK signaling pathway
# → MAPK 信号通路
# Lysosome biogenesis
# → 溶酶体生物发生
# Ovarian steroidogenesis
# → 卵巢类固醇生成
# Phosphatidylinositol signaling system
# → 磷脂酰肌醇信号系统
# Oxytocin signaling pathway
# → 催产素信号通路
# Sphingolipid metabolism
# → 鞘脂代谢
# VEGF signaling pathway
# → VEGF 信号通路
# Long-term depression
# → 长时程抑制
# Arachidonic acid metabolism
# → 花生四烯酸代谢
# Glutamatergic synapse
# → 谷氨酸能突触
# Long-term potentiation
# → 长时程增强
# alpha-Linolenic acid metabolism
# → α-亚麻酸代谢
# C_to_U_Plots | GO BP | UP（编辑增加）
# early endosome
# → 早期内体
# growth cone
# → 生长锥
# secretory granule membrane
# → 分泌颗粒膜
# nuclear membrane
# → 核膜
# site of polarized growth
# → 极性生长位点
# core mediator complex
# → 核心 Mediator 复合体
# C_to_U_Plots | KEGG | UP（编辑增加）
# Spinocerebellar ataxia
# → 脊髓小脑性共济失调
# Endocytosis
# → 内吞作用
# Fanconi anemia pathway
# → 范可尼贫血通路
# Oxytocin signaling pathway
# → 催产素信号通路
# Amoebiasis
# → 阿米巴病
# Human cytomegalovirus infection
# → 人巨细胞病毒感染
# Cholinergic synapse
# → 胆碱能突触
# Parathyroid hormone synthesis, secretion and action
# → 甲状旁腺激素的合成、分泌和作用
# Glutamatergic synapse
# → 谷氨酸能突触
# Folate biosynthesis
# → 叶酸生物合成
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Leishmaniasis
# → 利什曼病
# Circadian rhythm
# → 昼夜节律
# Vascular smooth muscle contraction
# → 血管平滑肌收缩
# Apoptosis
# → 细胞凋亡
# Score [1] 五大维度交叉富集看板
# 1_C2U_Edit_and_DE_Sig_GO_BP.csv
# regulation of apoptotic cell clearance
# → 凋亡细胞清除的调控
# complement activation, GZMK pathway
# → 补体激活，GZMK 通路
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# complement activation, lectin pathway
# → 补体激活，凝集素通路
# regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应调控
# activation of membrane attack complex
# → 膜攻击复合体激活
# MHC class II protein complex assembly
# → MHC II 类蛋白复合体组装
# peptide antigen assembly with MHC class II protein complex
# → 肽抗原与 MHC II 类蛋白复合体组装
# detection of bacterium
# → 细菌检测
# fibroblast activation
# → 成纤维细胞激活
# 1_C2U_Edit_and_DE_Sig_KEGG.csv
# Antigen processing and presentation
# → 抗原加工与呈递
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Phagosome
# → 吞噬体
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Viral myocarditis
# → 病毒性心肌炎
# 2_A2I_Edit_and_DE_Sig_GO_BP.csv
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽或多糖抗原加工与呈递
# monocyte differentiation
# → 单核细胞分化
# antigen processing and presentation of exogenous peptide antigen
# → 外源性肽抗原加工与呈递
# antigen processing and presentation of exogenous antigen
# → 外源性抗原加工与呈递
# positive regulation of monocyte differentiation
# → 单核细胞分化的正调控
# positive regulation of chemokine (C-X-C motif) ligand 2 production
# → 趋化因子 C-X-C 基序配体 2 产生的正调控
# interleukin-10 production
# → 白细胞介素-10 产生
# regulation of interleukin-10 production
# → 白细胞介素-10 产生的调控
# 2_A2I_Edit_and_DE_Sig_KEGG.csv
# Tuberculosis
# → 结核病
# Antigen processing and presentation
# → 抗原加工与呈递
# Asthma
# → 哮喘
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Intestinal immune network for IgA production
# → IgA 产生相关的肠道免疫网络
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# N-Glycan biosynthesis
# → N-糖链生物合成
# 3_Any_Edit_and_DE_Sig_GO_BP.csv
# regulation of B cell proliferation
# → B 细胞增殖的调控
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# regulation of lymphocyte proliferation
# → 淋巴细胞增殖的调控
# regulation of mononuclear cell proliferation
# → 单核细胞增殖的调控
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of endogenous antigen
# → 内源性抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应
# B cell proliferation
# → B 细胞增殖
# regulation of leukocyte proliferation
# → 白细胞增殖的调控
# 3_Any_Edit_and_DE_Sig_KEGG.csv
# Antigen processing and presentation
# → 抗原加工与呈递
# Tuberculosis
# → 结核病
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Alcoholic liver disease
# → 酒精性肝病
# 4_Any_Edit_and_DE_Upregulated_GO_BP.csv
# regulation of B cell proliferation
# → B 细胞增殖的调控
# B cell proliferation
# → B 细胞增殖
# negative regulation of B cell proliferation
# → B 细胞增殖的负调控
# cristae formation
# → 嵴形成
# MyD88-dependent toll-like receptor signaling pathway
# → MyD88 依赖性 Toll 样受体信号通路
# regulation of B cell activation
# → B 细胞活化的调控
# immune response-activating cell surface receptor signaling pathway
# → 激活免疫反应的细胞表面受体信号通路
# T cell receptor signaling pathway
# → T 细胞受体信号通路
# adrenergic receptor signaling pathway
# → 肾上腺素能受体信号通路
# negative regulation of B cell activation
# → B 细胞活化的负调控
# 5_Any_Edit_and_DE_Downregulated_GO_BP.csv
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of endogenous antigen
# → 内源性抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽或多糖抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen
# → 外源性肽抗原加工与呈递
# lymphocyte mediated immunity
# → 淋巴细胞介导的免疫
# antigen processing and presentation
# → 抗原加工与呈递
# adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应
# antigen processing and presentation of exogenous antigen
# → 外源性抗原加工与呈递
# 5_Any_Edit_and_DE_Downregulated_KEGG.csv
# Antigen processing and presentation
# → 抗原加工与呈递
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Tuberculosis
# → 结核病
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Phagosome
# → 吞噬体
# Score [0.5]
# A_to_I_Plots | GO BP | DOWN（编辑减少）
# phosphate ion transport
# → 磷酸根离子转运
# phosphate ion transmembrane transport
# → 磷酸根离子跨膜转运
# lysosomal membrane
# → 溶酶体膜
# A_to_I_Plots | KEGG | DOWN（编辑减少）
# Terpenoid backbone biosynthesis
# → 萜类骨架生物合成
# Influenza A
# → 甲型流感
# Measles
# → 麻疹
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Virion - Ebolavirus, Lyssavirus and Morbillivirus
# → 病毒颗粒：埃博拉病毒、狂犬病毒属和麻疹病毒属
# NOD-like receptor signaling pathway
# → NOD 样受体信号通路
# Coronavirus disease - COVID-19
# → 冠状病毒病 COVID-19
# Huntington disease
# → 亨廷顿病
# Bacterial invasion of epithelial cells
# → 细菌侵入上皮细胞
# Leishmaniasis
# → 利什曼病
# Taurine and hypotaurine metabolism
# → 牛磺酸和亚牛磺酸代谢
# Fatty acid biosynthesis
# → 脂肪酸生物合成
# Endocrine and other factor-regulated calcium reabsorption
# → 内分泌及其他因子调控的钙重吸收
# Fatty acid metabolism
# → 脂肪酸代谢
# Phospholipase D signaling pathway
# → 磷脂酶 D 信号通路
# A_to_I_Plots | GO BP | UP（编辑增加）
# ribonucleoprotein complex biogenesis
# → 核糖核蛋白复合体生成
# regulation of viral life cycle
# → 病毒生命周期的调控
# regulation of type I interferon production
# → I 型干扰素产生的调控
# type I interferon production
# → I 型干扰素产生
# regulation of apoptotic DNA fragmentation
# → 凋亡性 DNA 片段化的调控
# cytoplasmic pattern recognition receptor signaling pathway
# → 细胞质模式识别受体信号通路
# mRNA splicing, via spliceosome
# → 通过剪接体进行的 mRNA 剪接
# protein folding
# → 蛋白质折叠
# positive regulation of innate immune response
# → 先天免疫反应的正调控
# positive regulation of response to biotic stimulus
# → 对生物性刺激反应的正调控
# transcription initiation at RNA polymerase II promoter
# → RNA 聚合酶 II 启动子处的转录起始
# establishment of protein localization to organelle
# → 蛋白质定位至细胞器的建立
# maturation of SSU-rRNA
# → 小亚基 rRNA 成熟
# tRNA processing
# → tRNA 加工
# mitochondrion organization
# → 线粒体组织
# A_to_I_Plots | KEGG | UP（编辑增加）
# Tuberculosis
# → 结核病
# Apoptosis
# → 细胞凋亡
# Lysosome biogenesis
# → 溶酶体生物发生
# Protein processing in endoplasmic reticulum
# → 内质网中的蛋白质加工
# N-Glycan biosynthesis
# → N-糖链生物合成
# Fatty acid metabolism
# → 脂肪酸代谢
# Prion disease
# → 朊病毒病
# Amyotrophic lateral sclerosis
# → 肌萎缩侧索硬化症
# Terpenoid backbone biosynthesis
# → 萜类骨架生物合成
# Shigellosis
# → 志贺菌病
# Various types of N-glycan biosynthesis
# → 多种类型的 N-糖链生物合成
# Parkinson disease
# → 帕金森病
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Ubiquinone and other terpenoid-quinone biosynthesis
# → 泛醌及其他萜醌类生物合成
# Huntington disease
# → 亨廷顿病
# C_to_U_Plots | GO BP | DOWN（编辑减少）
# phosphate ion transmembrane transport
# → 磷酸根离子跨膜转运
# phosphate transmembrane transporter activity
# → 磷酸盐跨膜转运蛋白活性
# GTPase regulator activity
# → GTP 酶调节因子活性
# intracellularly gated calcium channel activity
# → 细胞内门控钙通道活性
# phosphotransferase activity, phosphate group as acceptor
# → 磷酸转移酶活性，以磷酸基团为受体
# small GTPase binding
# → 小 GTP 酶结合
# kinesin binding
# → 驱动蛋白结合
# transition metal ion transmembrane transporter activity
# → 过渡金属离子跨膜转运蛋白活性
# early endosome
# → 早期内体
# nuclear speck
# → 核斑
# lysosomal membrane
# → 溶酶体膜
# C_to_U_Plots | KEGG | DOWN（编辑减少）
# Phosphatidylinositol signaling system
# → 磷脂酰肌醇信号系统
# Long-term depression
# → 长时程抑制
# MAPK signaling pathway
# → MAPK 信号通路
# Glycerophospholipid metabolism
# → 甘油磷脂代谢
# Oxytocin signaling pathway
# → 催产素信号通路
# Long-term potentiation
# → 长时程增强
# Biosynthesis of nucleotide sugars
# → 核苷酸糖生物合成
# Glutamatergic synapse
# → 谷氨酸能突触
# Amino sugar and nucleotide sugar metabolism
# → 氨基糖和核苷酸糖代谢
# Platelet activation
# → 血小板活化
# Lysosome biogenesis
# → 溶酶体生物发生
# GnRH signaling pathway
# → 促性腺激素释放激素信号通路
# Ovarian steroidogenesis
# → 卵巢类固醇生成
# Inflammatory mediator regulation of TRP channels
# → 炎症介质对 TRP 通道的调控
# Sphingolipid metabolism
# → 鞘脂代谢
# C_to_U_Plots | GO BP | UP（编辑增加）
# early endosome
# → 早期内体
# growth cone
# → 生长锥
# site of polarized growth
# → 极性生长位点
# secretory granule membrane
# → 分泌颗粒膜
# C_to_U_Plots | KEGG | UP（编辑增加）
# Spinocerebellar ataxia
# → 脊髓小脑性共济失调
# Fanconi anemia pathway
# → 范可尼贫血通路
# Endocytosis
# → 内吞作用
# Amoebiasis
# → 阿米巴病
# Oxytocin signaling pathway
# → 催产素信号通路
# Folate biosynthesis
# → 叶酸生物合成
# Cholinergic synapse
# → 胆碱能突触
# Parathyroid hormone synthesis, secretion and action
# → 甲状旁腺激素的合成、分泌和作用
# Leukocyte transendothelial migration
# → 白细胞跨内皮迁移
# Glutamatergic synapse
# → 谷氨酸能突触
# Autophagy - animal
# → 自噬 - 动物
# Human cytomegalovirus infection
# → 人巨细胞病毒感染
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Circadian rhythm
# → 昼夜节律
# Leishmaniasis
# → 利什曼病
# Score [0.5] 五大维度交叉富集看板
# 1_C2U_Edit_and_DE_Sig_GO_BP.csv
# regulation of apoptotic cell clearance
# → 凋亡细胞清除的调控
# complement activation, GZMK pathway
# → 补体激活，GZMK 通路
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# complement activation, lectin pathway
# → 补体激活，凝集素通路
# activation of membrane attack complex
# → 膜攻击复合体激活
# MHC class II protein complex assembly
# → MHC II 类蛋白复合体组装
# peptide antigen assembly with MHC class II protein complex
# → 肽抗原与 MHC II 类蛋白复合体组装
# detection of bacterium
# → 细菌检测
# fibroblast activation
# → 成纤维细胞激活
# regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应调控
# 1_C2U_Edit_and_DE_Sig_KEGG.csv
# Antigen processing and presentation
# → 抗原加工与呈递
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Phagosome
# → 吞噬体
# Viral myocarditis
# → 病毒性心肌炎
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# 2_A2I_Edit_and_DE_Sig_GO_BP.csv
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽或多糖抗原加工与呈递
# monocyte differentiation
# → 单核细胞分化
# antigen processing and presentation of exogenous peptide antigen
# → 外源性肽抗原加工与呈递
# antigen processing and presentation of exogenous antigen
# → 外源性抗原加工与呈递
# positive regulation of monocyte differentiation
# → 单核细胞分化的正调控
# positive regulation of chemokine (C-X-C motif) ligand 2 production
# → 趋化因子 C-X-C 基序配体 2 产生的正调控
# interleukin-10 production
# → 白细胞介素-10 产生
# regulation of interleukin-10 production
# → 白细胞介素-10 产生的调控
# 2_A2I_Edit_and_DE_Sig_KEGG.csv
# Tuberculosis
# → 结核病
# Antigen processing and presentation
# → 抗原加工与呈递
# Asthma
# → 哮喘
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Intestinal immune network for IgA production
# → IgA 产生相关的肠道免疫网络
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# N-Glycan biosynthesis
# → N-糖链生物合成
# 3_Any_Edit_and_DE_Sig_GO_BP.csv
# regulation of B cell proliferation
# → B 细胞增殖的调控
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# regulation of lymphocyte proliferation
# → 淋巴细胞增殖的调控
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# regulation of mononuclear cell proliferation
# → 单核细胞增殖的调控
# antigen processing and presentation of endogenous antigen
# → 内源性抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应
# B cell proliferation
# → B 细胞增殖
# negative regulation of lymphocyte proliferation
# → 淋巴细胞增殖的负调控
# 3_Any_Edit_and_DE_Sig_KEGG.csv
# Antigen processing and presentation
# → 抗原加工与呈递
# Tuberculosis
# → 结核病
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Allograft rejection
# → 同种异体移植排斥
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Alcoholic liver disease
# → 酒精性肝病
# 4_Any_Edit_and_DE_Upregulated_GO_BP.csv
# regulation of B cell proliferation
# → B 细胞增殖的调控
# B cell proliferation
# → B 细胞增殖
# negative regulation of B cell proliferation
# → B 细胞增殖的负调控
# cristae formation
# → 嵴形成
# MyD88-dependent toll-like receptor signaling pathway
# → MyD88 依赖性 Toll 样受体信号通路
# 5_Any_Edit_and_DE_Downregulated_GO_BP.csv
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of endogenous antigen
# → 内源性抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽或多糖抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen
# → 外源性肽抗原加工与呈递
# lymphocyte mediated immunity
# → 淋巴细胞介导的免疫
# antigen processing and presentation
# → 抗原加工与呈递
# adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应
# antigen processing and presentation of exogenous antigen
# → 外源性抗原加工与呈递
# 5_Any_Edit_and_DE_Downregulated_KEGG.csv
# Antigen processing and presentation
# → 抗原加工与呈递
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Tuberculosis
# → 结核病
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Phagosome
# → 吞噬体
# Score [0]
# A_to_I_Plots | GO BP | DOWN（编辑减少）
# lysosomal membrane
# → 溶酶体膜
# phosphate ion transport
# → 磷酸根离子转运
# phosphate ion transmembrane transport
# → 磷酸根离子跨膜转运
# A_to_I_Plots | KEGG | DOWN（编辑减少）
# Huntington disease
# → 亨廷顿病
# Citrate cycle (TCA cycle)
# → 柠檬酸循环（三羧酸循环，TCA 循环）
# Influenza A
# → 甲型流感
# Measles
# → 麻疹
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# NOD-like receptor signaling pathway
# → NOD 样受体信号通路
# Virion - Ebolavirus, Lyssavirus and Morbillivirus
# → 病毒颗粒：埃博拉病毒、狂犬病毒属和麻疹病毒属
# Leishmaniasis
# → 利什曼病
# Taurine and hypotaurine metabolism
# → 牛磺酸和亚牛磺酸代谢
# Fatty acid biosynthesis
# → 脂肪酸生物合成
# Endocrine and other factor-regulated calcium reabsorption
# → 内分泌及其他因子调控的钙重吸收
# Coronavirus disease - COVID-19
# → 冠状病毒病 COVID-19
# Terpenoid backbone biosynthesis
# → 萜类骨架生物合成
# Epstein-Barr virus infection
# → EB 病毒感染
# mTOR signaling pathway
# → mTOR 信号通路
# A_to_I_Plots | GO BP | UP（编辑增加）
# ribonucleoprotein complex biogenesis
# → 核糖核蛋白复合体生成
# regulation of viral life cycle
# → 病毒生命周期的调控
# regulation of type I interferon production
# → I 型干扰素产生的调控
# type I interferon production
# → I 型干扰素产生
# establishment of protein localization to organelle
# → 蛋白质定位至细胞器的建立
# pattern recognition receptor signaling pathway
# → 模式识别受体信号通路
# regulation of apoptotic DNA fragmentation
# → 凋亡性 DNA 片段化的调控
# innate immune response-activating signaling pathway
# → 激活先天免疫反应的信号通路
# tRNA processing
# → tRNA 加工
# protein folding
# → 蛋白质折叠
# positive regulation of response to biotic stimulus
# → 对生物性刺激反应的正调控
# mitochondrion organization
# → 线粒体组织
# maturation of SSU-rRNA
# → 小亚基 rRNA 成熟
# mRNA splicing, via spliceosome
# → 通过剪接体进行的 mRNA 剪接
# endoplasmic reticulum protein-containing complex
# → 含蛋白质的内质网复合体
# A_to_I_Plots | KEGG | UP（编辑增加）
# Tuberculosis
# → 结核病
# Ubiquinone and other terpenoid-quinone biosynthesis
# → 泛醌及其他萜醌类生物合成
# Lysosome biogenesis
# → 溶酶体生物发生
# Protein processing in endoplasmic reticulum
# → 内质网中的蛋白质加工
# Apoptosis
# → 细胞凋亡
# N-Glycan biosynthesis
# → N-糖链生物合成
# Fatty acid metabolism
# → 脂肪酸代谢
# Coronavirus disease - COVID-19
# → 冠状病毒病 COVID-19
# Peroxisome
# → 过氧化物酶体
# Terpenoid backbone biosynthesis
# → 萜类骨架生物合成
# Prion disease
# → 朊病毒病
# Amyotrophic lateral sclerosis
# → 肌萎缩侧索硬化症
# Various types of N-glycan biosynthesis
# → 多种类型的 N-糖链生物合成
# Shigellosis
# → 志贺菌病
# Parkinson disease
# → 帕金森病
# C_to_U_Plots | GO BP | DOWN（编辑减少）
# phosphate ion transmembrane transport
# → 磷酸根离子跨膜转运
# phosphate transmembrane transporter activity
# → 磷酸盐跨膜转运蛋白活性
# solute:proton symporter activity
# → 溶质:质子同向转运体活性
# GTPase regulator activity
# → GTP 酶调节因子活性
# intracellularly gated calcium channel activity
# → 细胞内门控钙通道活性
# phosphotransferase activity, phosphate group as acceptor
# → 磷酸转移酶活性，以磷酸基团为受体
# small GTPase binding
# → 小 GTP 酶结合
# kinesin binding
# → 驱动蛋白结合
# transition metal ion transmembrane transporter activity
# → 过渡金属离子跨膜转运蛋白活性
# C_to_U_Plots | KEGG | DOWN（编辑减少）
# Phosphatidylinositol signaling system
# → 磷脂酰肌醇信号系统
# Long-term depression
# → 长时程抑制
# Glycerophospholipid metabolism
# → 甘油磷脂代谢
# Oxytocin signaling pathway
# → 催产素信号通路
# Long-term potentiation
# → 长时程增强
# Glutamatergic synapse
# → 谷氨酸能突触
# Biosynthesis of nucleotide sugars
# → 核苷酸糖生物合成
# Amino sugar and nucleotide sugar metabolism
# → 氨基糖和核苷酸糖代谢
# Platelet activation
# → 血小板活化
# MAPK signaling pathway
# → MAPK 信号通路
# Lysosome biogenesis
# → 溶酶体生物发生
# GnRH signaling pathway
# → 促性腺激素释放激素信号通路
# Ovarian steroidogenesis
# → 卵巢类固醇生成
# Inflammatory mediator regulation of TRP channels
# → 炎症介质对 TRP 通道的调控
# Sphingolipid metabolism
# → 鞘脂代谢
# C_to_U_Plots | GO BP | UP（编辑增加）
# early endosome
# → 早期内体
# growth cone
# → 生长锥
# site of polarized growth
# → 极性生长位点
# secretory granule membrane
# → 分泌颗粒膜
# C_to_U_Plots | KEGG | UP（编辑增加）
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Fanconi anemia pathway
# → 范可尼贫血通路
# Endocytosis
# → 内吞作用
# Spinocerebellar ataxia
# → 脊髓小脑性共济失调
# Amoebiasis
# → 阿米巴病
# Oxytocin signaling pathway
# → 催产素信号通路
# Folate biosynthesis
# → 叶酸生物合成
# Cholinergic synapse
# → 胆碱能突触
# Parathyroid hormone synthesis, secretion and action
# → 甲状旁腺激素的合成、分泌和作用
# Leukocyte transendothelial migration
# → 白细胞跨内皮迁移
# Glutamatergic synapse
# → 谷氨酸能突触
# Autophagy - animal
# → 自噬 - 动物
# Human cytomegalovirus infection
# → 人巨细胞病毒感染
# Coronavirus disease - COVID-19
# → 冠状病毒病 COVID-19
# Circadian rhythm
# → 昼夜节律
# Score [0] 五大维度交叉富集看板
# 1_C2U_Edit_and_DE_Sig_GO_BP.csv
# regulation of apoptotic cell clearance
# → 凋亡细胞清除的调控
# complement activation, GZMK pathway
# → 补体激活，GZMK 通路
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# complement activation, lectin pathway
# → 补体激活，凝集素通路
# activation of membrane attack complex
# → 膜攻击复合体激活
# MHC class II protein complex assembly
# → MHC II 类蛋白复合体组装
# peptide antigen assembly with MHC class II protein complex
# → 肽抗原与 MHC II 类蛋白复合体组装
# detection of bacterium
# → 细菌检测
# fibroblast activation
# → 成纤维细胞激活
# regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应调控
# 1_C2U_Edit_and_DE_Sig_KEGG.csv
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Phagosome
# → 吞噬体
# Viral myocarditis
# → 病毒性心肌炎
# Leishmaniasis
# → 利什曼病
# 2_A2I_Edit_and_DE_Sig_GO_BP.csv
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽或多糖抗原加工与呈递
# monocyte differentiation
# → 单核细胞分化
# antigen processing and presentation of exogenous peptide antigen
# → 外源性肽抗原加工与呈递
# antigen processing and presentation of exogenous antigen
# → 外源性抗原加工与呈递
# positive regulation of monocyte differentiation
# → 单核细胞分化的正调控
# positive regulation of chemokine (C-X-C motif) ligand 2 production
# → 趋化因子 C-X-C 基序配体 2 产生的正调控
# interleukin-10 production
# → 白细胞介素-10 产生
# regulation of interleukin-10 production
# → 白细胞介素-10 产生的调控
# 2_A2I_Edit_and_DE_Sig_KEGG.csv
# Tuberculosis
# → 结核病
# Antigen processing and presentation
# → 抗原加工与呈递
# Asthma
# → 哮喘
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Intestinal immune network for IgA production
# → IgA 产生相关的肠道免疫网络
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# N-Glycan biosynthesis
# → N-糖链生物合成
# 3_Any_Edit_and_DE_Sig_GO_BP.csv
# chemokine (C-X-C motif) ligand 2 production
# → 趋化因子 C-X-C 基序配体 2 产生
# regulation of chemokine (C-X-C motif) ligand 2 production
# → 趋化因子 C-X-C 基序配体 2 产生的调控
# regulation of B cell proliferation
# → B 细胞增殖的调控
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# regulation of lymphocyte proliferation
# → 淋巴细胞增殖的调控
# antigen processing and presentation of endogenous antigen
# → 内源性抗原加工与呈递
# regulation of mononuclear cell proliferation
# → 单核细胞增殖的调控
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# B cell proliferation
# → B 细胞增殖
# 3_Any_Edit_and_DE_Sig_KEGG.csv
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Tuberculosis
# → 结核病
# Antigen processing and presentation
# → 抗原加工与呈递
# Allograft rejection
# → 同种异体移植排斥
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Alcoholic liver disease
# → 酒精性肝病
# 4_Any_Edit_and_DE_Upregulated_GO_BP.csv
# regulation of B cell proliferation
# → B 细胞增殖的调控
# B cell proliferation
# → B 细胞增殖
# negative regulation of B cell proliferation
# → B 细胞增殖的负调控
# cristae formation
# → 嵴形成
# MyD88-dependent toll-like receptor signaling pathway
# → MyD88 依赖性 Toll 样受体信号通路
# 5_Any_Edit_and_DE_Downregulated_GO_BP.csv
# antigen processing and presentation of peptide antigen
# → 肽抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的外源性肽抗原加工与呈递
# antigen processing and presentation of endogenous antigen
# → 内源性抗原加工与呈递
# antigen processing and presentation of peptide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽抗原加工与呈递
# antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# → 通过 MHC II 类分子进行的肽或多糖抗原加工与呈递
# antigen processing and presentation of exogenous peptide antigen
# → 外源性肽抗原加工与呈递
# antigen processing and presentation
# → 抗原加工与呈递
# lymphocyte mediated immunity
# → 淋巴细胞介导的免疫
# adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# → 基于免疫球蛋白超家族结构域构成的免疫受体体细胞重组的适应性免疫反应
# antigen processing and presentation of exogenous antigen
# → 外源性抗原加工与呈递
# 5_Any_Edit_and_DE_Downregulated_KEGG.csv
# Herpes simplex virus 1 infection
# → 单纯疱疹病毒 1 型感染
# Antigen processing and presentation
# → 抗原加工与呈递
# Tuberculosis
# → 结核病
# Staphylococcus aureus infection
# → 金黄色葡萄球菌感染
# Allograft rejection
# → 同种异体移植排斥
# Type I diabetes mellitus
# → 1 型糖尿病
# Graft-versus-host disease
# → 移植物抗宿主病
# Systemic lupus erythematosus
# → 系统性红斑狼疮
# Autoimmune thyroid disease
# → 自身免疫性甲状腺病
# Phagosome
# → 吞噬体
