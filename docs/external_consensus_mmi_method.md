# External-consensus macrophage maturation index

## Rationale

The original exploratory maturation score used the top markers of Classical Mono identified in the carotid dataset itself. That design can create circular evidence when the same dataset is used both to select the signature and to demonstrate a maturation trend. The publication pipeline therefore replaces it with a fixed external-consensus macrophage maturation index (MMI). No gene is selected by differential expression in the carotid study data.

## External signatures

The score combines three independent external sources:

- **GSE5099:** MSigDB C7 day-7 M-CSF macrophage-high and unstimulated monocyte-high signatures.
- **GSE11864:** MSigDB C7 CSF1-cultured macrophage-high and untreated monocyte-high signatures.
- **HPCA:** the top 100 Macrophage-versus-Monocyte and Monocyte-versus-Macrophage markers obtained from the external Human Primary Cell Atlas reference with `SingleR::getClassicMarkers()`.

The exploratory 20-gene C1Q tissue-macrophage module is not included as a source because it has no matched monocyte-negative counterpart and primarily represents tissue adaptation/lipid handling rather than maturation alone. Individual C1Q-family genes are retained when they are independently selected by GSE5099, GSE11864, or the HPCA marker comparison; in the frozen list, `C1QA`, `C1QB`, and `C1QC` originate from HPCA rather than from the excluded exploratory module.

Genes are deduplicated within each direction. Six genes reported in opposite directions by different external sources (`DPY19L4`, `FMNL2`, `HGS`, `RDX`, `TRAM1`, and `ZBTB33`) are removed from both directions. The frozen final lists contain 476 mature-positive and 476 monocyte/immature-negative genes and are stored in `data/gene_sets/external_consensus_mmi_no_c1q.csv`. Source-specific membership is retained in `data/gene_sets/external_consensus_mmi_no_c1q_sources.csv`.

## AUCell calculation

AUCell is applied to the unintegrated Seurat `RNA` assay `counts` layer. Genes are ranked within each cell with `AUCell_buildRankings()`. Normalized AUC values are calculated with `AUCell_calcAUC(normAUC = TRUE)` and `aucMaxRank` set to 5% of ranked genes. Missing signature genes are excluded by intersecting each fixed list with RNA-assay features before scoring.

For cell *i*, the maturation index is:

```text
MMI_i = AUC_i(mature-positive genes) - AUC_i(monocyte/immature-negative genes)
```

Higher values indicate greater enrichment of the external macrophage-maturation program relative to the external monocyte/immaturity program. Scores are relative transcriptional indices and should not be interpreted as direct measurements of developmental time or lineage origin.

The original macrophage polarization index is retained:

```text
MPI_i = AUC_i(M1 genes) - AUC_i(M2 genes)
```

The fixed M1 and M2 lists remain defined in `scripts/04b_auc_macrophage_signatures.R` to preserve the original analysis.

## Statistical analysis

UMAP and violin plots display cell-level score distributions for visualization. Statistical inference is performed on the median score of each biological sample, identified by `Source_GSE`, `Patient_ID`, and `Sample_Type`.

- Core-versus-adjacent comparisons use a paired Wilcoxon signed-rank test when at least three complete patient pairs are available; otherwise a biological-sample-level Wilcoxon rank-sum test is used.
- Foam cells1, Foam cells2, and LAM comparisons use paired biological-sample medians when both subtypes are present in the same sample.
- Benjamini-Hochberg correction is applied within each displayed family of tests.
- The analysis unit, sample counts, raw P values, and adjusted P values are exported with the figures.

Because lesion region and source cohort may be partially confounded in pooled public data, region-level results should be interpreted as descriptive unless supported by within-patient pairs or an appropriately specified cohort-aware model.

## Reproducibility outputs

`scripts/04c_publication_mpi_mmi_plots.R` exports:

- the clean and source-specific gene-set tables;
- per-cell AUC and MMI values;
- biological-sample and sample-by-subtype medians;
- raw and multiplicity-adjusted statistical tests;
- AUCell parameters and gene-set sizes;
- `sessionInfo()`;
- PDF and 300-dpi PNG versions of all figures.

The plotting geometry follows the original `single_cell.R` figures, while statistical annotations use biological replicates to avoid cell-level pseudoreplication.
