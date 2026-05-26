# Final Coverage Audit

Audit source: `/public3/DSC/single_cell/GSE159677_Carotid_MainProject`

Audit date: 2026-05-26

## Source R Files

| Source file | Repository status |
|---|---|
| `single_cell.R` | Archived and used as the primary reference for dataset loading, integration, Mo/Ma analysis, custom AUCell signatures, and Monocle3 logic. |
| `hdWGCNA/Mo_Ma/MM_WGCNA.R` | Archived and modularized into `scripts/03_hdWGCNA_mo_ma.R`. |
| `monocle3/monocle_MM/MM_monocle.R` | Archived and modularized into `scripts/05_monocle_pseudotime.R`. |
| `GSEA/GSEA_GSVA.R` | Archived and partially modularized into `scripts/06_gsea_gsva.R`. |
| `macSpectrum.R` | Archived and retained as optional `scripts/04_macSpectrum_scores.R`; excluded from the default pipeline. |
| `hdWGCNA/hdWGCNA.R` | Superseded by the newer Mo/Ma-specific `MM_WGCNA.R`. |
| `untitled9.R` | Supplemental plotting scratch script; archived status only. |

## Dataset Loading Coverage

All single-cell GSE datasets referenced by the source `single_cell.R` loader are represented in `configs/config.yaml`, `metadata/sample_manifest.csv`, and `scripts/01_prepare_seurat_objects.R`.

| Dataset | Status |
|---|---|
| `GSE260657` | Covered |
| `GSE247238` | Covered |
| `GSE131778` | Covered with source path `/public3/DSC/single_cell/GSE131778_Coronary_AC/GSE131778_human_coronary_scRNAseq.txt` |
| `GSE210152` | Covered with source path `/public3/DSC/single_cell/GSE210152_Carotid_AC/GSE210152_raw.RDS` |
| `GSE155468` | Covered |
| `GSE159677` | Covered |
| `GSE213740` | Covered; source notes overlap with `GSE216860` |
| `GSE234077` | Covered |
| `GSE224273` | Covered |
| `GSE253903` | Covered |
| `GSE216860` | Covered |

## Default Pipeline Coverage

| Analysis block | Status |
|---|---|
| Multi-dataset Seurat object construction | Covered |
| QC, integration, clustering, UMAP, cluster markers | Covered |
| Mo/Ma hdWGCNA modules and hub genes | Covered |
| Custom M1/M2/Mono immaturity AUCell scoring from `single_cell.R` | Covered and set as default |
| Monocle3 trajectory, pseudotime, APOBEC3A trajectory plot, graph-test genes | Covered |
| Foam cell vs Macrophage GSEA/GSVA | Covered |
| APOBEC3A-KO clone13 and clone37 GSEA/GSVA | Covered |
| `macSpectrum::macspec()` scoring | Optional auxiliary script only |

## Archived But Not Fully Automated

These source blocks remain preserved in `archive/original_scripts/` but are not yet default, noninteractive pipeline steps:

| Source block | Reason |
|---|---|
| scMetabolism KEGG heatmaps from `GSEA_GSVA.R` | Requires `scMetabolism` and generates selected metabolism heatmaps; not included in default GSEA/GSVA script yet. |
| Combined `clone13+clone37` APOBEC3A-KO comparison | Source code exists but current default script covers clone13 and clone37 separately. |
| Manual/interactive Monocle3 branch selection using `choose_graph_segments()` | Requires interactive branch selection; not safe for unattended pipeline execution. |
| Advanced branch-dynamics GLM/AUC tables: `T1_Global_Trajectory_Dependent_Genes.csv`, `T2_Model_Based_Hetero_Genes_Final.csv`, `T3_T5_Bifurcation_Key_Genes.xlsx` | Depends on branch-specific objects and modeling setup from late exploratory sections. |
| Late-stage figure polish: `MainFigure_Top20_Complete_Story.pdf`, `F4.1_Final_Clean_Heatmap.pdf`, `F5_Bifurcation_Key_TFs.pdf`, `F2.6_Pseudotime_Density_Celltype_Tracks.*`, `F3.1_APOBEC3A_*` | These are publication-figure construction blocks rather than core reproducible analysis steps. |
| `untitled9.R` selected metabolic pathway heatmap | Scratch plotting script depends on objects already in memory. |

## Validation Performed

- Parsed all non-archive R scripts with `Rscript scripts/validate_repo.R`.
- Checked that all configured source paths exist.
- Compared `GSE[0-9]+` identifiers in `single_cell.R` against repository config, metadata, scripts, and docs.
