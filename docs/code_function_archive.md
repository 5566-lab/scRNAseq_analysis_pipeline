# Code Function Archive

Source directory: `/public3/DSC/single_cell/GSE159677_Carotid_MainProject`

## Retained Implementation Sources

| Source | Status | Reason |
|---|---|---|
| `single_cell.R` | Reference only | Latest whole-workflow script; too broad for direct pipeline use, split into modular scripts. |
| `hdWGCNA/Mo_Ma/MM_WGCNA.R` | Retained for hdWGCNA logic | Newer and more focused than `hdWGCNA/hdWGCNA.R` for monocyte/macrophage analysis. |
| `monocle3/monocle_MM/MM_monocle.R` | Retained for pseudotime logic | Newer focused Monocle3 workflow; replaces duplicated blocks in `single_cell.R`. |
| `GSEA/GSEA_GSVA.R` | Retained for pathway logic | Main GSEA/GSVA and APOBEC3A-KO comparison script. |
| `macSpectrum.R` | Retained for MPI/AMDI scoring | Focused macrophage state scoring. |
| `untitled9.R` | Supplemental only | Plotting scratch code; not included as a primary pipeline step. |
| `hdWGCNA/hdWGCNA.R` | Superseded | Older generic Mono/Macro hdWGCNA script. |

## Normalized Pipeline Mapping

| New script | Source logic |
|---|---|
| `scripts/01_prepare_seurat_objects.R` | Dataset loading blocks from `single_cell.R`. |
| `scripts/02_integrate_cluster.R` | QC, integration, clustering, marker output blocks from `single_cell.R`. |
| `scripts/03_hdWGCNA_mo_ma.R` | Latest Mo/Ma hdWGCNA workflow from `hdWGCNA/Mo_Ma/MM_WGCNA.R`. |
| `scripts/04_macSpectrum_scores.R` | MPI/AMDI scoring from `macSpectrum.R`. |
| `scripts/05_monocle_pseudotime.R` | Trajectory workflow from `monocle3/monocle_MM/MM_monocle.R`. |
| `scripts/06_gsea_gsva.R` | scMetabolism, GSEA, GSVA, APOBEC3A-KO comparison from `GSEA/GSEA_GSVA.R`. |

## Loaded Datasets

`scripts/01_prepare_seurat_objects.R` now covers the source loading blocks for:

- `GSE260657`
- `GSE247238`
- `GSE131778`
- `GSE210152`
- `GSE155468`
- `GSE159677`
- `GSE213740`
- `GSE216860`

The exploratory script notes that `GSE213740` overlaps with samples in `GSE216860`; this repository keeps both loaders and records that relationship in `metadata/sample_manifest.csv` so the decision to exclude one can be made explicitly during analysis.
