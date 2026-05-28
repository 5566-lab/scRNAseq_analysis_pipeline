# scRNA-seq Analysis Pipeline

Configuration-driven R pipeline for public scRNA-seq analysis, including Seurat integration, monocyte/macrophage hdWGCNA, custom AUCell M1/M2 macrophage signature scoring, Monocle3 pseudotime analysis, and GSEA/GSVA pathway comparison with APOBEC3A-KO RNA-seq results.

## Pipeline

1. `scripts/01_prepare_seurat_objects.R` builds per-study Seurat objects from public count matrices.
2. `scripts/02_integrate_cluster.R` performs QC, integration, clustering, marker detection, and annotation-ready outputs.
3. `scripts/03_hdWGCNA_mo_ma.R` runs hdWGCNA on monocyte/macrophage subsets.
4. `scripts/04b_auc_macrophage_signatures.R` computes custom AUCell M1/M2/Mono immaturity scores from the curated gene sets in the source script.
5. `scripts/05_monocle_pseudotime.R` performs Monocle3 trajectory and branch analysis.
6. `scripts/06_gsea_gsva.R` performs GSEA, GSVA, and APOBEC3A-KO pathway comparisons.

The AUCell step also regenerates the S2.6 MPI/AMDI violin plots, LAM/Foam-cell density plot, and Foam marker violin panels with the updated font sizing.

Optional auxiliary scoring:

```bash
Rscript scripts/04_macSpectrum_scores.R --config configs/config.yaml
```

`04_macSpectrum_scores.R` uses the `macSpectrum` package model and is not part of the default pipeline. The default macrophage polarization and maturation indices are the AUCell scores from `single_cell.R`.

Run the full pipeline:

```bash
Rscript scripts/run_pipeline.R --config configs/config.yaml
```

Run one step:

```bash
Rscript scripts/03_hdWGCNA_mo_ma.R --config configs/config.yaml
```

## Configuration

Edit `configs/config.yaml` for input paths, output directories, filtering thresholds, analysis parameters, and APOBEC3A-KO count/DEG files. No script requires project-specific absolute paths outside the config file.

The current dataset loader covers `GSE260657`, `GSE247238`, `GSE131778`, `GSE210152`, `GSE155468`, `GSE159677`, `GSE213740`, `GSE234077`, `GSE224273`, `GSE253903`, and `GSE216860`.

See `docs/final_coverage_audit.md` for the final coverage check against the exploratory scripts, including modules that remain archived but are not automated in the default pipeline.

## Source Code Policy

Original exploratory scripts are summarized in `docs/code_function_archive.md`. When multiple scripts had overlapping functions, the latest or most focused version was retained as the implementation source:

- `single_cell.R` is treated as the latest whole-workflow reference.
- `hdWGCNA/Mo_Ma/MM_WGCNA.R` supersedes the older generic `hdWGCNA/hdWGCNA.R` for Mo/Ma module analysis.
- `monocle3/monocle_MM/MM_monocle.R` supersedes duplicated Monocle3 blocks in `single_cell.R`.
- `untitled9.R` is treated as a supplemental plotting scratch script, not a pipeline step.
