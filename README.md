# APOBEC3A in Human Atherosclerosis: Multi-omics Analysis Pipeline

Publication code for the study **Single-Cell and Spatial Transcriptomics Identify APOBEC3A as a Myeloid Regulator in Human Atherosclerosis**. The repository links public single-cell RNA-seq, macrophage-state scoring, Monocle3 trajectories, Xenium/Visium/GeoMx spatial analyses, APOBEC3A-knockout bulk RNA-seq, RNA editing, and supplementary-table generation in one ordered workflow.

## Analysis graph

```text
Public scRNA-seq matrices
  -> QC / integration / annotation
  -> monocyte-macrophage hdWGCNA
  -> MPI and external MMI scoring
  -> Monocle3 pseudotime and divergent fates
  -> spatial projection and core-like-region statistics

APOBEC3A-KO BAM files
  -> featureCounts -> clone13 / clone37 / combined DESeq2
  -> JACUSA2 call-2 -> editing-ratio limma
  -> editing-expression and cross-clone consistency

All machine-readable outputs
  -> Supplementary Tables S2-S10
```

## Repository layout

| Path | Purpose |
|---|---|
| `configs/config.yaml` | All machine-specific inputs and analysis thresholds |
| `metadata/` | Public-dataset manifest, frozen annotation map, and 6 WT + 6 KO manifest |
| `data/gene_sets/` | Frozen MPI, Top200 MMI, and APOBEC3A-positive signatures |
| `workflow/01_scrna/` | scRNA-seq preparation, integration, hdWGCNA, and exact publication source |
| `workflow/02_scoring/` | External GEO ranking, deterministic Top200 selection, AUCell MPI/MMI figures |
| `workflow/03_trajectory/` | Monocle3 publication figures |
| `workflow/04_spatial/` | Xenium, Visium, GeoMx, and core-like-region analyses |
| `workflow/05_bulk_rnaseq/` | Unified clone13, clone37, and combined DESeq2/GSEA/GSVA analysis |
| `workflow/06_rna_editing/` | JACUSA2 calling, editing statistics, and contrast consistency |
| `workflow/07_supplementary_tables/` | Supplementary Tables S2-S10 |
| `scripts/run_pipeline.R` | Ordered stage runner |
| `scripts/validate_repo.R` | Static and scientific-invariant validation |

## Run

Edit `configs/config.yaml`, then validate without launching expensive analyses:

```bash
Rscript scripts/validate_repo.R --config configs/config.yaml
Rscript scripts/run_pipeline.R --config configs/config.yaml --stage all --dry-run
```

Run selected stages:

```bash
Rscript scripts/run_pipeline.R --config configs/config.yaml \
  --stage scoring,trajectory

Rscript scripts/run_pipeline.R --config configs/config.yaml \
  --stage bulk,rna-editing
```

Available stages are `scrna`, `scoring`, `trajectory`, `pathway`, `spatial`, `bulk`, `rna-editing`, and `supplementary`. Outputs are written below `results/` and are ignored by Git.

## Reproducibility notes

The final annotated Seurat object, Monocle3 CDS, spatial Seurat objects, BAM files, reference genome, and raw public matrices are too large for Git. Their paths are declared in `configs/config.yaml`; accession-level provenance is retained in `metadata/`. `workflow/01_scrna/single_cell_publication.R` is the exact current publication source and is retained for audit, but it contains the original interactive branch-selection block. Automated final figures instead use the frozen non-interactive scripts in `workflow/03_trajectory/`.

Software requirements are listed in `DESCRIPTION`, `environment.yml`, and `requirements.txt`. JACUSA2 v2.0.4, Java, featureCounts, and a GRCh38 reference/annotation are external command-line requirements.

## License

MIT. Public datasets remain subject to their original repository terms and citation requirements.
