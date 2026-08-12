# Reproducibility and Statistical Design

## Single-cell RNA-seq

Public human atherosclerosis scRNA-seq datasets are loaded as study-specific Seurat objects, filtered using thresholds in `configs/config.yaml`, normalized, integrated, clustered, and annotated. The exact current whole-project source is retained as `workflow/01_scrna/single_cell_publication.R`; non-interactive preparation and checkpoint generation are implemented in the numbered scripts in the same directory. Sample provenance is recorded in `metadata/sample_manifest.csv`.

## MPI and MMI

MPI is the normalized AUCell activity of the 145-gene M1-like signature minus that of the 165-gene M2-like signature. Signature membership and literature provenance are frozen in `data/gene_sets/mpi_signatures.csv`.

MMI is independent of the carotid discovery data. GSE5099 and GSE11864 are re-analysed using donor-aware limma models for macrophage versus monocyte differentiation. Genes must be measured in both datasets and have the same effect direction. Within each direction, source-specific ranks are converted to percentiles and averaged; the best 200 concordant macrophage-associated and 200 concordant monocyte-associated genes are retained. MMI is `AUCell(macrophage maturation) - AUCell(monocyte-associated)`. The complete selection table and all excluded direction conflicts are versioned under `data/gene_sets/`.

## Pseudotime

Monocyte/macrophage trajectories are inferred with Monocle3. Final publication branch membership is read from `data_pseudotime.rds`, preserving the original interactive selections. The frozen checkpoint contains 6,712 Fate 1 cells and 9,450 Fate 2 cells with a shared trunk; scripts stop if the required `subset_1` and `subset_2` columns are absent. Later checkpoints with altered branch membership are not valid substitutes.

## Spatial transcriptomics

Xenium, Visium, and GeoMx analyses are ordered under `workflow/04_spatial/`. The workflow covers dataset QC, APOBEC3A detection and co-expression, myeloid label transfer, severity associations, neighborhood composition, pseudobulk differential expression, and algorithmically defined necrotic-core-like regions. GeoMx results are interpreted at ROI level, not as single-cell subtype calls. Core-like regions remain algorithmic spatial classes unless independently pathologically validated.

## Bulk RNA-seq

featureCounts integer counts are matched to `metadata/rna_editing_sample_manifest.tsv`. THP-1 samples are excluded. Clone13 and clone37 are analysed separately using `design = ~ condition`; the combined 12-sample analysis uses `design = ~ clone + condition`. The condition coefficient is KO versus WT. Benjamini-Hochberg adjusted P values are used for DEG classification. GO, GSEA, and GSVA use the same fitted contrasts.

## RNA editing

JACUSA2 v2.0.4 `call-2` receives WT BAMs as `cond1` and KO BAMs as `cond2`; the exact order is written to a manifest. Clone13, clone37, and combined calls are produced by the same Python implementation. Sites are filtered by JACUSA2 score, minimum coverage in the configured number of replicates, and the JACUSA2helper robust filter. C-to-U candidates include C>T and complementary G>A events; A-to-I candidates include A>G and complementary T>C events.

Editing ratios are modelled with limma. The combined model includes clone before condition; clone-specific models contain condition only. `delta_ratio` is always KO minus WT. Existing nominal `P < 0.10` and absolute delta at least 0.05 labels are retained to reproduce the current tables; BH-FDR and an FDR-based status are exported alongside them. This distinction must be stated when reporting exploratory versus FDR-supported sites.

## Large data and checkpoints

Large RDS, BAM, FASTQ, public count matrices, JACUSA2 output, and reference files are not committed. Their configured paths are validated with:

```bash
Rscript scripts/validate_repo.R --config configs/config.yaml --check-inputs
```

Every result directory is generated and ignored by Git. The repository contains no credentials, raw human sequencing reads, or patient-level protected data.
