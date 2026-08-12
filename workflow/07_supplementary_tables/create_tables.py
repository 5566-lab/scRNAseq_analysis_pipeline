#!/usr/bin/env python3
"""Create publication-style supplementary Excel tables from local analysis outputs."""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Iterable

import pandas as pd
from scipy.stats import ttest_ind


REPO_ROOT = Path(
    os.environ.get(
        "PIPELINE_REPO_ROOT",
        Path(__file__).resolve().parents[2],
    )
).resolve()
ROOT = Path(
    os.environ.get("SUPPLEMENTARY_SOURCE_ROOT", "/public3/DSC/single_cell")
).resolve()
OUT = Path(
    os.environ.get("SUPPLEMENTARY_OUTPUT_DIR", REPO_ROOT / "results" / "supplementary_tables")
).resolve()
OUT.mkdir(parents=True, exist_ok=True)
GENCODE_GFF = Path(
    os.environ.get(
        "GENCODE_GFF",
        "/dell_1/dsc/RNA_seq/data/GFF/gencode.v46.annotation.gff3",
    )
)


def p(rel: str) -> Path:
    path = Path(rel)
    return path if path.is_absolute() else ROOT / rel


def read_csv(rel: str, **kwargs) -> pd.DataFrame:
    path = p(rel)
    if not path.exists():
        return pd.DataFrame({"note": [f"Missing source file: {path}"]})
    df = pd.read_csv(path, **kwargs)
    df["source_file"] = str(path)
    return df


def read_excel_sheet(rel: str, sheet: str) -> pd.DataFrame:
    path = p(rel)
    if not path.exists():
        return pd.DataFrame({"note": [f"Missing source file: {path}"]})
    df = pd.read_excel(path, sheet_name=sheet)
    df["source_file"] = str(path)
    return df


def clean_sheet_name(name: str) -> str:
    bad = "[]:*?/\\"
    for ch in bad:
        name = name.replace(ch, "_")
    return name[:31]


def write_xlsx(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            if df is None:
                df = pd.DataFrame()
            sheet = clean_sheet_name(name)
            df.to_excel(writer, sheet_name=sheet, index=False)
            ws = writer.book[sheet]
            ws.freeze_panes = "A2"
            if ws.max_row > 1 and ws.max_column > 0:
                ws.auto_filter.ref = ws.dimensions
            for column_cells in ws.columns:
                header = str(column_cells[0].value or "")
                width = min(max(len(header) + 2, 12), 45)
                ws.column_dimensions[column_cells[0].column_letter].width = width


def add_source(df: pd.DataFrame, src: str) -> pd.DataFrame:
    out = df.copy()
    out["source_file"] = str(p(src))
    return out


def first_existing(paths: Iterable[str]) -> str | None:
    for rel in paths:
        if p(rel).exists():
            return rel
    return None


_GENE_INTERVALS: dict[str, list[tuple[int, int, str]]] | None = None


def parse_gff_gene_name(attr: str) -> str:
    for item in attr.split(";"):
        if item.startswith("gene_name="):
            return item.split("=", 1)[1]
        if item.startswith("Name="):
            return item.split("=", 1)[1]
    return ""


def gene_intervals() -> dict[str, list[tuple[int, int, str]]]:
    global _GENE_INTERVALS
    if _GENE_INTERVALS is not None:
        return _GENE_INTERVALS
    intervals: dict[str, list[tuple[int, int, str]]] = {}
    if not GENCODE_GFF.exists():
        _GENE_INTERVALS = intervals
        return intervals
    with GENCODE_GFF.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0].removeprefix("chr")
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            name = parse_gff_gene_name(parts[8])
            if name and not name.startswith("ENSG"):
                intervals.setdefault(chrom, []).append((start, end, name))
    for chrom in intervals:
        intervals[chrom].sort()
    _GENE_INTERVALS = intervals
    return intervals


def annotate_sites_with_genes(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "chr" not in df.columns or "pos" not in df.columns:
        return df
    intervals = gene_intervals()
    if not intervals:
        return df
    genes = []
    regions = []
    for chrom, pos in zip(df["chr"], df["pos"]):
        chrom = str(chrom).removeprefix("chr")
        try:
            pos_int = int(pos)
        except (TypeError, ValueError):
            genes.append(pd.NA)
            regions.append("not_annotated")
            continue
        hit = None
        for start, end, gene in intervals.get(chrom, []):
            if start > pos_int:
                break
            if start <= pos_int <= end:
                hit = gene
                break
        if hit:
            genes.append(hit)
            regions.append("gene_overlap")
        else:
            genes.append(pd.NA)
            regions.append("intergenic_or_unannotated")
    out = df.copy()
    out["gene_name"] = genes
    out["region_annotation"] = regions
    return out


def dataset_overview() -> pd.DataFrame:
    qc = read_csv("Result/Sample_Cell_Counts_QC_Summary.csv")
    sc_sum = pd.DataFrame()
    if "Source_GSE" in qc.columns:
        sc_sum = (
            qc.groupby("Source_GSE", dropna=False)
            .agg(
                n_samples=("Patient_ID", "nunique"),
                n_GSM=("GSM_ID", lambda x: x.dropna().nunique()),
                cells_before_QC=("Count_Before_QC", "sum"),
                cells_after_QC=("Count_After_QC", "sum"),
            )
            .reset_index()
        )
    code_info = pd.DataFrame(
        [
            ["GSE131778", "scRNA-seq", "single count matrix", "coronary atherosclerotic core", "atherosclerotic core"],
            ["GSE155468", "scRNA-seq", "file-based count matrices", "ascending aorta control/proximal adjacent", "control/proximal adjacent"],
            ["GSE159677", "scRNA-seq", "10x-format carotid samples", "carotid proximal adjacent and atherosclerotic core", "PA vs AC"],
            ["GSE216860", "scRNA-seq", "10x-format public dataset", "proximal adjacent/control aorta in code", "proximal adjacent/control"],
            ["GSE224273", "scRNA-seq", "10x-format public dataset", "carotid/coronary atherosclerotic core in code", "atherosclerotic core"],
            ["GSE234077", "scRNA-seq", "10x-format public dataset", "carotid/coronary atherosclerotic core in code", "atherosclerotic core"],
            ["GSE247238", "scRNA-seq", "10x matrix files plus GEO metadata", "carotid plaque", "stable vs unstable"],
            ["GSE253903", "scRNA-seq", "10x-format public dataset", "carotid atherosclerotic core", "atherosclerotic core"],
            ["GSE260657", "scRNA-seq", "multiple carotid sample files", "carotid plaque", "stable vs unstable"],
        ],
        columns=["dataset", "platform", "platform_detail_from_code", "sample_type_or_region_from_code", "plaque_status_or_severity_from_code"],
    )
    overview = code_info.merge(sc_sum, left_on="dataset", right_on="Source_GSE", how="left").drop(columns=["Source_GSE"], errors="ignore")

    spatial_qc = read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_dataset_qc_and_inclusion_Xenium_GeoMx_only.csv")
    if "dataset" in spatial_qc.columns:
        spatial = spatial_qc.rename(
            columns={
                "n_units": "included_cells_ROI_or_spots",
                "n_samples_or_fovs": "n_samples_ROI_or_FOV",
            }
        )
        spatial["sample_type_or_region_from_code"] = spatial.get("analysis_level", "")
        spatial["plaque_status_or_severity_from_code"] = "Mild/Moderate/Severe where encoded"
        spatial["platform_detail_from_code"] = spatial["analysis_level"].astype(str)
        spatial = spatial[
            [
                "dataset",
                "platform",
                "platform_detail_from_code",
                "sample_type_or_region_from_code",
                "plaque_status_or_severity_from_code",
                "included_cells_ROI_or_spots",
                "n_samples_ROI_or_FOV",
                "n_features",
                "n_myeloid_units",
                "usable_for_A3A",
                "note",
                "source_file",
            ]
        ]
        overview["included_cells_ROI_or_spots"] = overview["cells_after_QC"]
        overview["n_samples_ROI_or_FOV"] = overview["n_samples"]
        overview = pd.concat([overview, spatial], ignore_index=True, sort=False)
    overview["publication_table"] = "Supplementary Table 2"
    return overview


def table2() -> None:
    sample_qc = read_csv("Result/Sample_Cell_Counts_QC_Summary.csv")
    mps_sample = read_csv("Supplementary_Tables_Publication/_intermediate/mps_sample_metadata_summary.csv")
    xenium_overview = read_csv("spatial/atherosclerosis_spatial_publication/results/tables/xenium_cell_overview.csv")
    xenium_fov = read_csv("spatial/atherosclerosis_spatial_publication/results/tables/xenium_representative_fovs_for_spatial_maps.csv")
    geomx_roi = read_csv("spatial/atherosclerosis_spatial_publication/results/figures_1/GeoMx_ROI_raw_image_check/GeoMx_ROI_table_used_for_overview.csv")
    geomx_qc = read_csv("spatial/atherosclerosis_spatial_publication/results/tables/geomx_roi_qc_by_grade_localisation.csv")
    notes = pd.DataFrame(
        {
            "note": [
                "Experimental oligonucleotide information is intentionally excluded per user request.",
                "scRNA-seq platform details not explicitly encoded in local scripts are reported as code-derived input format rather than inferred protocol.",
                "GSE315246 Xenium and GSE277170 GeoMx metadata are taken from spatial publication tables under spatial/atherosclerosis_spatial_publication/results/tables.",
            ]
        }
    )
    write_xlsx(
        OUT / "Supplementary_Table_2_Public_datasets_and_sample_metadata.xlsx",
        {
            "Table S2A Dataset overview": dataset_overview(),
            "Table S2B scRNA sample QC": sample_qc,
            "Table S2C MPS sample meta": mps_sample,
            "Table S2D Xenium severity": xenium_overview,
            "Table S2E Xenium FOVs": xenium_fov,
            "Table S2F GeoMx ROI metadata": geomx_roi,
            "Table S2G GeoMx ROI QC": geomx_qc,
            "Table S2H Notes": notes,
        },
    )


def table3() -> None:
    all_markers = read_csv("Result/RPCA_All_Cell_Markers.csv")
    mps_main_markers = read_csv("Result/Celltype_raw1_all_markers.csv")
    mps_markers = read_csv("Result/figer_new/ST2_Ma_Mo_Marker.csv")
    mono_markers = read_csv("Result/figer_new/ST2_Mono_cluster_markers.csv")
    cols = ["gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"]
    def reorder(df: pd.DataFrame) -> pd.DataFrame:
        keep = [c for c in cols if c in df.columns] + [c for c in df.columns if c not in cols and c != "Unnamed: 0"]
        return df[keep]
    write_xlsx(
        OUT / "Supplementary_Table_3_Cell_type_annotation_markers.xlsx",
        {
            "Table S3A Main cell markers": reorder(all_markers),
            "Table S3B MPS subtype markers": reorder(mps_main_markers),
            "Table S3C MPS marker backup": reorder(mps_markers),
            "Table S3D Mono cluster markers": reorder(mono_markers),
        },
    )


def add_foam_lam(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in ["Foam cells1", "Foam cells2", "LAM"]:
        if col not in out.columns:
            out[col] = 0
    subtype_cols = [
        "Classical Mono",
        "CX3CR1+ TRM",
        "Foam cells1",
        "Foam cells2",
        "Inflammatory Mono",
        "ISG+ Mono",
        "LAM",
        "LYVE1+ TRM",
        "Non-classical Mono",
        "Transitional Mac",
    ]
    present_subtypes = [c for c in subtype_cols if c in out.columns]
    is_proportion_sheet = False
    if present_subtypes:
        numeric = out[present_subtypes].apply(pd.to_numeric, errors="coerce")
        is_proportion_sheet = bool((numeric.max(numeric_only=True).max() <= 1.000001))
    out["Foam_cells1_plus_Foam_cells2"] = out["Foam cells1"] + out["Foam cells2"]
    out["Foam_cells1_plus_Foam_cells2_plus_LAM"] = out["Foam cells1"] + out["Foam cells2"] + out["LAM"]
    if is_proportion_sheet:
        out["Foam_cells1_plus_Foam_cells2_fraction"] = out["Foam_cells1_plus_Foam_cells2"]
        out["Foam_cells1_plus_Foam_cells2_plus_LAM_fraction"] = out["Foam_cells1_plus_Foam_cells2_plus_LAM"]
    elif "total_MPS_cells" in out.columns:
        denom = out["total_MPS_cells"].replace(0, pd.NA)
        out["Foam_cells1_plus_Foam_cells2_fraction"] = out["Foam_cells1_plus_Foam_cells2"] / denom
        out["Foam_cells1_plus_Foam_cells2_plus_LAM_fraction"] = out["Foam_cells1_plus_Foam_cells2_plus_LAM"] / denom
    group_cols = [c for c in ["AC_PA", "Sample_Type", "Source_GSE", "Patient_ID"] if c in out.columns]
    if group_cols:
        group_col = group_cols[0]
        out.insert(0, "comparison_group", out[group_col])
        order = ["comparison_group", group_col] + [c for c in out.columns if c not in {"comparison_group", group_col}]
        out = out[order]
    return out


def table4() -> None:
    sheets = {}
    for grp in ["AC_PA", "Sample_Type", "Source_GSE", "Patient_ID"]:
        for kind in ["counts", "proportions"]:
            rel = f"Supplementary_Tables_Publication/_intermediate/composition_{kind}_by_{grp}.csv"
            sheets[f"S4 {kind} by {grp}"] = add_foam_lam(read_csv(rel))
    notes = pd.DataFrame(
        {
            "note": [
                "MPS composition was computed from Result/figer_new/group_result.rds metadata.",
                "Requested proximal adjacent vs atherosclerotic core is represented by Sample_Type/AC_PA sheets.",
                "Stable vs vulnerable is represented where encoded in AC_PA as Carotid Stable Plaque and Carotid Unstable Plaque.",
            ]
        }
    )
    sheets["Table S4I Notes"] = notes
    write_xlsx(OUT / "Supplementary_Table_4_MPS_composition_statistics.xlsx", sheets)


def wgcna_summary() -> pd.DataFrame:
    genes = read_csv("Result/figer_new/hdWGCNA/Mo_Ma/data_genes.csv")
    modules = read_csv("Result/figer_new/hdWGCNA/Mo_Ma/data_modules.csv")
    go = read_csv("Result/figer_new/hdWGCNA/GO/GO_all_pathways.csv")
    if "module" in genes.columns:
        genes_sorted = genes.sort_values(["module", "kME"], ascending=[True, False])
        hub = genes_sorted.groupby("module").head(20).groupby("module")["gene_name"].apply(lambda x: ";".join(map(str, x))).reset_index(name="top20_hub_genes_by_kME")
    else:
        hub = pd.DataFrame()
    if "Module" in go.columns:
        go_sorted = go.sort_values(["Module", "p.adjust"], ascending=[True, True])
        topgo = (
            go_sorted.groupby("Module")
            .head(5)
            .groupby("Module")
            .agg(
                top5_enriched_terms=("Description", lambda x: ";".join(map(str, x))),
                min_top5_adjusted_p=("p.adjust", "min"),
            )
            .reset_index()
        )
        topgo["module"] = topgo["Module"].astype(str).str.lower()
    else:
        topgo = pd.DataFrame()
    mod_counts = pd.DataFrame()
    if "module" in modules.columns:
        mod_counts = modules.groupby(["module", "color"], dropna=False).size().reset_index(name="n_module_genes")
    out = mod_counts.merge(hub, on="module", how="left").merge(topgo.drop(columns=["Module"], errors="ignore"), on="module", how="left")
    out["publication_table"] = "Supplementary Table 5"
    return out


def wgcna_notes() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "note": [
                "Table S5A is limited to machine-readable module outputs: module size, hub genes by kME, and top enriched GO/KEGG terms.",
                "The previous top_categories column was derived from the Category column in Result/figer_new/hdWGCNA/GO/GO_all_pathways.csv. That source Category field is a manually curated Chinese pathway class and is retained only in Table S5D as source provenance.",
                "The previous cell_state_annotation_from_enrichment column was removed because it was not an independently computed cell-state annotation; it duplicated the source pathway Category labels.",
                "The previous region_bias_from_code placeholder column was removed from Table S5A because no saved machine-readable region-bias result table was found among the current hdWGCNA outputs.",
                "The script GSE159677_Carotid_MainProject/hdWGCNA/Mo_Ma/MM_WGCNA.R contains code to compare module scores between atherosclerotic core and proximal adjacent regions by Wilcoxon test, but the current available outputs do not include a saved CSV/XLSX result for that comparison.",
            ]
        }
    )


def table5() -> None:
    write_xlsx(
        OUT / "Supplementary_Table_5_WGCNA_module_annotation.xlsx",
        {
            "Table S5A Module summary": wgcna_summary(),
            "Table S5B Gene kME": read_csv("Result/figer_new/hdWGCNA/Mo_Ma/data_genes.csv"),
            "Table S5C Module assignment": read_csv("Result/figer_new/hdWGCNA/Mo_Ma/data_modules.csv"),
            "Table S5D GO all pathways": read_csv("Result/figer_new/hdWGCNA/GO/GO_all_pathways.csv"),
            "Table S5E Module GO summary": read_csv("Result/figer_new/hdWGCNA/GO/module_summary.csv"),
            "Table S5F Notes": wgcna_notes(),
        },
    )


def long_from_signature_excel(rel: str) -> pd.DataFrame:
    path = p(rel)
    frames = []
    if not path.exists():
        return pd.DataFrame({"note": [f"Missing source file: {path}"]})
    xls = pd.ExcelFile(path)
    for sheet in xls.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        df["signature_sheet"] = sheet
        df["source_file"] = str(path)
        frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False)


def split_marker_genes(df: pd.DataFrame) -> pd.DataFrame:
    if "marker_genes" not in df.columns:
        return df
    rows = []
    for _, row in df.iterrows():
        genes = str(row["marker_genes"]).split(";")
        for gene in genes:
            r = row.drop(labels=["marker_genes"]).to_dict()
            r["gene"] = gene.strip()
            rows.append(r)
    return pd.DataFrame(rows)


def scoring_gene_set(
    df: pd.DataFrame,
    gene_col: str,
    gene_set: str,
    score_or_index: str,
    role: str,
    source_file: str,
    source_sheet: str = "",
) -> pd.DataFrame:
    out = df.copy()
    out = out.rename(columns={gene_col: "gene", "Source": "source"})
    if "source" not in out.columns:
        out["source"] = ""
    out = out[["gene", "source"]].dropna(subset=["gene"]).drop_duplicates()
    out.insert(0, "gene_set", gene_set)
    out.insert(1, "score_or_index", score_or_index)
    out["role_in_scoring"] = role
    out["source_file"] = str(p(source_file))
    out["source_sheet"] = source_sheet
    return out


def table6() -> None:
    gene_set_dir = REPO_ROOT / "data" / "gene_sets"
    s6a = pd.read_csv(gene_set_dir / "mpi_signatures.csv")[
        ["gene_set", "source", "gene"]
    ]
    s6b = pd.read_csv(gene_set_dir / "mmi_top200_signatures.csv")[
        ["gene_set", "source", "gene"]
    ]
    s6c = pd.read_csv(gene_set_dir / "apobec3a_positive_myeloid_signature.csv")[
        ["gene_set", "source", "gene"]
    ]
    write_xlsx(
        OUT / "Supplementary_Table_6_Signature_gene_sets_used_in_scoring.xlsx",
        {
            "Table S6A MPI signatures": s6a,
            "Table S6B MMI signatures": s6b,
            "Table S6C A3Apos myeloid sig": s6c,
        },
    )


def table7() -> None:
    global_genes = read_csv("Result/figer_new/monocle3/T1_Global_Trajectory_Dependent_Genes.csv")
    branch_final = read_csv("Result/figer_new/monocle3/T2_Model_Based_Hetero_Genes_Final.csv")
    branch_full = read_csv("Result/figer_new/monocle3/T2_Heterogeneous_Genes_Fate1_vs_Fate2.csv")
    xls = p("Result/figer_new/monocle3/T3_T5_Bifurcation_Key_Genes.xlsx")
    flags = []
    if xls.exists():
        for sheet, flag in [
            ("Transient_Drivers", "transient_candidate"),
            ("Early_Switch_Fate2", "early_switch_Fate2"),
            ("Early_Switch_Fate1", "early_switch_Fate1"),
        ]:
            df = pd.read_excel(xls, sheet_name=sheet)
            df["candidate_flag"] = flag
            df["source_sheet"] = sheet
            df["source_file"] = str(xls)
            flags.append(df)
    flags_df = pd.concat(flags, ignore_index=True, sort=False) if flags else pd.DataFrame()
    if not branch_final.empty and "gene_short_name" in branch_final.columns:
        branch_final["Fate_bias"] = branch_final["AUC_Diff"].apply(lambda x: "Fate2-biased" if x > 0 else ("Fate1-biased" if x < 0 else "No bias"))
        branch_final["pseudo_time_effect"] = branch_final["estimate"]
        branch_final["Branch_Fate_interaction_p"] = branch_final["p_value"]
        branch_final["Branch_Fate_interaction_q"] = branch_final["q_value"]
        branch_final["cluster_C1_C4"] = "Not assigned in available output"
    write_xlsx(
        OUT / "Supplementary_Table_7_Pseudotime_and_fate_branch_dynamic_genes.xlsx",
        {
            "Table S7A Global trajectory": global_genes,
            "Table S7B Fate branch final": branch_final,
            "Table S7C Fate branch full": branch_full,
            "Table S7D Switch candidates": flags_df,
            "Table S7E Notes": pd.DataFrame({"note": ["C1-C4 heatmap cluster assignment was not found as a machine-readable column in current outputs."]}),
        },
    )


def add_dataset_label(df: pd.DataFrame, label: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "analysis_label", label)
    return out


def table8() -> None:
    sheets = {
        "Table S8A clone13 DEG": read_csv("/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/ko/3A_ALL.csv"),
        "Table S8B clone37 DEG": read_csv("DEG_ALL.csv"),
        "Table S8C combined DEG": read_csv("GSE159677_Carotid_MainProject/GSVA/Differential_Expression_Results_13+37_vs_WT.csv"),
        "Table S8D clone13 GO up": read_csv("/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/ko/APOBEC3A_up_GO.BP.csv"),
        "Table S8E clone13 GO down": read_csv("/dsk2/data/C-to-U/APOBEC3A/00.mergeRawFq/ko/APOBEC3A_down_GO.BP.csv"),
        "Table S8F clone13 GSEA": read_csv("GSE159677_Carotid_MainProject/GSEA/A3A_KO/gsea_results.csv"),
        "Table S8G clone13 GSVA": read_csv("GSE159677_Carotid_MainProject/GSEA/A3A_KO/GSVA_diff_pathways_all.csv"),
        "Table S8H clone37 GSEA": read_csv("GSE159677_Carotid_MainProject/GSEA/A3A_KO/clone37/gsea_results_37.csv"),
        "Table S8I clone37 GSVA": read_csv("GSE159677_Carotid_MainProject/GSEA/A3A_KO/clone37/GSVA_diff_pathways_37.csv"),
        "Table S8J combined GSEA": read_csv("GSE159677_Carotid_MainProject/GSEA/13+37_gsea_results.csv"),
        "Table S8K Notes": pd.DataFrame(
            {
                "note": [
                    "The previous Table S8B clone13 DEG sig sheet was removed because 3A_total.csv is a significant DEG subset of 3A_ALL.csv; Table S8A retains the full clone13 DEG table with the sig column.",
                    "Single-cell Foam vs Macrophage GSEA/GSVA outputs were removed from this bulk RNA-seq supplementary table because they belong to the single-cell pathway analysis, not APOBEC3A-KO bulk RNA-seq.",
                    "Combined clone13+37 GSEA is included as Table S8J from 13+37_gsea_results.csv. A separate, clearly named combined clone13+37 GSVA result file was not found in the current outputs, so no combined GSVA sheet is included.",
                ]
            }
        ),
    }
    write_xlsx(OUT / "Supplementary_Table_8_Bulk_RNA_seq_differential_expression_and_enrichment.xlsx", sheets)


def parse_counts(field: str) -> list[int]:
    vals = []
    for x in str(field).split(","):
        try:
            vals.append(int(x))
        except ValueError:
            vals.append(0)
    vals = (vals + [0, 0, 0, 0])[:4]
    return vals


def ratio_for_ref(counts: list[int], ref: str) -> float:
    total = sum(counts)
    if total <= 0:
        return math.nan
    # JACUSA base fields are A,C,G,T. A-to-I appears as A>G on the
    # reference strand and T>C on the complementary strand.
    if ref == "A":
        return counts[2] / total
    if ref == "T":
        return counts[1] / total
    return math.nan


def jacusa_raw_to_atoi(raw_rel: str, label: str, out_rel: str) -> pd.DataFrame:
    out_path = p(out_rel)
    if out_path.exists() and out_path.stat().st_size > 0:
        cached = pd.read_csv(out_path)
        if "gene_name" not in cached.columns:
            cached["gene_name"] = "not_annotated_from_raw_JACUSA2"
        if "region_annotation" not in cached.columns:
            cached["region_annotation"] = "not_annotated_from_raw_JACUSA2"
        if cached["gene_name"].astype(str).str.contains("not_annotated_from_raw_JACUSA2", na=False).any():
            cached = annotate_sites_with_genes(cached)
            cached.to_csv(out_path, index=False)
        return cached
    raw_path = p(raw_rel)
    if not raw_path.exists():
        return pd.DataFrame({"note": [f"Missing source file: {raw_path}"]})
    rows = []
    headers = None
    with raw_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#contig"):
                headers = line.strip().lstrip("#").split("\t")
                continue
            if headers is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(headers):
                continue
            rec = dict(zip(headers, parts))
            ref = rec.get("ref", "")
            if ref not in {"A", "T"}:
                continue
            try:
                score = float(rec.get("score", "nan"))
            except ValueError:
                continue
            if not math.isfinite(score) or score < 2:
                continue
            wt_counts = [parse_counts(rec.get(f"bases1{i}", "")) for i in range(1, 4)]
            ko_counts = [parse_counts(rec.get(f"bases2{i}", "")) for i in range(1, 4)]
            if any(sum(x) < 10 for x in wt_counts) or any(sum(x) <= 10 for x in ko_counts):
                continue
            wt_ratios = [ratio_for_ref(x, ref) for x in wt_counts]
            ko_ratios = [ratio_for_ref(x, ref) for x in ko_counts]
            if any(math.isnan(x) for x in wt_ratios + ko_ratios):
                continue
            if len(set(wt_ratios)) == 1 or len(set(ko_ratios)) == 1:
                p_value = math.nan
            else:
                p_value = float(ttest_ind(wt_ratios, ko_ratios, equal_var=False, nan_policy="omit").pvalue)
            mean_wt = sum(wt_ratios) / len(wt_ratios)
            mean_ko = sum(ko_ratios) / len(ko_ratios)
            fold_change = mean_ko / mean_wt if mean_wt != 0 else math.inf
            rows.append(
                {
                    "editing_analysis": label,
                    "chr": rec.get("contig"),
                    "pos": int(rec.get("start", 0)),
                    "end_pos": int(rec.get("end", 0)),
                    "strand": rec.get("strand", "."),
                    "ref": ref,
                    "alt": "G" if ref == "A" else "C",
                    "editing_type": "A_to_I_candidate_A>G" if ref == "A" else "A_to_I_candidate_T>C_complement",
                    "gene_name": "not_annotated_from_raw_JACUSA2",
                    "region_annotation": "not_annotated_from_raw_JACUSA2",
                    "JACUSA2_score": score,
                    "wt_rep1_ratio": wt_ratios[0],
                    "wt_rep2_ratio": wt_ratios[1],
                    "wt_rep3_ratio": wt_ratios[2],
                    "ko_rep1_ratio": ko_ratios[0],
                    "ko_rep2_ratio": ko_ratios[1],
                    "ko_rep3_ratio": ko_ratios[2],
                    "mean_wt_ratio": mean_wt,
                    "mean_ko_ratio": mean_ko,
                    "delta_KO_minus_WT": mean_ko - mean_wt,
                    "fold_change_KO_vs_WT": fold_change,
                    "log2_fold_change": math.log2(fold_change) if fold_change > 0 and math.isfinite(fold_change) else math.nan,
                    "p_value": p_value,
                    "source_file": str(raw_path),
                }
            )
    df = pd.DataFrame(rows)
    df = annotate_sites_with_genes(df)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return df


def editing_table(rel: str, label: str) -> pd.DataFrame:
    df = read_csv(rel)
    if "seqnames" in df.columns:
        df = df.rename(columns={"seqnames": "chr", "start": "pos", "end": "end_pos", "score": "JACUSA2_score", "p_values": "p_value"})
        df.insert(0, "editing_analysis", label)
        df["strand"] = "not_encoded"
        df["alt_for_C_to_U_screen"] = df["ref"].map(lambda x: "T" if x == "C" else ("A" if x == "G" else pd.NA))
        df["region_annotation"] = "gene_overlap"
        df["delta_T_KO_minus_WT"] = df.get("mean_ko_T", pd.Series(dtype=float)) - df.get("mean_wt_T", pd.Series(dtype=float))
        df["delta_A_KO_minus_WT"] = df.get("mean_ko_A", pd.Series(dtype=float)) - df.get("mean_wt_A", pd.Series(dtype=float))
    return df


PUBLIC7_NGS_ROOT = Path(
    os.environ.get(
        "PUBLIC7_NGS_RESULT_ROOT",
        "/public7/DSC_Public7/DSC/DSC/20251004_NGS/result",
    )
)


def read_public7_excel_first_sheet(path: Path, label: str, extra: dict[str, object] | None = None) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame({"analysis_label": [label], "note": [f"Missing source file: {path}"]})
    df = pd.read_excel(path)
    df.insert(0, "analysis_label", label)
    if extra:
        for key, value in extra.items():
            df[key] = value
    df["source_file"] = str(path)
    df["source_sheet"] = "Sheet 1"
    return df


def read_public7_integrated(path: Path, label: str, score_cutoff: float, algorithm: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame({"analysis_label": [label], "note": [f"Missing source file: {path}"]})
    frames = []
    xls = pd.ExcelFile(path)
    for sheet in xls.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        df.insert(0, "analysis_label", label)
        df.insert(1, "editing_type", "C_to_U" if sheet.startswith("C_to_U") else ("A_to_I" if sheet.startswith("A_to_I") else "unknown"))
        df.insert(2, "statistical_test", "limma" if "limma" in sheet else ("ttest" if "ttest" in sheet else algorithm))
        df.insert(3, "score_cutoff", score_cutoff)
        df["source_file"] = str(path)
        df["source_sheet"] = sheet
        frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False)


def public7_s9_sources() -> pd.DataFrame:
    rows = [
        {
            "analysis_label": "M013 clone13",
            "root": str(PUBLIC7_NGS_ROOT / "M013_result_Score1_without_THP1"),
            "included_files": "Integrated_Results_Score1_Hybrid.xlsx; Target_Atherosclerosis_Lipid_SiteLevel_Report_Hybrid.xlsx; Final_True_Dimensions_Targets_Summary.xlsx; ALL_Pathways_Gene_Summary_Hybrid.xlsx",
        },
        {
            "analysis_label": "M037 clone37",
            "root": str(PUBLIC7_NGS_ROOT / "M037_result_Score1"),
            "included_files": "Integrated_Results_Score1_Hybrid.xlsx; Target_Atherosclerosis_Lipid_SiteLevel_Report_Hybrid.xlsx; Final_True_Dimensions_Targets_Summary.xlsx; ALL_Pathways_Gene_Summary_Hybrid.xlsx",
        },
        {
            "analysis_label": "combined all samples",
            "root": str(PUBLIC7_NGS_ROOT / "all_sample_result_PureLimma_no_THP1_WTcond1"),
            "included_files": "Score_1/Integrated_Results_Score1_PureLimma.xlsx; Score_1/Target_Atherosclerosis_Lipid_SiteLevel_Report_Limma_Score1.xlsx; Score_1/ALL_Pathways_Gene_Summary_Limma_Score1.xlsx; sample_metadata_from_manifest.csv",
        },
    ]
    return pd.DataFrame(rows)


def table9() -> None:
    m013_root = PUBLIC7_NGS_ROOT / "M013_result_Score1_without_THP1"
    m037_root = PUBLIC7_NGS_ROOT / "M037_result_Score1"
    combined_root = PUBLIC7_NGS_ROOT / "all_sample_result_PureLimma_no_THP1_WTcond1"

    m013_integrated = read_public7_integrated(
        m013_root / "Final_Results_Score1_Hybrid/Integrated_Results_Score1_Hybrid.xlsx",
        "M013 clone13",
        1.0,
        "Hybrid",
    )
    m037_integrated = read_public7_integrated(
        m037_root / "Final_Results_Score1_Hybrid/Integrated_Results_Score1_Hybrid.xlsx",
        "M037 clone37",
        1.0,
        "Hybrid",
    )
    combined_integrated = read_public7_integrated(
        combined_root / "Score_1/Final_Results_Score1_PureLimma/Integrated_Results_Score1_PureLimma.xlsx",
        "combined all samples",
        1.0,
        "PureLimma",
    )

    target_sites = pd.concat(
        [
            read_public7_excel_first_sheet(m013_root / "Target_Atherosclerosis_Lipid_SiteLevel_Report_Hybrid.xlsx", "M013 clone13", {"score_cutoff": 1.0, "algorithm": "Hybrid"}),
            read_public7_excel_first_sheet(m037_root / "Target_Atherosclerosis_Lipid_SiteLevel_Report_Hybrid.xlsx", "M037 clone37", {"score_cutoff": 1.0, "algorithm": "Hybrid"}),
            read_public7_excel_first_sheet(combined_root / "Score_1/Target_Atherosclerosis_Lipid_SiteLevel_Report_Limma_Score1.xlsx", "combined all samples", {"score_cutoff": 1.0, "algorithm": "PureLimma"}),
        ],
        ignore_index=True,
        sort=False,
    )
    target_summary = pd.concat(
        [
            read_public7_excel_first_sheet(m013_root / "Final_True_Dimensions_Targets_Summary.xlsx", "M013 clone13", {"score_cutoff": 1.0, "algorithm": "Hybrid"}),
            read_public7_excel_first_sheet(m037_root / "Final_True_Dimensions_Targets_Summary.xlsx", "M037 clone37", {"score_cutoff": 1.0, "algorithm": "Hybrid"}),
        ],
        ignore_index=True,
        sort=False,
    )
    pathway_summary = pd.concat(
        [
            read_public7_excel_first_sheet(m013_root / "ALL_Pathways_Gene_Summary_Hybrid.xlsx", "M013 clone13"),
            read_public7_excel_first_sheet(m037_root / "ALL_Pathways_Gene_Summary_Hybrid.xlsx", "M037 clone37"),
            read_public7_excel_first_sheet(combined_root / "Score_1/ALL_Pathways_Gene_Summary_Limma_Score1.xlsx", "combined all samples"),
        ],
        ignore_index=True,
        sort=False,
    )
    metadata = read_csv(str(combined_root / "sample_metadata_from_manifest.csv"))
    notes = pd.DataFrame(
        {
            "note": [
                "Supplementary Table 9 was regenerated from the public7 20251004_NGS result directories supplied by the user, replacing the older /public3/DSC/single_cell/jacusa/JACUSA2_output-derived table.",
                "M013 and M037 sheets use Score1 Hybrid integrated results, including both C_to_U and A_to_I results from limma and t-test source sheets.",
                "The combined all-sample sheets use the Score_1 PureLimma outputs under all_sample_result_PureLimma_no_THP1_WTcond1.",
                "Target site sheets collect the atherosclerosis/lipid site-level reports from the three analyses.",
                "Pathway summary collects ALL_Pathways_Gene_Summary outputs from M013, M037, and the combined Score_1 analysis.",
            ]
        }
    )
    write_xlsx(
        OUT / "Supplementary_Table_9_RNA_editing_results.xlsx",
        {
            "Table S9A M013 integrated": m013_integrated,
            "Table S9B M037 integrated": m037_integrated,
            "Table S9C combined integrated": combined_integrated,
            "Table S9D target lipid sites": target_sites,
            "Table S9E target summary": target_summary,
            "Table S9F pathway summary": pathway_summary,
            "Table S9G sample metadata": metadata,
            "Table S9H source files": public7_s9_sources(),
            "Table S9I Notes": notes,
        },
    )


def table10() -> None:
    """Spatial statistics for Xenium, GeoMx, and core-like region analyses."""
    notes = pd.DataFrame(
        {
            "note": [
                "Supplementary Table 10 collects spatial statistics from spatial/atherosclerosis_spatial_publication/results/tables.",
                "Xenium analyses are cell-level statistics from GSE315246, including myeloid composition, APOBEC3A detection, neighborhood enrichment, and myeloid pseudobulk severity DE.",
                "GeoMx analyses are ROI-level WTA statistics from GSE277170; GeoMx marker/program projections should not be interpreted as single-cell subtype calls.",
                "Core-like regions are algorithmically identified necrotic-core-like regions based on local LAM/Foam Cell enrichment; they should be described as core-like regions unless independently pathologically validated.",
                "ISG/core-like region sheets compare APOBEC3A-related or ISG-like myeloid program scores across Core-colocalized, Peri-core adjacent, and Other region classes.",
            ]
        }
    )
    sheets = {
        "Table S10A Dataset overview": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_dataset_qc_and_inclusion_Xenium_GeoMx_only.csv"),
        "Table S10B A3A detection": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_cross_platform_a3a_detection_summary_Xenium_GeoMx_only.csv"),
        "Table S10C Severity tests": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_severity_preference_tests_Xenium_GeoMx_only.csv"),
        "Table S10D Xenium MoMa comp": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_Xenium_MoMa_composition_by_severity_publication_counts_and_proportions.csv"),
        "Table S10E Xenium A3A celltype": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_xenium_a3a_by_celltype_severity.csv"),
        "Table S10F Xenium A3A substate": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_xenium_myeloid_substate_a3a_by_severity.csv"),
        "Table S10G Xenium programs": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_spatial_myeloid_program_a3a_detection_Xenium_GeoMx_only.csv"),
        "Table S10H Xenium neighbors": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_xenium_neighborhood_composition_k20.csv"),
        "Table S10I Xenium myeloid DE": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_xenium_myeloid_pseudobulk_severity_de.csv"),
        "Table S10J GeoMx A3A summary": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_geomx_Plaque_Adventitia_ROI_A3A_summary.csv"),
        "Table S10K GeoMx A3A tests": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_geomx_Plaque_Adventitia_ROI_A3A_tests.csv"),
        "Table S10L GeoMx programs": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_geomx_ROI_myeloid_program_summary.csv"),
        "Table S10M GeoMx program tests": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_geomx_A3A_high_vs_negative_program_tests.csv"),
        "Table S10N GeoMx marker corr": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_geomx_APOBEC3A_marker_spearman.csv"),
        "Table S10O Core FOV summary": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_necrotic_core_like_FOV_summary_unified_area.csv"),
        "Table S10P Core cluster summary": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_necrotic_core_like_cluster_summary_unified_area.csv"),
        "Table S10Q Core area summary": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_necrotic_core_like_area_distribution_summary_unified_area.csv"),
        "Table S10R Core thresholds": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_necrotic_core_like_threshold_diagnostics.csv"),
        "Table S10S A3A core-like": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_APOBEC3A_by_necrotic_core_like_region_summary.csv"),
        "Table S10T ISG region summary": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_core_region3_ISG_from_scratch_ALL_FOV_FINAL_myeloid_ISG_score_spatial_region_bar_summary_ALL_FOVS.csv"),
        "Table S10U ISG region Kruskal": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_core_region3_ISG_from_scratch_ALL_FOV_FINAL_myeloid_ISG_score_spatial_region_kruskal_ALL_FOVS.csv"),
        "Table S10V ISG region pairwise": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_core_region3_ISG_from_scratch_ALL_FOV_FINAL_myeloid_ISG_score_spatial_region_pairwise_wilcox_BH_ALL_FOVS.csv"),
        "Table S10W ISG region per FOV": read_csv("spatial/atherosclerosis_spatial_publication/results/tables/myeloid_story_core_region3_ISG_from_scratch_ALL_FOV_FINAL_myeloid_spatial_region_summary_PER_FOV.csv"),
        "Table S10X Notes": notes,
    }
    write_xlsx(OUT / "Supplementary_Table_10_Spatial_statistics.xlsx", sheets)


def main() -> None:
    table2()
    table3()
    table4()
    table5()
    table6()
    table7()
    table8()
    table9()
    table10()
    manifest = []
    for xlsx in sorted(OUT.glob("Supplementary_Table_*.xlsx")):
        xl = pd.ExcelFile(xlsx)
        for sheet in xl.sheet_names:
            header = pd.read_excel(xlsx, sheet_name=sheet, nrows=0)
            n_rows = len(pd.read_excel(xlsx, sheet_name=sheet, usecols=[0]))
            manifest.append({"file": str(xlsx), "sheet": sheet, "n_rows": n_rows, "n_columns": len(header.columns)})
    write_xlsx(OUT / "Supplementary_Tables_manifest.xlsx", {"Manifest": pd.DataFrame(manifest)})
    print(f"Wrote supplementary tables to {OUT}")


if __name__ == "__main__":
    main()
