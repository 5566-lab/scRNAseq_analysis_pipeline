#!/usr/bin/env python3
"""Resolve GSE277170 GeoMx WTA DCC counts through the Bruker WTA PKC."""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import tarfile
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(os.environ.get("AST_ROOT", "workflow/04_spatial")).resolve()
SPATIAL_ROOT = Path(
    os.environ.get("AST_SPATIAL_ROOT", "/public3/DSC/single_cell/spatial")
).resolve()
OUTPUT_ROOT = Path(os.environ.get("AST_OUTPUT_ROOT", ROOT / "results")).resolve()
TABLE_DIR = OUTPUT_ROOT / "tables"
REF_DIR = ROOT / "reference" / "geomx"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
REF_DIR.mkdir(parents=True, exist_ok=True)

PKC_URL = "https://brukerspatialbiology.com/wp-content/uploads/Hs_R_NGS_WTA_v1.0.pkc_.zip"
PKC_ZIP = REF_DIR / "Hs_R_NGS_WTA_v1.0.pkc_.zip"
PKC_PATH = REF_DIR / "Hs_R_NGS_WTA_v1.0.pkc"
DCC_TAR = SPATIAL_ROOT / "GSE277170_RAW.tar"
META_CSV = TABLE_DIR / "raw_gse_sample_metadata.csv"

APOBEC_TARGETS = [
    "APOBEC3A",
    "APOBEC3A_B",
    "APOBEC3B",
    "APOBEC1",
    "APOBEC2",
    "APOBEC3C",
    "APOBEC3D",
    "APOBEC3F",
    "APOBEC3G",
    "APOBEC3H",
    "APOBEC4",
]

MARKER_TARGETS = [
    "CD45",
    "PTPRC",
    "CD4",
    "CD14",
    "FCGR3A",
    "ITGAM",
    "CSF1R",
    "LST1",
    "AIF1",
    "CD68",
    "CD163",
    "MRC1",
    "C1QA",
    "C1QB",
    "C1QC",
    "SPP1",
    "TREM2",
    "APOE",
    "LIPA",
    "MSR1",
    "MARCO",
    "HLA-DRA",
    "CXCL9",
    "CXCL10",
    "ISG15",
    "IFI6",
    "IFI27",
    "IFIT1",
    "IFIT2",
    "IFIT3",
]


def download_pkc() -> None:
    if PKC_PATH.exists() and PKC_PATH.stat().st_size > 0:
        return
    req = urllib.request.Request(PKC_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        PKC_ZIP.write_bytes(resp.read())
    with zipfile.ZipFile(PKC_ZIP) as zf:
        zf.extract("Hs_R_NGS_WTA_v1.0.pkc", REF_DIR)


def load_pkc() -> tuple[dict[str, str], list[dict[str, str]], set[str]]:
    download_pkc()
    data = json.loads(PKC_PATH.read_text())
    rts_to_target: dict[str, str] = {}
    probe_rows: list[dict[str, str]] = []
    target_names = set()
    for target in data["Targets"]:
        target_name = target["DisplayName"]
        target_names.add(target_name)
        for probe in target["Probes"]:
            rts_id = probe["RTS_ID"]
            rts_to_target[rts_id] = target_name
            probe_rows.append(
                {
                    "target": target_name,
                    "rts_id": rts_id,
                    "probe_id": str(probe.get("ProbeID", "")),
                    "probe_display_name": probe.get("DisplayName", ""),
                    "code_class": target.get("CodeClass", ""),
                    "gene_id": ";".join(map(str, probe.get("GeneID", []))),
                    "accession": ";".join(probe.get("Accession", [])),
                    "genome_coordinates": ";".join(probe.get("GenomeCoordinates", [])),
                }
            )
    return rts_to_target, probe_rows, target_names


def read_metadata() -> dict[str, dict[str, str]]:
    with META_CSV.open(newline="") as fh:
        rows = [row for row in csv.DictReader(fh) if row.get("gse") == "GSE277170"]
    return {row["gsm"]: row for row in rows}


def parse_dcc_text(text: str, rts_to_target: dict[str, str]) -> tuple[dict[str, str], dict[str, int], int, int]:
    ngs_metrics: dict[str, str] = {}
    target_counts: defaultdict[str, int] = defaultdict(int)
    total_probe_counts = 0
    n_detected_probes = 0
    in_ngs = False
    in_code = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line == "<NGS_Processing_Attributes>":
            in_ngs = True
            continue
        if line == "</NGS_Processing_Attributes>":
            in_ngs = False
            continue
        if line == "<Code_Summary>":
            in_code = True
            continue
        if line == "</Code_Summary>":
            in_code = False
            continue
        if in_ngs and "," in line:
            key, value = line.split(",", 1)
            if key in {"Raw", "Trimmed", "Stitched", "Aligned", "umiQ30", "rtsQ30"}:
                ngs_metrics[key] = value.strip('"')
        if in_code and "," in line:
            rts_id, value = line.split(",", 1)
            count = int(float(value))
            target = rts_to_target.get(rts_id)
            if target is None:
                raise ValueError(f"RTS_ID {rts_id} was not found in {PKC_PATH}")
            target_counts[target] += count
            total_probe_counts += count
            n_detected_probes += 1
    return ngs_metrics, dict(target_counts), total_probe_counts, n_detected_probes


def main() -> None:
    rts_to_target, probe_rows, target_names = load_pkc()
    metadata = read_metadata()
    requested_targets = sorted(set(APOBEC_TARGETS + MARKER_TARGETS) & target_names)

    with (TABLE_DIR / "geomx_wta_probe_map.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(probe_rows[0].keys()))
        writer.writeheader()
        writer.writerows(probe_rows)

    roi_rows = []
    marker_rows = []
    with tarfile.open(DCC_TAR, "r") as tf:
        dcc_members = sorted(member.name for member in tf.getmembers() if member.isfile() and member.name.endswith(".dcc.gz"))
        for member in dcc_members:
            gsm_match = re.match(r"(GSM\d+)", Path(member).name)
            if gsm_match is None:
                continue
            gsm = gsm_match.group(1)
            fh = tf.extractfile(member)
            if fh is None:
                continue
            text = gzip.decompress(fh.read()).decode("utf-8", errors="replace")
            ngs_metrics, target_counts, total_probe_counts, n_detected_probes = parse_dcc_text(text, rts_to_target)
            meta = metadata.get(gsm, {})
            is_ntc = meta.get("title", "").lower() == "no template control"
            is_negative_control = meta.get("localisation", "").lower() == "negative control"

            row = {
                "gsm": gsm,
                "dcc_file": Path(member).name,
                "title": meta.get("title", ""),
                "grade": meta.get("grade", ""),
                "localisation": meta.get("localisation", ""),
                "subset": meta.get("subset", ""),
                "roi": meta.get("roi", ""),
                "is_no_template_control": is_ntc,
                "is_negative_control_roi": is_negative_control,
                "total_probe_counts": total_probe_counts,
                "n_detected_probes": n_detected_probes,
                "n_detected_targets": len(target_counts),
            }
            row.update(ngs_metrics)
            for target in APOBEC_TARGETS:
                row[f"{target}_count"] = target_counts.get(target, 0)
            row["APOBEC3A_positive"] = target_counts.get("APOBEC3A", 0) > 0
            roi_rows.append(row)

            for target in requested_targets:
                marker_rows.append(
                    {
                        "gsm": gsm,
                        "target": target,
                        "count": target_counts.get(target, 0),
                        "grade": meta.get("grade", ""),
                        "localisation": meta.get("localisation", ""),
                        "subset": meta.get("subset", ""),
                        "roi": meta.get("roi", ""),
                        "is_no_template_control": is_ntc,
                        "is_negative_control_roi": is_negative_control,
                    }
                )

    roi_path = TABLE_DIR / "geomx_wta_apobec_roi_counts.csv"
    with roi_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(roi_rows[0].keys()))
        writer.writeheader()
        writer.writerows(roi_rows)

    marker_path = TABLE_DIR / "geomx_wta_marker_roi_counts.csv"
    with marker_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(marker_rows[0].keys()))
        writer.writeheader()
        writer.writerows(marker_rows)

    n_total = len(roi_rows)
    n_ntc = sum(row["is_no_template_control"] for row in roi_rows)
    n_tissue = n_total - n_ntc
    n_pos = sum(row["APOBEC3A_positive"] and not row["is_no_template_control"] for row in roi_rows)
    print(f"Wrote {roi_path}")
    print(f"Mapped {len(rts_to_target)} RTS_ID probes to {len(target_names)} GeoMx WTA targets")
    print(f"APOBEC3A detected in {n_pos}/{n_tissue} non-NTC GeoMx ROI/DCC files")


if __name__ == "__main__":
    main()
