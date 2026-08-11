#!/usr/bin/env python3
"""Build metadata and QC manifests for raw spatial archives.

The script is intentionally lightweight: it does not unpack full archives to
disk and does not create large count matrices. It extracts sample-level GEO
metadata, file inventories, and compact QC summaries from the local tar files.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import re
import subprocess
import tarfile
import time
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(os.environ.get("AST_ROOT", "workflow/04_spatial")).resolve()
SPATIAL_ROOT = Path(
    os.environ.get("AST_SPATIAL_ROOT", "/public3/DSC/single_cell/spatial")
).resolve()
OUTPUT_ROOT = Path(os.environ.get("AST_OUTPUT_ROOT", ROOT / "results")).resolve()
TABLE_DIR = OUTPUT_ROOT / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
SOFT_CACHE_DIR = TABLE_DIR / "raw_geo_soft_cache"
SOFT_CACHE_DIR.mkdir(parents=True, exist_ok=True)

GSE_INFO = {
    "GSE277441": {
        "archive": SPATIAL_ROOT / "GSE277441_RAW.tar",
        "platform_hint": "NanoString CosMx SMI",
        "role": "cell-level validation of atherosclerosis severity",
    },
    "GSE277170": {
        "archive": SPATIAL_ROOT / "GSE277170_RAW.tar",
        "platform_hint": "NanoString GeoMx DSP",
        "role": "ROI-level validation of the same coronary artery progression study",
    },
    "GSE283269": {
        "archive": SPATIAL_ROOT / "GSE283269_RAW.tar",
        "platform_hint": "10x Visium FFPE",
        "role": "early diffuse intimal thickening coronary artery spatial reference",
    },
}


def cache_path_for_url(url: str) -> Path:
    acc_match = re.search(r"[?&]acc=([^&]+)", url)
    targ_match = re.search(r"[?&]targ=([^&]+)", url)
    acc = acc_match.group(1) if acc_match else "geo"
    targ = targ_match.group(1) if targ_match else "self"
    return SOFT_CACHE_DIR / f"{acc}.{targ}.soft"


def fetch_text(url: str) -> str:
    cache_path = cache_path_for_url(url)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return cache_path.read_text(errors="replace")

    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=45) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            cache_path.write_text(text)
            return text
        except Exception as exc:  # pragma: no cover - network-dependent fallback
            last_error = exc
            time.sleep(1 + attempt)

    try:
        proc = subprocess.run(
            ["curl", "-L", "--retry", "3", "--retry-delay", "2", "-sS", url],
            check=True,
            capture_output=True,
            text=True,
        )
        cache_path.write_text(proc.stdout)
        return proc.stdout
    except Exception as exc:  # pragma: no cover - only reached when network fails
        raise RuntimeError(f"Failed to fetch {url}") from (last_error or exc)


def parse_soft_samples(gse: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    series_url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        f"acc={gse}&targ=self&form=text&view=full"
    )
    sample_url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        f"acc={gse}&targ=gsm&form=text&view=full"
    )
    series_text = fetch_text(series_url)
    sample_text = fetch_text(sample_url)

    series = {"gse": gse}
    for line in series_text.splitlines():
        if line.startswith("!Series_title = "):
            series["series_title"] = line.split(" = ", 1)[1]
        elif line.startswith("!Series_summary = "):
            series.setdefault("series_summary", "")
            series["series_summary"] += (" " if series["series_summary"] else "") + line.split(" = ", 1)[1]
        elif line.startswith("!Series_overall_design = "):
            series["overall_design"] = line.split(" = ", 1)[1]
        elif line.startswith("!Series_platform_id = "):
            series.setdefault("platform_ids", [])
            series["platform_ids"].append(line.split(" = ", 1)[1])
    if isinstance(series.get("platform_ids"), list):
        series["platform_ids"] = ";".join(series["platform_ids"])

    samples: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    for line in sample_text.splitlines():
        if line.startswith("^SAMPLE = "):
            if current:
                samples.append(current)
            current = {"gse": gse, "gsm": line.split(" = ", 1)[1]}
        elif current is None:
            continue
        elif line.startswith("!Sample_title = "):
            current["title"] = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_source_name_ch1 = "):
            current["source_name"] = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = "):
            payload = line.split(" = ", 1)[1]
            if ":" in payload:
                key, value = payload.split(":", 1)
                current[key.strip().lower().replace(" ", "_")] = value.strip()
            else:
                current.setdefault("characteristics", "")
                current["characteristics"] += ("; " if current["characteristics"] else "") + payload
    if current:
        samples.append(current)
    return series, samples


def archive_file_manifest(gse: str, archive: Path) -> list[dict[str, str]]:
    rows = []
    with tarfile.open(archive, "r") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            filename = Path(member.name).name
            gsm_match = re.match(r"(GSM\d+)", filename)
            rows.append(
                {
                    "gse": gse,
                    "gsm": gsm_match.group(1) if gsm_match else "",
                    "filename": filename,
                    "size_bytes": str(member.size),
                    "suffix": filename.split("_", 1)[1] if "_" in filename else filename,
                }
            )
    return rows


def read_tar_gzip_text(tf: tarfile.TarFile, member_name: str) -> str:
    fh = tf.extractfile(member_name)
    if fh is None:
        return ""
    data = fh.read()
    if member_name.endswith(".gz"):
        return gzip.decompress(data).decode("utf-8", errors="replace")
    return data.decode("utf-8", errors="replace")


def cosmx_qc(gse: str, archive: Path) -> list[dict[str, str]]:
    rows = []
    with tarfile.open(archive, "r") as tf:
        metadata_members = [m.name for m in tf.getmembers() if m.isfile() and m.name.endswith("metadata_file.csv.gz")]
        for member in metadata_members:
            text = read_tar_gzip_text(tf, member)
            reader = csv.DictReader(io.StringIO(text))
            n_cells = 0
            fovs = set()
            total_counts = 0.0
            total_features = 0.0
            for row in reader:
                n_cells += 1
                fovs.add(row.get("fov", ""))
                total_counts += float(row.get("nCount_RNA", 0) or 0)
                total_features += float(row.get("nFeature_RNA", 0) or 0)
            gsm = re.match(r"(GSM\d+)", Path(member).name).group(1)
            rows.append(
                {
                    "gse": gse,
                    "gsm": gsm,
                    "platform": "CosMx",
                    "n_cells": n_cells,
                    "n_fovs": len([f for f in fovs if f]),
                    "mean_nCount_RNA": round(total_counts / n_cells, 3) if n_cells else "",
                    "mean_nFeature_RNA": round(total_features / n_cells, 3) if n_cells else "",
                }
            )
    return rows


def geomx_dcc_qc(gse: str, archive: Path) -> list[dict[str, str]]:
    rows = []
    with tarfile.open(archive, "r") as tf:
        dcc_members = [m.name for m in tf.getmembers() if m.isfile() and m.name.endswith(".dcc.gz")]
        for member in dcc_members:
            text = read_tar_gzip_text(tf, member)
            gsm = re.match(r"(GSM\d+)", Path(member).name).group(1)
            metrics: dict[str, str] = {"gse": gse, "gsm": gsm, "platform": "GeoMx"}
            in_ngs = False
            in_code = False
            n_probes = 0
            total_probe_counts = 0
            for line in text.splitlines():
                line = line.strip()
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
                        metrics[key] = value.strip('"')
                if in_code and "," in line:
                    _, value = line.split(",", 1)
                    try:
                        total_probe_counts += int(float(value))
                        n_probes += 1
                    except ValueError:
                        pass
            metrics["n_probes"] = n_probes
            metrics["total_probe_counts"] = total_probe_counts
            rows.append(metrics)
    return rows


def visium_metrics_qc(gse: str, archive: Path) -> list[dict[str, str]]:
    rows = []
    with tarfile.open(archive, "r") as tf:
        metric_members = [m.name for m in tf.getmembers() if m.isfile() and m.name.endswith("metrics_summary.csv.gz")]
        for member in metric_members:
            text = read_tar_gzip_text(tf, member)
            reader = csv.DictReader(io.StringIO(text))
            row = next(reader)
            gsm = re.match(r"(GSM\d+)", Path(member).name).group(1)
            out = {
                "gse": gse,
                "gsm": gsm,
                "platform": "Visium",
                "sample_id": row.get("Sample ID", ""),
                "spots_under_tissue": row.get("Number of Spots Under Tissue", ""),
                "number_of_reads": row.get("Number of Reads", ""),
                "mean_reads_per_spot": row.get("Mean Reads per Spot", ""),
                "median_genes_per_spot": row.get("Median Genes per Spot", ""),
                "median_umi_counts_per_spot": row.get("Median UMI Counts per Spot", ""),
                "genes_detected": row.get("Genes Detected", ""),
                "sequencing_saturation": row.get("Sequencing Saturation", ""),
            }
            if "sample2d1" in member.lower():
                out["published_qc_note"] = "removed from downstream analysis according to GEO series design"
            else:
                out["published_qc_note"] = ""
            rows.append(out)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = sorted({k for row in rows for k in row.keys()})
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    series_rows = []
    sample_rows = []
    file_rows = []
    qc_rows = []

    for gse, info in GSE_INFO.items():
        series, samples = parse_soft_samples(gse)
        series.update(
            {
                "platform_hint": info["platform_hint"],
                "role": info["role"],
                "archive": str(info["archive"]),
                "archive_size_gb": round(info["archive"].stat().st_size / 1024**3, 3),
            }
        )
        series_rows.append(series)
        sample_rows.extend(samples)
        file_rows.extend(archive_file_manifest(gse, info["archive"]))

        if gse == "GSE277441":
            qc_rows.extend(cosmx_qc(gse, info["archive"]))
        elif gse == "GSE277170":
            qc_rows.extend(geomx_dcc_qc(gse, info["archive"]))
        elif gse == "GSE283269":
            qc_rows.extend(visium_metrics_qc(gse, info["archive"]))

    write_csv(TABLE_DIR / "raw_gse_series_metadata.csv", series_rows)
    write_csv(TABLE_DIR / "raw_gse_sample_metadata.csv", sample_rows)
    write_csv(TABLE_DIR / "raw_archive_file_manifest.csv", file_rows)
    write_csv(TABLE_DIR / "raw_archive_qc_summary.csv", qc_rows)

    sample_by_gse = defaultdict(int)
    files_by_gse = defaultdict(int)
    for row in sample_rows:
        sample_by_gse[row["gse"]] += 1
    for row in file_rows:
        files_by_gse[row["gse"]] += 1
    summary = []
    for gse, info in GSE_INFO.items():
        summary.append(
            {
                "gse": gse,
                "platform_hint": info["platform_hint"],
                "n_geo_samples": sample_by_gse[gse],
                "n_raw_files": files_by_gse[gse],
                "archive_size_gb": round(info["archive"].stat().st_size / 1024**3, 3),
            }
        )
    write_csv(TABLE_DIR / "raw_archive_sample_summary.csv", summary)


if __name__ == "__main__":
    main()
