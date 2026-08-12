#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JACUSA2_ALL_no_THP1_ordered.py

目的：
1. 排除所有文件夹名或 BAM 路径中含有 THP_1 的样品，即去掉 clone13#2 批次。
2. 固定 JACUSA2 输入顺序：cond1 = WT，cond2 = KO。
3. 将 Python 实际传给 JACUSA2 的样品顺序写入 jacusa_sample_manifest.tsv。
4. R 下游必须读取该 manifest，不再手写/猜测样品顺序。
"""

import csv
import glob
import os
import re
import subprocess
import time
import multiprocessing
from dataclasses import dataclass
from typing import List, Tuple

# ====================================================================
# 1. 软件与参考基因组/注释文件路径
# ====================================================================
JACUSA_JAR_PATH = os.environ.get(
    "RNA_EDITING_JACUSA_JAR",
    "/dell_1/dsc/RNA_seq/data/jacusa/JACUSA2-2.0.4/JACUSA_v2.0.4.jar",
)
REF_FASTA_PATH = os.environ.get(
    "RNA_EDITING_REFERENCE_FASTA",
    "/dell_1/dsc/RNA_seq/data/index/GRCh38.dna.Ensembl_release_104.fa",
)
GTF_PATH = os.environ.get(
    "RNA_EDITING_ANNOTATION_GTF",
    "/dell_1/dsc/RNA_seq/data/index/Homo_sapiens.GRCh38.104.gtf",
)
FEATURECOUNTS_PATH = os.environ.get("RNA_EDITING_FEATURECOUNTS", "featureCounts")

# ====================================================================
# 2. 输入输出路径
# ====================================================================
ALIGNMENT_BASE_DIR = os.environ.get(
    "RNA_EDITING_ALIGNMENT_DIR",
    "/public7/DSC_Public7/DSC/DSC/20251004_NGS/Aligned",
)

# 新目录：避免混用旧的 7vs7 / cond1=KO 缓存和结果
JACUSA_OUTPUT_DIR = os.environ.get(
    "RNA_EDITING_JACUSA_OUTPUT_DIR",
    "results/rna_editing/jacusa2_combined",
)
JACUSA_OUT_FILE = os.path.join(JACUSA_OUTPUT_DIR, "Aggregated_WTcond1_KOcond2_no_THP1.out")
SAMPLE_MANIFEST_FILE = os.path.join(JACUSA_OUTPUT_DIR, "jacusa_sample_manifest.tsv")

EXPRESSION_OUTPUT_DIR = os.environ.get(
    "RNA_EDITING_EXPRESSION_OUTPUT_DIR",
    "results/bulk_rnaseq/counts",
)
EXPRESSION_OUT_FILE = os.path.join(EXPRESSION_OUTPUT_DIR, "gene_counts_all_no_THP1.txt")

# ====================================================================
# 3. 分组规则
# ====================================================================
# 先排除，再分组。只要路径或文件夹名中含 THP_1，就不进入 WT/KO。
EXCLUDE_KEYWORDS = ["THP_1"]

# WT：clone13 的 M0-WT-* 和 clone37 的 WT.*
WT_PREFIXES = ("M0-WT", "WT.")

# KO：clone13 的 M0-13-* 和 clone37 的 KO.*
KO_PREFIXES = ("M0-13", "KO.")

THREADS = os.environ.get("RNA_EDITING_THREADS", "8")
MIN_COVERAGE = os.environ.get("RNA_EDITING_MIN_COVERAGE", "10")
CONTRAST = os.environ.get("RNA_EDITING_CONTRAST", "combined")
VALIDATE_ONLY = os.environ.get("RNA_EDITING_VALIDATE_ONLY", "0").lower() in {
    "1", "true", "yes", "y"
}


@dataclass(frozen=True)
class SampleRecord:
    sample_name: str
    condition: str      # WT or KO
    batch: str          # clone13 or clone37
    bam_path: str


def natural_key(text: str):
    """自然排序：WT.2 排在 WT.10 前。"""
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", text)]


def infer_batch(sample_name: str) -> str:
    if sample_name.startswith(("M0-WT", "M0-13")):
        return "clone13"
    if sample_name.startswith(("WT.", "KO.")):
        return "clone37"
    return "Unknown_Batch"


def classify_sample(folder_name: str, bam_path: str) -> str:
    text = f"{folder_name} {bam_path}"
    if any(key in text for key in EXCLUDE_KEYWORDS):
        return "EXCLUDE"
    if folder_name.startswith(WT_PREFIXES):
        return "WT"
    if folder_name.startswith(KO_PREFIXES):
        return "KO"
    return "UNKNOWN"


def collect_samples(base_dir: str) -> Tuple[List[SampleRecord], List[SampleRecord], List[Tuple[str, str]]]:
    if not os.path.exists(base_dir):
        raise FileNotFoundError(f"比对目录不存在: {base_dir}")

    wt_samples: List[SampleRecord] = []
    ko_samples: List[SampleRecord] = []
    skipped: List[Tuple[str, str]] = []

    for folder_name in sorted(os.listdir(base_dir), key=natural_key):
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        bam_files = sorted(glob.glob(os.path.join(folder_path, "*.sorted.bam")), key=natural_key)
        if not bam_files:
            skipped.append((folder_name, "未找到 *.sorted.bam"))
            continue

        bam_path = bam_files[0]
        group = classify_sample(folder_name, bam_path)

        if group == "EXCLUDE":
            skipped.append((folder_name, "排除：含 THP_1，去掉 clone13#2 批次"))
            continue
        if group == "UNKNOWN":
            skipped.append((folder_name, "跳过：无法按规则识别为 WT 或 KO"))
            continue

        rec = SampleRecord(
            sample_name=folder_name,
            condition=group,
            batch=infer_batch(folder_name),
            bam_path=bam_path,
        )
        if group == "WT":
            wt_samples.append(rec)
        else:
            ko_samples.append(rec)

    # 这里的排序就是最终传给 JACUSA2 的真实顺序；R 读取 manifest 后完全照这个顺序对应。
    wt_samples = sorted(wt_samples, key=lambda x: natural_key(x.sample_name))
    ko_samples = sorted(ko_samples, key=lambda x: natural_key(x.sample_name))
    return wt_samples, ko_samples, skipped


def validate_samples(wt_samples: List[SampleRecord], ko_samples: List[SampleRecord]) -> None:
    if not wt_samples or not ko_samples:
        raise RuntimeError(f"WT 或 KO 样本为空：WT={len(wt_samples)}, KO={len(ko_samples)}")

    bad = [s for s in wt_samples + ko_samples if "THP_1" in s.sample_name or "THP_1" in s.bam_path]
    if bad:
        raise RuntimeError("仍检测到 THP_1 样本未被排除：\n" + "\n".join(str(x) for x in bad))

    expected = 6 if CONTRAST == "combined" else 3
    if len(wt_samples) != expected or len(ko_samples) != expected:
        print(f"⚠️ 警告：{CONTRAST} 不是预期的 {expected} WT vs {expected} KO：WT={len(wt_samples)}, KO={len(ko_samples)}")
        print("   如果你的目录中样本数确实改变，可忽略；否则请检查分组规则。")


def write_manifest(wt_samples: List[SampleRecord], ko_samples: List[SampleRecord], out_file: str) -> None:
    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    rows = []
    for i, rec in enumerate(wt_samples, start=1):
        rows.append({
            "JACUSA_Condition": "cond1",
            "JACUSA_Replicate": f"rep{i}",
            "Condition": "WT",
            "Sample_Name": rec.sample_name,
            "Batch": rec.batch,
            "BAM": rec.bam_path,
            "Order_In_JACUSA_Command": i,
        })
    for i, rec in enumerate(ko_samples, start=1):
        rows.append({
            "JACUSA_Condition": "cond2",
            "JACUSA_Replicate": f"rep{i}",
            "Condition": "KO",
            "Sample_Name": rec.sample_name,
            "Batch": rec.batch,
            "BAM": rec.bam_path,
            "Order_In_JACUSA_Command": i,
        })

    with open(out_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def print_order(wt_samples: List[SampleRecord], ko_samples: List[SampleRecord], skipped: List[Tuple[str, str]]) -> None:
    print("\n========== JACUSA2 样本顺序确认 ==========")
    print("cond1 = WT，传给 JACUSA2 的第一个组：")
    for i, rec in enumerate(wt_samples, start=1):
        print(f"  cond1_rep{i}\tWT\t{rec.sample_name}\t{rec.batch}\t{rec.bam_path}")

    print("\ncond2 = KO，传给 JACUSA2 的第二个组：")
    for i, rec in enumerate(ko_samples, start=1):
        print(f"  cond2_rep{i}\tKO\t{rec.sample_name}\t{rec.batch}\t{rec.bam_path}")

    if skipped:
        print("\n被排除/跳过的目录：")
        for name, reason in skipped:
            print(f"  {name}\t{reason}")
    print("==========================================\n")


def check_required_files(paths: List[str]) -> None:
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError("以下必需文件不存在：\n" + "\n".join(missing))


def run_jacusa_aggregated(wt_samples: List[SampleRecord], ko_samples: List[SampleRecord], threads: str, min_coverage: str) -> None:
    print("\n[进程 A] >>> 启动 JACUSA2 联合 RNA 编辑分析")
    os.makedirs(JACUSA_OUTPUT_DIR, exist_ok=True)

    wt_bams = [x.bam_path for x in wt_samples]
    ko_bams = [x.bam_path for x in ko_samples]

    check_required_files([JACUSA_JAR_PATH, REF_FASTA_PATH] + wt_bams + ko_bams)

    jacusa_cmd = [
        "java", "-Xmx50g", "-jar", JACUSA_JAR_PATH, "call-2",
        "-c", min_coverage,
        "-q", "25",
        "-m", "30",
        "-p", threads,
        "-r", JACUSA_OUT_FILE,
        "-R", REF_FASTA_PATH,
        ",".join(wt_bams),   # cond1 = WT
        ",".join(ko_bams),   # cond2 = KO
    ]

    print("[进程 A] JACUSA2 输出:", JACUSA_OUT_FILE)
    print("[进程 A] manifest:", SAMPLE_MANIFEST_FILE)

    try:
        subprocess.run(jacusa_cmd, check=True)
        print("\n[进程 A] >>> JACUSA2 执行成功！")
    except subprocess.CalledProcessError as e:
        print(f"\n[进程 A] JACUSA2 执行失败: {e}")
        raise


def run_featurecounts(all_samples: List[SampleRecord], threads: str) -> None:
    print("\n[进程 B] >>> 启动 featureCounts 基因表达定量分析")
    os.makedirs(EXPRESSION_OUTPUT_DIR, exist_ok=True)

    all_bams = [x.bam_path for x in all_samples]
    check_required_files([GTF_PATH] + all_bams)

    fc_cmd = [
        FEATURECOUNTS_PATH,
        "-T", threads,
        "-p",
        "-t", "exon",
        "-g", "gene_id",
        "-a", GTF_PATH,
        "-o", EXPRESSION_OUT_FILE,
    ] + all_bams

    try:
        subprocess.run(fc_cmd, check=True)
        print("\n[进程 B] >>> featureCounts 执行成功！")
    except FileNotFoundError:
        print(f"\n[进程 B] 错误：找不到命令 {FEATURECOUNTS_PATH}")
        raise
    except subprocess.CalledProcessError as e:
        print(f"\n[进程 B] featureCounts 执行失败: {e}")
        raise


def main() -> None:
    start_time = time.time()
    print("========== 联合分析任务开始 ==========")
    print(f"当前时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    wt_samples, ko_samples, skipped = collect_samples(ALIGNMENT_BASE_DIR)
    if CONTRAST not in {"combined", "clone13", "clone37"}:
        raise ValueError(f"不支持的 RNA_EDITING_CONTRAST: {CONTRAST}")
    if CONTRAST != "combined":
        wt_samples = [sample for sample in wt_samples if sample.batch == CONTRAST]
        ko_samples = [sample for sample in ko_samples if sample.batch == CONTRAST]
    validate_samples(wt_samples, ko_samples)
    write_manifest(wt_samples, ko_samples, SAMPLE_MANIFEST_FILE)
    print_order(wt_samples, ko_samples, skipped)

    all_samples = wt_samples + ko_samples

    print("汇总信息:")
    print(f" - WT 样本数量 / cond1: {len(wt_samples)}")
    print(f" - KO 样本数量 / cond2: {len(ko_samples)}")
    print(f" - 总计 BAM 数量: {len(all_samples)}")
    print(f" - JACUSA2 输出文件: {JACUSA_OUT_FILE}")
    print(f" - 样本顺序 manifest: {SAMPLE_MANIFEST_FILE}")
    print("=================================================\n")

    if VALIDATE_ONLY:
        print("验证完成：未启动 JACUSA2 或 featureCounts。")
        return

    run_fc_flag = not os.path.exists(EXPRESSION_OUT_FILE)
    if not run_fc_flag:
        print(f"🚀 [状态检测] 发现 readcount 文件已存在: {EXPRESSION_OUT_FILE}")
        print("🚀 [动作拦截] 将跳过 featureCounts，仅运行 JACUSA2！\n")

    process_jacusa = multiprocessing.Process(
        target=run_jacusa_aggregated,
        args=(wt_samples, ko_samples, THREADS, MIN_COVERAGE),
    )

    if run_fc_flag:
        process_fc = multiprocessing.Process(
            target=run_featurecounts,
            args=(all_samples, THREADS),
        )

    process_jacusa.start()
    if run_fc_flag:
        process_fc.start()

    process_jacusa.join()
    if run_fc_flag:
        process_fc.join()

    if process_jacusa.exitcode != 0:
        raise RuntimeError(f"JACUSA2 进程失败，exitcode={process_jacusa.exitcode}")
    if run_fc_flag and process_fc.exitcode != 0:
        raise RuntimeError(f"featureCounts 进程失败，exitcode={process_fc.exitcode}")

    elapsed = time.time() - start_time
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)

    print("\n========== 所有任务结束 ==========")
    print(f"总耗时: {int(hours)}小时 {int(minutes)}分钟 {int(seconds)}秒")
    print(f"JACUSA2 输出: {JACUSA_OUT_FILE}")
    print(f"样本 manifest: {SAMPLE_MANIFEST_FILE}")
    if run_fc_flag:
        print(f"featureCounts 输出: {EXPRESSION_OUT_FILE}")
    else:
        print("featureCounts 输出: [文件已存在，已跳过执行]")


if __name__ == "__main__":
    main()
