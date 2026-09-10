#!/usr/bin/env python3
"""
09_tad_boundary_indel.py — 重建 TAD 边界 vs 内部的序列分歧与 InDel 密度.

输入:
  - MUMmer 1-to-1 delta:     02.collinearity/03_mummer/{sT_vs_Ath,sA_vs_Aar}.1to1.delta
  - MUMmer coords:           同目录 .coords (show-coords -rclTH)
  - 共线 10 kb bin:          02.collinearity/04_collinear_pairs/Asu_{Ath,Aar}_orth_bins_10000.tsv
  - Asu TAD BED:             01.single/Asu/TADs/HiCExplorer/10000/Asu_domains.bed
  - TAD conservation:        results/tad/10000/tad_conservation_Asu_{Ath,Aar}.tsv

输出 (results/tad_features/):
  tad_boundary_divergence_perTAD.tsv
    pair, tad_class, boundary_div, interior_div,
    boundary_indel_per_kb, interior_indel_per_kb
  tad_boundary_divergence_summary.tsv
    pair, tad_class, n, boundary_div_median, interior_div_median, div_wilcoxon_p,
    boundary_indel_median, interior_indel_median, indel_wilcoxon_p

定义:
  - 边界 = TAD 起止位置两侧各 250 bp 的 Asu 区间 (共 4 段, 合计 ≤1000 bp);
    内部 = TAD 区间去除两端各 250 bp 后的 Asu 区间.
  - Asu 区间经共线 10 kb bin 线性插值投影到亲本 (Ath/Aar) 坐标, 合并重叠.
  - InDel = MUMmer alignment 中 ref base 为 '.' (qry 插入) 或 qry base 为 '.'
    (ref 插入, 即 qry 缺失) 的 SNP 行, 每行计 1 bp, 以 ref 坐标 P1 锚定;
    投影区间内全部 indel 碱基数之和 / 投影有效长度 (kb) → indel/kb.
    indel 同时包含 ref 插入与 qry 插入 (两方向).
  - div = 100 - %idy, 以 alignment 在投影区间内的 ref 覆盖碱基数为权重的长度加权平均.
  - tad_class 按 conservation_score 分: ≥0.8 conserved, 0.4–0.8 moderate, <0.4 diverged.
"""
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
from bisect import bisect_left, bisect_right
from scipy.stats import wilcoxon

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
# Locate MUMmer v4.x binaries: prefer PATH, then $MUMMER_HOME, then ~/software/mummer-4.0.1.
def _mummer_bin():
    exe = shutil.which("show-snps")
    if exe:
        return os.path.dirname(exe)
    env = os.environ.get("MUMMER_HOME")
    if env and os.path.exists(os.path.join(env, "show-snps")):
        return env
    home = os.path.join(os.environ.get("HOME", ""), "software", "mummer-4.0.1")
    return home if os.path.exists(os.path.join(home, "show-snps")) else ""


MUMMER_BIN = _mummer_bin()
DELTA = {
    "Asu_Ath": f"{BASE}/02.collinearity/03_mummer/sT_vs_Ath.1to1.delta",
    "Asu_Aar": f"{BASE}/02.collinearity/03_mummer/sA_vs_Aar.1to1.delta",
}
COORDS = {
    "Asu_Ath": f"{BASE}/02.collinearity/03_mummer/sT_vs_Ath.coords",
    "Asu_Aar": f"{BASE}/02.collinearity/03_mummer/sA_vs_Aar.coords",
}
ORTH_BINS = {
    "Asu_Ath": f"{BASE}/02.collinearity/04_collinear_pairs/Asu_Ath_orth_bins_10000.tsv",
    "Asu_Aar": f"{BASE}/02.collinearity/04_collinear_pairs/Asu_Aar_orth_bins_10000.tsv",
}
ASU_TAD = f"{BASE}/01.single/Asu/TADs/HiCExplorer/10000/Asu_domains.bed"
TAD_CONS = {
    "Asu_Ath": f"{BASE}/03.pairwise_subgenome_conservation/results/tad/10000/tad_conservation_Asu_Ath.tsv",
    "Asu_Aar": f"{BASE}/03.pairwise_subgenome_conservation/results/tad/10000/tad_conservation_Asu_Aar.tsv",
}
OUT_DIR = f"{BASE}/03.pairwise_subgenome_conservation/results/tad_features"
os.makedirs(OUT_DIR, exist_ok=True)

BOUNDARY_EXT = 250


def log(msg):
    print(msg, flush=True)


def build_indel_index(delta_path):
    """调用 show-snps -rTH, 返回 dict[ref_chr] -> sorted ndarray of indel ref positions (每 base 1 计数)."""
    out = subprocess.run([f"{MUMMER_BIN}/show-snps", "-r", "-T", "-H", delta_path],
                         capture_output=True, text=True, check=True)
    idx = {}
    for line in out.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 12:
            continue
        ref_base, qry_base = p[1], p[2]
        if ref_base != "." and qry_base != ".":
            continue  # 纯 SNP, 非 indel
        chrom = p[10]
        try:
            pos = int(p[0])
        except ValueError:
            continue
        idx.setdefault(chrom, []).append(pos)
    for ch in idx:
        idx[ch] = np.sort(np.array(idx[ch], dtype=np.int64))
    return idx


def parse_coords(path):
    """show-coords -rclTH 13 列."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("=") or not line.strip() or line.startswith("NUCMER") or line.startswith("/"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            try:
                rows.append(dict(
                    r_s=int(p[0]), r_e=int(p[1]), q_s=int(p[2]), q_e=int(p[3]),
                    idy=float(p[6]), ref_chr=p[11], qry_chr=p[12]))
            except (ValueError, IndexError):
                continue
    return rows


def build_orth_bins(path):
    df = pd.read_csv(path, sep="\t", header=None, usecols=range(6),
                     names=["asu_chr", "asu_s", "asu_e", "ref_chr", "ref_s", "ref_e"])
    out = {}
    for ch, g in df.groupby("asu_chr"):
        g = g.sort_values("asu_s")
        out[ch] = list(zip(g["asu_s"].values, g["asu_e"].values,
                           g["ref_chr"].values, g["ref_s"].values, g["ref_e"].values))
    return out


def project_interval(bins, asu_s, asu_e):
    out = []
    for bs, be, rch, rs, re in bins:
        if be <= asu_s or bs >= asu_e:
            continue
        ov_s = max(asu_s, bs)
        ov_e = min(asu_e, be)
        if ov_e <= ov_s:
            continue
        span = be - bs
        fs = (ov_s - bs) / span
        fe = (ov_e - bs) / span
        if re >= rs:
            r_s2 = int(round(rs + fs * (re - rs)))
            r_e2 = int(round(rs + fe * (re - rs)))
        else:
            r_s2 = int(round(rs - fs * (rs - re)))
            r_e2 = int(round(rs - fe * (rs - re)))
        if r_e2 < r_s2:
            r_s2, r_e2 = r_e2, r_s2
        out.append((rch, r_s2, r_e2))
    return out


def merge_intervals(ivs):
    by_chr = {}
    for ch, s, e in ivs:
        by_chr.setdefault(ch, []).append((min(s, e), max(s, e)))
    merged = {}
    for ch, lst in by_chr.items():
        lst.sort()
        out = []
        for s, e in lst:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[ch] = out
    return merged


def total_len(merged):
    return sum(e - s for lst in merged.values() for s, e in lst)


def count_indel(merged, indel_idx):
    total = 0
    for ch, ivs in merged.items():
        pos = indel_idx.get(ch)
        if pos is None or len(pos) == 0:
            continue
        for s, e in ivs:
            lo = bisect_left(pos, s)
            hi = bisect_right(pos, e)
            total += int(hi - lo)
    return total


def weighted_div(merged, coords):
    by_chr = {}
    for c in coords:
        by_chr.setdefault(c["ref_chr"], []).append(c)
    for ch in by_chr:
        by_chr[ch].sort(key=lambda x: x["r_s"])
    tot_len = 0.0
    tot_div = 0.0
    for ch, ivs in merged.items():
        for c in by_chr.get(ch, []):
            c_s, c_e = c["r_s"], c["r_e"]
            for s, e in ivs:
                ov_s = max(c_s, s)
                ov_e = min(c_e, e)
                if ov_e > ov_s:
                    L = ov_e - ov_s
                    tot_len += L
                    tot_div += L * (100.0 - c["idy"])
    return tot_div / tot_len if tot_len else np.nan


def tad_class(score):
    if pd.isna(score):
        return np.nan
    if score >= 0.8:
        return "conserved"
    if score >= 0.4:
        return "moderate"
    return "diverged"


def main():
    asu_tads = pd.read_csv(ASU_TAD, sep="\t", header=None,
                           names=["chr", "start", "end", "name", "score",
                                  "strand", "thickStart", "thickEnd", "rgb"])
    rows = []
    for pair in ("Asu_Ath", "Asu_Aar"):
        log(f"--- {pair}: show-snps ---")
        indel_idx = build_indel_index(DELTA[pair])
        n_indel = sum(len(v) for v in indel_idx.values())
        log(f"    indel bases in ref coords: {n_indel:,}")
        coords = parse_coords(COORDS[pair])
        bins = build_orth_bins(ORTH_BINS[pair])
        cons = pd.read_csv(TAD_CONS[pair], sep="\t")
        score_map = dict(zip(cons["asu_tad_name"], cons["conservation_score"]))

        for ch, g in asu_tads.groupby("chr"):
            if ch not in bins:
                continue
            for _, t in g.iterrows():
                ts, te = int(t["start"]), int(t["end"])
                if te - ts < 2 * BOUNDARY_EXT + 1000:
                    continue
                bnd_ivs = []
                for a, b in [(ts - BOUNDARY_EXT, ts + BOUNDARY_EXT),
                             (te - BOUNDARY_EXT, te + BOUNDARY_EXT)]:
                    bnd_ivs.extend(project_interval(bins[ch], max(0, a), b))
                bnd_merged = merge_intervals(bnd_ivs)
                int_ivs = project_interval(bins[ch], ts + BOUNDARY_EXT, te - BOUNDARY_EXT)
                int_merged = merge_intervals(int_ivs)

                bnd_len = total_len(bnd_merged)
                int_len = total_len(int_merged)
                if bnd_len < 100 or int_len < 1000:
                    continue
                rows.append({
                    "pair": pair,
                    "tad_class": tad_class(score_map.get(t["name"], np.nan)),
                    "boundary_div": weighted_div(bnd_merged, coords),
                    "interior_div": weighted_div(int_merged, coords),
                    "boundary_indel_per_kb": count_indel(bnd_merged, indel_idx) / (bnd_len / 1000.0),
                    "interior_indel_per_kb": count_indel(int_merged, indel_idx) / (int_len / 1000.0),
                })
    df = pd.DataFrame(rows).dropna(subset=["tad_class"])
    df.to_csv(f"{OUT_DIR}/tad_boundary_divergence_perTAD.tsv", sep="\t", index=False)
    log(f"perTAD rows: {len(df):,}")

    summ = []
    for (pair, cls), sub in df.groupby(["pair", "tad_class"]):
        d_div = sub[["boundary_div", "interior_div"]].dropna()
        d_ind = sub[["boundary_indel_per_kb", "interior_indel_per_kb"]].dropna()
        p_div = wilcoxon(d_div["boundary_div"], d_div["interior_div"]).pvalue if len(d_div) >= 5 else np.nan
        p_ind = wilcoxon(d_ind["boundary_indel_per_kb"], d_ind["interior_indel_per_kb"]).pvalue if len(d_ind) >= 5 else np.nan
        summ.append({
            "pair": pair, "tad_class": cls, "n": len(sub),
            "boundary_div_median": round(sub["boundary_div"].median(), 4),
            "interior_div_median": round(sub["interior_div"].median(), 4),
            "div_wilcoxon_p": p_div,
            "boundary_indel_median": round(sub["boundary_indel_per_kb"].median(), 3),
            "interior_indel_median": round(sub["interior_indel_per_kb"].median(), 3),
            "indel_wilcoxon_p": p_ind,
        })
    sdf = pd.DataFrame(summ)
    sdf.to_csv(f"{OUT_DIR}/tad_boundary_divergence_summary.tsv", sep="\t", index=False)
    log(sdf.to_string(index=False))
    log("=== DONE ===")


if __name__ == "__main__":
    main()
