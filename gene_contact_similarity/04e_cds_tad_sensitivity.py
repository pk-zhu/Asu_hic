#!/usr/bin/env python3
"""
04e. CDS 密度 ↔ TAD conservation score 的分辨率 sensitivity:

  TAD 分辨率 5 / 10 / 25 kb × 4 个比较 (Asu_Ath, Asu_Aar, syn_Asu_vs_Ath, syn_Asu_vs_Aar)
  = 12 行 Spearman ρ / p / n.

TAD source 已经含亚基因组过滤 (syn_Asu_vs_Asu 同名但已分别过 sT / sA),
不需要再行过滤. 对每个 TAD: 取 bin 段 CDS coverage, 与 conservation_score
算 Spearman ρ.

输出 -> results/density_link/cds_density_tad_sensitivity.tsv
"""
import os
import csv
import bisect
from collections import defaultdict
from scipy.stats import spearmanr

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
CDS_BED = os.path.join(BASE, "04.cds_synteny_3d_conservation/results/cds/Asu_cds.bed")
GENE_BED = os.path.join(BASE, "04.cds_synteny_3d_conservation/results/cds/Asu_genes.bed")
TAD_BASE = os.path.join(BASE, "03.pairwise_subgenome_conservation/results/tad")
OUT = os.path.join(BASE, "04.cds_synteny_3d_conservation/results/density_link/cds_density_tad_sensitivity.tsv")

COMPARISONS = ["Asu_Ath", "Asu_Aar", "syn_Asu_vs_Ath", "syn_Asu_vs_Aar"]
SOURCE = {c: "syn" if c.startswith("syn_") else "nat" for c in COMPARISONS}
TAD_RES = [5000, 10000, 25000]


def _load_bed(path):
    feats = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            feats[p[0]].append((int(p[1]), int(p[2])))
    for c in feats:
        feats[c].sort()
    return feats


def load_cds():
    return _load_bed(CDS_BED)


def load_gene_bed():
    return _load_bed(GENE_BED)


def coverage(feats_chr, start, end):
    """区间内特征碱基覆盖比例 (旧 CDS coverage 口径, 保留作对照)."""
    if end <= start or not feats_chr:
        return 0.0
    starts = [s for s, _ in feats_chr]
    j = bisect.bisect_left(starts, start)
    if j > 0:
        j -= 1
    total = 0
    for k in range(j, len(feats_chr)):
        s, e = feats_chr[k]
        if s >= end:
            break
        if e <= start:
            continue
        total += min(e, end) - max(s, start)
    return total / (end - start)


def count_features(feats_chr, start, end):
    """与区间有任意重叠的特征个数."""
    if end <= start or not feats_chr:
        return 0
    starts = [s for s, _ in feats_chr]
    j = bisect.bisect_left(starts, start)
    if j > 0:
        j -= 1
    n = 0
    for k in range(j, len(feats_chr)):
        s, e = feats_chr[k]
        if s >= end:
            break
        if e > start:
            n += 1
    return n


def tad_cds_corr(genes_bed, cds_bed, tad_rows):
    """tad_rows: list of (chr, start, end, score).
    返回三口径 (gene dens ρ/p, cds dens ρ/p, coverage ρ/p, n)."""
    gd, cd, covs, scores = [], [], [], []
    for ac, s, e, score in tad_rows:
        mb = max((e - s) / 1e6, 1e-9)
        gd.append(count_features(genes_bed.get(ac, []), s, e) / mb)
        cd.append(count_features(cds_bed.get(ac, []), s, e) / mb)
        covs.append(coverage(cds_bed.get(ac, []), s, e))
        scores.append(score)
    n = len(scores)
    if n < 5:
        return (float("nan"),) * 6 + (n,)
    rg, pg = spearmanr(gd, scores)
    rc, pc = spearmanr(cd, scores)
    rv, pv = spearmanr(covs, scores)
    return (float(rg), float(pg), float(rc), float(pc), float(rv), float(pv), n)


def load_tads(path, four_col):
    """解析 tad_conservation 文件. nat 为 11 列, syn 为 4 列."""
    rows = []
    with open(path) as f:
        header = next(csv.reader(f, delimiter="\t"))
        for row in csv.reader(f, delimiter="\t"):
            if four_col:
                # syn 4-col: asu_chr, asu_start, asu_end, conservation_score
                ac, s, e, score = row[0], int(row[1]), int(row[2]), float(row[3])
            else:
                # nat 11-col: asu_chr, asu_tad_start, asu_tad_end, ...
                ac, s, e = row[0], int(row[1]), int(row[2])
                score = float(row[8])
            rows.append((ac, s, e, score))
    return rows


def _r(x):
    return None if x != x else round(x, 6)


def main():
    cds = load_cds()
    genes_bed = load_gene_bed()
    rows = []
    for res in TAD_RES:
        for comp in COMPARISONS:
            path = os.path.join(TAD_BASE, str(res), f"tad_conservation_{comp}.tsv")
            if not os.path.exists(path):
                print(f"[skip] missing {path}")
                continue
            four_col = SOURCE[comp] == "syn"
            tad_rows = load_tads(path, four_col)
            rg, pg, rc, pc, rv, pv, n = tad_cds_corr(genes_bed, cds, tad_rows)
            rows.append({"tad_resolution_bp": res, "comparison": comp,
                         "source": SOURCE[comp],
                         "spearman_rho": _r(rg), "spearman_p": _r(pg),
                         "spearman_rho_cds_density": _r(rc), "spearman_p_cds_density": _r(pc),
                         "spearman_rho_cds_coverage": _r(rv), "spearman_p_cds_coverage": _r(pv),
                         "n_tads": n})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cols = ["tad_resolution_bp", "comparison", "source",
            "spearman_rho", "spearman_p",
            "spearman_rho_cds_density", "spearman_p_cds_density",
            "spearman_rho_cds_coverage", "spearman_p_cds_coverage", "n_tads"]
    with open(OUT, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"-> {OUT}")
    for r in rows:
        rho = r["spearman_rho"]
        rho_s = f"{rho:+.4f}" if rho is not None else "NA    "
        tag = f"TAD {r['tad_resolution_bp']/1000:.0f}kb"
        print(f"  {tag:<10}{r['comparison']:<18}{r['source']:<4}geneDens ρ={rho_s}  n={r['n_tads']}")


if __name__ == "__main__":
    main()
