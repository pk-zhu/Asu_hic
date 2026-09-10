#!/usr/bin/env python3
"""
04. CDS 密度 ↔ 3D 保守性 (Part B2): 所有 4 层 + 自然 + 合成 + MWU.

把 Asu 侧区域的 CDS 覆盖率 (密度) 与 03 的四层 per-region 保守性关联:
  compartment (100k): 每 bin CDS 密度  -> 比较 conserved vs A_to_B/B_to_A (MWU)
  TAD (10k):          每 TAD CDS 密度  -> 与 conservation_score 相关 (Spearman)
  loop (5k):          loop anchor CDS 密度 -> 比较各 conservation 类 (MWU)
  matrix:             per-gene matrix_corr vs 基因区 CDS 密度 -> 相关 (Spearman)

比较 (4 路):
  Asu_Ath              自然 sT  → Ath
  Asu_Aar              自然 sA  → Aar
  syn_Asu_vs_Ath       合成 sT  → Ath
  syn_Asu_vs_Aar       合成 sA  → Aar

输出 -> results/density_link/
  cds_density_vs_{compartment,tad,loop,matrix}.tsv         (per-bin/region)
  density_link_summary.tsv                                   (汇总, 含 MWU p 值)
"""
import os
import csv
from collections import defaultdict
from scipy.stats import spearmanr, mannwhitneyu

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
CONS_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results")
CDS_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/cds")
GENE3D_NAT_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/gene_3d/10000")
GENE3D_SYN_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/gene_3d_syn_Asu/10000")
OUTPUT_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/density_link")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COMPARISONS = ["Asu_Ath", "Asu_Aar", "syn_Asu_vs_Ath", "syn_Asu_vs_Aar"]
SUB = {
    "Asu_Ath": "sT", "Asu_Aar": "sA",
    "syn_Asu_vs_Ath": "sT", "syn_Asu_vs_Aar": "sA",
}
SOURCE = {"Asu_Ath": "nat", "Asu_Aar": "nat",
          "syn_Asu_vs_Ath": "syn", "syn_Asu_vs_Aar": "syn"}
WINDOW = 50000  # 与 03 一致, 基因区扩展


def _load_bed3(fname):
    feats = defaultdict(list)
    with open(os.path.join(CDS_DIR, fname)) as f:
        for line in f:
            p = line.split("\t")[:3]
            feats[p[0]].append((int(p[1]), int(p[2])))
    for c in feats:
        feats[c].sort()
    return feats


def load_cds():
    # CDS/外显子片段级 bed (284k 行, 一个基因多个 CDS 片段)
    return _load_bed3("Asu_cds.bed")


def load_genes():
    # 基因级 bed (58k 行, 一个基因一行), 用于基因数密度
    return _load_bed3("Asu_genes.bed")


def coverage(cds_chr, start, end):
    """区间内 CDS 碱基覆盖比例 (fraction of bp), 旧口径, 保留供对照与下游列."""
    if end <= start or not cds_chr:
        return 0.0
    import bisect
    starts = [s for s, _ in cds_chr]
    j = bisect.bisect_left(starts, start)
    if j > 0:
        j -= 1
    total = 0
    for k in range(j, len(cds_chr)):
        s, e = cds_chr[k]
        if s >= end:
            break
        if e <= start:
            continue
        total += min(e, end) - max(s, start)
    return total / (end - start)


def count_features(feats_chr, start, end):
    """与区间有任意重叠的特征个数 (bed 已按 start 排序)."""
    if end <= start or not feats_chr:
        return 0
    import bisect
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


def densities(genes_chr, cds_chr, start, end):
    """返回 (gene_density_per_mb, cds_density_per_mb): 单位长度计数密度.

    主口径 gene density = 区间内重叠基因数 / 区间长度(Mb);
    对照 cds density = 区间内重叠 CDS(外显子)片段数 / 区间长度(Mb).
    """
    mb = max((end - start) / 1e6, 1e-9)
    return (count_features(genes_chr, start, end) / mb,
            count_features(cds_chr, start, end) / mb)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def group_stats(values_by_group):
    out = {}
    for g, vals in values_by_group.items():
        if not vals:
            out[g] = (0, float("nan"), float("nan"))
            continue
        import numpy as np
        arr = np.array(vals, dtype=float)
        out[g] = (len(arr), float(np.mean(arr)), float(np.median(arr)))
    return out


def write_rows(rows, fname):
    path = os.path.join(OUTPUT_DIR, fname)
    if not rows:
        open(path, "w").close()
        print(f"  -> {path} (empty)")
        return
    with open(path, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  -> {path} ({len(rows)} rows)")


def write_summary(summary):
    sout = os.path.join(OUTPUT_DIR, "density_link_summary.tsv")
    with open(sout, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["layer", "comparison", "metric", "value", "n"], delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    print(f"\n汇总 -> {sout}")
    for s in summary:
        print(f"  {s['layer']:<12}{s['comparison']:<18}{s['metric']:<35}{s['value']}  (n={s['n']})")


def main():
    cds = load_cds()
    genes = load_genes()
    summary = []

    # ---------- compartment ----------
    comp_rows = []
    # 主统计用 gene density (genes/Mb); cds fragment density 作对照; coverage 保留
    comp_groups = defaultdict(lambda: defaultdict(lambda: {"gene": [], "cds": [], "cov": []}))
    for comp in COMPARISONS:
        path = os.path.join(CONS_DIR, "compartment", f"compartment_conservation_{comp}.tsv")
        if not os.path.exists(path):
            print(f"  [skip compartment/{comp}] missing {path}")
            continue
        for row in read_tsv(path):
            ac = row["asu_chr"]
            s, e = int(row["asu_start"]), int(row["asu_end"])
            cov = coverage(cds.get(ac, []), s, e)
            gdens, cdens = densities(genes.get(ac, []), cds.get(ac, []), s, e)
            ctype = row["conservation_type"]
            comp_rows.append({"comparison": comp, "subgenome": SUB[comp],
                              "source": SOURCE[comp], "asu_chr": ac,
                              "asu_start": s, "asu_end": e, "conservation_type": ctype,
                              "gene_density_per_mb": round(gdens, 3),
                              "cds_density_per_mb": round(cdens, 3),
                              "cds_coverage": round(cov, 6)})
            grp = ctype  # conserved | A_to_B | B_to_A, 三类细分
            comp_groups[comp][grp]["gene"].append(gdens)
            comp_groups[comp][grp]["cds"].append(cdens)
            comp_groups[comp][grp]["cov"].append(cov)
    write_rows(comp_rows, "cds_density_vs_compartment.tsv")

    def _mean(vals):
        import numpy as np
        return float(np.mean(vals)) if vals else float("nan")

    for comp in COMPARISONS:
        for grp in ("conserved", "A_to_B", "B_to_A"):
            d = comp_groups[comp].get(grp)
            if not d:
                continue
            n = len(d["gene"])
            summary.append({"layer": "compartment", "comparison": comp,
                            "metric": f"mean_gene_density[{grp}]",
                            "value": round(_mean(d["gene"]), 3), "n": n})
            summary.append({"layer": "compartment", "comparison": comp,
                            "metric": f"mean_cds_density[{grp}]",
                            "value": round(_mean(d["cds"]), 3), "n": n})
            summary.append({"layer": "compartment", "comparison": comp,
                            "metric": f"mean_cds_cov[{grp}]",
                            "value": round(_mean(d["cov"]), 6), "n": n})
        # MWU 三对, 主口径 gene density, 另出 cds density 与 coverage 对照
        cons = comp_groups[comp].get("conserved", {"gene": [], "cds": [], "cov": []})
        atob = comp_groups[comp].get("A_to_B", {"gene": [], "cds": [], "cov": []})
        btoa = comp_groups[comp].get("B_to_A", {"gene": [], "cds": [], "cov": []})
        for pa, pb, label in [(cons, atob, "conserved_vs_A_to_B"),
                              (cons, btoa, "conserved_vs_B_to_A"),
                              (atob, btoa, "A_to_B_vs_B_to_A")]:
            for metric_key, tag in [("gene", "geneDens"), ("cds", "cdsDens"), ("cov", "cov")]:
                va, vb = pa[metric_key], pb[metric_key]
                if len(va) > 5 and len(vb) > 5:
                    u, p = mannwhitneyu(va, vb, alternative="two-sided")
                    summary.append({"layer": "compartment", "comparison": comp,
                                    "metric": f"mwu_{label}_{tag}_stat",
                                    "value": round(float(u), 2),
                                    "n": min(len(va), len(vb))})
                    summary.append({"layer": "compartment", "comparison": comp,
                                    "metric": f"mwu_{label}_{tag}_p",
                                    "value": f"{p:.4e}",
                                    "n": min(len(va), len(vb))})

    # ---------- TAD ----------
    tad_rows = []
    for comp in COMPARISONS:
        path = os.path.join(CONS_DIR, "tad", "10000", f"tad_conservation_{comp}.tsv")
        if not os.path.exists(path):
            print(f"  [skip tad/{comp}] missing {path}")
            continue
        gd_list, cd_list, cov_list, scores = [], [], [], []
        for row in read_tsv(path):
            ac = row["asu_chr"]
            s = int(row.get("asu_tad_start") or row.get("asu_start"))
            e = int(row.get("asu_tad_end") or row.get("asu_end"))
            cov = coverage(cds.get(ac, []), s, e)
            gdens, cdens = densities(genes.get(ac, []), cds.get(ac, []), s, e)
            score = float(row["conservation_score"])
            tad_rows.append({"comparison": comp, "subgenome": SUB[comp],
                             "source": SOURCE[comp], "asu_chr": ac,
                             "asu_tad_start": s, "asu_tad_end": e,
                             "conservation_score": round(score, 6),
                             "gene_density_per_mb": round(gdens, 3),
                             "cds_density_per_mb": round(cdens, 3),
                             "cds_coverage": round(cov, 6)})
            gd_list.append(gdens); cd_list.append(cdens)
            cov_list.append(cov); scores.append(score)
        if len(scores) > 2:
            for vals, tag in [(gd_list, "geneDens"), (cd_list, "cdsDens"), (cov_list, "cov")]:
                rho, p = spearmanr(vals, scores)
                summary.append({"layer": "tad", "comparison": comp,
                                "metric": f"spearman_{tag}_vs_score",
                                "value": round(float(rho), 6), "n": len(scores)})
                summary.append({"layer": "tad", "comparison": comp,
                                "metric": f"spearman_{tag}_p",
                                "value": f"{p:.4e}", "n": len(scores)})
    write_rows(tad_rows, "cds_density_vs_tad.tsv")

    # ---------- loop (仅自然 Asu, syn 当前 03 pipeline 未输出 loop 文件) ----------
    loop_rows = []
    loop_groups = defaultdict(lambda: defaultdict(lambda: {"gene": [], "cds": [], "cov": []}))
    for comp in ["Asu_Ath", "Asu_Aar"]:
        path = os.path.join(CONS_DIR, "loop", f"loop_conservation_{comp}.tsv")
        if not os.path.exists(path):
            continue
        for row in read_tsv(path):
            c1, s1, e1 = row["asu_chr1"], int(row["asu_start1"]), int(row["asu_end1"])
            c2, s2, e2 = row["asu_chr2"], int(row["asu_start2"]), int(row["asu_end2"])
            cov = (coverage(cds.get(c1, []), s1, e1) + coverage(cds.get(c2, []), s2, e2)) / 2
            g1, c1d = densities(genes.get(c1, []), cds.get(c1, []), s1, e1)
            g2, c2d = densities(genes.get(c2, []), cds.get(c2, []), s2, e2)
            gdens, cdens = (g1 + g2) / 2, (c1d + c2d) / 2
            cls = row["conservation"]
            loop_rows.append({"comparison": comp, "subgenome": SUB[comp],
                              "source": "nat", "asu_chr1": c1,
                              "asu_start1": s1, "asu_chr2": c2, "asu_start2": s2,
                              "conservation": cls,
                              "anchor_gene_density_per_mb": round(gdens, 3),
                              "anchor_cds_density_per_mb": round(cdens, 3),
                              "anchor_cds_coverage": round(cov, 6)})
            grp = "not_conserved" if cls in ("not_conserved", "unmappable") else "conserved"
            loop_groups[comp][grp]["gene"].append(gdens)
            loop_groups[comp][grp]["cds"].append(cdens)
            loop_groups[comp][grp]["cov"].append(cov)
    write_rows(loop_rows, "cds_density_vs_loop.tsv")
    for comp in ["Asu_Ath", "Asu_Aar"]:
        for grp in ("conserved", "not_conserved"):
            d = loop_groups[comp].get(grp)
            if not d:
                continue
            n = len(d["gene"])
            summary.append({"layer": "loop", "comparison": comp,
                            "metric": f"mean_anchor_gene_density[{grp}]",
                            "value": round(_mean(d["gene"]), 3), "n": n})
            summary.append({"layer": "loop", "comparison": comp,
                            "metric": f"mean_anchor_cds_density[{grp}]",
                            "value": round(_mean(d["cds"]), 3), "n": n})
            summary.append({"layer": "loop", "comparison": comp,
                            "metric": f"mean_anchor_cds_cov[{grp}]",
                            "value": round(_mean(d["cov"]), 6), "n": n})
        cons = loop_groups[comp].get("conserved", {"gene": [], "cds": [], "cov": []})
        nc = loop_groups[comp].get("not_conserved", {"gene": [], "cds": [], "cov": []})
        for metric_key, tag in [("gene", "geneDens"), ("cds", "cdsDens"), ("cov", "cov")]:
            va, vb = cons[metric_key], nc[metric_key]
            if len(va) > 5 and len(vb) > 5:
                u, p = mannwhitneyu(va, vb, alternative="two-sided")
                summary.append({"layer": "loop", "comparison": comp,
                                "metric": f"mwu_conserved_vs_not_conserved_{tag}_stat",
                                "value": round(float(u), 2),
                                "n": min(len(va), len(vb))})
                summary.append({"layer": "loop", "comparison": comp,
                                "metric": f"mwu_conserved_vs_not_conserved_{tag}_p",
                                "value": f"{p:.4e}",
                                "n": min(len(va), len(vb))})

    # ---------- matrix (per-gene) ----------
    mat_rows = []
    for comp in COMPARISONS:
        gene_dir = GENE3D_SYN_DIR if SOURCE[comp] == "syn" else GENE3D_NAT_DIR
        path = os.path.join(gene_dir, f"gene_3d_conservation_{comp}.tsv")
        if not os.path.exists(path):
            print(f"  [skip matrix/{comp}] missing {path}")
            continue
        gd_list, cd_list, cov_list, corrs = [], [], [], []
        for row in read_tsv(path):
            if row["matrix_corr"] == "NA":
                continue
            ac = row["asu_chr"]
            s = max(0, int(row["asu_start"]) - WINDOW)
            e = int(row["asu_end"]) + WINDOW
            cov = coverage(cds.get(ac, []), s, e)
            gdens, cdens = densities(genes.get(ac, []), cds.get(ac, []), s, e)
            mc = float(row["matrix_corr"])
            mat_rows.append({"comparison": comp, "subgenome": SUB[comp],
                             "source": SOURCE[comp], "gene_id": row["gene_id"],
                             "asu_chr": ac, "matrix_corr": round(mc, 6),
                             "gene_density_per_mb": round(gdens, 3),
                             "cds_density_per_mb": round(cdens, 3),
                             "region_cds_coverage": round(cov, 6)})
            gd_list.append(gdens); cd_list.append(cdens)
            cov_list.append(cov); corrs.append(mc)
        if len(corrs) > 2:
            for vals, tag in [(gd_list, "geneDens"), (cd_list, "cdsDens"), (cov_list, "cov")]:
                rho, p = spearmanr(vals, corrs)
                summary.append({"layer": "matrix", "comparison": comp,
                                "metric": f"spearman_{tag}_vs_corr",
                                "value": round(float(rho), 6), "n": len(corrs)})
                summary.append({"layer": "matrix", "comparison": comp,
                                "metric": f"spearman_{tag}_p",
                                "value": f"{p:.4e}", "n": len(corrs)})
    write_rows(mat_rows, "cds_density_vs_matrix.tsv")

    write_summary(summary)
    print("04 done.")


if __name__ == "__main__":
    main()
