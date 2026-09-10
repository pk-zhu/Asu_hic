#!/usr/bin/env python3
"""
02_bin_3d_vs_expr.py — bin 级 3D 保守性 (03/) × bin 内基因表达 (08).

把 Asu 每个 bin 内的基因 TPM 聚合 (mean / median), 与 bin 的 3D 类型 join:
  compartment_conservation : 100kb bin, conservation_type ∈ {conserved, A_to_B, B_to_A}
  tad_conservation         : TAD 区域, conservation_score
  loop_conservation        : 5kb anchor 对, conservation ∈ {highly/moderately/weakly/not_conserved, unmappable}

输出 results/bin_3d_expr/
  - bin_expr_compartment_Asu_X.tsv   (一行一个 bin + 聚合后的 TPM)
  - bin_expr_tad_Asu_X.tsv           (一行一个 TAD)
  - bin_expr_loop_Asu_X.tsv          (一行一个 loop, anchor1+anchor2 的表达均值)
  - summary_bin_expr.tsv             (各分类下 expr 均值, MWU)
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from common import (
    load_tpm, tpm_summary, load_genes_bed, PWCONS_DIR, OUT_BASE, asu_subgenome,
)

OUT_DIR = os.path.join(OUT_BASE, "bin_3d_expr")
os.makedirs(OUT_DIR, exist_ok=True)


def build_gene_index(species):
    """gene_id × chr/start/end + mean_tpm/log_mean (取 Asu_genes.bed + Asu.tpm)"""
    g = load_genes_bed(species)
    t = tpm_summary(load_tpm(species))
    df = g.merge(t, left_on="gene_id", right_index=True, how="left")
    df["subgenome"] = df["chr"].map(asu_subgenome) if species == "Asu" else species
    return df


def agg_expr_in_region(genes, chrom, start, end):
    """返回 region 内基因数 + 平均 TPM + 平均 log_mean"""
    sub = genes[(genes["chr"] == chrom)
                & (genes["start"] < end)
                & (genes["end"] > start)
                & genes["mean_tpm"].notna()]
    if len(sub) == 0:
        return 0, np.nan, np.nan
    return len(sub), sub["mean_tpm"].mean(), sub["log_mean"].mean()


def do_compartment(genes):
    """compartment 100kb bin -> bin 内基因表达均值"""
    print("  [compartment]")
    rows_all = []
    for comp in ("Asu_Ath", "Asu_Aar"):
        path = os.path.join(PWCONS_DIR, "compartment", f"compartment_conservation_{comp}.tsv")
        if not os.path.exists(path): continue
        df = pd.read_csv(path, sep="\t")
        recs = []
        for r in df.itertuples():
            n, mt, lm = agg_expr_in_region(genes, r.asu_chr, r.asu_start, r.asu_end)
            recs.append((r.asu_chr, r.asu_start, r.asu_end,
                         r.asu_compartment, r.ref_compartment, r.conservation_type,
                         n, mt, lm))
        out = pd.DataFrame(recs, columns=["asu_chr","asu_start","asu_end",
                                          "asu_compartment","ref_compartment","conservation_type",
                                          "n_genes","mean_tpm","mean_log"])
        out["comparison"] = comp
        out.to_csv(os.path.join(OUT_DIR, f"bin_expr_compartment_{comp}.tsv"), sep="\t", index=False)
        rows_all.append(out)
        sub = out[out.n_genes > 0]
        for ctype, g in sub.groupby("conservation_type"):
            print(f"    {comp} {ctype}: n_bins={len(g)} mean_log={g['mean_log'].mean():.3f}")
    return pd.concat(rows_all, ignore_index=True) if rows_all else None


def do_tad(genes):
    """每个 TAD 区域内基因表达 (用 10kb 分辨率的 TAD)"""
    print("  [tad @ 10000bp]")
    rows_all = []
    for comp in ("Asu_Ath", "Asu_Aar"):
        path = os.path.join(PWCONS_DIR, "tad", "10000", f"tad_conservation_{comp}.tsv")
        if not os.path.exists(path): continue
        df = pd.read_csv(path, sep="\t")
        recs = []
        for r in df.itertuples():
            n, mt, lm = agg_expr_in_region(genes, r.asu_chr, r.asu_tad_start, r.asu_tad_end)
            recs.append((r.asu_chr, r.asu_tad_start, r.asu_tad_end, r.asu_tad_name,
                         r.conservation_score, n, mt, lm))
        out = pd.DataFrame(recs, columns=["asu_chr","asu_tad_start","asu_tad_end","asu_tad_name",
                                          "conservation_score","n_genes","mean_tpm","mean_log"])
        out["comparison"] = comp
        out.to_csv(os.path.join(OUT_DIR, f"bin_expr_tad_{comp}.tsv"), sep="\t", index=False)
        rows_all.append(out)
        sub = out[out.n_genes > 0]
        if len(sub) > 10:
            rho, p = spearmanr(sub["conservation_score"], sub["mean_log"])
            print(f"    {comp} TAD: n={len(sub)} Spearman(conservation × mean_log)={rho:.3f}, p={p:.3g}")
    return pd.concat(rows_all, ignore_index=True) if rows_all else None


def do_loop(genes):
    """loop anchor1+anchor2 的 mean_log 平均"""
    print("  [loop]")
    rows_all = []
    for comp in ("Asu_Ath", "Asu_Aar"):
        path = os.path.join(PWCONS_DIR, "loop", f"loop_conservation_{comp}.tsv")
        if not os.path.exists(path): continue
        df = pd.read_csv(path, sep="\t")
        recs = []
        for r in df.itertuples():
            n1, mt1, lm1 = agg_expr_in_region(genes, r.asu_chr1, int(r.asu_start1), int(r.asu_end1))
            n2, mt2, lm2 = agg_expr_in_region(genes, r.asu_chr2, int(r.asu_start2), int(r.asu_end2))
            recs.append((r.asu_chr1, r.asu_start1, r.asu_end1,
                         r.asu_chr2, r.asu_start2, r.asu_end2,
                         r.conservation, n1, n2,
                         np.nanmean([lm1, lm2])))
        out = pd.DataFrame(recs, columns=["asu_chr1","asu_start1","asu_end1",
                                          "asu_chr2","asu_start2","asu_end2",
                                          "conservation","n_genes_a1","n_genes_a2","mean_log_anchors"])
        out["comparison"] = comp
        out.to_csv(os.path.join(OUT_DIR, f"bin_expr_loop_{comp}.tsv"), sep="\t", index=False)
        rows_all.append(out)
        sub = out[out["mean_log_anchors"].notna()]
        for cls, g in sub.groupby("conservation"):
            print(f"    {comp} loop {cls}: n={len(g)} mean_log_anchors={g['mean_log_anchors'].mean():.3f}")
    return pd.concat(rows_all, ignore_index=True) if rows_all else None


def write_summary(comp_all, tad_all, loop_all):
    rows = []
    if comp_all is not None:
        for comp in ("Asu_Ath", "Asu_Aar"):
            s = comp_all[(comp_all.comparison == comp) & (comp_all.n_genes > 0)]
            for ct, g in s.groupby("conservation_type"):
                rows.append({"comparison": comp, "layer": "compartment", "category": ct,
                             "n_bins": len(g),
                             "mean_mean_log": round(g.mean_log.mean(), 4),
                             "mean_genes_per_bin": round(g.n_genes.mean(), 2)})
    if tad_all is not None:
        for comp in ("Asu_Ath", "Asu_Aar"):
            s = tad_all[(tad_all.comparison == comp) & (tad_all.n_genes > 0)].copy()
            s["bucket"] = pd.cut(s["conservation_score"], [-0.01, 0.25, 0.5, 0.75, 1.01],
                                 labels=["0-0.25","0.25-0.5","0.5-0.75","0.75-1"])
            for b, g in s.groupby("bucket", observed=True):
                rows.append({"comparison": comp, "layer": "tad", "category": f"score_{b}",
                             "n_bins": len(g),
                             "mean_mean_log": round(g.mean_log.mean(), 4),
                             "mean_genes_per_bin": round(g.n_genes.mean(), 2)})
    if loop_all is not None:
        for comp in ("Asu_Ath", "Asu_Aar"):
            s = loop_all[(loop_all.comparison == comp) & loop_all.mean_log_anchors.notna()]
            for cls, g in s.groupby("conservation"):
                rows.append({"comparison": comp, "layer": "loop", "category": cls,
                             "n_bins": len(g),
                             "mean_mean_log": round(g.mean_log_anchors.mean(), 4),
                             "mean_genes_per_bin": round((g.n_genes_a1 + g.n_genes_a2).mean()/2, 2)})
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "summary_bin_expr.tsv"), sep="\t", index=False)


def main():
    print("=== 02_bin_3d_vs_expr (bin/TAD/loop × Asu TPM) ===")
    genes = build_gene_index("Asu")
    print(f"  Asu 基因: {len(genes)}, 有 TPM: {genes['mean_tpm'].notna().sum()}")
    comp_all = do_compartment(genes)
    tad_all  = do_tad(genes)
    loop_all = do_loop(genes)
    write_summary(comp_all, tad_all, loop_all)
    print(f"  保存: {OUT_DIR}/")


if __name__ == "__main__":
    main()
