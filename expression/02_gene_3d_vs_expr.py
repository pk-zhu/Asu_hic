#!/usr/bin/env python3
"""
01_gene_3d_vs_expr.py — Asu 基因 3D 保守性 (04/gene_3d) × Asu 表达 (08).

对 Asu 的每个 gene_id:
  - 关联到 04/gene_3d 表里的 matrix_corr / compartment_type / tad_conserved_10000
  - 关联到 08/Asu.tpm.tsv 的 mean_tpm / log_mean / expressed
按 (1) sT vs sA 子基因组分组 (2) Asu_Ath vs Asu_Aar 比较
输出: results/gene_3d_expr/
  - per_gene_expr_3d_Asu_Ath.tsv  / per_gene_expr_3d_Asu_Aar.tsv
  - summary_by_compartment.tsv
  - summary_by_tad.tsv
  - summary_matrix_corr_bins.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from common import (
    load_tpm, tpm_summary, GENE3D_DIR, OUT_BASE, asu_subgenome,
)

RES = 10000  # gene_3d 默认分辨率
OUT_DIR = os.path.join(OUT_BASE, "gene_3d_expr")
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    print(f"=== 01_gene_3d_vs_expr (gene_3d @ {RES}bp × Asu TPM) ===")

    tpm = tpm_summary(load_tpm("Asu"))
    print(f"  Asu TPM 矩阵基因数: {len(tpm)}")

    rows_all = []
    for comp in ("Asu_Ath", "Asu_Aar"):
        g3d_path = os.path.join(GENE3D_DIR, str(RES), f"gene_3d_conservation_{comp}.tsv")
        g3d = pd.read_csv(g3d_path, sep="\t")
        g3d["subgenome"] = g3d["asu_chr"].map(asu_subgenome)
        # 仅保留与 comparison 对应的亚基因组 (Asu_Ath 只看 sT, Asu_Aar 只看 sA)
        want = "sT" if comp == "Asu_Ath" else "sA"
        g3d = g3d[g3d["subgenome"] == want].copy()
        # join 表达
        m = g3d.merge(tpm, left_on="gene_id", right_index=True, how="left")
        m["comparison"] = comp
        out_path = os.path.join(OUT_DIR, f"per_gene_expr_3d_{comp}.tsv")
        m.to_csv(out_path, sep="\t", index=False)
        print(f"  {comp}: n={len(m)}, expr-joined={m['mean_tpm'].notna().sum()}  →  {out_path}")
        rows_all.append(m)

    merged = pd.concat(rows_all, ignore_index=True)

    # ---- summary 1: by compartment_type ----
    rows = []
    for comp in ("Asu_Ath", "Asu_Aar"):
        sub = merged[(merged.comparison == comp) & merged.mean_tpm.notna()]
        for ctype, g in sub.groupby("compartment_type"):
            rows.append({
                "comparison": comp,
                "compartment_type": ctype,
                "n_genes": len(g),
                "mean_tpm":   round(g["mean_tpm"].mean(), 4),
                "median_tpm": round(g["mean_tpm"].median(), 4),
                "mean_log":   round(g["log_mean"].mean(), 4),
                "expressed_rate": round(g["expressed"].mean(), 4),
            })
        # MWU: conserved vs A_to_B
        cons = sub[sub.compartment_type == "conserved"]["log_mean"].dropna()
        a2b  = sub[sub.compartment_type == "A_to_B"]["log_mean"].dropna()
        if len(cons) > 10 and len(a2b) > 10:
            stat, p = mannwhitneyu(cons, a2b, alternative="two-sided")
            print(f"  [{comp}] MWU log_mean(conserved vs A_to_B): U={stat:.0f}, p={p:.3g}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "summary_by_compartment.tsv"),
                              sep="\t", index=False)

    # ---- summary 2: by tad_conserved_10000 ----
    rows = []
    for comp in ("Asu_Ath", "Asu_Aar"):
        sub = merged[(merged.comparison == comp) & merged.mean_tpm.notna()]
        for tcons, g in sub.groupby("tad_conserved_10000"):
            rows.append({
                "comparison": comp,
                "tad_conserved_10000": int(tcons),
                "n_genes": len(g),
                "mean_tpm":   round(g["mean_tpm"].mean(), 4),
                "median_tpm": round(g["mean_tpm"].median(), 4),
                "mean_log":   round(g["log_mean"].mean(), 4),
                "expressed_rate": round(g["expressed"].mean(), 4),
            })
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "summary_by_tad.tsv"),
                              sep="\t", index=False)

    # ---- summary 3: matrix_corr 分箱 vs 表达 ----
    rows = []
    bins = [-1, 0.0, 0.3, 0.6, 0.8, 1.01]
    labels = ["<0", "0-0.3", "0.3-0.6", "0.6-0.8", ">0.8"]
    for comp in ("Asu_Ath", "Asu_Aar"):
        sub = merged[(merged.comparison == comp)
                     & merged.mean_tpm.notna()
                     & merged.matrix_corr.notna()].copy()
        sub["mc_bin"] = pd.cut(sub["matrix_corr"], bins=bins, labels=labels)
        for lab, g in sub.groupby("mc_bin", observed=True):
            rows.append({
                "comparison": comp,
                "matrix_corr_bin": str(lab),
                "n_genes": len(g),
                "mean_tpm":   round(g["mean_tpm"].mean(), 4),
                "median_tpm": round(g["mean_tpm"].median(), 4),
                "mean_log":   round(g["log_mean"].mean(), 4),
                "expressed_rate": round(g["expressed"].mean(), 4),
            })
        rho, p = spearmanr(sub["matrix_corr"], sub["log_mean"], nan_policy="omit")
        print(f"  [{comp}] Spearman matrix_corr × log_mean: rho={rho:.4f}, p={p:.3g}, n={len(sub)}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "summary_matrix_corr_bins.tsv"),
                              sep="\t", index=False)

    print(f"  保存: {OUT_DIR}/")


if __name__ == "__main__":
    main()
