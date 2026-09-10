#!/usr/bin/env python3
"""04b_2kb_sensitivity.py — T6 检验: 2kb 分辨率下 R2/R3/R4 matrix_corr 结论是否成立.

审稿人指出 10kb ±50kb matrix_corr 退化为 compartment 尺度 (median 0.87-0.90,
方差压缩), 2kb (median 0.66) 才有区分度. 稿件 R2 (TE density vs matrix rho=-0.65),
R3(vi) (matrix_corr vs |log2FC| |rho|<0.02), R4(iii) (env vs non-env matrix_corr
p=3e-6) 全部基于 10kb. 若 2kb 下结论改变, 这些需重写.

本脚本对 Asu_Ath / Asu_Aar 各做四项检验, 每项同时算 2kb (04/2000/) 与 10kb 基线:

  R2  : TE density (te_coverage_frac) vs matrix_corr            (Spearman)
  R4  : env vs non-env matrix_corr                               (ttest + MWU)
  R3vi: matrix_corr vs |log2FC|, 分 env / non-env                (Spearman)
  CDS : region_cds_coverage vs matrix_corr                       (Spearman)

输出: results/gene_3d/2000/sensitivity_2kb_summary.tsv
  列: comparison, test, resolution, n, statistic, p,
      env_mean_or_na, nonenv_mean_or_na, verdict

verdict 判定 (仅 2kb 行; 10kb 行标 "baseline"):
  R2   : stable  if rho_2kb < 0 and |rho_2kb| > 0.3  (强负相关保持)
  R4   : stable  if delta_2kb > 0 and p_2kb < 0.05   (env 仍显著高于 non-env)
  R3vi : stable  if |rho_2kb| < 0.05                  (弱耦合保持)
  CDS  : stable  if rho_2kb > 0 and rho_2kb > 0.2     (强正相关保持)

禁止修改输入文件. 脚本可重复运行.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind, mannwhitneyu

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
P_2KB = {
    "Asu_Ath": f"{BASE}/04.cds_synteny_3d_conservation/results/gene_3d/2000/gene_3d_conservation_Asu_Ath.tsv",
    "Asu_Aar": f"{BASE}/04.cds_synteny_3d_conservation/results/gene_3d/2000/gene_3d_conservation_Asu_Aar.tsv",
}
P_TE = f"{BASE}/06.env_genes_3d_conservation/results/gene_3d/env_gene_te_density.tsv"
P_EXPR = {
    "Asu_Ath": f"{BASE}/10.env_3d_expression_integration/results/env_3d_expr/per_env_gene_3d_expr_Asu_Ath.tsv",
    "Asu_Aar": f"{BASE}/10.env_3d_expression_integration/results/env_3d_expr/per_env_gene_3d_expr_Asu_Aar.tsv",
}
P_CDS = f"{BASE}/04.cds_synteny_3d_conservation/results/density_link/cds_density_vs_matrix.tsv"
OUT = f"{BASE}/04.cds_synteny_3d_conservation/results/gene_3d/2000/sensitivity_2kb_summary.tsv"

# Manuscript quoted 10kb baselines (for cross-check; computed values may differ slightly)
MANUSCRIPT_R2 = {"Asu_Ath": -0.507, "Asu_Aar": -0.626}
MANUSCRIPT_R4 = {
    "Asu_Ath": {"p": 0.006, "delta": None},
    "Asu_Aar": {"p": 3e-6, "delta": 0.041},
}
MANUSCRIPT_CDS = {"Asu_Ath": 0.386, "Asu_Aar": 0.451}


def safe_spearman(x, y):
    """Spearman on non-null pairs. Returns (rho, p, n)."""
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 3:
        return np.nan, np.nan, n
    rho, p = spearmanr(x[mask], y[mask])
    if np.isnan(rho):  # zero variance
        return np.nan, np.nan, n
    return float(rho), float(p), n


def safe_ttest(env_vals, nonenv_vals):
    """Welch t-test + Mann-Whitney U. Returns (env_mean, nonenv_mean, delta, p_t, p_mwu, n_env, n_nonenv)."""
    env_vals = env_vals.dropna()
    nonenv_vals = nonenv_vals.dropna()
    n_env = len(env_vals)
    n_nonenv = len(nonenv_vals)
    if n_env < 2 or n_nonenv < 2:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, n_env, n_nonenv)
    t, p_t = ttest_ind(env_vals, nonenv_vals, equal_var=False)
    u, p_mwu = mannwhitneyu(env_vals, nonenv_vals, alternative="two-sided")
    return (float(env_vals.mean()), float(nonenv_vals.mean()),
            float(env_vals.mean() - nonenv_vals.mean()),
            float(p_t), float(p_mwu), n_env, n_nonenv)


def fmt(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    te = pd.read_csv(P_TE, sep="\t")
    cds_all = pd.read_csv(P_CDS, sep="\t")
    rows = []  # output rows

    for cmp in ["Asu_Ath", "Asu_Aar"]:
        print(f"\n{'='*70}")
        print(f"=== {cmp} ===")
        print(f"{'='*70}")

        # ---- 1. Merge ----
        m2kb = pd.read_csv(P_2KB[cmp], sep="\t")[["gene_id", "matrix_corr"]].rename(
            columns={"matrix_corr": "matrix_corr_2kb"})
        expr = pd.read_csv(P_EXPR[cmp], sep="\t")[
            ["gene_id", "matrix_corr", "is_env", "abs_log2fc"]].rename(
            columns={"matrix_corr": "matrix_corr_10kb"})
        df = (m2kb.merge(expr, on="gene_id", how="inner")
                  .merge(te[["gene_id", "te_coverage_frac"]], on="gene_id", how="left"))
        print(f"Merged n (gene_id intersection 04/2000 + 10/env_3d_expr): {len(df)}")
        print(f"  matrix_corr_2kb  non-null: {df['matrix_corr_2kb'].notna().sum()}")
        print(f"  matrix_corr_10kb non-null: {df['matrix_corr_10kb'].notna().sum()}")
        print(f"  te_coverage_frac non-null: {df['te_coverage_frac'].notna().sum()}")
        print(f"  abs_log2fc       non-null: {df['abs_log2fc'].notna().sum()}")
        print(f"  is_env=1: {int((df['is_env']==1).sum())}  is_env=0: {int((df['is_env']==0).sum())}")
        print(f"  median matrix_corr_2kb : {df['matrix_corr_2kb'].median():.4f}")
        print(f"  median matrix_corr_10kb: {df['matrix_corr_10kb'].median():.4f}")
        print(f"  IQR   matrix_corr_2kb : [{df['matrix_corr_2kb'].quantile(.25):.4f}, {df['matrix_corr_2kb'].quantile(.75):.4f}]")
        print(f"  IQR   matrix_corr_10kb: [{df['matrix_corr_10kb'].quantile(.25):.4f}, {df['matrix_corr_10kb'].quantile(.75):.4f}]")

        # ---- 2. R2: TE density vs matrix_corr (Spearman) ----
        print(f"\n--- R2: TE density vs matrix_corr (Spearman) ---")
        print(f"  manuscript 10kb baseline: rho={MANUSCRIPT_R2[cmp]}")
        rho_10, p_10, n_10 = safe_spearman(df["te_coverage_frac"], df["matrix_corr_10kb"])
        rho_2, p_2, n_2 = safe_spearman(df["te_coverage_frac"], df["matrix_corr_2kb"])
        print(f"  10kb: rho={rho_10:.4f}  p={p_10:.2e}  n={n_10}")
        print(f"  2kb : rho={rho_2:.4f}  p={p_2:.2e}  n={n_2}")
        verdict_r2 = "stable" if (rho_2 < 0 and abs(rho_2) > 0.3) else "changed"
        print(f"  verdict(2kb): {verdict_r2}  (criterion: rho<0 & |rho|>0.3)")
        rows.append(dict(comparison=cmp, test="R2_TE_vs_matrix", resolution="10kb",
                         n=n_10, statistic=rho_10, p=p_10,
                         env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict="baseline"))
        rows.append(dict(comparison=cmp, test="R2_TE_vs_matrix", resolution="2kb",
                         n=n_2, statistic=rho_2, p=p_2,
                         env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict=verdict_r2))

        # ---- 3. R4: env vs non-env matrix_corr (ttest + MWU) ----
        print(f"\n--- R4: env vs non-env matrix_corr (ttest + MWU) ---")
        print(f"  manuscript 10kb baseline: p={MANUSCRIPT_R4[cmp]['p']}, delta={MANUSCRIPT_R4[cmp]['delta']}")
        for res, col in [("10kb", "matrix_corr_10kb"), ("2kb", "matrix_corr_2kb")]:
            env_v = df.loc[df["is_env"] == 1, col]
            nonenv_v = df.loc[df["is_env"] == 0, col]
            em, nm, delta, p_t, p_mwu, n_e, n_ne = safe_ttest(env_v, nonenv_v)
            print(f"  {res}: env_mean={em:.4f}  nonenv_mean={nm:.4f}  delta={delta:.4f}  "
                  f"p_ttest={p_t:.2e}  p_mwu={p_mwu:.2e}  n_env={n_e}  n_nonenv={n_ne}")
            if res == "10kb":
                verdict_r4 = "baseline"
            else:
                verdict_r4 = "stable" if (delta > 0 and p_t < 0.05) else "changed"
                print(f"  verdict(2kb): {verdict_r4}  (criterion: delta>0 & p<0.05)")
            rows.append(dict(comparison=cmp, test="R4_env_vs_nonenv_ttest", resolution=res,
                             n=n_e + n_ne, statistic=delta, p=p_t,
                             env_mean_or_na=em, nonenv_mean_or_na=nm, verdict=verdict_r4))
            rows.append(dict(comparison=cmp, test="R4_env_vs_nonenv_mwu", resolution=res,
                             n=n_e + n_ne, statistic=delta, p=p_mwu,
                             env_mean_or_na=em, nonenv_mean_or_na=nm, verdict=verdict_r4))

        # ---- 4. R3vi: matrix_corr vs |log2FC| (Spearman, split env/non-env) ----
        print(f"\n--- R3vi: matrix_corr vs |log2FC| (Spearman, split env/non-env) ---")
        print(f"  manuscript 10kb baseline: |rho|<0.02")
        for grp, mask in [("env", df["is_env"] == 1), ("nonenv", df["is_env"] == 0)]:
            sub = df.loc[mask]
            for res, col in [("10kb", "matrix_corr_10kb"), ("2kb", "matrix_corr_2kb")]:
                rho, p, n = safe_spearman(sub[col], sub["abs_log2fc"])
                print(f"  {grp} {res}: rho={rho:.4f}  p={p:.2e}  n={n}")
                if res == "10kb":
                    verdict_r3 = "baseline"
                else:
                    verdict_r3 = "stable" if abs(rho) < 0.05 else "changed"
                    if res == "2kb":
                        print(f"    verdict(2kb): {verdict_r3}  (criterion: |rho|<0.05)")
                rows.append(dict(comparison=cmp, test=f"R3vi_matrix_vs_log2fc_{grp}",
                                 resolution=res, n=n, statistic=rho, p=p,
                                 env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict=verdict_r3))

        # ---- 5. CDS density vs matrix_corr (Spearman) ----
        print(f"\n--- CDS: region_cds_coverage vs matrix_corr (Spearman) ---")
        print(f"  manuscript 10kb baseline: rho={MANUSCRIPT_CDS[cmp]}")
        cds_cmp = cds_all[cds_all["comparison"] == cmp][
            ["gene_id", "matrix_corr", "region_cds_coverage"]].rename(
            columns={"matrix_corr": "matrix_corr_10kb_cds"})
        cdf = m2kb.merge(cds_cmp, on="gene_id", how="inner")
        print(f"  CDS merged n: {len(cdf)}  (cds_density_vs_matrix has its own 10kb matrix_corr)")
        rho_10c, p_10c, n_10c = safe_spearman(cdf["region_cds_coverage"], cdf["matrix_corr_10kb_cds"])
        rho_2c, p_2c, n_2c = safe_spearman(cdf["region_cds_coverage"], cdf["matrix_corr_2kb"])
        print(f"  10kb: rho={rho_10c:.4f}  p={p_10c:.2e}  n={n_10c}")
        print(f"  2kb : rho={rho_2c:.4f}  p={p_2c:.2e}  n={n_2c}")
        verdict_cds = "stable" if (rho_2c > 0 and rho_2c > 0.2) else "changed"
        print(f"  verdict(2kb): {verdict_cds}  (criterion: rho>0 & rho>0.2)")
        rows.append(dict(comparison=cmp, test="CDS_density_vs_matrix", resolution="10kb",
                         n=n_10c, statistic=rho_10c, p=p_10c,
                         env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict="baseline"))
        rows.append(dict(comparison=cmp, test="CDS_density_vs_matrix", resolution="2kb",
                         n=n_2c, statistic=rho_2c, p=p_2c,
                         env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict=verdict_cds))

    # ---- Output TSV ----
    out = pd.DataFrame(rows, columns=[
        "comparison", "test", "resolution", "n", "statistic", "p",
        "env_mean_or_na", "nonenv_mean_or_na", "verdict"])
    out.to_csv(OUT, sep="\t", index=False, na_rep="NA")
    print(f"\n{'='*70}")
    print(f"Wrote {OUT}  ({len(out)} rows)")
    print(f"{'='*70}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
