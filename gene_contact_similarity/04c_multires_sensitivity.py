#!/usr/bin/env python3
"""04c_multires_sensitivity.py — 4 分辨率 (2/5/10/25kb) per-gene matrix_corr 敏感性分析.

补强 R2/R3(vi)/R4(iii) 的 matrix 关联结论. 04b 只有 2kb vs 10kb, 本脚本新增
5kb 与 25kb, 得 4 分辨率曲线. per-gene matrix_corr 是基因区 ±50kb 子矩阵对齐
Pearson (非 SCC).

Step 1: import 03_gene_region_3d_conservation.py 的 run_matrix_resolution,
        跑 (5000, 5000) 与 (25000, 25000). 不改原脚本 MATRIX_CONFIGS.
        产出 04/results/gene_3d/{5000,25000}/gene_3d_conservation_Asu_{Ath,Aar}.tsv
        若已存在则跳过.

Step 2: 对 2/5/10/25kb 四个分辨率 × Asu_Ath/Asu_Aar 两 comparison, 各做 4 项检验:
  R2  : TE density (te_coverage_frac) vs matrix_corr       (Spearman rho+p+n)
  CDS : region_cds_coverage vs matrix_corr                 (Spearman rho+p+n)
  R4  : env vs non-env matrix_corr                         (Welch ttest delta+p+n_env+n_nonenv)
  R3vi: matrix_corr vs |log2FC| 分 env/non-env             (Spearman rho+p+n)

Step 3: 输出 sensitivity_multires_summary.tsv
  列: resolution, comparison, test, n, statistic, p,
      env_mean_or_na, nonenv_mean_or_na, verdict

verdict 判定 (10kb = baseline; 2/5/25kb = stable/changed):
  R2   : stable if rho < 0 and |rho| > 0.3
  CDS  : stable if rho > 0 and rho > 0.2
  R4   : stable if delta > 0 and p < 0.05
  R3vi : stable if |rho| < 0.05

不修改原 2kb/10kb TSV, 只新增 5kb/25kb + 汇总 TSV. 脚本可重复运行.
"""
import os
import sys
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
SCRIPTS = f"{BASE}/04.cds_synteny_3d_conservation/scripts"
GENE3D_RESULTS = f"{BASE}/04.cds_synteny_3d_conservation/results/gene_3d"

# Per-gene matrix_corr TSVs for each resolution
P_RES = {
    2000: {
        "Asu_Ath": f"{GENE3D_RESULTS}/2000/gene_3d_conservation_Asu_Ath.tsv",
        "Asu_Aar": f"{GENE3D_RESULTS}/2000/gene_3d_conservation_Asu_Aar.tsv",
    },
    5000: {
        "Asu_Ath": f"{GENE3D_RESULTS}/5000/gene_3d_conservation_Asu_Ath.tsv",
        "Asu_Aar": f"{GENE3D_RESULTS}/5000/gene_3d_conservation_Asu_Aar.tsv",
    },
    10000: {
        "Asu_Ath": f"{GENE3D_RESULTS}/10000/gene_3d_conservation_Asu_Ath.tsv",
        "Asu_Aar": f"{GENE3D_RESULTS}/10000/gene_3d_conservation_Asu_Aar.tsv",
    },
    25000: {
        "Asu_Ath": f"{GENE3D_RESULTS}/25000/gene_3d_conservation_Asu_Ath.tsv",
        "Asu_Aar": f"{GENE3D_RESULTS}/25000/gene_3d_conservation_Asu_Aar.tsv",
    },
}

# Covariate files
P_TE = f"{BASE}/06.env_genes_3d_conservation/results/gene_3d/env_gene_te_density.tsv"
P_EXPR = {
    "Asu_Ath": f"{BASE}/10.env_3d_expression_integration/results/env_3d_expr/per_env_gene_3d_expr_Asu_Ath.tsv",
    "Asu_Aar": f"{BASE}/10.env_3d_expression_integration/results/env_3d_expr/per_env_gene_3d_expr_Asu_Aar.tsv",
}
P_CDS = f"{BASE}/04.cds_synteny_3d_conservation/results/density_link/cds_density_vs_matrix.tsv"

OUT = f"{GENE3D_RESULTS}/sensitivity_multires_summary.tsv"

RESOLUTIONS = [2000, 5000, 10000, 25000]
RES_LABEL = {2000: "2kb", 5000: "5kb", 10000: "10kb", 25000: "25kb"}
BASELINE_RES = 10000
COMPARISONS = ["Asu_Ath", "Asu_Aar"]


def safe_spearman(x, y):
    """Spearman on non-null pairs. Returns (rho, p, n)."""
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 3:
        return np.nan, np.nan, n
    rho, p = spearmanr(x[mask], y[mask])
    if np.isnan(rho):
        return np.nan, np.nan, n
    return float(rho), float(p), n


def safe_ttest(env_vals, nonenv_vals):
    """Welch t-test. Returns (env_mean, nonenv_mean, delta, p, n_env, n_nonenv)."""
    env_vals = env_vals.dropna()
    nonenv_vals = nonenv_vals.dropna()
    n_env = len(env_vals)
    n_nonenv = len(nonenv_vals)
    if n_env < 2 or n_nonenv < 2:
        return (np.nan, np.nan, np.nan, np.nan, n_env, n_nonenv)
    t, p = ttest_ind(env_vals, nonenv_vals, equal_var=False)
    return (float(env_vals.mean()), float(nonenv_vals.mean()),
            float(env_vals.mean() - nonenv_vals.mean()),
            float(p), n_env, n_nonenv)


def fmt(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def step1_run_new_resolutions():
    """Run 5kb and 25kb per-gene matrix_corr if outputs missing.

    Uses importlib to load 03_gene_region_3d_conservation.py without modifying
    its MATRIX_CONFIGS. Calls run_matrix_resolution(5000,5000) and (25000,25000).
    """
    spec = importlib.util.spec_from_file_location(
        "gene3d_03", f"{SCRIPTS}/03_gene_region_3d_conservation.py")
    gene3d = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gene3d)

    for matrix_res, orth_res in [(5000, 5000), (25000, 25000)]:
        paths = [P_RES[matrix_res][c] for c in COMPARISONS]
        if all(os.path.exists(p) for p in paths):
            print(f"[Step 1] matrix_res={matrix_res} 已存在, 跳过")
            continue
        print(f"[Step 1] 运行 matrix_res={matrix_res} orth_res={orth_res} ...")
        gene3d.run_matrix_resolution(matrix_res, orth_res)


def load_covariates(cmp):
    """Load TE density, expr (is_env + abs_log2fc), CDS density for a comparison."""
    te = pd.read_csv(P_TE, sep="\t")[["gene_id", "te_coverage_frac"]]
    expr = pd.read_csv(P_EXPR[cmp], sep="\t")[["gene_id", "is_env", "abs_log2fc"]]
    cds_all = pd.read_csv(P_CDS, sep="\t")
    cds = cds_all[cds_all["comparison"] == cmp][["gene_id", "region_cds_coverage"]]
    return te, expr, cds


def run_tests_for_resolution(cmp, resolution, te, expr, cds):
    """Run 4 tests for one (comparison, resolution). Returns list of row dicts."""
    res_label = RES_LABEL[resolution]
    is_baseline = (resolution == BASELINE_RES)
    col = f"matrix_corr_{res_label}"

    # Per-gene matrix_corr at this resolution
    df_res = pd.read_csv(P_RES[resolution][cmp], sep="\t")[["gene_id", "matrix_corr"]]
    df_res = df_res.rename(columns={"matrix_corr": col})

    # Merge with expr (is_env, abs_log2fc) and TE density
    df = (df_res.merge(expr, on="gene_id", how="inner")
                .merge(te, on="gene_id", how="left"))

    rows = []

    # ---- 1. R2: TE density vs matrix_corr (Spearman) ----
    rho, p, n = safe_spearman(df["te_coverage_frac"], df[col])
    verdict = "baseline" if is_baseline else (
        "stable" if (rho < 0 and abs(rho) > 0.3) else "changed")
    rows.append(dict(
        resolution=res_label, comparison=cmp, test="R2_TE_vs_matrix",
        n=n, statistic=rho, p=p,
        env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict=verdict))

    # ---- 2. CDS: region_cds_coverage vs matrix_corr (Spearman) ----
    cdf = df_res.merge(cds, on="gene_id", how="inner")
    rho_c, p_c, n_c = safe_spearman(cdf["region_cds_coverage"], cdf[col])
    verdict_c = "baseline" if is_baseline else (
        "stable" if (rho_c > 0 and rho_c > 0.2) else "changed")
    rows.append(dict(
        resolution=res_label, comparison=cmp, test="CDS_density_vs_matrix",
        n=n_c, statistic=rho_c, p=p_c,
        env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict=verdict_c))

    # ---- 3. R4: env vs non-env matrix_corr (Welch ttest) ----
    env_v = df.loc[df["is_env"] == 1, col]
    nonenv_v = df.loc[df["is_env"] == 0, col]
    em, nm, delta, p_t, n_e, n_ne = safe_ttest(env_v, nonenv_v)
    verdict_r4 = "baseline" if is_baseline else (
        "stable" if (delta > 0 and p_t < 0.05) else "changed")
    rows.append(dict(
        resolution=res_label, comparison=cmp, test="R4_env_vs_nonenv_ttest",
        n=n_e + n_ne, statistic=delta, p=p_t,
        env_mean_or_na=em, nonenv_mean_or_na=nm, verdict=verdict_r4))

    # ---- 4. R3vi: matrix_corr vs |log2FC| (Spearman, split env/non-env) ----
    for grp, mask in [("env", df["is_env"] == 1), ("nonenv", df["is_env"] == 0)]:
        sub = df.loc[mask]
        rho_r3, p_r3, n_r3 = safe_spearman(sub[col], sub["abs_log2fc"])
        verdict_r3 = "baseline" if is_baseline else (
            "stable" if abs(rho_r3) < 0.05 else "changed")
        rows.append(dict(
            resolution=res_label, comparison=cmp,
            test=f"R3vi_matrix_vs_log2fc_{grp}",
            n=n_r3, statistic=rho_r3, p=p_r3,
            env_mean_or_na="NA", nonenv_mean_or_na="NA", verdict=verdict_r3))

    return rows


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    # Step 1: produce 5kb/25kb per-gene TSVs if missing
    print(f"\n{'='*70}")
    print("Step 1: 5kb/25kb per-gene matrix_corr")
    print(f"{'='*70}")
    step1_run_new_resolutions()

    # Step 2-3: 4 resolutions x 2 comparisons x 4 tests
    print(f"\n{'='*70}")
    print("Step 2-3: 4 分辨率敏感性检验")
    print(f"{'='*70}")
    all_rows = []
    for cmp in COMPARISONS:
        te, expr, cds = load_covariates(cmp)
        for res in RESOLUTIONS:
            print(f"\n=== {cmp} @ {RES_LABEL[res]} ===")
            rows = run_tests_for_resolution(cmp, res, te, expr, cds)
            for r in rows:
                print(f"  {r['test']:<35} n={r['n']:<6} "
                      f"stat={fmt(r['statistic']):<14} p={fmt(r['p']):<14} "
                      f"verdict={r['verdict']}")
            all_rows.extend(rows)

    out = pd.DataFrame(all_rows, columns=[
        "resolution", "comparison", "test", "n", "statistic", "p",
        "env_mean_or_na", "nonenv_mean_or_na", "verdict"])
    out.to_csv(OUT, sep="\t", index=False, na_rep="NA")
    print(f"\n{'='*70}")
    print(f"Wrote {OUT}  ({len(out)} rows)")
    print(f"{'='*70}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
