#!/usr/bin/env python3
"""
07b_expression_trajectory_parental.py — per-side vs-parent 表达分析 (共线落点统一版, Z-score 口径).

同源对应一律采用与 matrix_corr 一致的共线落点映射 (见
04.cds_synteny_3d_conservation): 据 nat gene_3d 表中 Asu 基因的亲本共线坐标
(par_chr/par_start/par_end) 取中点, 在亲本基因 BED 中找覆盖该中点的基因作为
表达侧亲本对应. 不使用任何序列 ortholog/RBH; 亚基因组 dominance (sT vs sA)
是另一分析, 保留 sT-sA 配对, 不在此脚本.

每侧单一有效基因池 (nat Asu + syn Asu + 亲本三边 mean TPM >= 0.01),
在该池内对每阶段 log2fc 跨基因做 robust Z-score 标准化:
  Delta_stage = mean log1p(TPM)_Asu_stage - mean log1p(TPM)_par
  z_stage = (Delta_stage - median(Delta_stage)) / std(Delta_stage)

(1) per-side vs-parent bias (图 5E/F):
    |z_stage| > 1 -> biased, 否则 close-to-parent. 两阶段同源判据, 只是分别取边际.
    阈值扫描 |z| = 0.5 / 1.0 / 1.5 / 2.0.

(2) 表达轨迹 4 类 (图 5G):
    (nat, syn) 二分组合:
      Maintained (双 |z|<=1), Remodeled (双 |z|>1),
      Synthetic-shift (|z_syn|>1, |z_nat|<=1),
      Natural-shift (|z_nat|>1, |z_syn|<=1).

输出 results/ortholog_expr_allo738/:
  expression_bias_per_side.tsv        (长格式: gene, side, prefix, par_gene, log2fc, z, bias_sig)
  expression_bias_per_side_summary.tsv
  expression_bias_threshold_sensitivity.tsv
  expression_trajectory_parental.tsv  (side, gene, par_gene, delta_nat, delta_syn,
                                       z_nat, z_syn, expr_class)
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT_BASE  # noqa: E402

OUT_DIR = os.path.join(OUT_BASE, "ortholog_expr_allo738")
# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
GENE3D = f"{BASE}/04.cds_synteny_3d_conservation/results/gene_3d/10000"
BED = f"{BASE}/04.cds_synteny_3d_conservation/results/cds"
TPM = f"{BASE}/08.expression/04.count/matrix"

Z_THRESH = 1.0
TPM_MIN = 0.01
Z_SCAN = (0.5, 1.0, 1.5, 2.0)

EXPR_TRAJ_ORDER = ["Maintained", "Remodeled", "Synthetic-shift", "Natural-shift"]


def log(msg):
    print(msg, flush=True)


def load_tpm(species):
    """返回 (gene 列表, log1p(tpm) 矩阵 [gene x rep], mean tpm 向量, {gene: row idx})."""
    df = pd.read_csv(f"{TPM}/{species}.tpm.tsv", sep="\t")
    gid = df.columns[0]
    genes = df[gid].values
    num = [c for c in df.columns if c != gid]
    v = df[num].values.astype(float)
    return genes, np.log1p(v), v.mean(axis=1), {g: i for i, g in enumerate(genes)}


def gene_bed(species):
    bed = pd.read_csv(f"{BED}/{species}_genes.bed", sep="\t", header=None,
                      usecols=[0, 1, 2, 3])
    bed.columns = ["chrom", "start", "end", "gid"]
    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)
    return {ch: g for ch, g in bed.groupby("chrom")}


def parent_gene_at(bed_idx, chrom, mid):
    if chrom not in bed_idx:
        return None
    g = bed_idx[chrom]
    starts = g["start"].values
    ends = g["end"].values
    pos = np.searchsorted(starts, mid, side="right") - 1
    if pos < 0:
        return None
    if starts[pos] <= mid <= ends[pos]:
        return g["gid"].values[pos]
    return None


def robust_z(x):
    x = np.asarray(x, dtype=float)
    s = x.std(ddof=0)
    return (x - np.median(x)) / s if s > 0 else np.zeros_like(x)


def classify(z_nat, z_syn, thresh=Z_THRESH):
    if abs(z_nat) <= thresh and abs(z_syn) <= thresh:
        return "Maintained"
    if abs(z_nat) > thresh and abs(z_syn) > thresh:
        return "Remodeled"
    if abs(z_nat) <= thresh < abs(z_syn):
        return "Synthetic-shift"
    if abs(z_syn) <= thresh < abs(z_nat):
        return "Natural-shift"
    return "Remodeled"


def build_side(side, par_species, want_prefix,
               nat_log, nat_mean, nat_g2r,
               syn_log, syn_mean, syn_g2r,
               par_log, par_mean, par_g2r, bed_idx):
    """从 nat gene_3d 取某侧 Asu 基因 + 共线落点亲本基因, 三边 TPM>=0.01 池."""
    comp = f"Asu_{par_species}"
    g3d = pd.read_csv(f"{GENE3D}/gene_3d_conservation_{comp}.tsv", sep="\t")
    g3d = g3d[g3d["asu_chr"].str.startswith(want_prefix)].copy()
    recs = []
    for _, r in g3d.iterrows():
        asu_gene = r["gene_id"]
        if asu_gene not in nat_g2r or asu_gene not in syn_g2r:
            continue
        mid = (r["par_start"] + r["par_end"]) / 2.0
        par_gene = parent_gene_at(bed_idx, r["par_chr"], mid)
        if par_gene is None or par_gene not in par_g2r:
            continue
        recs.append((asu_gene, par_gene))
    df = pd.DataFrame(recs, columns=["gene", "par_gene"]).drop_duplicates("gene")
    ai_nat = np.array([nat_g2r[g] for g in df["gene"]])
    ai_syn = np.array([syn_g2r[g] for g in df["gene"]])
    pi = np.array([par_g2r[pg] for pg in df["par_gene"]])
    keep = ((nat_mean[ai_nat] >= TPM_MIN) &
            (syn_mean[ai_syn] >= TPM_MIN) &
            (par_mean[pi] >= TPM_MIN))
    return df.iloc[np.where(keep)[0]].reset_index(drop=True)


def main():
    log("=== 07b per-side vs-parent (共线落点, TPM log1p, Z-score bias/trajectory) ===")
    os.makedirs(OUT_DIR, exist_ok=True)

    _, nat_log, nat_mean, nat_g2r = load_tpm("Asu")
    _, syn_log, syn_mean, syn_g2r = load_tpm("Allo738")
    _, ath_log, ath_mean, ath_g2r = load_tpm("Ath")
    _, aar_log, aar_mean, aar_g2r = load_tpm("Aar")
    ath_bed = gene_bed("Ath")
    aar_bed = gene_bed("Aar")

    sides = [
        ("sT", "Ath", "sT", ath_bed, ath_log, ath_mean, ath_g2r),
        ("sA", "Aar", "sA", aar_bed, aar_log, aar_mean, aar_g2r),
    ]

    bias_rows = []
    summary_rows = []
    traj_rows = []
    sens_rows = []
    chi_rows = []

    for side, par_sp, want, bed_idx, par_log, par_mean, par_g2r in sides:
        pool = build_side(side, par_sp, want,
                          nat_log, nat_mean, nat_g2r,
                          syn_log, syn_mean, syn_g2r,
                          par_log, par_mean, par_g2r, bed_idx)
        n_pool = len(pool)
        log(f"\n  [{side}] 有效池 (nat+syn+par 三边 TPM>={TPM_MIN}): {n_pool:,}")

        ai_nat = np.array([nat_g2r[g] for g in pool["gene"]])
        ai_syn = np.array([syn_g2r[g] for g in pool["gene"]])
        pi = np.array([par_g2r[pg] for pg in pool["par_gene"]])
        d_nat = nat_log[ai_nat].mean(axis=1) - par_log[pi].mean(axis=1)
        d_syn = syn_log[ai_syn].mean(axis=1) - par_log[pi].mean(axis=1)
        z_nat = robust_z(d_nat)
        z_syn = robust_z(d_syn)
        l2fc_nat = d_nat / np.log(2)
        l2fc_syn = d_syn / np.log(2)

        # (1) per-side bias: |z|>1
        for prefix, l2fc, z in [("nat", l2fc_nat, z_nat), ("syn", l2fc_syn, z_syn)]:
            sig = np.abs(z) > Z_THRESH
            n_bias = int(sig.sum())
            summary_rows.append({
                "side": side, "sample": prefix, "n_total": n_pool,
                "n_biased": n_bias, "pct_biased": round(n_bias / n_pool * 100, 2),
                "n_balanced": n_pool - n_bias,
                "pct_balanced": round((n_pool - n_bias) / n_pool * 100, 2),
            })
            log(f"    {prefix} {side}: n={n_pool:,} biased(|z|>{Z_THRESH})={n_bias:,} "
                f"({n_bias/n_pool*100:.1f}%)")
            for i in range(n_pool):
                bias_rows.append({
                    "gene": pool["gene"].iloc[i], "side": side, "prefix": prefix,
                    "par_gene": pool["par_gene"].iloc[i],
                    "log2fc": l2fc[i], "z": z[i], "bias_sig": int(sig[i]),
                })

        # 阈值敏感性
        for thr in Z_SCAN:
            for prefix, z in [("nat", z_nat), ("syn", z_syn)]:
                nb = int((np.abs(z) > thr).sum())
                sens_rows.append({
                    "threshold": thr, "sample": prefix, "side": side,
                    "n_total": n_pool, "n_biased": nb,
                    "pct_biased": round(nb / n_pool * 100, 2),
                })

        # sT vs sA 在每阶段的 χ², 在外层第二轮有 sA 后计算, 此处先存 z
        # (2) trajectory 4 类, 主阈值 1.0
        cls = [classify(z_nat[i], z_syn[i]) for i in range(n_pool)]
        for i in range(n_pool):
            traj_rows.append({
                "side": side, "gene": pool["gene"].iloc[i],
                "par_gene": pool["par_gene"].iloc[i],
                "delta_nat": d_nat[i], "delta_syn": d_syn[i],
                "z_nat": z_nat[i], "z_syn": z_syn[i], "expr_class": cls[i],
            })
        vc = pd.Series(cls).value_counts().reindex(EXPR_TRAJ_ORDER, fill_value=0)
        log(f"    [{side}] 表达轨迹 4 类:")
        for cat in EXPR_TRAJ_ORDER:
            log(f"      {cat}: {int(vc[cat]):,} ({vc[cat]/n_pool*100:.1f}%)")

    df_bias = pd.DataFrame(bias_rows)
    df_bias.to_csv(f"{OUT_DIR}/expression_bias_per_side.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(
        f"{OUT_DIR}/expression_bias_per_side_summary.tsv", sep="\t", index=False)
    pd.DataFrame(sens_rows).to_csv(
        f"{OUT_DIR}/expression_bias_threshold_sensitivity.tsv", sep="\t", index=False)

    # χ²: sT vs sA biased/balanced per sample
    summary = pd.DataFrame(summary_rows)
    from scipy.stats import chi2_contingency
    for prefix in ("nat", "syn"):
        st = summary[(summary.side == "sT") & (summary["sample"] == prefix)].iloc[0]
        sa = summary[(summary.side == "sA") & (summary["sample"] == prefix)].iloc[0]
        cont = np.array([[st["n_biased"], st["n_balanced"]],
                         [sa["n_biased"], sa["n_balanced"]]])
        chi2, p, dof, _ = chi2_contingency(cont)
        N = int(cont.sum())
        V = float(np.sqrt(chi2 / N)) if N else float("nan")
        chi_rows.append({"test": f"chi2_{prefix}_sT_vs_sA",
                         "chi2": round(chi2, 2), "dof": int(dof), "p": p,
                         "cramers_v": round(V, 4), "n": N})
        log(f"\n  {prefix} sT vs sA bias: χ²={chi2:.2f} p={p:.2e} V={V:.4f}")
    pd.DataFrame(chi_rows).to_csv(
        f"{OUT_DIR}/expression_bias_per_side_chi2.tsv", sep="\t", index=False)

    df_traj = pd.DataFrame(traj_rows)
    df_traj.to_csv(f"{OUT_DIR}/expression_trajectory_parental.tsv",
                   sep="\t", index=False)

    log(f"\n  保存:")
    log(f"    {OUT_DIR}/expression_bias_per_side.tsv")
    log(f"    {OUT_DIR}/expression_bias_per_side_summary.tsv")
    log(f"    {OUT_DIR}/expression_bias_threshold_sensitivity.tsv")
    log(f"    {OUT_DIR}/expression_bias_per_side_chi2.tsv")
    log(f"    {OUT_DIR}/expression_trajectory_parental.tsv")
    log("=== DONE ===")


if __name__ == "__main__":
    main()
