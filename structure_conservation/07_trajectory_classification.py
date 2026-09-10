#!/usr/bin/env python3
"""
07_trajectory_classification.py — 3D基因组演化轨迹分类.

利用 syn_Asu (合成F3) + Asu (天然) + Parent (Ath/Aar) 的 per-gene matrix_corr,
将每个基因分为四类轨迹:

  Parental Retention:    syn≈Parent 且 Asu≈Parent  (亲本继承)
  Immediate Remodeling:  syn≠Parent 且 Asu≠Parent  (早期重塑, 长期保留)
  Transient Perturbation: syn≠Parent 且 Asu≈Parent (早期扰动后恢复)
  Gradual Remodeling:    syn≈Parent 且 Asu≠Parent (长期演化逐渐形成)

阈值: matrix_corr 的 33rd percentile (下三分位数), 低于此值为"diverged".

输出: results/trajectory/
  trajectory_{sT_Ath,sA_Aar}.tsv       per-gene 分类
  trajectory_summary.tsv                汇总统计
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
GENE_3D_NAT = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/gene_3d/10000")
GENE_3D_SYN = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/gene_3d_syn_Asu/10000")
OUT_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results/trajectory")
os.makedirs(OUT_DIR, exist_ok=True)

COMPARISONS = [
    ("sT_Ath", "Asu_Ath", "syn_Asu_vs_Ath", "sT"),
    ("sA_Aar", "Asu_Aar", "syn_Asu_vs_Aar", "sA"),
]


def log(msg):
    print(msg, flush=True)


def classify_trajectory(syn_corr, nat_corr, threshold):
    """四类轨迹分类."""
    syn_cons = syn_corr >= threshold
    nat_cons = nat_corr >= threshold
    if syn_cons and nat_cons:
        return "Parental_Retention"
    elif not syn_cons and not nat_cons:
        return "Immediate_Remodeling"
    elif not syn_cons and nat_cons:
        return "Transient_Perturbation"
    else:
        return "Gradual_Remodeling"


def main():
    log("=== 07_trajectory_classification ===")

    all_summaries = []

    for label, nat_comp, syn_comp, subgenome in COMPARISONS:
        log(f"\n{'='*60}")
        log(f"  {label} ({subgenome})")
        log(f"{'='*60}")

        # Load natural Asu per-gene data
        nat_file = os.path.join(GENE_3D_NAT, f"gene_3d_conservation_{nat_comp}.tsv")
        nat = pd.read_csv(nat_file, sep="\t")
        nat = nat[["gene_id", "matrix_corr"]].copy()
        nat = nat.dropna(subset=["matrix_corr"])
        nat = nat.rename(columns={"matrix_corr": "nat_corr"})
        log(f"  Natural Asu genes: {len(nat):,}")

        # Load syn_Asu per-gene data
        syn_file = os.path.join(GENE_3D_SYN, f"gene_3d_conservation_{syn_comp}.tsv")
        syn = pd.read_csv(syn_file, sep="\t")
        syn = syn[["gene_id", "matrix_corr"]].copy()
        syn = syn.dropna(subset=["matrix_corr"])
        syn = syn.rename(columns={"matrix_corr": "syn_corr"})
        log(f"  Synthetic syn_Asu genes: {len(syn):,}")

        # Merge
        merged = nat.merge(syn, on="gene_id", how="inner")
        log(f"  Merged genes: {len(merged):,}")

        # Compute threshold (33rd percentile of all corr values)
        all_corr = np.concatenate([merged["nat_corr"].values, merged["syn_corr"].values])
        threshold = np.percentile(all_corr, 33)
        log(f"  Threshold (33rd percentile): {threshold:.4f}")

        # Classify
        merged["trajectory"] = merged.apply(
            lambda r: classify_trajectory(r["syn_corr"], r["nat_corr"], threshold), axis=1)

        # Summary
        counts = merged["trajectory"].value_counts()
        total = len(merged)
        log(f"\n  Trajectory distribution:")
        for traj in ["Parental_Retention", "Immediate_Remodeling",
                      "Transient_Perturbation", "Gradual_Remodeling"]:
            n = counts.get(traj, 0)
            mean_nat = merged[merged["trajectory"] == traj]["nat_corr"].mean()
            mean_syn = merged[merged["trajectory"] == traj]["syn_corr"].mean()
            log(f"    {traj}: {n:6,} ({n/total*100:5.1f}%)  "
                f"nat_corr={mean_nat:.3f}  syn_corr={mean_syn:.3f}")

        # Save per-gene
        out_cols = ["gene_id", "nat_corr", "syn_corr", "trajectory"]
        merged[out_cols].to_csv(
            os.path.join(OUT_DIR, f"trajectory_{label}.tsv"), sep="\t", index=False)

        # Summary row (threshold 一并落盘: Figure 4A 要画这条 33 分位参考线,
        # 若绘图端自己重算, 一旦上下游的 dropna 口径或分位算法不同就会与散点
        # 分色错位。改由本脚本输出唯一权威值。)
        for traj in ["Parental_Retention", "Immediate_Remodeling",
                      "Transient_Perturbation", "Gradual_Remodeling"]:
            sub = merged[merged["trajectory"] == traj]
            all_summaries.append(dict(
                comparison=label, subgenome=subgenome, trajectory=traj,
                n=len(sub), fraction=round(len(sub)/total, 4),
                mean_nat_corr=round(sub["nat_corr"].mean(), 4),
                mean_syn_corr=round(sub["syn_corr"].mean(), 4),
                threshold_p33=round(float(threshold), 6),
            ))

        # Compare syn vs nat overall
        log(f"\n  Overall: nat_corr mean={merged['nat_corr'].mean():.4f}, "
            f"syn_corr mean={merged['syn_corr'].mean():.4f}")
        log(f"  Δ(syn - nat) = {merged['syn_corr'].mean() - merged['nat_corr'].mean():.4f}")

    # Save summary
    summary = pd.DataFrame(all_summaries)
    summary.to_csv(os.path.join(OUT_DIR, "trajectory_summary.tsv"), sep="\t", index=False)

    # sT vs sA comparison
    log(f"\n{'='*60}")
    log(f"  sT vs sA comparison")
    log(f"{'='*60}")
    for traj in ["Parental_Retention", "Immediate_Remodeling",
                  "Transient_Perturbation", "Gradual_Remodeling"]:
        st = summary[(summary["trajectory"] == traj) & (summary["subgenome"] == "sT")]
        sa = summary[(summary["trajectory"] == traj) & (summary["subgenome"] == "sA")]
        if len(st) > 0 and len(sa) > 0:
            log(f"  {traj}: sT={st.iloc[0]['fraction']:.3f}  sA={sa.iloc[0]['fraction']:.3f}")

    # Chi-square test for sT vs sA distribution
    st_data = summary[summary["subgenome"] == "sT"].set_index("trajectory")["n"]
    sa_data = summary[summary["subgenome"] == "sA"].set_index("trajectory")["n"]
    traj_order = ["Parental_Retention", "Immediate_Remodeling",
                   "Transient_Perturbation", "Gradual_Remodeling"]
    cont = pd.DataFrame({"sT": [st_data.get(t, 0) for t in traj_order],
                          "sA": [sa_data.get(t, 0) for t in traj_order]},
                         index=traj_order)
    if cont.values.sum() > 0:
        chi2, p, dof, _ = chi2_contingency(cont)
        log(f"\n  χ² sT vs sA: chi2={chi2:.2f}, dof={dof}, p={p:.4f}")

    log(f"\n 保存: {OUT_DIR}/")
    log("=== DONE ===")


if __name__ == "__main__":
    main()