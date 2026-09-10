#!/usr/bin/env python3
"""07d_expr_vs_3d_trajectory.py — 4×4 交叉表: 三维轨迹 (PR/IR/TP/GR) × 表达轨迹 (Maintained/Remodeled/Synthetic-shift/Natural-shift).

三维轨迹来自 03.pairwise_subgenome_conservation/results/trajectory/trajectory_{sT_Ath,sA_Aar}.tsv
(per-gene matrix_corr 第 33 百分位阈值, 见 Methods 4.6)。
表达轨迹来自 07b 输出 expression_trajectory_parental.tsv (Δ Z-score, |z|≤1, 见 Methods 4.9)。

输出 long-format 表 + χ²/Cramér's V。
"""
import os
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
BASE = ROOT
EXPR_TSV = f"{BASE}/09.expression_3d_integration/results/ortholog_expr_allo738/expression_trajectory_parental.tsv"
TRAJ_TSV = {
    "sT": f"{BASE}/03.pairwise_subgenome_conservation/results/trajectory/trajectory_sT_Ath.tsv",
    "sA": f"{BASE}/03.pairwise_subgenome_conservation/results/trajectory/trajectory_sA_Aar.tsv",
}
OUT_DIR = f"{BASE}/09.expression_3d_integration/results/expr_vs_3d_trajectory"
TABLE_OUT = f"{BASE}/11.writing/tables/Table_S15_expr_vs_3d_trajectory.tsv"
os.makedirs(OUT_DIR, exist_ok=True)

TRAJ3D_MAP = {"Parental_Retention": "PR",
              "Immediate_Remodeling": "IR",
              "Transient_Perturbation": "TP",
              "Gradual_Remodeling": "GR"}
TRAJ3D_ORDER = ["PR", "IR", "TP", "GR"]
EXPR_ORDER = ["Maintained", "Remodeled", "Synthetic-shift", "Natural-shift"]


def main():
    expr = pd.read_csv(EXPR_TSV, sep="\t")
    rows = []
    stats = []
    for side in ("sT", "sA"):
        e = expr[expr.side == side][["gene", "expr_class"]].rename(columns={"gene": "gene_id"})
        t = pd.read_csv(TRAJ_TSV[side], sep="\t")
        t["traj_3d"] = t["trajectory"].map(TRAJ3D_MAP)
        m = e.merge(t[["gene_id", "traj_3d"]], on="gene_id", how="inner")
        n_side = len(m)
        ct = pd.crosstab(m["traj_3d"], m["expr_class"]).reindex(
            index=TRAJ3D_ORDER, columns=EXPR_ORDER, fill_value=0)
        n_traj = ct.sum(axis=1)
        for t3 in TRAJ3D_ORDER:
            for ec in EXPR_ORDER:
                n_cell = int(ct.loc[t3, ec])
                rows.append({
                    "side": side,
                    "traj_3d": t3,
                    "expr_class": ec,
                    "n": n_cell,
                    "n_traj_3d": int(n_traj[t3]),
                    "pct_within_traj_3d": round(n_cell / n_traj[t3] * 100, 2),
                })
        chi2, p, dof, _ = chi2_contingency(ct.values)
        r, c = ct.shape
        V = float(np.sqrt(chi2 / (n_side * min(r - 1, c - 1))))
        stats.append({"side": side, "n": n_side, "chi2": round(chi2, 2),
                      "dof": int(dof), "p": float(p), "cramers_v": round(V, 4)})
        print(f"[{side}] n={n_side:,}  chi2={chi2:.2f}  p={p:.3e}  V={V:.4f}")
        print(ct.to_string())
        print((ct.div(n_traj, axis=0) * 100).round(1).to_string())
        print()

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT_DIR}/expr_vs_3d_trajectory.tsv", sep="\t", index=False)
    out.to_csv(TABLE_OUT, sep="\t", index=False)
    pd.DataFrame(stats).to_csv(f"{OUT_DIR}/expr_vs_3d_trajectory_stats.tsv",
                               sep="\t", index=False)
    print(f"saved: {TABLE_OUT}")
    print(f"saved: {OUT_DIR}/expr_vs_3d_trajectory_stats.tsv")


if __name__ == "__main__":
    main()
