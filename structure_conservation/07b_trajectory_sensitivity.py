#!/usr/bin/env python3
"""Sensitivity analysis: trajectory classification across percentile thresholds.

Redo classification at 25th/33rd/40th/50th percentile of matrix_corr and
recompute (Parental_Retention, Immediate_Remodeling, Transient_Perturbation,
Gradual_Remodeling) proportions plus chi-square between sT and sA.
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
BASE = os.path.join(ROOT, "03.pairwise_subgenome_conservation", "results", "trajectory")
OUT = os.path.join(BASE, "trajectory_sensitivity.tsv")

st = pd.read_csv(f"{BASE}/trajectory_sT_Ath.tsv", sep="\t").dropna(subset=["nat_corr", "syn_corr"])
sa = pd.read_csv(f"{BASE}/trajectory_sA_Aar.tsv", sep="\t").dropna(subset=["nat_corr", "syn_corr"])


def classify(df, q):
    all_vals = pd.concat([df["nat_corr"], df["syn_corr"]])
    thr = all_vals.quantile(q)
    nat_ok = df["nat_corr"] >= thr
    syn_ok = df["syn_corr"] >= thr
    trj = np.where(nat_ok & syn_ok, "PR",
          np.where(~nat_ok & ~syn_ok, "IR",
          np.where(~syn_ok & nat_ok, "TP", "GR")))
    return trj, thr


rows = []
for q in [0.25, 0.33, 0.40, 0.50]:
    st_trj, st_thr = classify(st, q)
    sa_trj, sa_thr = classify(sa, q)
    st_counts = pd.Series(st_trj).value_counts().to_dict()
    sa_counts = pd.Series(sa_trj).value_counts().to_dict()
    for cat in ["PR", "IR", "TP", "GR"]:
        st_counts.setdefault(cat, 0)
        sa_counts.setdefault(cat, 0)

    ct = np.array([[st_counts[c] for c in ["PR", "IR", "TP", "GR"]],
                   [sa_counts[c] for c in ["PR", "IR", "TP", "GR"]]])
    chi2, p, _, _ = chi2_contingency(ct)

    n_st = int(sum(st_counts.values()))
    n_sa = int(sum(sa_counts.values()))
    for cat in ["PR", "IR", "TP", "GR"]:
        rows.append({
            "percentile": f"{int(q*100)}th",
            "st_threshold": round(float(st_thr), 4),
            "sa_threshold": round(float(sa_thr), 4),
            "subgenome": "sT",
            "category": cat,
            "n": int(st_counts[cat]),
            "fraction": round(st_counts[cat] / n_st, 4),
        })
        rows.append({
            "percentile": f"{int(q*100)}th",
            "st_threshold": round(float(st_thr), 4),
            "sa_threshold": round(float(sa_thr), 4),
            "subgenome": "sA",
            "category": cat,
            "n": int(sa_counts[cat]),
            "fraction": round(sa_counts[cat] / n_sa, 4),
        })
    rows.append({
        "percentile": f"{int(q*100)}th",
        "st_threshold": round(float(st_thr), 4),
        "sa_threshold": round(float(sa_thr), 4),
        "subgenome": "chi2_sT_vs_sA",
        "category": "chi2",
        "n": round(float(chi2), 2),
        "fraction": float(p),
    })
    print(f"{int(q*100)}th percentile (sT thr={st_thr:.3f} sA thr={sa_thr:.3f}):")
    print(f"  sT: PR={st_counts['PR']/n_st:.3f} IR={st_counts['IR']/n_st:.3f} "
          f"TP={st_counts['TP']/n_st:.3f} GR={st_counts['GR']/n_st:.3f}")
    print(f"  sA: PR={sa_counts['PR']/n_sa:.3f} IR={sa_counts['IR']/n_sa:.3f} "
          f"TP={sa_counts['TP']/n_sa:.3f} GR={sa_counts['GR']/n_sa:.3f}")
    print(f"  chi2={chi2:.1f}  p={p:.2e}")

pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
print(f"saved: {OUT}")
