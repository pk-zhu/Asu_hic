#!/usr/bin/env python3
"""08_sample_spearman.py — Asu 品系样本间 Spearman 相关矩阵 (支撑 Figure 5A).

样本: 6 个天然 Asu 种子系 + 3 个合成 syn_Asu (Allo738) 生物学重复。
输入: 08.expression/04.count/matrix/{Asu,Allo738}.tpm.tsv 的原始 TPM, 取两者共有基因。
输出: 09.expression_3d_integration/results/ortholog_expr_allo738/sample_spearman.tsv

原先该矩阵在 11.writing/figures/plot_fig5.py 内即算即缓存, 并把 TSV 写进 figures/
目录, 违反 figures/common.py "绘图脚本只读 TSV" 的约定。统计计算移到本脚本,
绘图端改为纯读取。
"""
import os
import sys

import pandas as pd
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
EXPR_MATRIX = f"{BASE}/08.expression/04.count/matrix"
OUT_DIR = f"{BASE}/09.expression_3d_integration/results/ortholog_expr_allo738"
OUT_TSV = os.path.join(OUT_DIR, "sample_spearman.tsv")

ASU_NAT_COLS = ["SRR3676009", "SRR3676011", "SRR3676012",
                "SRR3676013", "SRR3676014", "SRR3676015"]
ASU_SYN_COLS = ["SRR12880912", "SRR12880913", "SRR12880914"]


def main():
    print("=== 08_sample_spearman ===", flush=True)
    asu_tpm = pd.read_csv(f"{EXPR_MATRIX}/Asu.tpm.tsv",
                          sep="\t", index_col=0)[ASU_NAT_COLS]
    allo_tpm = pd.read_csv(f"{EXPR_MATRIX}/Allo738.tpm.tsv",
                           sep="\t", index_col=0)[ASU_SYN_COLS]
    common = asu_tpm.index.intersection(allo_tpm.index)
    print(f"  nat samples: {len(ASU_NAT_COLS)}  syn samples: {len(ASU_SYN_COLS)}")
    print(f"  common genes: {len(common):,}")

    mat = pd.concat([asu_tpm.loc[common], allo_tpm.loc[common]], axis=1)
    rho, _ = spearmanr(mat.values, axis=0)
    corr = pd.DataFrame(rho, index=mat.columns, columns=mat.columns)

    os.makedirs(OUT_DIR, exist_ok=True)
    corr.to_csv(OUT_TSV, sep="\t", float_format="%.4f")
    print(f"  shape: {corr.shape}")

    nat = corr.loc[ASU_NAT_COLS, ASU_NAT_COLS].values
    syn = corr.loc[ASU_SYN_COLS, ASU_SYN_COLS].values
    off = corr.loc[ASU_NAT_COLS, ASU_SYN_COLS].values
    import numpy as np
    tri = lambda m: m[np.triu_indices_from(m, k=1)]
    print(f"  within nat: {tri(nat).min():.3f} - {tri(nat).max():.3f}")
    print(f"  within syn: {tri(syn).min():.3f} - {tri(syn).max():.3f}")
    print(f"  nat vs syn: {off.min():.3f} - {off.max():.3f}")
    print(f"  保存: {OUT_TSV}")
    print("=== DONE ===")


if __name__ == "__main__":
    main()
