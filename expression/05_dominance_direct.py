#!/usr/bin/env python3
"""
05_dominance_direct.py — sT-sA 直接 RBH 计算 subgenome dominance, 替代 Ath 桥接 triplet.
输出覆盖 09/results/ortholog_expr/subgenome_dominance.tsv (格式兼容下游).
"""
import os, sys, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT_BASE, load_tpm, tpm_summary, normalize_gene_id, asu_subgenome
from common import GENES_BED_DIR, TPM_MIN

OUT_DIR = os.path.join(OUT_BASE, "ortholog_expr")
RBH_DIR = os.path.join(OUT_BASE, "mmseqs_rbh")
os.makedirs(OUT_DIR, exist_ok=True)

RBH_COLS = ["query", "target", "fident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "tstart", "tend", "evalue", "bits"]


def log(msg):
    print(msg, flush=True)


def main():
    log("=== 05_dominance_direct (sT-sA direct RBH) ===")

    # Load expression
    asu_expr = tpm_summary(load_tpm("Asu"))[["log_mean", "mean_tpm", "expressed"]].reset_index()
    asu_expr = asu_expr[asu_expr["expressed"] == 1]
    log(f"  Expressed genes (mean TPM >= {TPM_MIN}): {len(asu_expr):,}")
    asu_idx = asu_expr.set_index("gene_id")["log_mean"]

    # Load sT-sA RBH
    rbh = pd.read_csv(os.path.join(RBH_DIR, "sT_sA_rbh"), sep="\t", header=None,
                      names=RBH_COLS, usecols=[0, 1, 2, 3, 10, 11])
    rbh["sT_gene"] = rbh["query"].apply(lambda x: normalize_gene_id(x, "Asu"))
    rbh["sA_gene"] = rbh["target"].apply(lambda x: normalize_gene_id(x, "Asu"))
    log(f"  RBH pairs: {len(rbh):,}  fident median={rbh['fident'].median():.3f}")

    # Verify subgenome assignment
    gene_chr = pd.read_csv(os.path.join(GENES_BED_DIR, "Asu_genes.bed"),
                           sep="\t", header=None,
                           names=["chr", "start", "end", "gene_id", "strand"])
    rbh = rbh.merge(gene_chr[["chr", "gene_id"]].rename(columns={"chr": "sT_chr", "gene_id": "sT_gene"}),
                    on="sT_gene", how="left")
    rbh = rbh.merge(gene_chr[["chr", "gene_id"]].rename(columns={"chr": "sA_chr", "gene_id": "sA_gene"}),
                    on="sA_gene", how="left")
    rbh["sT_sub"] = rbh["sT_chr"].map(asu_subgenome)
    rbh["sA_sub"] = rbh["sA_chr"].map(asu_subgenome)
    bad = rbh[(rbh["sT_sub"] != "sT") | (rbh["sA_sub"] != "sA")]
    if len(bad) > 0:
        log(f"  WARNING: {len(bad)} pairs with wrong subgenome assignment")
    rbh = rbh[(rbh["sT_sub"] == "sT") & (rbh["sA_sub"] == "sA")]
    log(f"  After subgenome filter: {len(rbh):,}")

    # Merge expression
    dom = rbh[["sT_gene", "sA_gene", "fident"]].copy()
    dom = dom.merge(asu_idx.rename("sT_log"), left_on="sT_gene", right_index=True, how="inner")
    dom = dom.merge(asu_idx.rename("sA_log"), left_on="sA_gene", right_index=True, how="inner")
    log(f"  With expression: {len(dom):,}")

    # Compute dominance
    dom["log2fc_sT_over_sA"] = (dom["sT_log"] - dom["sA_log"]) / np.log(2)
    log2 = np.log(2)
    cond_st = dom["sT_log"] > dom["sA_log"] + log2
    cond_sa = dom["sA_log"] > dom["sT_log"] + log2
    dom["dominant"] = np.where(cond_st, "sT", np.where(cond_sa, "sA", "balanced"))

    # Add ath_gene column for format compatibility (set to NA)
    dom["ath_gene"] = "NA"
    dom["ath_log"] = np.nan
    dom = dom[["sT_gene", "sA_gene", "ath_gene", "sT_log", "sA_log", "ath_log",
               "log2fc_sT_over_sA", "dominant"]]

    dom.to_csv(os.path.join(OUT_DIR, "subgenome_dominance.tsv"), sep="\t", index=False)

    v = dom.dominant.value_counts(normalize=True).to_dict()
    log(f"  Dominance: sT={v.get('sT',0):.3f}  sA={v.get('sA',0):.3f}  balanced={v.get('balanced',0):.3f}")
    log(f"  n={len(dom):,}")

    # Also compute dominance_sensitivity (1.5x/2x/3x thresholds)
    sens_rows = []
    for threshold in [1.5, 2.0, 3.0]:
        thresh_log = np.log2(threshold)
        st = dom["log2fc_sT_over_sA"] > thresh_log
        sa = dom["log2fc_sT_over_sA"] < -thresh_log
        bal = ~(st | sa)
        total = len(dom)
        sens_rows.append({
            "threshold": f"{threshold}x",
            "n_total": total,
            "n_sT_dom": int(st.sum()),
            "n_sA_dom": int(sa.sum()),
            "n_balanced": int(bal.sum()),
            "sT_dom_rate": float(st.mean()),
            "sA_dom_rate": float(sa.mean()),
            "balanced_rate": float(bal.mean()),
        })
    pd.DataFrame(sens_rows).to_csv(os.path.join(OUT_DIR, "dominance_sensitivity.tsv"),
                                    sep="\t", index=False)
    for r in sens_rows:
        log(f"  {r['threshold']}: sT={r['sT_dom_rate']:.3f} sA={r['sA_dom_rate']:.3f} balanced={r['balanced_rate']:.3f}")

    log(f"\n  保存: {OUT_DIR}/")
    log("=== DONE ===")


if __name__ == "__main__":
    main()