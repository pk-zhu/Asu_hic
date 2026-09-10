#!/usr/bin/env python3
"""
06_dominance_allo738.py — Allo738 (第二个 synthetic A. suecica) subgenome dominance.
用 sT-sA 直接 RBH 判定 pair, 用 Allo738 TPM 计算 dominance.
输出与 05_dominance_direct.py 平行, 供跨样本对比 (Asu vs Allo738).
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT_BASE, EXPR_DIR, tpm_summary, normalize_gene_id, asu_subgenome
from common import GENES_BED_DIR, TPM_MIN

OUT_DIR = os.path.join(OUT_BASE, "ortholog_expr_allo738")
RBH_DIR = os.path.join(OUT_BASE, "mmseqs_rbh")
os.makedirs(OUT_DIR, exist_ok=True)

RBH_COLS = ["query", "target", "fident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "tstart", "tend", "evalue", "bits"]


def log(msg):
    print(msg, flush=True)


def load_tpm_allo738():
    """读 Allo738 TPM 矩阵 (与 Asu 相同 gene_id 空间)."""
    path = os.path.join(EXPR_DIR, "Allo738.tpm.tsv")
    df = pd.read_csv(path, sep="\t").rename(columns={"Gene_Id": "gene_id"})
    return df.set_index("gene_id")


def main():
    log("=== 06_dominance_allo738 ===")

    asu_expr = tpm_summary(load_tpm_allo738())[["log_mean", "mean_tpm", "expressed"]].reset_index()
    asu_expr = asu_expr[asu_expr["expressed"] == 1]
    log(f"  Expressed genes (mean TPM >= {TPM_MIN}): {len(asu_expr):,}")
    asu_idx = asu_expr.set_index("gene_id")["log_mean"]

    rbh = pd.read_csv(os.path.join(RBH_DIR, "sT_sA_rbh"), sep="\t", header=None,
                     names=RBH_COLS, usecols=[0, 1, 2, 3, 10, 11])
    rbh["sT_gene"] = rbh["query"].apply(lambda x: normalize_gene_id(x, "Asu"))
    rbh["sA_gene"] = rbh["target"].apply(lambda x: normalize_gene_id(x, "Asu"))
    log(f"  RBH pairs: {len(rbh):,}")

    gene_chr = pd.read_csv(os.path.join(GENES_BED_DIR, "Asu_genes.bed"),
                          sep="\t", header=None,
                          names=["chr", "start", "end", "gene_id", "strand"])
    rbh = rbh.merge(gene_chr[["chr", "gene_id"]].rename(columns={"chr": "sT_chr", "gene_id": "sT_gene"}),
                   on="sT_gene", how="left")
    rbh = rbh.merge(gene_chr[["chr", "gene_id"]].rename(columns={"chr": "sA_chr", "gene_id": "sA_gene"}),
                   on="sA_gene", how="left")
    rbh["sT_sub"] = rbh["sT_chr"].map(asu_subgenome)
    rbh["sA_sub"] = rbh["sA_chr"].map(asu_subgenome)
    rbh = rbh[(rbh["sT_sub"] == "sT") & (rbh["sA_sub"] == "sA")]
    log(f"  After subgenome filter: {len(rbh):,}")

    dom = rbh[["sT_gene", "sA_gene", "fident"]].copy()
    dom = dom.merge(asu_idx.rename("sT_log"), left_on="sT_gene", right_index=True, how="inner")
    dom = dom.merge(asu_idx.rename("sA_log"), left_on="sA_gene", right_index=True, how="inner")
    log(f"  With Allo738 expression: {len(dom):,}")

    dom["log2fc_sT_over_sA"] = (dom["sT_log"] - dom["sA_log"]) / np.log(2)
    log2 = np.log(2)
    cond_st = dom["sT_log"] > dom["sA_log"] + log2
    cond_sa = dom["sA_log"] > dom["sT_log"] + log2
    dom["dominant"] = np.where(cond_st, "sT", np.where(cond_sa, "sA", "balanced"))

    dom["ath_gene"] = "NA"
    dom["ath_log"] = np.nan
    dom = dom[["sT_gene", "sA_gene", "ath_gene", "sT_log", "sA_log", "ath_log",
              "log2fc_sT_over_sA", "dominant"]]

    dom.to_csv(os.path.join(OUT_DIR, "subgenome_dominance.tsv"), sep="\t", index=False)

    v = dom.dominant.value_counts(normalize=True).to_dict()
    log(f"  Allo738 Dominance: sT={v.get('sT',0):.3f}  sA={v.get('sA',0):.3f}  balanced={v.get('balanced',0):.3f}")
    log(f"  n={len(dom):,}")

    sens_rows = []
    for threshold in [1.5, 2.0, 3.0]:
        thresh_log = np.log2(threshold)
        st = dom["log2fc_sT_over_sA"] > thresh_log
        sa = dom["log2fc_sT_over_sA"] < -thresh_log
        bal = ~(st | sa)
        sens_rows.append({
            "threshold": f"{threshold}x",
            "n_total": len(dom),
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

    # Cross-sample comparison: Asu vs Allo738
    asu_dom_path = os.path.join(OUT_BASE, "ortholog_expr", "subgenome_dominance.tsv")
    if os.path.exists(asu_dom_path):
        asu_dom = pd.read_csv(asu_dom_path, sep="\t")
        log(f"\n  === Asu vs Allo738 cross-sample comparison ===")
        av = asu_dom.dominant.value_counts(normalize=True).to_dict()
        log(f"  Asu:     sT={av.get('sT',0):.3f}  sA={av.get('sA',0):.3f}  balanced={av.get('balanced',0):.3f}")
        log(f"  Allo738: sT={v.get('sT',0):.3f}  sA={v.get('sA',0):.3f}  balanced={v.get('balanced',0):.3f}")

        # Merge on (sT_gene, sA_gene) and compute correlation
        merged = asu_dom[["sT_gene", "sA_gene", "log2fc_sT_over_sA", "dominant"]].rename(
            columns={"log2fc_sT_over_sA": "log2fc_asu", "dominant": "dom_asu"}
        ).merge(
            dom[["sT_gene", "sA_gene", "log2fc_sT_over_sA", "dominant"]].rename(
                columns={"log2fc_sT_over_sA": "log2fc_allo", "dominant": "dom_allo"}
            ), on=["sT_gene", "sA_gene"], how="inner"
        )
        from scipy.stats import spearmanr, pearsonr
        rho, p = spearmanr(merged["log2fc_asu"], merged["log2fc_allo"])
        r, pp = pearsonr(merged["log2fc_asu"], merged["log2fc_allo"])
        log(f"  n_common={len(merged):,}")
        log(f"  log2fc Spearman rho={rho:.3f} (p={p:.2e})")
        log(f"  log2fc Pearson  r  ={r:.3f} (p={pp:.2e})")

        # Confusion matrix
        ct = pd.crosstab(merged["dom_asu"], merged["dom_allo"], margins=True)
        log(f"  Confusion (rows=Asu, cols=Allo738):")
        log(ct.to_string())
        agree = (merged["dom_asu"] == merged["dom_allo"]).mean()
        log(f"  Categorical agreement: {agree:.3f}")

        merged.to_csv(os.path.join(OUT_DIR, "dominance_asu_vs_allo738.tsv"), sep="\t", index=False)

        # 把上面只 log 出来的统计量落盘, 供 Figure 5C/5D 与正文 2.5 引用,
        # 避免绘图脚本自己重算 spearmanr (分析与绘图解耦)。
        pd.DataFrame([
            dict(metric="n_common", value=len(merged), p=np.nan),
            dict(metric="log2fc_spearman_rho", value=float(rho), p=float(p)),
            dict(metric="log2fc_pearson_r", value=float(r), p=float(pp)),
            dict(metric="categorical_agreement", value=float(agree), p=np.nan),
        ]).to_csv(os.path.join(OUT_DIR, "dominance_asu_vs_allo738_stats.tsv"),
                  sep="\t", index=False)

    log("=== DONE ===")


if __name__ == "__main__":
    main()
