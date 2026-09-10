#!/usr/bin/env python3
"""09_gene_region_expression_joint.py — 基因级 3D region (comp × TAD) × 表达 联合分类.

对每个 Asu gene (取自 09 的 expression_triplets, 12,578 个 triplet × 2 sides):
  1. 用 Asu_genes.bed 找基因所在染色体 + 起止.
  2. 把基因投到对应 (comparison) 的 100kb compartment bin:
       sT-side → nat_sT_Ath 或 syn_sT_Ath
       sA-side → nat_sA_Aar 或 syn_sA_Aar
     命中 bin → 取该 bin 的 comp_conservation_type (per-bin dedup).
     comp 4 类: conserved / A_to_B / B_to_A / no_bin.
  3. 在 (comparison) 的 TAD 表里找 overlap 最长的 TAD → conservation_score.
     TAD 4 类: perfect (>=1.0) / high (>=0.7) / partial (>=0.4) / low (<0.4).
  4. 与 expression 的 4 类 (sT_expr_class / sA_expr_class) 联合计数.

输出:
  09/results/expr_3d_joint/compartment_expression_joint_summary.tsv (4 expr × 4 comp × 4 case)
  09/results/expr_3d_joint/tad_expression_joint_summary.tsv (4 expr × 4 tad × 4 case)
  09/results/expr_3d_joint/gene_region_class_{label}.tsv (per-gene 详情, 4 文件)
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT_BASE, PWCONS_DIR, load_genes_bed  # noqa: E402

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
COMP_DIR = os.path.join(PWCONS_DIR, "compartment")
TAD_DIR  = os.path.join(PWCONS_DIR, "tad/10000")
GENE_BED_TSV = os.path.join(BASE, "04.cds_synteny_3d_conservation/results/cds/Asu_genes.bed")
EXPR_TSV = os.path.join(OUT_BASE, "ortholog_expr_allo738/expression_trajectory_parental.tsv")

OUT_DIR = os.path.join(OUT_BASE, "expr_3d_joint")
os.makedirs(OUT_DIR, exist_ok=True)

CASES = [
    # (label, subgenome_side, comp_file, tad_file)
    ("nat_sT_Ath", "sT", "compartment_conservation_Asu_Ath.tsv",
     "tad_conservation_Asu_Ath.tsv"),
    ("nat_sA_Aar", "sA", "compartment_conservation_Asu_Aar.tsv",
     "tad_conservation_Asu_Aar.tsv"),
    ("syn_sT_Ath", "sT", "compartment_conservation_syn_Asu_vs_Ath.tsv",
     "tad_conservation_syn_Asu_vs_Ath.tsv"),
    ("syn_sA_Aar", "sA", "compartment_conservation_syn_Asu_vs_Aar.tsv",
     "tad_conservation_syn_Asu_vs_Aar.tsv"),
]

EXPR_ORDER = ["Maintained", "Remodeled", "Synthetic-shift", "Natural-shift"]

COMP_4_ORDER = ["conserved", "A_to_B", "B_to_A", "no_bin"]
TAD_4_ORDER = ["perfect", "high", "partial", "low"]
TAD_4_THRESH = [1.0, 0.7, 0.4]


def tad_4class(score):
    if pd.isna(score):
        return "no_tad"
    if score >= TAD_4_THRESH[0]:
        return "perfect"
    if score >= TAD_4_THRESH[1]:
        return "high"
    if score >= TAD_4_THRESH[2]:
        return "partial"
    return "low"


def load_comp_per_bin(path):
    """per-bin dedup: 同一 bin 多个 ref 行, 优先级 conserved > A_to_B > B_to_A."""
    df = pd.read_csv(path, sep="\t")

    def merge_type(x):
        x = set(x)
        if "conserved" in x:
            return "conserved"
        if "A_to_B" in x:
            return "A_to_B"
        if "B_to_A" in x:
            return "B_to_A"
        return "mixed"

    per_bin = df.groupby(["asu_chr", "asu_start", "asu_end"]).agg(
        comp_type=("conservation_type", merge_type)
    ).reset_index()
    return per_bin


def load_tad(path):
    df = pd.read_csv(path, sep="\t")
    if "asu_tad_start" in df.columns:
        df = df.rename(columns={"asu_tad_start": "asu_start",
                                "asu_tad_end": "asu_end"})
    return df[["asu_chr", "asu_start", "asu_end", "conservation_score"]]


def assign_gene_to_bin(genes, comp_per_bin, tad):
    """对每个基因找其所在 comp bin + overlap 最长的 TAD."""
    comp_dict = {}
    for idx, row in comp_per_bin.iterrows():
        comp_dict.setdefault(row["asu_chr"], []).append((row["asu_start"],
                                                          row["asu_end"],
                                                          idx,
                                                          row["comp_type"]))
    tad_by_chr = {c: g for c, g in tad.groupby("asu_chr")}

    out_rows = []
    for gene, (chr_, start, end) in genes.items():
        comp_type = "no_bin"
        for (bs, be, idx, ct) in comp_dict.get(chr_, []):
            ov_s = max(bs, start)
            ov_e = min(be, end)
            if ov_e > ov_s:
                comp_type = ct
                break
        tad_score = np.nan
        ts = tad_by_chr.get(chr_)
        if ts is not None and len(ts):
            ts_arr = ts["asu_start"].values
            te_arr = ts["asu_end"].values
            ov_s = np.maximum(ts_arr, start)
            ov_e = np.minimum(te_arr, end)
            ov_len = np.maximum(ov_e - ov_s, 0)
            if ov_len.max() > 0:
                best = int(np.argmax(ov_len))
                tad_score = float(ts["conservation_score"].values[best])
        out_rows.append({
            "gene": gene,
            "asu_chr": chr_,
            "start": start,
            "end": end,
            "comp_type": comp_type,
            "tad_score": tad_score,
        })
    return pd.DataFrame(out_rows)


def main():
    expr = pd.read_csv(EXPR_TSV, sep="\t")
    triplets = expr[["sT_gene", "sA_gene", "sT_expr_class", "sA_expr_class"]].drop_duplicates(
        subset=["sT_gene", "sA_gene"])
    sT_genes = triplets.set_index("sT_gene")["sT_expr_class"].to_dict()
    sA_genes = triplets.set_index("sA_gene")["sA_expr_class"].to_dict()
    print(f"Triplets: {len(triplets):,}  sT genes: {len(sT_genes):,}  sA genes: {len(sA_genes):,}")

    gene_bed = load_genes_bed("Asu").drop_duplicates(subset=["gene_id"], keep="first")
    gene_bed = gene_bed.rename(columns={"chr": "asu_chr"})
    gene_dict = gene_bed.set_index("gene_id")[["asu_chr", "start", "end"]]
    print(f"Loaded gene bed: {len(gene_dict):,} Asu genes")

    rows_comp = []
    rows_tad = []
    for label, side, comp_file, tad_file in CASES:
        comp_per_bin = load_comp_per_bin(os.path.join(COMP_DIR, comp_file))
        tad = load_tad(os.path.join(TAD_DIR, tad_file))
        gene_map = sT_genes if side == "sT" else sA_genes
        target_genes = {g: gene_dict.loc[g].values.tolist() for g in gene_map
                        if g in gene_dict.index}
        print(f"\n[{label}]  comp bins: {len(comp_per_bin):,}  "
              f"tads: {len(tad):,}  target genes with coords: {len(target_genes):,}")
        assigned = assign_gene_to_bin(target_genes, comp_per_bin, tad)
        assigned["tad_class"] = assigned["tad_score"].apply(tad_4class)
        assigned["expr_class"] = assigned["gene"].map(gene_map)
        assigned["comparison"] = label
        assigned["side"] = side

        n_total = len(assigned)
        print(f"  assigned genes: {n_total:,}")
        assigned[["gene", "asu_chr", "start", "end",
                  "comp_type", "tad_score", "tad_class",
                  "expr_class"]].to_csv(
            os.path.join(OUT_DIR, f"gene_region_class_{label}.tsv"),
            sep="\t", index=False)

        for ec in EXPR_ORDER:
            sub_ec = assigned[assigned["expr_class"] == ec]
            n_ec = len(sub_ec)
            for rc in COMP_4_ORDER:
                n_rc = int((sub_ec["comp_type"] == rc).sum())
                pct = n_rc / n_ec * 100 if n_ec else 0.0
                rows_comp.append({
                    "comparison": label, "side": side,
                    "expr_class": ec, "comp_class": rc,
                    "n": n_rc, "n_expr_class": n_ec,
                    "pct_of_expr": round(pct, 2),
                })
        for ec in EXPR_ORDER:
            sub_ec = assigned[assigned["expr_class"] == ec]
            n_ec = len(sub_ec)
            for rc in TAD_4_ORDER + ["no_tad"]:
                n_rc = int((sub_ec["tad_class"] == rc).sum())
                pct = n_rc / n_ec * 100 if n_ec else 0.0
                rows_tad.append({
                    "comparison": label, "side": side,
                    "expr_class": ec, "tad_class": rc,
                    "n": n_rc, "n_expr_class": n_ec,
                    "pct_of_expr": round(pct, 2),
                })

    pc = os.path.join(OUT_DIR, "compartment_expression_joint_summary.tsv")
    pd.DataFrame(rows_comp).to_csv(pc, sep="\t", index=False)
    pt = os.path.join(OUT_DIR, "tad_expression_joint_summary.tsv")
    pd.DataFrame(rows_tad).to_csv(pt, sep="\t", index=False)
    print(f"\nSaved: {pc}")
    print(f"Saved: {pt}")


if __name__ == "__main__":
    main()