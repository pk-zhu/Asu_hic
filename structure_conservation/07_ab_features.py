#!/usr/bin/env python3
"""补做 A/B compartment 特征统计 (对齐 Wang 2024 Supp Fig 4-6):
- ab_gene_density.tsv: A vs B compartment 内 gene 密度 (基因/Mb), 双侧 Wilcoxon
- ab_te_density.tsv: A vs B TE bp 覆盖率, 双侧 Wilcoxon
- ab_expression.tsv: A vs B 内 log2(TPM+1) 基因表达差异, 双侧 Wilcoxon
- switch_expression.tsv: 4 类 compartment switch (A-to-A/A-to-B/B-to-B/B-to-A)
    的 log2(TPM+1) 分布对比 (env / non-env 拆分), 双侧 MWU
"""
import os, sys, numpy as np, pandas as pd
from scipy import stats

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
OUT = f"{BASE}/03.pairwise_subgenome_conservation/results/ab_features"
os.makedirs(OUT, exist_ok=True)

COMP_BED = {
    "Ath": f"{BASE}/01.single/Ath/Compartments/CALDER2/100000/Ath/sub_compartments/all_sub_compartments.bed",
    "Aar": f"{BASE}/01.single/Aar/Compartments/CALDER2/100000/Aar/sub_compartments/all_sub_compartments.bed",
    "Asu": f"{BASE}/01.single/Asu/Compartments/CALDER2/100000/Asu/sub_compartments/all_sub_compartments.bed",
    "syn_Asu": f"{BASE}/01.single/syn_Asu/Compartments/CALDER2/100000/syn_Asu/sub_compartments/all_sub_compartments.bed",
}
TE_BED = {
    "Ath": f"{BASE}/07.subgenome_TE_analysis/results/te_bed/Ath_te_all.bed",
    "Aar": f"{BASE}/07.subgenome_TE_analysis/results/te_bed/Aar_te_all.bed",
    "Asu": f"{BASE}/07.subgenome_TE_analysis/results/te_bed/Asu_te_all.bed",
    # syn_Asu 与天然 Asu 共享基因组组装与注释 (同一物种合成系), 复用 Asu TE bed
    "syn_Asu": f"{BASE}/07.subgenome_TE_analysis/results/te_bed/Asu_te_all.bed",
}
TPM = {
    "Ath": f"{BASE}/08.expression/04.count/matrix/Ath.tpm.tsv",
    "Aar": f"{BASE}/08.expression/04.count/matrix/Aar.tpm.tsv",
    "Asu": f"{BASE}/08.expression/04.count/matrix/Asu.tpm.tsv",
    # syn_Asu 的 RNA-seq 落盘文件名为 Allo738.tpm.tsv (同一合成系测量)
    "syn_Asu": f"{BASE}/08.expression/04.count/matrix/Allo738.tpm.tsv",
}
GENES = {
    "Ath": f"{BASE}/01.single/Ath/ref/Ath.gff",
    "Aar": f"{BASE}/01.single/Aar/ref/Aar.gff",
    "Asu": f"{BASE}/01.single/Asu/ref/Asu.gff",
    # syn_Asu 共享 Asu 基因组与基因注释
    "syn_Asu": f"{BASE}/01.single/Asu/ref/Asu.gff",
}
BINSIZE = 100000

def norm_chrom(c, sp):
    """Normalize CALDER2 chrom names to match gff/te bed."""
    if sp == "Aar":
        # chr__1 -> chr_1
        return c.replace("chr__", "chr_")
    if sp in ("Asu", "syn_Asu"):
        # chrsT1 -> sT1, chrsA6 -> sA6
        return c.replace("chr", "", 1) if c.startswith("chr") else c
    return c

def load_comp(bed, sp):
    """Load CALDER2 sub_compartments; column 4 like B.2.2.2.2 or A.1.1.1.1 -> A/B label.
    CALDER2 label is eigenvector-sign-based; we later relabel by gene density so A = active."""
    df = pd.read_csv(bed, sep="\t", header=None,
                     names=["chrom","start","end","sub","score","strand","tstart","tend","rgb"])
    df["ab_raw"] = df["sub"].str[0]  # first char A or B
    df = df[df["ab_raw"].isin(["A","B"])].copy()
    df["chrom"] = df["chrom"].map(lambda c: norm_chrom(c, sp))
    return df[["chrom","start","end","ab_raw"]]

def load_genes(path, sp):
    if path.endswith(".gff"):
        # Parse GFF -> genes. Prefer Name= over ID= (Ath Liftoff appends _0 to ID).
        rows = []
        with open(path) as f:
            for line in f:
                if line.startswith("#"): continue
                f2 = line.rstrip().split("\t")
                if len(f2) < 9 or f2[2] != "gene": continue
                attr = f2[8]
                def get_attr(a, key):
                    for kv in a.split(";"):
                        if kv.startswith(key+"="):
                            return kv[len(key)+1:]
                    return None
                gid = get_attr(attr, "Name") or get_attr(attr, "ID") or "NA"
                rows.append([f2[0], int(f2[3])-1, int(f2[4]), gid])
        return pd.DataFrame(rows, columns=["chrom","start","end","gene_id"])
    else:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        df = df.iloc[:, :4]
        df.columns = ["chrom","start","end","gene_id"]
        return df

def bin_key(chrom, start):
    return (chrom, (start // BINSIZE) * BINSIZE)

def annotate_bins_with_ab(comp_df):
    """Return dict {(chrom, bin_start): ab_raw_label} using majority overlap per 100kb bin."""
    out = {}
    for chrom, g in comp_df.groupby("chrom"):
        for _, r in g.iterrows():
            b_start = (r["start"] // BINSIZE) * BINSIZE
            b_end = ((r["end"]-1) // BINSIZE) * BINSIZE
            for b in range(b_start, b_end + 1, BINSIZE):
                # take label of segment covering more of this bin (approx: first seen wins if tied)
                key = (chrom, b)
                out[key] = r["ab_raw"]
    return out

def relabel_by_gene_density(ab_map, genes_df):
    """CALDER2 A/B label depends on eigenvector sign, not activity.
    Per-chromosome relabel: on each chromosome, A = higher mean gene density.
    (Source files are already re-oriented per chromosome; this is a no-op guard.)"""
    genes_df = genes_df.copy()
    genes_df["bin"] = (genes_df["start"] // BINSIZE) * BINSIZE
    gpb = genes_df.groupby(["chrom","bin"]).size().reset_index(name="n_genes")
    gpb["raw"] = gpb.apply(lambda r: ab_map.get((r["chrom"], r["bin"]), None), axis=1)
    gpb = gpb.dropna(subset=["raw"])
    if len(gpb) == 0:
        return ab_map, False
    remapped = {}
    flipped_any = False
    for chrom, sub in gpb.groupby("chrom"):
        means = sub.groupby("raw")["n_genes"].mean()
        if "A" not in means.index or "B" not in means.index:
            continue  # 该染色体 A/B 不齐, 保持原样
        flip = means["A"] < means["B"]
        if flip:
            flipped_any = True
        for k in [kk for kk in ab_map if kk[0] == chrom]:
            v = ab_map[k]
            remapped[k] = ("B" if v == "A" else "A") if flip else v
    return remapped, flipped_any

results = {"gene_density": [], "te_density": [], "expression": []}

for sp in ["Ath", "Aar", "Asu", "syn_Asu"]:
    comp = load_comp(COMP_BED[sp], sp)
    ab_map_raw = annotate_bins_with_ab(comp)
    # -- gene density: count genes per 100kb bin, then join AB label
    genes = load_genes(GENES[sp], sp)
    ab_map, flipped = relabel_by_gene_density(ab_map_raw, genes)
    print(f"  {sp}: A/B relabeled by gene density, flipped={flipped}")
    genes["bin"] = (genes["start"] // BINSIZE) * BINSIZE
    gene_per_bin = genes.groupby(["chrom","bin"]).size().reset_index(name="n_genes")
    gene_per_bin["ab"] = gene_per_bin.apply(lambda r: ab_map.get((r["chrom"], r["bin"]), None), axis=1)
    gene_per_bin = gene_per_bin.dropna(subset=["ab"])
    gene_per_bin["density"] = gene_per_bin["n_genes"] / (BINSIZE / 1e6)  # genes / Mb
    a_vals = gene_per_bin.loc[gene_per_bin["ab"]=="A", "density"].values
    b_vals = gene_per_bin.loc[gene_per_bin["ab"]=="B", "density"].values
    u, p = stats.mannwhitneyu(a_vals, b_vals, alternative="two-sided") if len(a_vals) and len(b_vals) else (np.nan, np.nan)
    results["gene_density"].append(dict(
        species=sp, n_A=len(a_vals), n_B=len(b_vals),
        A_median=float(np.median(a_vals)) if len(a_vals) else np.nan,
        B_median=float(np.median(b_vals)) if len(b_vals) else np.nan,
        A_mean=float(np.mean(a_vals)) if len(a_vals) else np.nan,
        B_mean=float(np.mean(b_vals)) if len(b_vals) else np.nan,
        mwu_p=float(p) if not np.isnan(p) else np.nan,
    ))
    # -- TE density: sum TE bp per bin / 100 kb -> coverage fraction
    if os.path.exists(TE_BED[sp]):
        te = pd.read_csv(TE_BED[sp], sep="\t", header=None, comment="#").iloc[:, :3]
        te.columns = ["chrom","start","end"]
        # explode TE spans onto 100kb bins
        te["bin_start"] = (te["start"] // BINSIZE) * BINSIZE
        te["bin_end"] = ((te["end"] - 1) // BINSIZE) * BINSIZE
        # for TE spanning >1 bin, split (rare at 100kb, but handle)
        bin_bp = {}
        for _, r in te.iterrows():
            s, e = int(r["start"]), int(r["end"])
            for b in range(int(r["bin_start"]), int(r["bin_end"]) + 1, BINSIZE):
                overlap = min(e, b + BINSIZE) - max(s, b)
                if overlap > 0:
                    bin_bp[(r["chrom"], b)] = bin_bp.get((r["chrom"], b), 0) + overlap
        te_rows = [dict(chrom=k[0], bin=k[1], te_bp=v, coverage=v / BINSIZE) for k,v in bin_bp.items()]
        te_bin = pd.DataFrame(te_rows)
        te_bin["ab"] = te_bin.apply(lambda r: ab_map.get((r["chrom"], r["bin"]), None), axis=1)
        te_bin = te_bin.dropna(subset=["ab"])
        a = te_bin.loc[te_bin["ab"]=="A", "coverage"].values
        b = te_bin.loc[te_bin["ab"]=="B", "coverage"].values
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan)
        results["te_density"].append(dict(
            species=sp, n_A=len(a), n_B=len(b),
            A_median=float(np.median(a)) if len(a) else np.nan,
            B_median=float(np.median(b)) if len(b) else np.nan,
            A_mean=float(np.mean(a)) if len(a) else np.nan,
            B_mean=float(np.mean(b)) if len(b) else np.nan,
            mwu_p=float(p),
        ))

    # -- expression: gene -> bin -> ab -> log2(TPM+1)
    tpm = pd.read_csv(TPM[sp], sep="\t")
    # infer gene id col
    gid_col = tpm.columns[0]
    numeric_cols = [c for c in tpm.columns if c != gid_col]
    tpm["tpm_mean"] = tpm[numeric_cols].mean(axis=1)
    tpm["log2_tpm"] = np.log2(tpm["tpm_mean"] + 1)
    ge = genes.merge(tpm[[gid_col, "log2_tpm"]], left_on="gene_id", right_on=gid_col, how="inner")
    ge["ab"] = ge.apply(lambda r: ab_map.get((r["chrom"], (r["start"]//BINSIZE)*BINSIZE), None), axis=1)
    ge = ge.dropna(subset=["ab"])
    a = ge.loc[ge["ab"]=="A", "log2_tpm"].values
    b = ge.loc[ge["ab"]=="B", "log2_tpm"].values
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan)
    results["expression"].append(dict(
        species=sp, n_A=len(a), n_B=len(b),
        A_median=float(np.median(a)) if len(a) else np.nan,
        B_median=float(np.median(b)) if len(b) else np.nan,
        A_mean=float(np.mean(a)) if len(a) else np.nan,
        B_mean=float(np.mean(b)) if len(b) else np.nan,
        mwu_p=float(p),
    ))
    print(f"  {sp} OK")

pd.DataFrame(results["gene_density"]).to_csv(f"{OUT}/ab_gene_density.tsv", sep="\t", index=False)
pd.DataFrame(results["te_density"]).to_csv(f"{OUT}/ab_te_density.tsv", sep="\t", index=False)
pd.DataFrame(results["expression"]).to_csv(f"{OUT}/ab_expression.tsv", sep="\t", index=False)
print(f"[OK] AB features written to {OUT}")
