#!/usr/bin/env python3
"""
05. 统一可视化

Part A (coverage):   results/coverage/coverage_plots.pdf
Part B1 (gene_3d):   results/gene_3d/gene_3d_plots.pdf
Part B2 (density):   results/density_link/density_link_plots.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = os.path.join(ROOT, "04.cds_synteny_3d_conservation", "results")
COV = os.path.join(BASE, "coverage")
G3D = os.path.join(BASE, "gene_3d")
DL = os.path.join(BASE, "density_link")

COL = {"sT": "#3B4992FF", "Ath": "#3B4992FF", "sA": "#EE0000FF", "Aar": "#EE0000FF",
       "Asu_Ath": "#3B4992FF", "Asu_Aar": "#EE0000FF"}


def plot_coverage():
    cov = pd.read_csv(os.path.join(COV, "cds_coverage_in_synteny.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(COV, "coverage_summary.tsv"), sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1) 逐染色体 dir1 覆盖率
    ax = axes[0]
    sides = ["sT", "Ath", "sA", "Aar"]
    xpos = 0
    xticks, xlabels = [], []
    for side in sides:
        sub = cov[cov["side"] == side].sort_values("chrom")
        xs = np.arange(xpos, xpos + len(sub))
        ax.bar(xs, sub["cds_coverage_fraction"], color=COL[side], width=0.8,
               label=side, edgecolor="black", linewidth=0.3)
        xticks += list(xs); xlabels += list(sub["chrom"])
        xpos += len(sub) + 1
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels, rotation=90, fontsize=7)
    ax.set_ylabel("CDS coverage in collinear bins")
    ax.set_title("Dir1: CDS coverage within synteny (per chrom)")
    ax.legend()

    # 2) dir1 整体
    ax = axes[1]
    ax.bar(summ["side"], summ["dir1_cds_coverage_in_synteny"],
           color=[COL[s] for s in summ["side"]], edgecolor="black")
    ax.set_ylabel("CDS coverage fraction")
    ax.set_title("Dir1: overall CDS coverage in synteny")
    for i, v in enumerate(summ["dir1_cds_coverage_in_synteny"]):
        ax.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=9)

    # 3) dir2 共线占比
    ax = axes[2]
    ax.bar(summ["side"], summ["dir2_cds_synteny_fraction"],
           color=[COL[s] for s in summ["side"]], edgecolor="black")
    ax.set_ylabel("CDS synteny fraction")
    ax.set_title("Dir2: fraction of CDS that is syntenic")
    for i, v in enumerate(summ["dir2_cds_synteny_fraction"]):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)

    plt.tight_layout()
    out = os.path.join(COV, "coverage_plots.pdf")
    plt.savefig(out, dpi=200); plt.close()
    print(f"  -> {out}")


def plot_gene3d():
    # 03_gene_region_3d_conservation.py 按 matrix 分辨率写到 G3D/{res}/, 默认看 10kb
    g3d_src = os.path.join(G3D, "10000")
    summ = pd.read_csv(os.path.join(g3d_src, "gene_3d_summary.tsv"), sep="\t")
    d_ath = pd.read_csv(os.path.join(g3d_src, "gene_3d_conservation_Asu_Ath.tsv"), sep="\t")
    d_aar = pd.read_csv(os.path.join(g3d_src, "gene_3d_conservation_Asu_Aar.tsv"), sep="\t")

    def corrs(df):
        v = pd.to_numeric(df["matrix_corr"], errors="coerce").dropna()
        return v.values

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1) matrix_corr 分布
    ax = axes[0]
    ax.hist(corrs(d_ath), bins=60, color=COL["sT"], alpha=0.6, density=True, label="sT vs Ath")
    ax.hist(corrs(d_aar), bins=60, color=COL["sA"], alpha=0.6, density=True, label="sA vs Aar")
    ax.set_xlabel("per-gene matrix correlation"); ax.set_ylabel("density")
    ax.set_title("B1: matrix conservation distribution"); ax.legend()

    # 2) mean/median matrix corr
    ax = axes[1]
    x = np.arange(len(summ)); w = 0.35
    ax.bar(x - w/2, summ["mean_matrix_corr"], w, label="mean", color="#888")
    ax.bar(x + w/2, summ["median_matrix_corr"], w, label="median", color="#ccc", edgecolor="black")
    ax.set_xticks(x); ax.set_xticklabels(summ["subgenome"])
    ax.set_ylabel("matrix correlation"); ax.set_title("B1: mean/median matrix corr"); ax.legend()

    # 3) TAD 保守率 (默认看 TAD@10kb)
    ax = axes[2]
    tad_col = "tad_conserved_rate_10000" if "tad_conserved_rate_10000" in summ.columns else "tad_conserved_rate"
    ax.bar(summ["subgenome"], summ[tad_col],
           color=[COL[s] for s in summ["subgenome"]], edgecolor="black")
    ax.set_ylabel("TAD conserved rate (10kb)"); ax.set_title("B1: TAD conservation rate")
    for i, v in enumerate(summ[tad_col]):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)

    plt.tight_layout()
    out = os.path.join(G3D, "gene_3d_plots.pdf")
    plt.savefig(out, dpi=200); plt.close()
    print(f"  -> {out}")


def plot_density_link():
    comp = pd.read_csv(os.path.join(DL, "cds_density_vs_compartment.tsv"), sep="\t")
    tad = pd.read_csv(os.path.join(DL, "cds_density_vs_tad.tsv"), sep="\t")
    loop = pd.read_csv(os.path.join(DL, "cds_density_vs_loop.tsv"), sep="\t")
    mat = pd.read_csv(os.path.join(DL, "cds_density_vs_matrix.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(DL, "density_link_summary.tsv"), sep="\t")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # 1) compartment boxplot: conserved vs switched
    ax = axes[0, 0]
    data, labels = [], []
    for c in ["Asu_Ath", "Asu_Aar"]:
        sub = comp[comp["comparison"] == c]
        grp = sub.assign(grp=np.where(sub["conservation_type"] == "conserved", "conserved", "switched"))
        for g in ["conserved", "switched"]:
            data.append(grp[grp["grp"] == g]["cds_coverage"].values)
            labels.append(f"{c.split('_')[1]}\n{g}")
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True)
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(COL["Asu_Ath"] if i < 2 else COL["Asu_Aar"]); patch.set_alpha(0.6)
    ax.set_ylabel("CDS coverage (100k bin)"); ax.set_title("B2 compartment: CDS density by conservation")

    # 2) TAD scatter
    ax = axes[0, 1]
    for c in ["Asu_Ath", "Asu_Aar"]:
        sub = tad[tad["comparison"] == c]
        ax.scatter(sub["cds_coverage"], sub["conservation_score"], s=10, alpha=0.5,
                   color=COL[c], label=c)
    ax.set_xlabel("TAD CDS coverage"); ax.set_ylabel("TAD conservation score")
    ax.set_title("B2 TAD: CDS density vs conservation score"); ax.legend()

    # 3) loop bar
    ax = axes[1, 0]
    rows = summ[summ["layer"] == "loop"]
    xs = np.arange(len(rows))
    ax.bar(xs, pd.to_numeric(rows["value"], errors="coerce"),
           color=["#999"]*len(rows), edgecolor="black")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{r.comparison}\n{r.metric.split('[')[-1].rstrip(']')}"
                        for r in rows.itertuples()], fontsize=8)
    ax.set_ylabel("mean anchor CDS coverage"); ax.set_title("B2 loop: anchor CDS density by class")

    # 4) matrix hexbin
    ax = axes[1, 1]
    hb = ax.hexbin(mat["region_cds_coverage"], mat["matrix_corr"], gridsize=50,
                   cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=ax, label="count")
    rho_rows = summ[(summ["layer"] == "matrix") & (summ["metric"] == "spearman_cov_vs_corr")]
    txt = "  ".join(f"{r.comparison}: ρ={float(r.value):.3f}" for r in rho_rows.itertuples())
    ax.set_xlabel("gene-region CDS coverage"); ax.set_ylabel("matrix correlation")
    ax.set_title(f"B2 matrix: CDS density vs conservation\n{txt}")

    plt.tight_layout()
    out = os.path.join(DL, "density_link_plots.pdf")
    plt.savefig(out, dpi=200); plt.close()
    print(f"  -> {out}")


def main():
    print("绘图...")
    plot_coverage()
    plot_gene3d()
    plot_density_link()
    print("05 done.")


if __name__ == "__main__":
    main()
