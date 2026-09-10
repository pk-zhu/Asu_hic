#!/usr/bin/env python3
"""
统一可视化: matrix/SCC + compartment + TAD + loop 四层保守性

输出到 results/ 下的 PDF 文件
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
RESULTS_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results")
MATRIX_DIR = os.path.join(RESULTS_DIR, "matrix")
COMPARTMENT_DIR = os.path.join(RESULTS_DIR, "compartment")
TAD_DIR = os.path.join(RESULTS_DIR, "tad")
LOOP_DIR = os.path.join(RESULTS_DIR, "loop")

COLOR_ATH = '#3B4992FF'
COLOR_AAR = '#EE0000FF'
COLOR_A = '#00A087FF'
COLOR_B = '#7E6148FF'
COLOR_CONSERVED = '#00A087FF'
COLOR_A_TO_B = '#EFC000FF'
COLOR_B_TO_A = '#8DB600FF'


def plot_scc():
    """SCC 可视化: 3 分辨率 × 2 亚基因组 柱状图 + stratum r 曲线"""
    print("  绘制 SCC 图...")

    summary = pd.read_csv(os.path.join(MATRIX_DIR, "scc_summary.tsv"), sep='\t')

    # 图1: SCC 柱状图 (分辨率 × 亚基因组)
    fig, ax = plt.subplots(figsize=(10, 6))
    resolutions = sorted(summary['resolution'].unique())
    x = np.arange(len(resolutions))
    width = 0.35

    ath_vals = []
    aar_vals = []
    for res in resolutions:
        ath_row = summary[(summary['resolution'] == res) & (summary['comparison'] == 'Asu_Ath')]
        aar_row = summary[(summary['resolution'] == res) & (summary['comparison'] == 'Asu_Aar')]
        ath_vals.append(ath_row['scc'].values[0] if len(ath_row) > 0 else 0)
        aar_vals.append(aar_row['scc'].values[0] if len(aar_row) > 0 else 0)

    bars1 = ax.bar(x - width/2, ath_vals, width, label='Asu-Ath (sT)', color=COLOR_ATH, edgecolor='black')
    bars2 = ax.bar(x + width/2, aar_vals, width, label='Asu-Aar (sA)', color=COLOR_AAR, edgecolor='black')

    ax.set_xlabel('Resolution (bp)', fontsize=13)
    ax.set_ylabel('SCC', fontsize=13)
    ax.set_title('Hi-C Matrix Conservation (SCC)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([str(r) for r in resolutions])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(y=0, color='black', linewidth=0.5)

    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 0.01,
                    f'{h:.3f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(MATRIX_DIR, 'scc_comparison.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

    # 图2: stratum r 曲线
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, comp, color, title in [
        (axes[0], 'Asu_Ath', COLOR_ATH, 'Asu-Ath (sT)'),
        (axes[1], 'Asu_Aar', COLOR_AAR, 'Asu-Aar (sA)'),
    ]:
        for res in resolutions:
            strata_file = os.path.join(MATRIX_DIR, f"scc_{res}_{comp}_strata.tsv")
            if not os.path.exists(strata_file):
                continue
            strata_df = pd.read_csv(strata_file, sep='\t')
            valid = strata_df[np.isfinite(strata_df['pearson_r'])]
            if len(valid) > 0:
                ax.plot(valid['dist_bins'] * res / 1e6, valid['pearson_r'],
                        lw=1.2, alpha=0.8, label=f'{res//1000}kb')

        ax.set_xlabel('Genomic distance (Mb)', fontsize=12)
        ax.set_ylabel('Stratum Pearson r', fontsize=12)
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(MATRIX_DIR, 'scc_strata_curves.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    scc_comparison.pdf, scc_strata_curves.pdf")


def plot_compartment():
    """Compartment 可视化: 堆积柱状图 + 按染色体"""
    print("  绘制 Compartment 图...")

    summary = pd.read_csv(os.path.join(COMPARTMENT_DIR, 'compartment_conservation_summary.tsv'), sep='\t')
    ath_results = pd.read_csv(os.path.join(COMPARTMENT_DIR, 'compartment_conservation_Asu_Ath.tsv'), sep='\t')
    aar_results = pd.read_csv(os.path.join(COMPARTMENT_DIR, 'compartment_conservation_Asu_Aar.tsv'), sep='\t')

    # 堆积柱状图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for ax, stats_row, title in [
        (ax1, summary[summary['comparison'] == 'Asu_Ath'].iloc[0], 'Asu-Ath (sT)'),
        (ax2, summary[summary['comparison'] == 'Asu_Aar'].iloc[0], 'Asu-Aar (sA)'),
    ]:
        cons = stats_row['conserved_rate'] * 100
        ab = stats_row['A_to_B_rate'] * 100
        ba = stats_row['B_to_A_rate'] * 100
        n = int(stats_row['total_bins'])

        ax.bar([title], [cons], color=COLOR_CONSERVED, label='Conserved', edgecolor='black')
        ax.bar([title], [ab], bottom=[cons], color=COLOR_A_TO_B, label='A→B', edgecolor='black')
        ax.bar([title], [ba], bottom=[cons + ab], color=COLOR_B_TO_A, label='B→A', edgecolor='black')
        ax.set_ylabel('Percentage (%)', fontsize=12)
        ax.set_title(f'{title}\n(n={n})', fontsize=13, fontweight='bold')
        ax.set_ylim(0, 105)
        ax.legend(loc='upper right', fontsize=9)
        ax.text(0, cons/2, f'{cons:.1f}%', ha='center', va='center', fontsize=11, fontweight='bold', color='white')
        ax.text(0, cons + ab/2, f'{ab:.1f}%', ha='center', va='center', fontsize=10)
        ax.text(0, cons + ab + ba/2, f'{ba:.1f}%', ha='center', va='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(COMPARTMENT_DIR, 'compartment_conservation_stacked.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

    # 按染色体
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    for ax, results, color, title in [
        (ax1, ath_results, COLOR_ATH, 'Asu-Ath (sT)'),
        (ax2, aar_results, COLOR_AAR, 'Asu-Aar (sA)'),
    ]:
        chrom_stats = results.groupby('asu_chr').apply(
            lambda x: pd.Series({
                'total': len(x),
                'conserved_rate': (x['conservation_type'] == 'conserved').mean() * 100
            })
        ).reset_index()

        x = np.arange(len(chrom_stats))
        bars = ax.bar(x, chrom_stats['conserved_rate'], 0.6, color=color, edgecolor='black', alpha=0.8)
        ax.set_xlabel('Chromosome', fontsize=12)
        ax.set_ylabel('Conservation Rate (%)', fontsize=12)
        ax.set_title(f'{title} - Compartment Conservation by Chromosome', fontsize=13, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(chrom_stats['asu_chr'], rotation=0)
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.3, axis='y')
        for i, bar in enumerate(bars):
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 1,
                    f'{h:.0f}%\n(n={int(chrom_stats["total"].iloc[i])})',
                    ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(COMPARTMENT_DIR, 'compartment_by_chromosome.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    compartment_conservation_stacked.pdf, compartment_by_chromosome.pdf")


def plot_tad():
    """TAD 可视化: 保守性分数分布 + 对比柱状图"""
    print("  绘制 TAD 图...")

    # 03_tad_conservation.py 把结果按分辨率写到 TAD_DIR/{res}/, 默认看 10kb
    tad_src = os.path.join(TAD_DIR, '10000')
    summary = pd.read_csv(os.path.join(tad_src, 'tad_conservation_summary.tsv'), sep='\t')
    ath_df = pd.read_csv(os.path.join(tad_src, 'tad_conservation_Asu_Ath.tsv'), sep='\t')
    aar_df = pd.read_csv(os.path.join(tad_src, 'tad_conservation_Asu_Aar.tsv'), sep='\t')

    # 图1: 分数分布直方图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for ax, df, color, title in [
        (ax1, ath_df, COLOR_ATH, 'Asu-Ath (sT)'),
        (ax2, aar_df, COLOR_AAR, 'Asu-Aar (sA)'),
    ]:
        if len(df) > 0:
            ax.hist(df['conservation_score'], bins=20, color=color, edgecolor='black', alpha=0.8)
            ax.axvline(df['conservation_score'].mean(), color='red', linestyle='--',
                       label=f'Mean={df["conservation_score"].mean():.3f}')
        ax.set_xlabel('Conservation Score', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{title}\n(n={len(df)})', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.set_xlim(0, 1.05)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(TAD_DIR, 'tad_score_distribution.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

    # 图2: 高保守率对比
    fig, ax = plt.subplots(figsize=(8, 6))
    labels = summary['comparison'].tolist()
    rates = [r * 100 for r in summary['high_conservation_rate'].tolist()] if 'high_conservation_rate' in summary.columns else [0]*len(summary)
    colors = [COLOR_ATH if 'Ath' in l else COLOR_AAR for l in labels]

    bars = ax.bar(labels, rates, color=colors, edgecolor='black', linewidth=1.2)
    ax.set_ylabel('High Conservation Rate (%)', fontsize=13)
    ax.set_title('TAD High Conservation (score ≥ 0.8)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3, axis='y')
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., h + 1,
                f'{h:.1f}%', ha='center', va='bottom', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(TAD_DIR, 'tad_high_conservation.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    tad_score_distribution.pdf, tad_high_conservation.pdf")


def plot_loop():
    """Loop 可视化: 分类柱状图"""
    print("  绘制 Loop 图...")

    summary = pd.read_csv(os.path.join(LOOP_DIR, 'loop_conservation_summary.tsv'), sep='\t')

    categories = ['highly_conserved', 'moderately_conserved', 'weakly_conserved', 'not_conserved', 'unmappable']
    labels = ['Highly\n(≤10kb)', 'Moderately\n(≤30kb)', 'Weakly\n(≤50kb)', 'Not\nConserved', 'Unmappable']
    colors = ['#00A087FF', '#3C5488FF', '#EFC000FF', '#E64B35FF', '#7E6148FF']

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, row, title in [
        (axes[0], summary[summary['comparison'] == 'Asu_Ath'].iloc[0], 'Asu-Ath (sT)'),
        (axes[1], summary[summary['comparison'] == 'Asu_Aar'].iloc[0], 'Asu-Aar (sA)'),
    ]:
        total = row['total']
        vals = [row.get(cat, 0) for cat in categories]
        pcts = [v / total * 100 if total > 0 else 0 for v in vals]

        bars = ax.bar(labels, pcts, color=colors, edgecolor='black')
        ax.set_ylabel('Percentage (%)', fontsize=12)
        ax.set_title(f'{title}\n(total={total})', fontsize=13, fontweight='bold')
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.3, axis='y')
        for bar, v in zip(bars, vals):
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 1,
                    f'{v}', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(LOOP_DIR, 'loop_conservation_classification.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("    loop_conservation_classification.pdf")


def main():
    print("=" * 70)
    print("统一可视化")
    print("=" * 70)

    plot_scc()
    plot_compartment()
    plot_tad()
    plot_loop()

    print(f"\n所有图已保存到: {RESULTS_DIR}/")
    print("=" * 70)


if __name__ == '__main__':
    main()
