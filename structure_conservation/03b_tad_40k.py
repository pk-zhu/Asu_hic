#!/usr/bin/env python3
"""
03b_tad_40k.py — 40k bp 分辨率 TAD 保守性分析.

复用 03_tad_conservation.py 的 calculate_tad_conservation 函数,
但只对有 40k TAD 的物种 (Aar, Asu) 计算.
Ath 在 HiCExplorer 下没有 40k TAD 输出, 因此跳过 Asu_Ath.

注: 共线性 bin 用现成的 10k bin (collinearity/04_collinear_pairs/Asu_Aar_orth_bins_10000.tsv),
因为 02 模块未生成 40k bin. bin 尺寸与 TAD 尺寸不严格匹配, 但仍是有效评估
(40k TAD 内的 10k bin 重叠率).

输出: results/tad/40000/
  tad_conservation_Asu_Aar.tsv
  tad_conservation_summary.tsv
"""
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chr_mapping import get_chr_mapping
from importlib import import_module
main_mod = import_module("03_tad_conservation")
calculate_tad_conservation = main_mod.calculate_tad_conservation
get_tad_paths = main_mod.get_tad_paths
get_collinear_paths = main_mod.get_collinear_paths
OUTPUT_BASE = main_mod.OUTPUT_BASE

RESOLUTION = 40000
COLLIN_RES = 10000  # 复用 10k bin


def get_paths_adapted(resolution, collin_res):
    """对没有 40k TAD 的 Ath 跳过, bin 用 10k."""
    return {
        'Aar': os.path.join(main_mod.SINGLE_DIR, f"Aar/TADs/HiCExplorer/{resolution}/Aar_domains.bed"),
        'Asu': os.path.join(main_mod.SINGLE_DIR, f"Asu/TADs/HiCExplorer/{resolution}/Asu_domains.bed"),
        'Asu_Aar': os.path.join(main_mod.COLLINEARITY_DIR, f"Asu_Aar_orth_bins_{collin_res}.tsv"),
        'Aar_Asu': os.path.join(main_mod.COLLINEARITY_DIR, f"Aar_Asu_orth_bins_{collin_res}.tsv"),
    }


def load_tad(path):
    """读 TAD BED: chr start end name score strand start2 end2 rgb."""
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'rgb'])
    return df[['chr', 'start', 'end', 'name']].to_dict('records')


def load_collinear(path):
    """读共线性 bin TSV: asu_chr asu_start asu_end other_chr other_start other_end."""
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['asu_chr', 'asu_start', 'asu_end', 'other_chr', 'other_start', 'other_end'])
    return df


def main():
    print(f"=== TAD 保守性分析 ({RESOLUTION}bp, bin={COLLIN_RES}bp) ===")

    paths = get_paths_adapted(RESOLUTION, COLLIN_RES)
    out_dir = os.path.join(OUTPUT_BASE, str(RESOLUTION))
    os.makedirs(out_dir, exist_ok=True)

    # 加载
    asu_tads = load_tad(paths['Asu'])
    aar_tads = load_tad(paths['Aar'])
    asu_aar_pairs = load_collinear(paths['Asu_Aar'])

    print(f"  Asu TADs: {len(asu_tads)}, Aar TADs: {len(aar_tads)}, bins: {len(asu_aar_pairs)}")

    # DataFrame → list of tuples (calculate_tad_conservation 期望 6 元组)
    asu_aar_list = list(asu_aar_pairs.itertuples(index=False, name=None))

    # 只算 Asu_Aar
    aar_df, aar_total, aar_unmapped = calculate_tad_conservation(
        asu_tads, aar_tads, asu_aar_list, 'Asu-Aar')
    aar_df.to_csv(os.path.join(out_dir, 'tad_conservation_Asu_Aar.tsv'),
                  sep='\t', index=False)

    summary = [{
        'comparison': 'Asu_Aar',
        'n_total_tads': aar_total,
        'n_tads': len(aar_df),
        'n_unmapped': aar_unmapped,
        'unmapped_rate': round(aar_unmapped / aar_total, 4) if aar_total > 0 else 0,
        'mean_score': round(aar_df['conservation_score'].mean(), 4) if len(aar_df) > 0 else 0,
        'median_score': round(aar_df['conservation_score'].median(), 4) if len(aar_df) > 0 else 0,
        'n_high_conservation': int((aar_df['conservation_score'] >= 0.8).sum()) if len(aar_df) > 0 else 0,
        'high_conservation_rate': round((aar_df['conservation_score'] >= 0.8).mean(), 4) if len(aar_df) > 0 else 0,
    }]
    pd.DataFrame(summary).to_csv(os.path.join(out_dir, 'tad_conservation_summary.tsv'),
                                   sep='\t', index=False)
    print(f"\n  保存: {out_dir}/tad_conservation_summary.tsv")
    print("=== DONE ===")


if __name__ == "__main__":
    main()