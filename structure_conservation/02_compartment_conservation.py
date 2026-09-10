#!/usr/bin/env python3
"""
Compartment 水平保守性分析 (100kb)

使用 CALDER2 100kb compartment 结果 + 100kb 共线性对
比较 Asu 亚基因组与亲本的 A/B compartment 状态
"""
import os
import pandas as pd
import numpy as np
from collections import defaultdict

from chr_mapping import get_chr_mapping

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
COLLINEARITY_DIR = os.path.join(BASE_DIR, "02.collinearity/04_collinear_pairs")
SINGLE_DIR = os.path.join(BASE_DIR, "01.single")
OUTPUT_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results/compartment")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COMPARTMENT_PATHS = {
    'Ath': os.path.join(SINGLE_DIR, "Ath/Compartments/CALDER2/100000/Ath/sub_compartments/all_sub_compartments.tsv"),
    'Aar': os.path.join(SINGLE_DIR, "Aar/Compartments/CALDER2/100000/Aar/sub_compartments/all_sub_compartments.tsv"),
    'Asu': os.path.join(SINGLE_DIR, "Asu/Compartments/CALDER2/100000/Asu/sub_compartments/all_sub_compartments.tsv"),
}

COLLINEAR_PAIRS_PATHS = {
    'Asu_Ath': os.path.join(COLLINEARITY_DIR, "Asu_Ath_orth_bins_100000.tsv"),
    'Asu_Aar': os.path.join(COLLINEARITY_DIR, "Asu_Aar_orth_bins_100000.tsv"),
}


def normalize_chr_name(chr_name, species):
    """标准化染色体名称以匹配共线性文件命名"""
    chr_name = str(chr_name)
    if species == 'Aar':
        # chr__1 -> chr_1
        if chr_name.startswith('chr__'):
            chr_name = 'chr_' + chr_name[5:]
    elif species == 'Asu':
        # chrsT1 -> sT1, chrsA6 -> sA6
        if chr_name.startswith('chrs'):
            chr_name = chr_name[3:]
    return chr_name


def load_compartments(filepath, species):
    """读取 compartment 数据, 返回 {chr: [(start, end, comp_state), ...]}"""
    compartments = defaultdict(list)
    df = pd.read_csv(filepath, sep='\t')
    for _, row in df.iterrows():
        chr_name = normalize_chr_name(str(row.iloc[0]), species)
        start = int(row['pos_start'])
        end = int(row['pos_end'])
        comp_name = str(row['comp_name'])
        comp_state = 'A' if comp_name.startswith('A') else 'B'
        compartments[chr_name].append((start, end, comp_state))
    return compartments


def load_collinear_pairs(filepath):
    """读取共线性 bin 对 (7列)
    过滤: 1) 固定染色体映射 2) bin去重(同asu_bin保留最高identity)"""
    all_pairs = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                identity = float(parts[6]) if len(parts) >= 7 else 0.0
                all_pairs.append((parts[0], int(parts[1]), int(parts[2]),
                                  parts[3], int(parts[4]), int(parts[5]), identity))

    chr_mapping = get_chr_mapping(filepath)
    after_chr = [p for p in all_pairs if p[3] == chr_mapping.get(p[0])]

    best_by_bin = {}
    for p in after_chr:
        key = (p[0], p[1])
        if key not in best_by_bin or p[6] > best_by_bin[key][6]:
            best_by_bin[key] = p
    filtered = list(best_by_bin.values())
    print(f"  映射+去重: {len(all_pairs)} -> {len(filtered)} "
          f"(chr={len(all_pairs)-len(after_chr)}, dup={len(after_chr)-len(filtered)})")
    return [(p[0], p[1], p[2], p[3], p[4], p[5]) for p in filtered]


def find_matching_compartment(bin_chr, bin_start, bin_end, compartment_dict):
    """为一个 bin 找到对应的 compartment 状态 (多数投票)"""
    if bin_chr not in compartment_dict:
        return None
    overlapping = []
    for comp_start, comp_end, comp_state in compartment_dict[bin_chr]:
        if comp_start < bin_end and comp_end > bin_start:
            overlapping.append(comp_state)
    if not overlapping:
        return None
    a_count = overlapping.count('A')
    b_count = overlapping.count('B')
    return 'A' if a_count >= b_count else 'B'


def analyze_compartment_conservation(asu_comps, ref_comps, collinear_pairs, comparison_label):
    """分析 compartment 保守性"""
    results = []
    conserved = 0
    a_to_b = 0
    b_to_a = 0
    total = 0
    unmapped = 0

    for asu_chr, asu_start, asu_end, ref_chr, ref_start, ref_end in collinear_pairs:
        asu_comp = find_matching_compartment(asu_chr, asu_start, asu_end, asu_comps)
        ref_comp = find_matching_compartment(ref_chr, ref_start, ref_end, ref_comps)

        if asu_comp is None or ref_comp is None:
            unmapped += 1
            continue

        total += 1

        if asu_comp == ref_comp:
            conserved += 1
            conservation_type = 'conserved'
        elif ref_comp == 'A' and asu_comp == 'B':
            a_to_b += 1
            conservation_type = 'A_to_B'
        else:
            b_to_a += 1
            conservation_type = 'B_to_A'

        results.append({
            'asu_chr': asu_chr,
            'asu_start': asu_start,
            'asu_end': asu_end,
            'ref_chr': ref_chr,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'asu_compartment': asu_comp,
            'ref_compartment': ref_comp,
            'conservation_type': conservation_type,
        })

    stats = {
        'total_collinear_bins': len(collinear_pairs),
        'total': total,
        'unmapped': unmapped,
        'unmapped_rate': unmapped / len(collinear_pairs) if len(collinear_pairs) > 0 else 0,
        'conserved': conserved,
        'conserved_rate': conserved / total if total > 0 else 0,
        'a_to_b': a_to_b,
        'a_to_b_rate': a_to_b / total if total > 0 else 0,
        'b_to_a': b_to_a,
        'b_to_a_rate': b_to_a / total if total > 0 else 0,
    }
    return pd.DataFrame(results), stats


def main():
    print("=" * 70)
    print("Compartment 保守性分析 (100kb)")
    print("=" * 70)

    # 1. 读取 compartment 数据
    print("\n1. 读取 compartment 数据...")
    asu_comps = load_compartments(COMPARTMENT_PATHS['Asu'], 'Asu')
    ath_comps = load_compartments(COMPARTMENT_PATHS['Ath'], 'Ath')
    aar_comps = load_compartments(COMPARTMENT_PATHS['Aar'], 'Aar')
    print(f"  Asu: {sum(len(v) for v in asu_comps.values())} bins")
    print(f"  Ath: {sum(len(v) for v in ath_comps.values())} bins")
    print(f"  Aar: {sum(len(v) for v in aar_comps.values())} bins")

    # 2. 读取共线性数据
    print("\n2. 读取共线性数据 (100kb)...")
    asu_ath_pairs = load_collinear_pairs(COLLINEAR_PAIRS_PATHS['Asu_Ath'])
    asu_aar_pairs = load_collinear_pairs(COLLINEAR_PAIRS_PATHS['Asu_Aar'])
    print(f"  Asu-Ath: {len(asu_ath_pairs)} 对")
    print(f"  Asu-Aar: {len(asu_aar_pairs)} 对")

    # 3. 分析
    print("\n3. 分析 Asu-Ath 亚基因组...")
    ath_results, ath_stats = analyze_compartment_conservation(
        asu_comps, ath_comps, asu_ath_pairs, 'Asu_Ath')

    print("\n4. 分析 Asu-Aar 亚基因组...")
    aar_results, aar_stats = analyze_compartment_conservation(
        asu_comps, aar_comps, asu_aar_pairs, 'Asu_Aar')

    # 5. 保存结果
    print("\n5. 保存结果...")
    ath_results.to_csv(os.path.join(OUTPUT_DIR, 'compartment_conservation_Asu_Ath.tsv'), sep='\t', index=False)
    aar_results.to_csv(os.path.join(OUTPUT_DIR, 'compartment_conservation_Asu_Aar.tsv'), sep='\t', index=False)

    # 6. 统计摘要
    print("\n" + "=" * 70)
    print("统计摘要")
    print("=" * 70)
    for label, stats in [("Asu-Ath", ath_stats), ("Asu-Aar", aar_stats)]:
        print(f"\n--- {label} ---")
        print(f"  Collinear bins: {stats['total_collinear_bins']}")
        print(f"  Mapped: {stats['total']}, Unmapped: {stats['unmapped']} ({stats['unmapped_rate']:.2%})")
        print(f"  Conserved: {stats['conserved']} ({stats['conserved_rate']:.2%})")
        print(f"  A→B: {stats['a_to_b']} ({stats['a_to_b_rate']:.2%})")
        print(f"  B→A: {stats['b_to_a']} ({stats['b_to_a_rate']:.2%})")

    summary_df = pd.DataFrame({
        'comparison': ['Asu_Ath', 'Asu_Aar'],
        'total_collinear_bins': [ath_stats['total_collinear_bins'], aar_stats['total_collinear_bins']],
        'total_bins': [ath_stats['total'], aar_stats['total']],
        'unmapped': [ath_stats['unmapped'], aar_stats['unmapped']],
        'unmapped_rate': [ath_stats['unmapped_rate'], aar_stats['unmapped_rate']],
        'conserved': [ath_stats['conserved'], aar_stats['conserved']],
        'conserved_rate': [ath_stats['conserved_rate'], aar_stats['conserved_rate']],
        'A_to_B': [ath_stats['a_to_b'], aar_stats['a_to_b']],
        'A_to_B_rate': [ath_stats['a_to_b_rate'], aar_stats['a_to_b_rate']],
        'B_to_A': [ath_stats['b_to_a'], aar_stats['b_to_a']],
        'B_to_A_rate': [ath_stats['b_to_a_rate'], aar_stats['b_to_a_rate']],
    })
    summary_df.to_csv(os.path.join(OUTPUT_DIR, 'compartment_conservation_summary.tsv'), sep='\t', index=False)
    print(f"\n汇总已保存: {OUTPUT_DIR}/compartment_conservation_summary.tsv")


if __name__ == '__main__':
    main()
