#!/usr/bin/env python3
"""
TAD 水平保守性分析 (multi-resolution: 5k, 10k, 25k)

使用 HiCExplorer TAD 结果 + 对应分辨率共线性对
方法: 为每个 Asu TAD 找到亲本中重叠 bin 最多的 TAD, 计算 conservation_score

Output: results/tad/{resolution}/
  tad_conservation_Asu_Ath.tsv
  tad_conservation_Asu_Aar.tsv
  tad_conservation_summary.tsv
"""
import os
import pandas as pd
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
OUTPUT_BASE = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results/tad")

RESOLUTIONS = [5000, 10000, 25000, 40000]

ASU_ATH_SUBGENOME = {'sT1', 'sT2', 'sT3', 'sT4', 'sT5'}
ASU_AAR_SUBGENOME = {'sA6', 'sA7', 'sA8', 'sA9', 'sA10', 'sA11', 'sA12', 'sA13'}


def get_tad_paths(resolution):
    res_str = str(resolution)
    return {
        'Ath': os.path.join(SINGLE_DIR, f"Ath/TADs/HiCExplorer/{res_str}/Ath_domains.bed"),
        'Aar': os.path.join(SINGLE_DIR, f"Aar/TADs/HiCExplorer/{res_str}/Aar_domains.bed"),
        'Asu': os.path.join(SINGLE_DIR, f"Asu/TADs/HiCExplorer/{res_str}/Asu_domains.bed"),
    }


def get_collinear_paths(resolution):
    res_str = str(resolution)
    return {
        'Asu_Ath': os.path.join(COLLINEARITY_DIR, f"Asu_Ath_orth_bins_{res_str}.tsv"),
        'Asu_Aar': os.path.join(COLLINEARITY_DIR, f"Asu_Aar_orth_bins_{res_str}.tsv"),
    }


def load_tads(bed_path):
    """加载 TAD domains.bed"""
    tads = []
    if not os.path.exists(bed_path):
        return tads
    with open(bed_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                tads.append({
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3]
                })
    return tads


def load_collinear_pairs(filepath):
    """读取共线性 bin 对 (7列)
    过滤: 1) 固定染色体映射 2) bin去重(同asu_bin保留最高identity)"""
    all_pairs = []
    if not os.path.exists(filepath):
        return []
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


def calculate_tad_conservation(asu_tads, other_tads, collinear_pairs, label):
    """计算 TAD 保守性: 为每个 Asu TAD 找到亲本中重叠 bin 最多的 TAD"""
    print(f"\n  计算 TAD 保守性: {label}")

    # 构建 Asu bin → 亲本 bin 映射
    asu_to_other = defaultdict(list)
    for asu_chr, asu_start, asu_end, ref_chr, ref_start, ref_end in collinear_pairs:
        asu_to_other[(asu_chr, asu_start, asu_end)].append((ref_chr, ref_start, ref_end))

    # 亲本 TAD 按染色体索引
    other_tads_by_chr = defaultdict(list)
    for tad in other_tads:
        other_tads_by_chr[tad['chr']].append(tad)

    results = []
    n_total_tads = len(asu_tads)
    n_unmapped = 0

    for asu_tad in asu_tads:
        asu_chr = asu_tad['chr']
        asu_start = asu_tad['start']
        asu_end = asu_tad['end']

        # 找到与这个 TAD 重叠的共线性 bin
        overlapping_bins = []
        target_bins = []
        for (bin_chr, bin_start, bin_end), tgt_list in asu_to_other.items():
            if bin_chr == asu_chr and not (bin_end <= asu_start or bin_start >= asu_end):
                overlapping_bins.append((bin_chr, bin_start, bin_end))
                target_bins.extend(tgt_list)

        if not overlapping_bins:
            n_unmapped += 1
            continue

        # 统计目标物种中与共线性 bin 重叠的 TAD
        target_chr_count = defaultdict(int)
        for tgt_bin in target_bins:
            target_chr_count[tgt_bin[0]] += 1

        if not target_chr_count:
            n_unmapped += 1
            continue

        best_target_chr = max(target_chr_count.items(), key=lambda x: x[1])[0]

        # 找到亲本中重叠 bin 最多的 TAD
        tads_in_target = other_tads_by_chr.get(best_target_chr, [])
        max_overlap = 0
        best_other_tad = None

        for other_tad in tads_in_target:
            bin_count = 0
            for tgt_bin in target_bins:
                if tgt_bin[0] == other_tad['chr']:
                    if not (tgt_bin[2] <= other_tad['start'] or tgt_bin[1] >= other_tad['end']):
                        bin_count += 1
            if bin_count > max_overlap:
                max_overlap = bin_count
                best_other_tad = other_tad

        if best_other_tad and max_overlap > 0:
            conservation_score = min(1.0, max_overlap / len(overlapping_bins))
            results.append({
                'asu_chr': asu_chr,
                'asu_tad_start': asu_start,
                'asu_tad_end': asu_end,
                'asu_tad_name': asu_tad['name'],
                'other_chr': best_other_tad['chr'],
                'other_tad_start': best_other_tad['start'],
                'other_tad_end': best_other_tad['end'],
                'other_tad_name': best_other_tad['name'],
                'conservation_score': conservation_score,
                'overlapping_bins_in_tad': max_overlap,
                'total_overlapping_bins': len(overlapping_bins),
            })

    df = pd.DataFrame(results)
    n_mapped = len(df)
    if n_mapped > 0:
        print(f"    TAD 保守性: mean={df['conservation_score'].mean():.4f}, "
              f"median={df['conservation_score'].median():.4f}, n={n_mapped}")
        print(f"    Unmapped: {n_unmapped}/{n_total_tads} ({n_unmapped/n_total_tads:.2%})")
    return df, n_total_tads, n_unmapped


def run_resolution(resolution):
    print(f"\n{'='*70}")
    print(f"TAD 保守性分析 ({resolution}bp)")
    print(f"{'='*70}")

    tad_paths = get_tad_paths(resolution)
    collinear_paths = get_collinear_paths(resolution)
    output_dir = os.path.join(OUTPUT_BASE, str(resolution))
    os.makedirs(output_dir, exist_ok=True)

    # 检查文件存在性
    missing = []
    for k, p in {**tad_paths, **collinear_paths}.items():
        if not os.path.exists(p):
            missing.append(f"{k}: {p}")
    if missing:
        print(f"  跳过: 缺少文件:\n    " + "\n    ".join(missing))
        return

    # 1. 加载 TAD
    print("\n1. 加载 TAD 数据...")
    asu_tads = load_tads(tad_paths['Asu'])
    ath_tads = load_tads(tad_paths['Ath'])
    aar_tads = load_tads(tad_paths['Aar'])
    print(f"  Asu TADs: {len(asu_tads)}")
    print(f"  Ath TADs: {len(ath_tads)}")
    print(f"  Aar TADs: {len(aar_tads)}")

    asu_tads_ath = [t for t in asu_tads if t['chr'] in ASU_ATH_SUBGENOME]
    asu_tads_aar = [t for t in asu_tads if t['chr'] in ASU_AAR_SUBGENOME]
    print(f"  Asu sT TADs: {len(asu_tads_ath)}")
    print(f"  Asu sA TADs: {len(asu_tads_aar)}")

    # 2. 加载共线性
    print(f"\n2. 加载共线性数据 ({resolution}bp)...")
    asu_ath_pairs = load_collinear_pairs(collinear_paths['Asu_Ath'])
    asu_aar_pairs = load_collinear_pairs(collinear_paths['Asu_Aar'])
    print(f"  Asu-Ath: {len(asu_ath_pairs)} 对")
    print(f"  Asu-Aar: {len(asu_aar_pairs)} 对")

    # 3. 计算 TAD 保守性
    print("\n3. 计算 TAD 保守性...")
    ath_df = ath_total = ath_unmapped = None
    if os.path.exists(tad_paths['Ath']) and os.path.exists(collinear_paths['Asu_Ath']):
        ath_df, ath_total, ath_unmapped = calculate_tad_conservation(asu_tads_ath, ath_tads, asu_ath_pairs, 'Asu-Ath')
    else:
        print(f"  跳过 Asu_Ath: 缺少 Ath TAD 或共线性 bin")
    aar_df, aar_total, aar_unmapped = calculate_tad_conservation(asu_tads_aar, aar_tads, asu_aar_pairs, 'Asu-Aar')

    # 4. 保存结果
    print("\n4. 保存结果...")
    if ath_df is not None:
        ath_df.to_csv(os.path.join(output_dir, 'tad_conservation_Asu_Ath.tsv'), sep='\t', index=False)
    aar_df.to_csv(os.path.join(output_dir, 'tad_conservation_Asu_Aar.tsv'), sep='\t', index=False)

    # 5. 汇总
    summary_data = []
    for label, df, n_total, n_unmapped in [
        ('Asu_Ath', ath_df, ath_total, ath_unmapped),
        ('Asu_Aar', aar_df, aar_total, aar_unmapped),
    ]:
        if df is None:
            continue
        if len(df) > 0:
            summary_data.append({
                'comparison': label,
                'n_total_tads': n_total,
                'n_tads': len(df),
                'n_unmapped': n_unmapped,
                'unmapped_rate': n_unmapped / n_total if n_total > 0 else 0,
                'mean_score': df['conservation_score'].mean(),
                'median_score': df['conservation_score'].median(),
                'n_high_conservation': int((df['conservation_score'] >= 0.8).sum()),
                'high_conservation_rate': float((df['conservation_score'] >= 0.8).mean()),
            })
        else:
            summary_data.append({'comparison': label, 'n_tads': 0, 'n_total_tads': n_total,
                                 'n_unmapped': n_unmapped, 'unmapped_rate': n_unmapped / n_total if n_total > 0 else 0})

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'tad_conservation_summary.tsv'), sep='\t', index=False)

    print(f"\n  汇总已保存: {output_dir}/tad_conservation_summary.tsv")
    print(summary_df.to_string(index=False))


def main():
    for res in RESOLUTIONS:
        run_resolution(res)
    print("\n" + "="*70)
    print("所有分辨率 TAD 保守性分析完成")
    print("="*70)


if __name__ == '__main__':
    main()
