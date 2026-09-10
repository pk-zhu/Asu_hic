#!/usr/bin/env python3
"""
Bin/matrix 水平保守性分析 — SCC (Stratum-adjusted Correlation Coefficient)

借鉴 Coolsecture (https://github.com/pk-zhu/Coolsecture) 的 SCC 实现
(见 src/coolsecture/cross_validate.py)。
核心改动: Coolsecture 对整条染色体算 SCC; 本脚本用共线性 bin 对齐后
对齐子矩阵算 SCC, 实现跨物种 (Asu 亚基因组 vs 亲本) 的 Hi-C 矩阵相关性比较。

三个分辨率: 10kb / 25kb / 100kb
两条路径: sT↔Ath, sA↔Aar
"""
import os
import math
from collections import defaultdict

import numpy as np
import pandas as pd
import cooler

from chr_mapping import get_chr_mapping

# ==================== 路径配置 ====================
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
OUTPUT_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results/matrix")

MCOOL_PATHS = {
    'Ath': os.path.join(SINGLE_DIR, "Ath/Matrix/Final/Ath.mcool"),
    'Aar': os.path.join(SINGLE_DIR, "Aar/Matrix/Final/Aar.mcool"),
    'Asu': os.path.join(SINGLE_DIR, "Asu/Matrix/Final/Asu.mcool"),
}

# 共线性路径配置: (collinearity_file, asu_species, ref_species, asu_subgenome_label)
COMPARISONS = [
    {
        'pairs': os.path.join(COLLINEARITY_DIR, "Asu_Ath_orth_bins_{res}.tsv"),
        'asu_species': 'Asu',
        'ref_species': 'Ath',
        'label': 'Asu_Ath',
    },
    {
        'pairs': os.path.join(COLLINEARITY_DIR, "Asu_Aar_orth_bins_{res}.tsv"),
        'asu_species': 'Asu',
        'ref_species': 'Aar',
        'label': 'Asu_Aar',
    },
]

RESOLUTIONS = [10000, 25000, 100000]
MAX_DIST_MB = 10.0
MIN_DIST_BINS = 1


# ==================== SCC 核心函数 (借鉴 Coolsecture) ====================

def _pearson_from_stats(n, sx, sy, sxx, syy, sxy):
    """从累积统计量计算 Pearson 相关系数"""
    if n < 2:
        return float("nan")
    vx = sxx - (sx * sx) / n
    vy = syy - (sy * sy) / n
    if vx <= 0 or vy <= 0:
        return float("nan")
    cov = sxy - (sx * sy) / n
    return cov / math.sqrt(vx * vy)


def _accumulate_stats(stats, d, x, y):
    """累积一个 (x, y) 对到距离 stratum d 的统计量中"""
    s = stats[d]
    s[0] += 1
    s[1] += x
    s[2] += y
    s[3] += x * x
    s[4] += y * y
    s[5] += x * y


def _scc_from_aligned_matrices(mat_a, mat_b, max_dist_bins, min_dist_bins):
    """
    对两个对齐的 N×N 稀疏矩阵计算 SCC。
    mat_a, mat_b: scipy sparse COO matrices, 相同 shape (N, N)
    每个矩阵的 bin i 在两个物种中对应共线性区域。
    """
    stats = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # 构建 mat_b 的字典: (i, j) -> value (仅上三角, j >= i)
    dict_b = {}
    for i, j, v in zip(mat_b.row, mat_b.col, mat_b.data):
        if j < i:
            i, j = j, i
        d = j - i
        if d < min_dist_bins or d > max_dist_bins:
            continue
        # 过滤 NaN
        if not np.isfinite(v):
            continue
        dict_b[(i, j)] = float(v)

    # 遍历 mat_a, 累积 (x, y) 对
    seen = set()
    for i, j, v in zip(mat_a.row, mat_a.col, mat_a.data):
        if j < i:
            i, j = j, i
        d = j - i
        if d < min_dist_bins or d > max_dist_bins:
            continue
        if not np.isfinite(v):
            continue
        y = dict_b.get((i, j), 0.0)
        _accumulate_stats(stats, d, float(v), y)
        seen.add((i, j))

    # 补充 mat_b 中有但 mat_a 中没有的 contact
    for (i, j), v in dict_b.items():
        if (i, j) in seen:
            continue
        d = j - i
        if d < min_dist_bins or d > max_dist_bins:
            continue
        _accumulate_stats(stats, d, 0.0, float(v))

    # 计算每层 r 和 SCC
    rows = []
    for d, (n, sx, sy, sxx, syy, sxy) in stats.items():
        r = _pearson_from_stats(n, sx, sy, sxx, syy, sxy)
        w = float(n) if np.isfinite(r) else float("nan")
        rows.append((d, n, r, w))
    rows.sort(key=lambda x: x[0])

    weights = [w for _, _, r, w in rows if np.isfinite(r) and np.isfinite(w) and w > 0]
    rvals = [r for _, _, r, w in rows if np.isfinite(r) and np.isfinite(w) and w > 0]
    if weights:
        scc = float(np.average(rvals, weights=weights))
    else:
        scc = float("nan")

    return rows, scc


# ==================== 辅助函数 ====================

def load_collinear_pairs(filepath):
    """
    加载共线性 bin 对 (7列: asu_chr asu_start asu_end ref_chr ref_start ref_end identity)
    过滤:
      1. 固定染色体映射 (基于 PAF matches, 不受分辨率影响)
      2. bin 去重 (同一 asu_bin 保留 identity 最高的)
    返回: list of (asu_chr, asu_start, asu_end, ref_chr, ref_start, ref_end)
    """
    all_pairs = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                identity = float(parts[6]) if len(parts) >= 7 else 0.0
                all_pairs.append((parts[0], int(parts[1]), int(parts[2]),
                                  parts[3], int(parts[4]), int(parts[5]), identity))

    # 1. 固定染色体映射
    chr_mapping = get_chr_mapping(filepath)
    after_chr = [p for p in all_pairs if p[3] == chr_mapping.get(p[0])]

    # 2. bin 去重: 同一 (asu_chr, asu_start) 保留 identity 最高的
    best_by_bin = {}
    for p in after_chr:
        key = (p[0], p[1])
        if key not in best_by_bin or p[6] > best_by_bin[key][6]:
            best_by_bin[key] = p

    filtered = list(best_by_bin.values())
    removed_chr = len(all_pairs) - len(after_chr)
    removed_dup = len(after_chr) - len(filtered)
    print(f"    染色体映射: {len(all_pairs)} -> {len(after_chr)} (移除 {removed_chr} 跨染色体)")
    print(f"    bin去重: {len(after_chr)} -> {len(filtered)} (移除 {removed_dup} 重复)")
    return [(p[0], p[1], p[2], p[3], p[4], p[5]) for p in filtered]


def extract_submatrix(clr, chrom, bin_starts, bin_size):
    """
    从 cooler 提取指定染色体上的子矩阵 (balanced)。
    fetch 整条染色体矩阵 (scipy sparse), 然后用 fancy indexing 取子矩阵。
    bin_starts: 该染色体上的 bin 起始坐标列表
    """
    n = len(bin_starts)
    if n < 2:
        return None
    # fetch 整条染色体矩阵 (转 CSR 以支持 fancy indexing)
    full_mat = clr.matrix(balance=True, sparse=True).fetch(chrom).tocsr()
    # bin 起始坐标 → 染色体内相对 bin index
    indices = [s // bin_size for s in bin_starts]
    # 子矩阵
    sub = full_mat[indices, :][:, indices].tocoo()
    # 过滤 NaN
    mask = np.isfinite(sub.data)
    sub.data = sub.data[mask]
    sub.row = sub.row[mask]
    sub.col = sub.col[mask]
    return sub


def compute_scc_for_comparison(comp, resolution):
    """
    对一条比较路径 (如 Asu_Ath) 在一个分辨率下计算 SCC。
    返回: (per_chrom_results, overall_scc, overall_strata)
    """
    pairs_file = comp['pairs'].format(res=resolution)
    if not os.path.exists(pairs_file):
        print(f"    共线性文件不存在: {pairs_file}")
        return [], float("nan"), []

    asu_species = comp['asu_species']
    ref_species = comp['ref_species']

    # 加载 mcool
    asu_mcool = f"{MCOOL_PATHS[asu_species]}::/resolutions/{resolution}"
    ref_mcool = f"{MCOOL_PATHS[ref_species]}::/resolutions/{resolution}"
    clr_asu = cooler.Cooler(asu_mcool)
    clr_ref = cooler.Cooler(ref_mcool)

    # 加载共线性对, 按 (asu_chr, ref_chr) 分组
    pairs = load_collinear_pairs(pairs_file)
    pairs_by_chrom = defaultdict(list)
    for asu_chr, asu_start, asu_end, ref_chr, ref_start, ref_end in pairs:
        pairs_by_chrom[(asu_chr, ref_chr)].append(
            (asu_start, ref_start))

    print(f"    {comp['label']}: {len(pairs)} 共线性对, {len(pairs_by_chrom)} 染色体对")

    max_dist_bins = int(MAX_DIST_MB * 1e6 / resolution)

    per_chrom_results = []
    all_strata = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])

    for (asu_chr, ref_chr) in sorted(pairs_by_chrom.keys()):
        chrom_pairs = pairs_by_chrom[(asu_chr, ref_chr)]
        asu_bin_starts = [p[0] for p in chrom_pairs]
        ref_bin_starts = [p[1] for p in chrom_pairs]
        n_bins = len(asu_bin_starts)

        if n_bins < 2:
            continue

        # 提取对齐子矩阵
        mat_asu = extract_submatrix(clr_asu, asu_chr, asu_bin_starts, resolution)
        mat_ref = extract_submatrix(clr_ref, ref_chr, ref_bin_starts, resolution)

        if mat_asu is None or mat_ref is None:
            continue

        # 计算 SCC
        strata, scc = _scc_from_aligned_matrices(mat_asu, mat_ref, max_dist_bins, MIN_DIST_BINS)

        n_strata = sum(1 for _, _, r, _ in strata if np.isfinite(r))
        total_weight = sum(w for _, _, r, w in strata if np.isfinite(r) and np.isfinite(w) and w > 0)

        per_chrom_results.append({
            'asu_chr': asu_chr,
            'ref_chr': ref_chr,
            'n_bins': n_bins,
            'n_strata': n_strata,
            'total_weight': total_weight,
            'scc': scc,
        })

        # 累积到全局 strata
        for d, n, r, w in strata:
            if np.isfinite(r) and w > 0:
                s = all_strata[d]
                s[0] += 1
                s[1] += r * w
                s[2] += w

        print(f"      {asu_chr}↔{ref_chr}: n_bins={n_bins}, scc={scc:.4f}, n_strata={n_strata}")

    # 计算 overall SCC (按 total_weight 加权平均 per-chromosome SCC)
    valid = [r for r in per_chrom_results if np.isfinite(r['scc']) and r['total_weight'] > 0]
    if valid:
        weights = [r['total_weight'] for r in valid]
        sccs = [r['scc'] for r in valid]
        overall_scc = float(np.average(sccs, weights=weights))
    else:
        overall_scc = float("nan")

    # 构建 overall strata rows (for plotting)
    overall_strata = []
    for d in sorted(all_strata.keys()):
        s = all_strata[d]
        if s[2] > 0:
            r = s[1] / s[2]
            overall_strata.append((d, int(s[0]), r, s[2]))

    return per_chrom_results, overall_scc, overall_strata


def main():
    print("=" * 70)
    print("Bin/matrix 水平保守性分析 (SCC)")
    print("=" * 70)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    summary_rows = []

    for resolution in RESOLUTIONS:
        print(f"\n{'=' * 70}")
        print(f"分辨率: {resolution} bp")
        print(f"{'=' * 70}")

        for comp in COMPARISONS:
            print(f"\n  处理 {comp['label']}...")

            per_chrom, overall_scc, overall_strata = compute_scc_for_comparison(comp, resolution)

            if not per_chrom:
                print(f"    无结果, 跳过")
                continue

            # 保存 per-chromosome 结果
            chrom_file = os.path.join(OUTPUT_DIR, f"scc_{resolution}_{comp['label']}.tsv")
            chrom_df = pd.DataFrame(per_chrom)
            chrom_df.to_csv(chrom_file, sep='\t', index=False)
            print(f"    保存: {os.path.basename(chrom_file)}")

            # 保存 per-stratum 结果
            strata_file = os.path.join(OUTPUT_DIR, f"scc_{resolution}_{comp['label']}_strata.tsv")
            with open(strata_file, 'w') as f:
                f.write("dist_bins\tn_chroms\tpearson_r\ttotal_weight\n")
                for d, n, r, w in overall_strata:
                    r_str = f"{r:.6f}" if np.isfinite(r) else "nan"
                    f.write(f"{d}\t{n}\t{r_str}\t{w:.6f}\n")
            print(f"    保存: {os.path.basename(strata_file)}")

            # 汇总
            n_total_bins = sum(r['n_bins'] for r in per_chrom)
            n_valid_bins = sum(r['n_bins'] for r in per_chrom if np.isfinite(r['scc']) and r['total_weight'] > 0)
            summary_rows.append({
                'resolution': resolution,
                'comparison': comp['label'],
                'scc': overall_scc,
                'n_chroms': len(per_chrom),
                'n_bins_total': n_total_bins,
                'n_bins_valid': n_valid_bins,
                'unmapped_bin_rate': 1 - n_valid_bins / n_total_bins if n_total_bins > 0 else 0,
            })

            print(f"    Overall SCC: {overall_scc:.4f}")
            print(f"    Valid bins: {n_valid_bins}/{n_total_bins} ({n_valid_bins/n_total_bins:.1%})")

    # 保存汇总
    summary_df = pd.DataFrame(summary_rows)
    summary_file = os.path.join(OUTPUT_DIR, "scc_summary.tsv")
    summary_df.to_csv(summary_file, sep='\t', index=False)

    print(f"\n{'=' * 70}")
    print("SCC 汇总:")
    print(summary_df.to_string(index=False))
    print(f"\n汇总已保存: {summary_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
