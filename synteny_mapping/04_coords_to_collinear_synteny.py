#!/usr/bin/env python3
"""
亚基因组感知的共线性 bin 对提取 (多分辨率, 多对多)
— MUMmer 版 (替换原 PAF 版 04_paf_to_collinear_synteny.py)

输入: 02.collinearity/03_mummer/{sT_vs_Ath,sA_vs_Aar}.coords
  show-coords -rclTH 输出, 13 列 tab 无表头:
    S1 E1 S2 E2 LEN1 LEN2 IDY LEN_R LEN_Q COV_R COV_Q TAG_R TAG_Q
  调用约定: nucmer <ref=亲本> <query=Asu亚基因组>
    → S1/E1=亲本 (1-based inclusive),
       S2/E2=Asu  (反链时 S2>E2),
       TAG_R=亲本染色体, TAG_Q=Asu 染色体,
       IDY=百分比 0-100.

策略 (与旧脚本一致, 保证下游 03/06 兼容):
  1. 解析比对, identity = IDY/100 (真实 base-level)
  2. 每个 Asu bin 找所有重叠比对, 按链向线性插值映射到亲本 bin (多对多)
  3. 邻近一致性过滤: 移除无相邻 bin 同 ref_chr 支持的孤立匹配

两条独立路径:
  sT 亚基因组 (sT1-sT5)   vs Ath  ->  sT_vs_Ath.coords
  sA 亚基因组 (sA6-sA13)  vs Aar  ->  sA_vs_Aar.coords

输出 (04_collinear_pairs/), 7 列 TSV:
  asu_chr  asu_start  asu_end  ref_chr  ref_start  ref_end  identity
  正向: Asu_Ath/Asu_Aar; 反向: Ath_Asu/Aar_Asu
  {res} ∈ {5000, 10000, 25000, 100000}
"""
import os
from collections import defaultdict

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
OUTPUT_DIR = os.path.join(BASE_DIR, "02.collinearity")
ALN_DIR = os.path.join(OUTPUT_DIR, "03_mummer")
BINS_DIR = os.path.join(OUTPUT_DIR, "02_genome_bins")
PAIRS_DIR = os.path.join(OUTPUT_DIR, "04_collinear_pairs")

BIN_SIZES = [5000, 10000, 25000, 40000, 100000]
# 防御性 identity 阈值; MUMmer 1-to-1 几乎不会触发 (一般 ≥0.7)
IDENTITY_THRESHOLD = 0.5

# 亚基因组定义
ST_SUBGENOME = {'sT1', 'sT2', 'sT3', 'sT4', 'sT5'}
SA_SUBGENOME = {'sA6', 'sA7', 'sA8', 'sA9', 'sA10', 'sA11', 'sA12', 'sA13'}

# 路径配置: (coords, query 亚基因组染色体集合, query 标签, 亲本标签)
SUBGENOME_PATHS = [
    {
        'coords': os.path.join(ALN_DIR, "sT_vs_Ath.coords"),
        'asu_chroms': ST_SUBGENOME,
        'asu_label': 'Asu',
        'ref_label': 'Ath',
    },
    {
        'coords': os.path.join(ALN_DIR, "sA_vs_Aar.coords"),
        'asu_chroms': SA_SUBGENOME,
        'asu_label': 'Asu',
        'ref_label': 'Aar',
    },
]


def load_genome_bins(bin_file, keep_chroms=None):
    """加载基因组 bin; keep_chroms 不为空时只保留指定染色体"""
    bins = []
    with open(bin_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0]
                if keep_chroms is not None and chrom not in keep_chroms:
                    continue
                bins.append((chrom, int(parts[1]), int(parts[2])))
    return bins


def parse_coords(coords_file):
    """解析 show-coords -rclTH 输出.

    列: S1 E1 S2 E2 LEN1 LEN2 IDY LEN_R LEN_Q COV_R COV_Q TAG_R TAG_Q
    约定: nucmer ref=亲本, query=Asu. 故 S1/E1=ref(亲本), S2/E2=query(Asu).
    Asu 反链时 S2>E2 (mummer 用方向编码 strand). 全部转 0-based 半开区间.
    输出 dict 字段与旧 parse_paf 完全一致, 让下游函数无感知差异.
    """
    alns = []
    n_filtered = 0
    with open(coords_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 13:
                continue
            ref_s, ref_e = int(parts[0]), int(parts[1])  # 1-based inclusive, ref s<=e
            qry_s, qry_e = int(parts[2]), int(parts[3])  # 反链时 qry_s>qry_e
            len_r, len_q = int(parts[4]), int(parts[5])
            idy = float(parts[6]) / 100.0  # 0-100 → 0-1
            ref_chr = parts[11]
            asu_chr = parts[12]

            # strand: query 端起止决定方向
            strand = '+' if qry_s <= qry_e else '-'
            if strand == '+':
                q_lo, q_hi = qry_s - 1, qry_e
            else:
                q_lo, q_hi = qry_e - 1, qry_s
            # ref 端 mummer 总是 s<=e
            r_lo, r_hi = ref_s - 1, ref_e

            if idy < IDENTITY_THRESHOLD:
                n_filtered += 1
                continue

            # 用于 find_all_matches_for_bin 的强度代理 (越大越优)
            matches = int(round(idy * max(len_r, len_q)))

            alns.append({
                'asu_chr': asu_chr,
                'asu_start': q_lo,
                'asu_end': q_hi,
                'ref_chr': ref_chr,
                'ref_start': r_lo,
                'ref_end': r_hi,
                'ref_mid': (r_lo + r_hi) // 2,
                'strand': strand,
                'matches': matches,
                'identity': idy,
            })
    print(f"    identity<{IDENTITY_THRESHOLD} 过滤: {n_filtered} 条比对被移除")
    return alns


def index_alignments_by_chr(alignments):
    by_chr = defaultdict(list)
    for aln in alignments:
        by_chr[aln['asu_chr']].append(aln)
    return by_chr


def find_all_matches_for_bin(bin_chr, bin_start, bin_end, alns_by_chr):
    """为一个 Asu bin 找到所有重叠的亲本比对 (按 matches 降序)"""
    matches = []
    for aln in alns_by_chr.get(bin_chr, []):
        if aln['asu_end'] <= bin_start or aln['asu_start'] >= bin_end:
            continue
        matches.append(aln)
    matches.sort(key=lambda x: x['matches'], reverse=True)
    return matches


def find_bin_containing_position(ref_chr, ref_pos, ref_bins_by_chr):
    """找到包含目标位置的亲本 bin (二分查找)"""
    bins = ref_bins_by_chr.get(ref_chr, [])
    if not bins:
        return None
    lo, hi = 0, len(bins) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if bins[mid][1] <= ref_pos < bins[mid][2]:
            return bins[mid]
        elif ref_pos < bins[mid][1]:
            hi = mid - 1
        else:
            lo = mid + 1
    return None


def map_asu_bin_to_ref_bin(asu_bin, aln, ref_bins_by_chr):
    """将 Asu bin 按比对中的相对位置映射到亲本 bin"""
    asu_chr, asu_start, asu_end = asu_bin
    asu_mid = (asu_start + asu_end) // 2

    aln_asu_start = aln['asu_start']
    aln_asu_end = aln['asu_end']
    aln_ref_start = aln['ref_start']
    aln_ref_end = aln['ref_end']

    if aln_asu_end == aln_asu_start:
        return None

    ratio = (asu_mid - aln_asu_start) / (aln_asu_end - aln_asu_start)
    if aln['strand'] == '+':
        ref_pos = aln_ref_start + ratio * (aln_ref_end - aln_ref_start)
    else:
        ref_pos = aln_ref_end - ratio * (aln_ref_end - aln_ref_start)

    ref_pos = int(ref_pos)
    return find_bin_containing_position(aln['ref_chr'], ref_pos, ref_bins_by_chr)


def filter_by_neighbor_consistency(asu_bin_to_refs, asu_bins_by_chr, bin_size):
    """
    邻近一致性过滤 (多对多版本):
    对每个 asu_bin 的每个 ref 匹配, 检查相邻 asu_bin 是否有同一 ref_chr 的匹配。
    如果没有相邻 bin 支持该 ref_chr, 则移除该匹配。
    asu_bin_to_refs: {asu_bin: [(ref_bin, identity), ...]}
    """
    threshold = bin_size * 3
    filtered = defaultdict(list)

    for asu_chr, bins in asu_bins_by_chr.items():
        sorted_bins = sorted(bins, key=lambda x: x[1])
        for i, asu_bin in enumerate(sorted_bins):
            if asu_bin not in asu_bin_to_refs:
                continue
            refs = asu_bin_to_refs[asu_bin]

            neighbor_ref_chrs = set()
            for neighbor_idx in [i - 1, i + 1]:
                if 0 <= neighbor_idx < len(sorted_bins):
                    neighbor = sorted_bins[neighbor_idx]
                    if neighbor in asu_bin_to_refs:
                        for nb_ref in asu_bin_to_refs[neighbor]:
                            neighbor_ref_chrs.add(nb_ref[0][0])

            kept = []
            for ref_bin, identity in refs:
                if ref_bin[0] in neighbor_ref_chrs:
                    kept.append((ref_bin, identity))

            if not kept and refs:
                # 无相邻支持但保留最佳匹配 (避免完全丢失孤立 bin)
                kept = [refs[0]]

            if kept:
                filtered[asu_bin] = kept

    return filtered


def write_pairs(filtered, out_asu_ref, out_ref_asu):
    """写正向 (Asu 在前) 与反向 (亲本在前) 两个文件; 含 identity 列
    filtered: {asu_bin: [(ref_bin, identity), ...]} (多对多)"""
    # 正向
    rows_asu_ref = []
    for asu_bin, refs in filtered.items():
        for ref_bin, identity in refs:
            rows_asu_ref.append((asu_bin, ref_bin, identity))
    rows_asu_ref.sort(key=lambda x: (x[0][0], x[0][1], x[1][0], x[1][1]))

    with open(out_asu_ref, 'w') as f:
        for asu_bin, ref_bin, identity in rows_asu_ref:
            f.write(f"{asu_bin[0]}\t{asu_bin[1]}\t{asu_bin[2]}\t"
                    f"{ref_bin[0]}\t{ref_bin[1]}\t{ref_bin[2]}\t{identity:.4f}\n")

    # 反向
    rows_ref_asu = [(ref_bin, asu_bin, identity) for asu_bin, ref_bin, identity in rows_asu_ref]
    rows_ref_asu.sort(key=lambda x: (x[0][0], x[0][1], x[1][0], x[1][1]))

    with open(out_ref_asu, 'w') as f:
        for ref_bin, asu_bin, identity in rows_ref_asu:
            f.write(f"{ref_bin[0]}\t{ref_bin[1]}\t{ref_bin[2]}\t"
                    f"{asu_bin[0]}\t{asu_bin[1]}\t{asu_bin[2]}\t{identity:.4f}\n")


def process_subgenome(cfg, bin_size):
    """处理一条亚基因组路径在一个分辨率下的共线性"""
    asu_label = cfg['asu_label']
    ref_label = cfg['ref_label']
    print(f"\n  [{asu_label}({ref_label}亚基因组) vs {ref_label}] @ {bin_size}bp")

    asu_bins = load_genome_bins(
        os.path.join(BINS_DIR, f"Asu_{bin_size}_bins.bed"),
        keep_chroms=cfg['asu_chroms'])
    ref_bins = load_genome_bins(
        os.path.join(BINS_DIR, f"{ref_label}_{bin_size}_bins.bed"))
    print(f"    Asu 亚基因组 bins: {len(asu_bins)}, {ref_label} bins: {len(ref_bins)}")

    ref_bins_by_chr = defaultdict(list)
    for b in ref_bins:
        ref_bins_by_chr[b[0]].append(b)
    for chrom in ref_bins_by_chr:
        ref_bins_by_chr[chrom].sort(key=lambda x: x[1])

    asu_bins_by_chr = defaultdict(list)
    for b in asu_bins:
        asu_bins_by_chr[b[0]].append(b)

    alignments = parse_coords(cfg['coords'])
    alns_by_chr = index_alignments_by_chr(alignments)
    print(f"    coords 比对数: {len(alignments)}")

    asu_bin_to_refs = defaultdict(list)
    for asu_bin in asu_bins:
        all_matches = find_all_matches_for_bin(asu_bin[0], asu_bin[1], asu_bin[2], alns_by_chr)
        for aln in all_matches:
            ref_bin = map_asu_bin_to_ref_bin(asu_bin, aln, ref_bins_by_chr)
            if ref_bin:
                asu_bin_to_refs[asu_bin].append((ref_bin, aln['identity']))

    n_pairs = sum(len(v) for v in asu_bin_to_refs.values())
    n_bins = len(asu_bin_to_refs)
    print(f"    初步匹配: {n_bins} bins, {n_pairs} 对 (多对多)")

    filtered = filter_by_neighbor_consistency(asu_bin_to_refs, asu_bins_by_chr, bin_size)
    n_filtered_pairs = sum(len(v) for v in filtered.values())
    n_filtered_bins = len(filtered)
    print(f"    过滤后: {n_filtered_bins} bins, {n_filtered_pairs} 对")

    out_asu_ref = os.path.join(PAIRS_DIR, f"{asu_label}_{ref_label}_orth_bins_{bin_size}.tsv")
    out_ref_asu = os.path.join(PAIRS_DIR, f"{ref_label}_{asu_label}_orth_bins_{bin_size}.tsv")
    write_pairs(filtered, out_asu_ref, out_ref_asu)
    print(f"    保存: {os.path.basename(out_asu_ref)} / {os.path.basename(out_ref_asu)}")


def main():
    print("=" * 60)
    print("亚基因组感知共线性 bin 对提取 (MUMmer coords 版)")
    print("=" * 60)

    os.makedirs(PAIRS_DIR, exist_ok=True)

    for bin_size in BIN_SIZES:
        print(f"\n{'-' * 60}\n分辨率: {bin_size} bp\n{'-' * 60}")
        for cfg in SUBGENOME_PATHS:
            process_subgenome(cfg, bin_size)

    print("\n完成!")


if __name__ == "__main__":
    main()
