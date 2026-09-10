#!/usr/bin/env python3
"""
06_syn_asu_conservation.py — syn_Asu (合成 Asu) vs Asu/Ath/Aar 三维结构保守性.

三条比较:
  syn_Asu_vs_Asu: 直接 bin 对齐 (同基因组同坐标), 量化合成 vs 天然差异
  syn_Asu_vs_Ath: 复用 Asu_Ath 共线性 bin 对, mcool 换为 syn_Asu
  syn_Asu_vs_Aar: 复用 Asu_Aar 共线性 bin 对, mcool 换为 syn_Asu

输出: results/matrix/  results/compartment/  results/tad/
"""
import os, sys, math, gzip
from collections import defaultdict
import numpy as np
import pandas as pd
import cooler

sys.path.insert(0, os.path.dirname(__file__))
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
OUTPUT_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results")

MCOOL_PATHS = {
    "Ath": os.path.join(SINGLE_DIR, "Ath/Matrix/Final/Ath.mcool"),
    "Aar": os.path.join(SINGLE_DIR, "Aar/Matrix/Final/Aar.mcool"),
    "Asu": os.path.join(SINGLE_DIR, "Asu/Matrix/Final/Asu.mcool"),
    "syn_Asu": os.path.join(SINGLE_DIR, "syn_Asu/Matrix/Final/syn_Asu.mcool"),
}

RESOLUTIONS = [10000, 25000, 100000]
MAX_DIST_MB = 10.0
MIN_DIST_BINS = 1
TAD_RESOLUTIONS = [5000, 10000, 25000, 40000]


def log(msg):
    print(msg, flush=True)


# ==================== SCC (复用 01_bin_matrix_conservation.py 核心) ====================

def _pearson_from_stats(n, sx, sy, sxx, syy, sxy):
    if n < 2: return float("nan")
    vx = sxx - (sx * sx) / n
    vy = syy - (sy * sy) / n
    if vx <= 0 or vy <= 0: return float("nan")
    return (sxy - (sx * sy) / n) / math.sqrt(vx * vy)


def _accumulate_stats(stats, d, x, y):
    s = stats[d]
    s[0] += 1; s[1] += x; s[2] += y
    s[3] += x * x; s[4] += y * y; s[5] += x * y


def scc_from_aligned_matrices(mat_a, mat_b, max_dist_bins, min_dist_bins):
    stats = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])
    dict_b = {}
    for i, j, v in zip(mat_b.row, mat_b.col, mat_b.data):
        if j < i: i, j = j, i
        d = j - i
        if d < min_dist_bins or d > max_dist_bins: continue
        if not np.isfinite(v): continue
        dict_b[(i, j)] = float(v)
    for i, j, v in zip(mat_a.row, mat_a.col, mat_a.data):
        if j < i: i, j = j, i
        d = j - i
        if d < min_dist_bins or d > max_dist_bins: continue
        if not np.isfinite(v): continue
        key = (i, j)
        if key in dict_b:
            _accumulate_stats(stats, d, float(v), dict_b[key])
    strata = []
    total_weight = 0
    weighted_sum = 0
    for d in sorted(stats.keys()):
        n, sx, sy, sxx, syy, sxy = stats[d]
        r = _pearson_from_stats(n, sx, sy, sxx, syy, sxy)
        if not np.isnan(r):
            strata.append((d, r, n))
            weighted_sum += r * n
            total_weight += n
    scc = weighted_sum / total_weight if total_weight > 0 else float("nan")
    return scc, strata


def compute_scc_for_comparison(label, cool_a, cool_b, resolution, collinear_pairs):
    """计算一对比较的 SCC."""
    res_str = str(resolution)
    max_dist_bins = int(MAX_DIST_MB * 1e6 / resolution)

    chr_results = []
    all_strata = []

    for chrom in cool_a.chromnames:
        # For identity pairs, both sides use same chromosome names
        # For collinear pairs, the mapping is in the collinear pairs themselves
        # Check if there are any collinear pairs for this chromosome
        chr_pairs = [p for p in collinear_pairs if p[0] == chrom]
        if len(chr_pairs) < 10: continue
        # Get the reference chromosome from the first pair
        ref_chrom = chr_pairs[0][3]
        # Get bins for this chromosome
        chr_bins_a = cool_a.bins()[:]
        chr_bins_a = chr_bins_a[chr_bins_a["chrom"] == chrom]
        chr_bins_b = cool_b.bins()[:]
        chr_bins_b = chr_bins_b[chr_bins_b["chrom"] == ref_chrom]

        if len(chr_bins_a) == 0 or len(chr_bins_b) == 0: continue

        # Build bin index mapping from collinear pairs
        bin_map = {}
        for pair in collinear_pairs:
            asu_chr, asu_start, asu_end, ref_chr, ref_start, ref_end = pair
            if asu_chr != chrom or ref_chr != ref_chrom: continue
            # Find local bin indices within this chromosome
            a_idx = chr_bins_a.index[(chr_bins_a["start"] == asu_start) & (chr_bins_a["end"] == asu_end)].tolist()
            b_idx = chr_bins_b.index[(chr_bins_b["start"] == ref_start) & (chr_bins_b["end"] == ref_end)].tolist()
            if len(a_idx) > 0 and len(b_idx) > 0:
                # Convert to local indices (0-based within chromosome)
                a_local = a_idx[0] - chr_bins_a.index[0]
                b_local = b_idx[0] - chr_bins_b.index[0]
                bin_map[a_local] = b_local

        if len(bin_map) < 10: continue

        # Get matrices as dense arrays for alignment
        mat_a = cool_a.matrix(balance=True, sparse=True).fetch(chrom).toarray()
        mat_b = cool_b.matrix(balance=True, sparse=True).fetch(ref_chrom).toarray()

        # Build aligned matrices using numpy fancy indexing
        a_idx = np.array(sorted(bin_map.keys()), dtype=int)
        b_idx = np.array([bin_map[k] for k in sorted(bin_map.keys())], dtype=int)
        n = len(a_idx)

        # Create meshgrid for upper triangle within max_dist
        rows, cols = np.triu_indices(n, k=1)
        dists = cols - rows
        valid = dists <= max_dist_bins
        rows, cols = rows[valid], cols[valid]

        if len(rows) < 100: continue

        # Extract values using fancy indexing
        a_i = np.minimum(a_idx[rows], a_idx[cols])
        a_j = np.maximum(a_idx[rows], a_idx[cols])
        b_i = np.minimum(b_idx[rows], b_idx[cols])
        b_j = np.maximum(b_idx[rows], b_idx[cols])

        va = mat_a[a_i, a_j]
        vb = mat_b[b_i, b_j]
        d_arr = a_j - a_i  # distance (in bins) shared with parent, per collinear pair
        mask = np.isfinite(va) & np.isfinite(vb)
        if mask.sum() < 100: continue

        # Stratum-adjusted SCC: per-distance Pearson r weighted by n, matching
        # 01_bin_matrix_conservation.py logic.
        from collections import defaultdict
        stats_by_d = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])
        for i_idx in np.where(mask)[0]:
            d = int(d_arr[i_idx])
            x = float(va[i_idx]); y = float(vb[i_idx])
            s = stats_by_d[d]
            s[0] += 1
            s[1] += x; s[2] += y
            s[3] += x * x; s[4] += y * y; s[5] += x * y

        strata = []
        weighted_sum = 0.0
        total_weight = 0
        for d in sorted(stats_by_d.keys()):
            n, sx, sy, sxx, syy, sxy = stats_by_d[d]
            if n < 2:
                continue
            mean_x = sx / n; mean_y = sy / n
            var_x = sxx / n - mean_x * mean_x
            var_y = syy / n - mean_y * mean_y
            cov_xy = sxy / n - mean_x * mean_y
            if var_x <= 0 or var_y <= 0:
                continue
            r = cov_xy / np.sqrt(var_x * var_y)
            if np.isnan(r):
                continue
            strata.append((d, r, n))
            weighted_sum += r * n
            total_weight += n

        if total_weight == 0:
            continue
        scc_val = weighted_sum / total_weight
        if not np.isnan(scc_val):
            chr_results.append((chrom, scc_val, len(bin_map)))
            for d, r, w in strata:
                all_strata.append((d, r, w, chrom))

    if chr_results:
        total_w = sum(c[2] for c in chr_results)
        overall = sum(c[1] * c[2] for c in chr_results) / total_w if total_w > 0 else float("nan")
        log(f"  {label} {res_str}: SCC={overall:.4f} ({len(chr_results)} chrs)")
        return overall, chr_results, all_strata
    return float("nan"), [], []


def generate_identity_pairs(cool, resolution):
    """为 syn_Asu vs Asu 生成 identity bin 对 (同染色体同坐标)."""
    pairs = []
    bins_df = cool.bins()[:]
    for chrom in cool.chromnames:
        chr_bins = bins_df[bins_df["chrom"] == chrom]
        for _, row in chr_bins.iterrows():
            pairs.append((chrom, int(row["start"]), int(row["end"]),
                         chrom, int(row["start"]), int(row["end"])))
    log(f"  Identity pairs: {len(pairs):,} bins across {len(cool.chromnames)} chrs")
    return pairs


def load_collinear_pairs(filepath):
    """加载共线性 bin 对."""
    pairs = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                pairs.append((parts[0], int(parts[1]), int(parts[2]),
                             parts[3], int(parts[4]), int(parts[5])))
    # Filter by chr_mapping
    chr_map = get_chr_mapping(filepath)
    filtered = [p for p in pairs if p[3] == chr_map.get(p[0])]
    log(f"  Collinear pairs: {len(pairs):,} -> {len(filtered):,} (chr filtered)")
    return filtered


def compare_compartment(label, asu_comp_file, ref_comp_file, resolution, collinear_pairs):
    """比较 compartment 保守性 (复用 02_compartment_conservation.py 逻辑).
    使用 majority vote: 对每个 100kb bin, 找所有重叠的 10kb sub-compartment, 多数 A/B 为 bin 状态."""
    def load_sub_compartments(path):
        """加载 CALDER2 sub_compartments.tsv, 返回 {chr: [(start, end, A/B)]}."""
        df = pd.read_csv(path, sep="\t")
        comps = defaultdict(list)
        for _, row in df.iterrows():
            chrom = str(row.iloc[0])
            # CALDER2 chr prefix: chrsT1->sT1, chr__1->chr_1
            if chrom.startswith("chr") and len(chrom) > 3 and chrom[3:].startswith("s"):
                chrom = chrom[3:]
            elif chrom.startswith("chr_") and "__" in chrom:
                chrom = chrom.replace("__", "_")
            comp_name = str(row.iloc[3])
            state = 'A' if comp_name.startswith('A') else 'B'
            comps[chrom].append((int(row.iloc[1]), int(row.iloc[2]), state))
        return comps

    def get_bin_state(bin_chr, bin_start, bin_end, comp_dict):
        """Majority vote: 找所有与 bin 重叠的 sub-compartment, 多数 A/B."""
        if bin_chr not in comp_dict:
            return None
        overlapping = []
        for cs, ce, state in comp_dict[bin_chr]:
            if cs < bin_end and ce > bin_start:
                overlapping.append(state)
        if not overlapping:
            return None
        return 'A' if overlapping.count('A') >= overlapping.count('B') else 'B'

    asu_comps = load_sub_compartments(asu_comp_file)
    ref_comps = load_sub_compartments(ref_comp_file)

    conserved = a_to_b = b_to_a = unmapped = 0
    for pair in collinear_pairs:
        a_state = get_bin_state(pair[0], pair[1], pair[2], asu_comps)
        r_state = get_bin_state(pair[3], pair[4], pair[5], ref_comps)
        if a_state is None or r_state is None:
            unmapped += 1; continue
        if a_state == r_state: conserved += 1
        elif a_state == 'A' and r_state == 'B': a_to_b += 1
        elif a_state == 'B' and r_state == 'A': b_to_a += 1
        else: unmapped += 1

    total = conserved + a_to_b + b_to_a
    if total > 0:
        log(f"  {label} compartment: conserved={conserved}/{total} ({conserved/total*100:.1f}%), "
            f"A->B={a_to_b}, B->A={b_to_a}, unmapped={unmapped}")

    # Also write per-bin compartment data
    bin_rows = []
    for pair in collinear_pairs:
        a_state = get_bin_state(pair[0], pair[1], pair[2], asu_comps)
        r_state = get_bin_state(pair[3], pair[4], pair[5], ref_comps)
        if a_state is None or r_state is None: continue
        cons_type = "conserved" if a_state == r_state else ("A_to_B" if a_state == 'A' else "B_to_A")
        bin_rows.append(dict(asu_chr=pair[0], asu_start=pair[1], asu_end=pair[2],
                            ref_chr=pair[3], ref_start=pair[4], ref_end=pair[5],
                            asu_compartment=a_state, ref_compartment=r_state,
                            conservation_type=cons_type))
    if bin_rows:
        pd.DataFrame(bin_rows).to_csv(
            os.path.join(OUTPUT_DIR, "compartment", f"compartment_conservation_{label}.tsv"),
            sep="\t", index=False)

    return dict(comparison=label, resolution=resolution, conserved=conserved,
                a_to_b=a_to_b, b_to_a=b_to_a, unmapped=unmapped,
                total=total, cons_rate=round(conserved/total, 4) if total > 0 else 0)


def compare_tad(label, asu_tad_file, ref_tad_file, resolution, collinear_pairs,
                subgenome_filter=None):
    """per-TAD 保守性: 为每个 Asu (或 syn_Asu) TAD 找亲本中重叠 bin 最多的 TAD.
    与 natural pipeline (03_tad_conservation.py) 完全一致.
    subgenome_filter: None (全部 TAD) 或 set of chrom names (只保留匹配亚基因组的 TAD).
    """
    def load_tads_list(path):
        tads = []
        if not os.path.exists(path):
            return tads
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    tads.append({
                        "chr": parts[0], "start": int(parts[1]),
                        "end": int(parts[2]), "name": parts[3],
                    })
        return tads

    asu_tads = load_tads_list(asu_tad_file)
    ref_tads = load_tads_list(ref_tad_file)
    if subgenome_filter is not None:
        asu_tads = [t for t in asu_tads if t["chr"] in subgenome_filter]
    if not asu_tads or not ref_tads:
        log(f"  {label} TAD: skip (missing bed or empty subgenome)")
        return None

    # Asu bin -> ref bin map
    asu_to_ref = defaultdict(list)
    for pair in collinear_pairs:
        asu_chr, asu_start, asu_end, ref_chr, ref_start, ref_end = pair[:6]
        asu_to_ref[(asu_chr, asu_start, asu_end)].append((ref_chr, ref_start, ref_end))

    ref_tads_by_chr = defaultdict(list)
    for t in ref_tads:
        ref_tads_by_chr[t["chr"]].append(t)

    scores = []
    n_unmapped = 0
    for atad in asu_tads:
        asu_chr, asu_s, asu_e = atad["chr"], atad["start"], atad["end"]
        overlapping_bins = []
        target_bins = []
        for (bin_chr, bin_s, bin_e), tgt_list in asu_to_ref.items():
            if bin_chr == asu_chr and not (bin_e <= asu_s or bin_s >= asu_e):
                overlapping_bins.append((bin_chr, bin_s, bin_e))
                target_bins.extend(tgt_list)
        if not overlapping_bins:
            n_unmapped += 1
            continue

        # majority target chrom
        target_chr_count = defaultdict(int)
        for tb in target_bins:
            target_chr_count[tb[0]] += 1
        if not target_chr_count:
            n_unmapped += 1
            continue
        best_ref_chr = max(target_chr_count.items(), key=lambda x: x[1])[0]

        tads_here = ref_tads_by_chr.get(best_ref_chr, [])
        max_overlap = 0
        for rt in tads_here:
            bc = 0
            for tb in target_bins:
                if tb[0] == rt["chr"] and not (tb[2] <= rt["start"] or tb[1] >= rt["end"]):
                    bc += 1
            if bc > max_overlap:
                max_overlap = bc

        if max_overlap > 0:
            scores.append(min(1.0, max_overlap / len(overlapping_bins)))

    if not scores:
        log(f"  {label} TAD: 0 scored / {len(asu_tads)} TADs")
        return None

    mean_score = float(np.mean(scores))
    high = sum(1 for s in scores if s >= 0.8)
    total = len(asu_tads)
    high_rate = high / len(scores)
    log(f"  {label} TAD: mean_score={mean_score:.3f}, "
        f"high_cons={high}/{len(scores)} ({high_rate*100:.1f}%), "
        f"n_total_tads={total}, unmapped={n_unmapped}")
    return dict(comparison=label, resolution=resolution,
                n_total_tads=total, n_tads=len(scores),
                n_unmapped=n_unmapped,
                mean_score=round(mean_score, 4),
                n_high=high,
                high_rate=round(high_rate, 4))


def main():
    log("=== 06_syn_asu_conservation ===")

    SUBGENOME = {
        "Ath": {"sT1", "sT2", "sT3", "sT4", "sT5"},
        "Aar": {"sA6", "sA7", "sA8", "sA9", "sA10", "sA11", "sA12", "sA13"},
        "Asu": None,  # compare all 13 chromosomes for syn_Asu vs Asu
    }

    # Three comparisons
    comps = [
        {"label": "syn_Asu_vs_Asu", "asu_sp": "syn_Asu", "ref_sp": "Asu",
         "collinear": "identity",  # direct bin alignment
         "asu_comp_res": 100000, "ref_comp_res": 100000,
         "tad_resolutions": TAD_RESOLUTIONS},
        {"label": "syn_Asu_vs_Ath", "asu_sp": "syn_Asu", "ref_sp": "Ath",
         "collinear": "Asu_Ath",
         "asu_comp_res": 100000, "ref_comp_res": 100000,
         "tad_resolutions": TAD_RESOLUTIONS},
        {"label": "syn_Asu_vs_Aar", "asu_sp": "syn_Asu", "ref_sp": "Aar",
         "collinear": "Asu_Aar",
         "asu_comp_res": 100000, "ref_comp_res": 100000,
         "tad_resolutions": TAD_RESOLUTIONS},
    ]

    # Output directories
    for sub in ["matrix", "compartment", "tad"]:
        os.makedirs(os.path.join(OUTPUT_DIR, sub), exist_ok=True)

    scc_summary = []
    comp_summary = []
    tad_summary = []

    for comp in comps:
        label = comp["label"]
        asu_sp = comp["asu_sp"]
        ref_sp = comp["ref_sp"]
        log(f"\n{'='*60}")
        log(f"  {label}")
        log(f"{'='*60}")

        # Load mcool
        asu_cool_uri = MCOOL_PATHS[asu_sp] + f"::/resolutions/10000"
        ref_cool_uri = MCOOL_PATHS[ref_sp] + f"::/resolutions/10000"

        asu_cool = cooler.Cooler(asu_cool_uri)
        ref_cool = cooler.Cooler(ref_cool_uri)

        # Get collinear pairs
        if comp["collinear"] == "identity":
            pairs = generate_identity_pairs(asu_cool, 10000)
        else:
            pair_file = os.path.join(COLLINEARITY_DIR, f"{comp['collinear']}_orth_bins_10000.tsv")
            if not os.path.exists(pair_file):
                log(f"  SKIP: collinear file not found: {pair_file}")
                continue
            pairs = load_collinear_pairs(pair_file)

        # SCC
        for res in RESOLUTIONS:
            asu_c = cooler.Cooler(MCOOL_PATHS[asu_sp] + f"::/resolutions/{res}")
            ref_c = cooler.Cooler(MCOOL_PATHS[ref_sp] + f"::/resolutions/{res}")

            # For identity comparison, regenerate pairs at correct resolution
            if comp["collinear"] == "identity":
                res_pairs = generate_identity_pairs(asu_c, res)
            else:
                res_pair_file = os.path.join(COLLINEARITY_DIR, f"{comp['collinear']}_orth_bins_{res}.tsv")
                res_pairs = load_collinear_pairs(res_pair_file)

            scc_val, chr_results, strata = compute_scc_for_comparison(label, asu_c, ref_c, res, res_pairs)
            if not np.isnan(scc_val):
                scc_summary.append(dict(comparison=label, resolution=res, SCC=round(scc_val, 4)))
                # Write strata
                if strata:
                    pd.DataFrame(strata, columns=["dist_bin", "pearson_r", "weight", "chrom"]).to_csv(
                        os.path.join(OUTPUT_DIR, "matrix", f"scc_{res}_{label}_strata.tsv"), sep="\t", index=False)

        # Compartment (use 100kb collinear pairs)
        if comp["collinear"] == "identity":
            comp_pairs = generate_identity_pairs(cooler.Cooler(MCOOL_PATHS[asu_sp] + "::/resolutions/100000"), 100000)
        else:
            comp_pair_file = os.path.join(COLLINEARITY_DIR, f"{comp['collinear']}_orth_bins_100000.tsv")
            comp_pairs = load_collinear_pairs(comp_pair_file) if os.path.exists(comp_pair_file) else []

        asu_comp_file = os.path.join(SINGLE_DIR, asu_sp,
            f"Compartments/CALDER2/{comp['asu_comp_res']}/{asu_sp}/sub_compartments/all_sub_compartments.tsv")
        ref_comp_file = os.path.join(SINGLE_DIR, ref_sp,
            f"Compartments/CALDER2/{comp['ref_comp_res']}/{ref_sp}/sub_compartments/all_sub_compartments.tsv")

        if os.path.exists(asu_comp_file) and os.path.exists(ref_comp_file) and comp_pairs:
            comp_res = compare_compartment(label, asu_comp_file, ref_comp_file,
                                          comp['asu_comp_res'], comp_pairs)
            if comp_res: comp_summary.append(comp_res)
        else:
            log(f"  SKIP compartment: missing files")

        # TAD (multi-resolution sweep)
        sub_filter = SUBGENOME.get(ref_sp)  # sT for Ath, sA for Aar, None for Asu
        for tad_res in comp["tad_resolutions"]:
            asu_tad_file = os.path.join(SINGLE_DIR, asu_sp,
                f"TADs/HiCExplorer/{tad_res}/{asu_sp}_domains.bed")
            ref_tad_file = os.path.join(SINGLE_DIR, ref_sp,
                f"TADs/HiCExplorer/{tad_res}/{ref_sp}_domains.bed")

            if os.path.exists(asu_tad_file) and os.path.exists(ref_tad_file):
                # Reuse pairs from 10k for TAD resolution >10k; for <10k need exact-res pairs
                # Since compare_tad only uses bin (chr,start,end) coords, the 10k pairs still work as long as bins are aligned at higher resolution boundaries
                tad_res_pair = pairs
                if comp["collinear"] != "identity":
                    res_pair_file = os.path.join(COLLINEARITY_DIR, f"{comp['collinear']}_orth_bins_{tad_res}.tsv")
                    if os.path.exists(res_pair_file):
                        tad_res_pair = load_collinear_pairs(res_pair_file)
                res_out = compare_tad(label, asu_tad_file, ref_tad_file,
                                      tad_res, tad_res_pair,
                                      subgenome_filter=sub_filter)
                if res_out: tad_summary.append(res_out)
            else:
                log(f"  SKIP TAD ({tad_res}bp): missing files")

    # Write summaries
    if scc_summary:
        pd.DataFrame(scc_summary).to_csv(os.path.join(OUTPUT_DIR, "matrix", "scc_summary_syn_asu.tsv"), sep="\t", index=False)
    if comp_summary:
        pd.DataFrame(comp_summary).to_csv(os.path.join(OUTPUT_DIR, "compartment", "compartment_conservation_syn_asu.tsv"), sep="\t", index=False)
    if tad_summary:
        pd.DataFrame(tad_summary).to_csv(os.path.join(OUTPUT_DIR, "tad", "tad_conservation_syn_asu.tsv"), sep="\t", index=False)

    log(f"\n=== DONE ===")
    log(f"  SCC: {len(scc_summary)} rows")
    log(f"  Compartment: {len(comp_summary)} rows")
    log(f"  TAD: {len(tad_summary)} rows")


if __name__ == "__main__":
    main()