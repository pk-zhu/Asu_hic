#!/usr/bin/env python3
"""
03. 基因区为中心的 3D 结构保守性 (multi-resolution)

对落在共线区的每个 Asu 基因 (sT->Ath, sA->Aar):
  1. 借共线性把 Asu 坐标线性插值映射到亲本坐标
  2. matrix 层: 取基因区 ±WINDOW 的 balanced 子矩阵, Asu vs 亲本对齐算 Pearson
     支持分辨率: 2k (用5k orth_bins映射), 10k
  3. TAD 层: 基因是否落在各自 TAD 内 (支持 5k, 10k, 25k)

输出 -> results/gene_3d_syn_Asu/{matrix_res}/
  gene_3d_conservation_syn_Asu_{Ath,Aar}.tsv   per-gene
  gene_3d_summary_syn_Asu.tsv                       sT vs sA 汇总
"""
import os
import sys
import numpy as np
import cooler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chr_mapping import ST_TO_ATH, SA_TO_AAR

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
COLLINEAR_DIR = os.path.join(BASE_DIR, "02.collinearity/04_collinear_pairs")
SINGLE_DIR = os.path.join(BASE_DIR, "01.single")
CDS_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/cds")
OUTPUT_BASE = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/gene_3d_syn_Asu")

# 支持的配置: (matrix_resolution, orth_bins_resolution)
MATRIX_CONFIGS = [
    (5000, 5000),   # 5k matrix, 5k orth_bins for mapping
    (10000, 10000), # 10k matrix, 10k orth_bins
    (25000, 25000), # 25k matrix, 25k orth_bins
]

TAD_RESOLUTIONS = [5000, 10000, 25000]
WINDOW = 50000          # 基因区两侧各扩展
MIN_SUBMATRIX = 4       # 子矩阵最小边长
COMP_RES = 100000


def load_collinear_map(fname, chr_map):
    """asu_chr -> sorted [(asu_mid, ref_mid)]"""
    raw = {}
    with open(os.path.join(COLLINEAR_DIR, fname)) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            ac, as_, ae = p[0], int(p[1]), int(p[2])
            rc, rs, re = p[3], int(p[4]), int(p[5])
            if chr_map.get(ac) != rc:
                continue
            raw.setdefault(ac, []).append(((as_ + ae) // 2, (rs + re) // 2))
    for ac in raw:
        raw[ac].sort()
    return raw


def map_coord(asu_chr, pos, cmap):
    arr = cmap.get(asu_chr)
    if not arr:
        return None
    amids = [a for a, _ in arr]
    rmids = [r for _, r in arr]
    n = len(arr)
    lo, hi = 0, n - 1
    if pos <= amids[0]:
        i0, i1 = 0, min(1, n - 1)
    elif pos >= amids[-1]:
        i0, i1 = max(0, n - 2), n - 1
    else:
        import bisect
        j = bisect.bisect_left(amids, pos)
        i0, i1 = j - 1, j
    a0, a1 = amids[i0], amids[i1]
    r0, r1 = rmids[i0], rmids[i1]
    if a1 == a0:
        return r0
    frac = (pos - a0) / (a1 - a0)
    return int(r0 + frac * (r1 - r0))


def load_genes(species, chrs):
    genes = []
    with open(os.path.join(CDS_DIR, f"{species}_genes.bed")) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] in chrs:
                genes.append((p[0], int(p[1]), int(p[2]), p[3]))
    return genes


def load_tads(species, resolution):
    tads = {}
    bed = os.path.join(SINGLE_DIR, species, "TADs/HiCExplorer", str(resolution),
                       f"{species}_domains.bed")
    if not os.path.exists(bed):
        return tads
    with open(bed) as f:
        for line in f:
            p = line.split("\t")
            tads.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    for c in tads:
        tads[c].sort()
    return tads


def load_compartments(comp):
    comp_path = os.path.join(
        BASE_DIR, "03.pairwise_subgenome_conservation/results/compartment",
        f"compartment_conservation_{comp}.tsv"
    )
    bins = {}
    with open(comp_path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ac, s, e, ctype = p[0], int(p[1]), int(p[2]), p[8]
            bins.setdefault(ac, []).append((s, e, ctype))
    for c in bins:
        bins[c].sort()
    return bins


def get_compartment(comp_bins, chrom, pos):
    arr = comp_bins.get(chrom, [])
    if not arr:
        return "NA"
    import bisect
    starts = [s for s, _, _ in arr]
    j = bisect.bisect_right(starts, pos) - 1
    if j < 0:
        return "NA"
    s, e, ctype = arr[j]
    if s <= pos < e:
        return ctype
    return "NA"


def in_tad(tads, chrom, pos):
    for s, e in tads.get(chrom, []):
        if s <= pos < e:
            return True
        if s > pos:
            break
    return False


def cache_chrom_matrices(clr, chroms):
    cache = {}
    for c in chroms:
        try:
            m = clr.matrix(balance=True).fetch(c)
            cache[c] = np.nan_to_num(m, nan=0.0)
        except Exception:
            cache[c] = None
    return cache


def submatrix(mat, start, end, resolution):
    if mat is None:
        return None
    n = mat.shape[0]
    i0 = max(0, (start - WINDOW) // resolution)
    i1 = min(n, (end + WINDOW) // resolution + 1)
    if i1 - i0 < MIN_SUBMATRIX:
        return None
    return mat[i0:i1, i0:i1]


def matrix_corr(asu_sub, par_sub):
    if asu_sub is None or par_sub is None:
        return np.nan
    k = min(asu_sub.shape[0], par_sub.shape[0])
    if k < MIN_SUBMATRIX:
        return np.nan
    a = asu_sub[:k, :k].flatten()
    b = par_sub[:k, :k].flatten()
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    r = np.corrcoef(a, b)[0, 1]
    return r if not np.isnan(r) else np.nan


def run_matrix_resolution(matrix_res, orth_res):
    print(f"\n{'='*70}")
    print(f"Matrix 分辨率: {matrix_res}bp (orth_bins: {orth_res}bp)")
    print(f"{'='*70}")

    comparisons = [
        ("syn_Asu_vs_Ath", f"Asu_Ath_orth_bins_{orth_res}.tsv", "Ath", ST_TO_ATH),
        ("syn_Asu_vs_Aar", f"Asu_Aar_orth_bins_{orth_res}.tsv", "Aar", SA_TO_AAR),
    ]

    output_dir = os.path.join(OUTPUT_BASE, str(matrix_res))
    os.makedirs(output_dir, exist_ok=True)

    summary = []
    for comp, fname, par_sp, cmap_def in comparisons:
        print(f"\n=== {comp} ===")
        cmap = load_collinear_map(fname, cmap_def)
        asu_chrs = set(cmap_def.keys())
        par_chrs = set(cmap_def.values())

        genes = load_genes("Asu", asu_chrs)
        print(f"  Asu 基因: {len(genes)}")

        # Load TADs for all resolutions
        asu_tads_all = {}
        par_tads_all = {}
        for tad_res in TAD_RESOLUTIONS:
            asu_tads_all[tad_res] = load_tads("Asu", tad_res)
            par_tads_all[tad_res] = load_tads(par_sp, tad_res)

        comp_bins = load_compartments(comp)
        print(f"  compartment bins: {sum(len(v) for v in comp_bins.values())}")

        print(f"  加载 Hi-C 矩阵 ({matrix_res}bp)...")
        asu_clr = cooler.Cooler(f"{SINGLE_DIR}/syn_Asu/Matrix/Final/syn_Asu.mcool::resolutions/{matrix_res}")
        par_clr = cooler.Cooler(f"{SINGLE_DIR}/{par_sp}/Matrix/Final/{par_sp}.mcool::resolutions/{matrix_res}")
        asu_mats = cache_chrom_matrices(asu_clr, asu_chrs)
        par_mats = cache_chrom_matrices(par_clr, par_chrs)

        rows = []
        n_done = 0
        for chrom, gs, ge, gid in genes:
            par_chr = cmap_def[chrom]
            gmid = (gs + ge) // 2
            pmid = map_coord(chrom, gmid, cmap)
            if pmid is None:
                continue
            ps = map_coord(chrom, gs, cmap)
            pe = map_coord(chrom, ge, cmap)
            if ps is None or pe is None:
                continue
            if ps > pe:
                ps, pe = pe, ps

            asu_sub = submatrix(asu_mats.get(chrom), gs, ge, matrix_res)
            par_sub = submatrix(par_mats.get(par_chr), ps, pe, matrix_res)
            corr = matrix_corr(asu_sub, par_sub)

            # TAD conservation for each resolution
            tad_results = {}
            for tad_res in TAD_RESOLUTIONS:
                asu_in = in_tad(asu_tads_all[tad_res], chrom, gmid)
                par_in = in_tad(par_tads_all[tad_res], par_chr, pmid)
                tad_results[f"tad_conserved_{tad_res}"] = int(asu_in and par_in)
                tad_results[f"asu_in_tad_{tad_res}"] = int(asu_in)
                tad_results[f"par_in_tad_{tad_res}"] = int(par_in)

            ctype = get_compartment(comp_bins, chrom, gmid)

            row = {
                "gene_id": gid, "asu_chr": chrom, "asu_start": gs, "asu_end": ge,
                "par_chr": par_chr, "par_start": ps, "par_end": pe,
                "matrix_corr": round(corr, 6) if not np.isnan(corr) else "NA",
                "compartment_type": ctype,
            }
            row.update(tad_results)
            rows.append(row)
            n_done += 1
            if n_done % 5000 == 0:
                print(f"    ...{n_done} genes")

        # Write per-gene
        import csv
        out = os.path.join(output_dir, f"gene_3d_conservation_{comp}.tsv")
        with open(out, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"  -> {out}  ({len(rows)} genes)")

        # Summary
        corrs = [r["matrix_corr"] for r in rows if r["matrix_corr"] != "NA"]
        corrs = np.array(corrs, dtype=float)
        comp_types = [r["compartment_type"] for r in rows if r["compartment_type"] != "NA"]
        comp_conserved_rate = sum(1 for c in comp_types if c == "conserved") / len(comp_types) if comp_types else 0.0
        sub_label = "sT" if "Ath" in comp else "sA"

        sum_row = {
            "comparison": comp, "subgenome": sub_label, "parent": par_sp,
            "n_genes": len(rows), "n_matrix_valid": len(corrs),
            "mean_matrix_corr": round(float(np.mean(corrs)), 6) if len(corrs) else "NA",
            "median_matrix_corr": round(float(np.median(corrs)), 6) if len(corrs) else "NA",
            "compartment_conserved_rate": round(float(comp_conserved_rate), 6),
        }
        for tad_res in TAD_RESOLUTIONS:
            tad_rate = np.mean([r[f"tad_conserved_{tad_res}"] for r in rows]) if rows else 0.0
            sum_row[f"tad_conserved_rate_{tad_res}"] = round(float(tad_rate), 6)

        summary.append(sum_row)
        print(f"  mean matrix_corr={sum_row['mean_matrix_corr']}  "
              f"TAD保守率(5k)={sum_row.get('tad_conserved_rate_5000','NA'):.4f}  "
              f"TAD保守率(10k)={sum_row.get('tad_conserved_rate_10000','NA'):.4f}  "
              f"TAD保守率(25k)={sum_row.get('tad_conserved_rate_25000','NA'):.4f}")

    import csv
    sout = os.path.join(output_dir, "gene_3d_summary_syn_Asu.tsv")
    with open(sout, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    print(f"\n汇总 -> {sout}")


def main():
    for matrix_res, orth_res in MATRIX_CONFIGS:
        run_matrix_resolution(matrix_res, orth_res)
    print("\n03 done.")


if __name__ == "__main__":
    main()
