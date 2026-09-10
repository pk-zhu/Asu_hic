#!/usr/bin/env python3
"""
03c_tad_threshold_sensitivity.py — TAD conservation_score 阈值敏感性分析.

对每个分辨率 (10k, 25k, 40k) 和每个保守阈值 sweep:
  阈值 ∈ [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

输出每 (resolution, threshold) 的 high_conservation 率 (4 个比较):
  - Asu_Ath (天然 sT vs Ath)
  - Asu_Aar (天然 sA vs Aar)
  - syn_Asu_vs_Ath (合成 sT vs Ath)
  - syn_Asu_vs_Aar (合成 sA vs Aar)

用于 fig3 panel F2 (阈值敏感性, syn_Asu 用虚线).

输出: results/tad/tad_threshold_sensitivity.tsv
"""
import os
import pandas as pd

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
TAD_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results/tad")
RESOLUTIONS = [10000, 25000, 40000]
THRESHOLDS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# (comparison key in TSV, filename pattern)
COMPARISONS = [
    ("Asu_Ath",         "Asu_Ath"),
    ("Asu_Aar",         "Asu_Aar"),
    ("syn_Asu_vs_Ath",  "syn_Asu_vs_Ath"),
    ("syn_Asu_vs_Aar",  "syn_Asu_vs_Aar"),
]


def main():
    rows = []
    for res in RESOLUTIONS:
        for comp_key, fname in COMPARISONS:
            tsv = os.path.join(TAD_DIR, str(res), f'tad_conservation_{fname}.tsv')
            if not os.path.exists(tsv):
                continue
            df = pd.read_csv(tsv, sep='\t')
            n_total = len(df)
            for thr in THRESHOLDS:
                n_high = int((df['conservation_score'] >= thr).sum())
                rate = n_high / n_total * 100 if n_total > 0 else 0
                rows.append({
                    'resolution': res,
                    'comparison': comp_key,
                    'threshold': thr,
                    'n_total': n_total,
                    'n_high': n_high,
                    'high_rate(%)': round(rate, 2),
                })

    df_out = pd.DataFrame(rows)
    out_path = os.path.join(TAD_DIR, 'tad_threshold_sensitivity.tsv')
    df_out.to_csv(out_path, sep='\t', index=False)
    print(f"保存: {out_path}\n")

    # 紧凑展示: 每分辨率 4 比较
    print("=== high_rate(%) per comparison per resolution ===")
    for res in RESOLUTIONS:
        print(f"\n  {res//1000} kb:")
        for thr in THRESHOLDS:
            line = [f"  thr={thr:.1f}"]
            for comp_key, _ in COMPARISONS:
                row = df_out[(df_out.resolution == res)
                             & (df_out.comparison == comp_key)
                             & (df_out.threshold == thr)]
                v = row['high_rate(%)'].values[0] if len(row) else float('nan')
                line.append(f"{comp_key}={v:5.1f}")
            print(' | '.join(line))


if __name__ == "__main__":
    main()