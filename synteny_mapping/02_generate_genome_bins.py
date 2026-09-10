#!/usr/bin/env python3
"""
将参考基因组切成多个分辨率的 bin
分辨率对应下游检测: 5kb(loops) / 10kb(TAD,compartment) / 25kb(TAD) / 100kb(compartment)
输出: 02_genome_bins/{species}_{size}_bins.bed
"""
import os

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE_DIR = ROOT
SINGLE_DIR = os.path.join(BASE_DIR, "01.single")
OUTPUT_DIR = os.path.join(BASE_DIR, "02.collinearity/02_genome_bins")

CHROM_SIZES = {
    'Ath': os.path.join(SINGLE_DIR, "Ath/ref/Ath.chrom.sizes"),
    'Aar': os.path.join(SINGLE_DIR, "Aar/ref/Aar.chrom.sizes"),
    'Asu': os.path.join(SINGLE_DIR, "Asu/ref/Asu.chrom.sizes")
}

BIN_SIZES = [5000, 10000, 25000, 40000, 100000]


def load_chrom_sizes(path):
    sizes = {}
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def main():
    print("=" * 60)
    print("生成多分辨率基因组 bin")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for species in ['Ath', 'Aar', 'Asu']:
        chrom_sizes = load_chrom_sizes(CHROM_SIZES[species])
        print(f"\n处理 {species} (染色体数: {len(chrom_sizes)})")

        for bin_size in BIN_SIZES:
            output_file = os.path.join(OUTPUT_DIR, f"{species}_{bin_size}_bins.bed")
            total_bins = 0
            with open(output_file, 'w') as f:
                for chrom, size in sorted(chrom_sizes.items()):
                    start = 0
                    while start < size:
                        end = min(start + bin_size, size)
                        f.write(f"{chrom}\t{start}\t{end}\n")
                        start += bin_size
                        total_bins += 1
            print(f"  {bin_size:>6} bp -> {total_bins:>7} bins  ({output_file})")

    print("\n完成!")


if __name__ == "__main__":
    main()
