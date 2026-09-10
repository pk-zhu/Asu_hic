#!/usr/bin/env python3
"""
04d. WINDOW 敏感性扫描 (matrix 层 window ∈ {25, 50, 75, 100, 200} kb, matrix_res=10kb).

设计:
  - 复用 03_gene_region_3d_conservation.py 与 04_syn_asu_gene_3d.py 的
    run_matrix_resolution(matrix_res, orth_res), 但在每次调用前通过
    importlib 重载模块并 patch 全局 WINDOW 与 OUTPUT_BASE, 实现窗口可参数化.
  - 矩阵分辨率固定 10 kb (10 kb orth_bins), 与正文 baseline 一致;
    变更的只有邻域窗口宽度 (WINDOW ∈ {25000, 50000, 75000, 100000, 200000}).
  - 不动现有 results/gene_3d/{2000,5000,10000,25000}/ 文件 (默认 50kb 窗口).
    敏感性输出位于 results/gene_3d_window_sensitivity/{window}/.

样本规模:
  - WINDOW=25 kb, matrix_res=10 kb, 基因长典型 2 至 4 kb:
    子矩阵边长 ≈ (2+50)/10 ≈ 5.2 bins (>=MIN_SUBMATRIX=4, 通过)
  - WINDOW≥50 kb: 子矩阵自然更宽, N 接近 100% 总基因数
  - WINDOW 越窄, NA 过滤后 N 越少, 需在最终表格里附带 (n_valid) 列.
"""
import os
import sys
import importlib.util
import argparse
import time

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
SCRIPTS = f"{BASE}/04.cds_synteny_3d_conservation/scripts"
OUT_BASE = f"{BASE}/04.cds_synteny_3d_conservation/results/gene_3d_window_sensitivity"

WINDOWS = [25000, 50000, 75000, 100000, 200000]
MATRIX_RES = 10000
ORTH_RES = 10000

NAT_SCRIPT = f"{SCRIPTS}/03_gene_region_3d_conservation.py"
SYN_SCRIPT = f"{SCRIPTS}/04_syn_asu_gene_3d.py"


def run_for_window(module_path: str, label: str, window: int):
    spec = importlib.util.spec_from_file_location(
        f"gene3d_{label}_{window}", module_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out_dir = os.path.join(OUT_BASE, f"win{window}")
    mod.OUTPUT_BASE = out_dir
    mod.WINDOW = window
    # 默认 WINDOW=50000 时 MIN_SUBMATRIX=4. 25kb 时 N 显著下降,
    # 我们仍保留 MIN_SUBMATRIX=4 不变, 让"门槛以上"的样本自然进入.
    mod.MIN_SUBMATRIX = 4
    print(f"\n[{label}] window={window}bp → out={out_dir}")
    t0 = time.time()
    mod.run_matrix_resolution(MATRIX_RES, ORTH_RES)
    print(f"  elapsed: {time.time()-t0:.1f}s")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--script", choices=["nat", "syn", "both"], default="both")
    args = parser.parse_args()

    for window in WINDOWS:
        if args.script in ("nat", "both"):
            run_for_window(NAT_SCRIPT, "nat", window)
        if args.script in ("syn", "both"):
            run_for_window(SYN_SCRIPT, "syn", window)
    print("\n04d done.")


if __name__ == "__main__":
    main()
