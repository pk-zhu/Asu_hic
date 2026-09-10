#!/usr/bin/env python3
"""
07d_trajectory_enrichment.py — 四类 3D 轨迹的基因功能富集 (sT / sA / ALL).

背景
----
07c 只对 Transient_Perturbation (TP) 做了富集。要判断 "某功能是否 TP 特有",
必须同时看另外三类轨迹 (Parental_Retention / Immediate_Remodeling /
Gradual_Remodeling), 否则无法排除该功能在全部轨迹中都富集的可能。本脚本
补齐四类轨迹, 供跨轨迹对照。

统计单元与检验方法与 07c 完全一致 (超几何上尾 + BH FDR), 直接复用其函数,
不重算注释。

输出
----
results/trajectory/traj_enrich/
  traj_enrich_{sT,sA,ALL}_{PR,IR,TP,GR}_GO.tsv     每轨迹富集结果
  traj_enrich_stress_FE.tsv                        胁迫相关 GO 的跨轨迹 FE 对照
  traj_enrich_summary.txt                          文字汇总

输入 (只读)
-----------
与 07c 相同: trajectory_sT_Ath.tsv / trajectory_sA_Aar.tsv / emapper / obo / gff
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module

_m = import_module("07c_tp_gene_enrichment")
parse_obo = _m.parse_obo
parse_emapper = _m.parse_emapper
parse_gff_genes = _m.parse_gff_genes
enrich = _m.enrich

# Project root. Override with the ARABIDOPSIS_HIC_ROOT environment variable;
# otherwise infer it from this script's location. This resolves correctly both
# in the original <root>/<module>/scripts layout and in the released
# code_release/<module> layout.
ROOT = os.environ.get(
    "ARABIDOPSIS_HIC_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
BASE = ROOT
TRJ_DIR = os.path.join(BASE, "03.pairwise_subgenome_conservation/results/trajectory")
EMAPPER = os.path.join(
    BASE, "05.function_annotation/results/emapper_out/Asu/Asu_emapper.emapper.annotations")
OBO = os.path.join(BASE, "05.function_annotation/resources/go-basic.obo")
GFF = os.path.join(BASE, "01.single/Asu/ref/Asu.gff")
OUT_DIR = os.path.join(TRJ_DIR, "traj_enrich")

TRAJ = {
    "PR": "Parental_Retention",
    "IR": "Immediate_Remodeling",
    "TP": "Transient_Perturbation",
    "GR": "Gradual_Remodeling",
}

# 与寒冷 / 干燥 / 渗透 / 氧化 / 解毒 相关的 GO, 用于跨轨迹 FE 对照。
# 生物胁迫条目一并列出, 以便区分 "非生物" 与 "生物" 两条主题。
STRESS_GO = {
    # 非生物: 温度
    "GO:0009409": "response to cold",
    "GO:0009266": "response to temperature stimulus",
    "GO:0009628": "response to abiotic stimulus",
    # 非生物: 水分 / 渗透 / 盐
    "GO:0009414": "response to water deprivation",
    "GO:0009415": "response to water",
    "GO:0006970": "response to osmotic stress",
    "GO:0009651": "response to salt stress",
    "GO:0009737": "response to abscisic acid",
    # 氧化 / 解毒 (脱水与低温的下游共通胁迫)
    "GO:0006979": "response to oxidative stress",
    "GO:0098754": "detoxification",
    "GO:0009636": "response to toxic substance",
    "GO:0006749": "glutathione metabolic process",
    # 生物胁迫 (对照)
    "GO:0006952": "defense response",
    "GO:0009607": "response to biotic stimulus",
    "GO:0050832": "defense response to fungus",
    "GO:0006032": "chitin catabolic process",
    "GO:0004568": "chitinase activity",
    # 细胞周期 (IR 主题对照)
    "GO:0007049": "cell cycle",
    "GO:0051301": "cell division",
    "GO:0051321": "meiotic cell cycle",
}


def log(msg):
    print(msg, flush=True)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    log("=== 07d_trajectory_enrichment ===")

    st = pd.read_csv(f"{TRJ_DIR}/trajectory_sT_Ath.tsv", sep="\t")
    sa = pd.read_csv(f"{TRJ_DIR}/trajectory_sA_Aar.tsv", sep="\t")
    all_traj = pd.concat([st, sa], ignore_index=True)

    log("  载入注释与坐标 ...")
    go_names, go_obsolete, alt2canon = parse_obo(OBO)
    gene2go, _kegg, _pfam, _desc, _name = parse_emapper(
        EMAPPER, alt2canon, go_obsolete)
    coord = parse_gff_genes(GFF)
    log(f"  emapper loci: {len(gene2go):,}  GFF coords: {len(coord):,}")

    fe_rows = []
    counts = []
    for label, dfg in [("sT", st), ("sA", sa), ("ALL", all_traj)]:
        bg_all = set(dfg.gene_id)
        universe = {g for g in bg_all if gene2go.get(g)}
        log(f"\n  === {label} === background={len(bg_all):,} "
            f"with-GO={len(universe):,}")
        for code, full in TRAJ.items():
            fg = set(dfg[dfg.trajectory == full].gene_id) & universe
            if not fg:
                log(f"  [{code}] 前景为空, 跳过")
                continue
            df = enrich(fg, universe, gene2go, coord, go_names)
            out = f"{OUT_DIR}/traj_enrich_{label}_{code}_GO.tsv"
            df.to_csv(out, sep="\t", index=False)
            n_sig = int((df.fdr_bh < 0.05).sum())
            frac = len(fg) / len(universe)
            log(f"  [{code}] fg={len(fg):,} ({frac:.1%} of universe) "
                f"terms={len(df):,} FDR<0.05={n_sig}")
            counts.append({"subgenome": label, "trajectory": code,
                           "n_fg": len(fg), "n_universe": len(universe),
                           "fg_fraction": round(frac, 4), "n_sig": n_sig})
            for _, r in df.head(8).iterrows():
                if r.fdr_bh >= 0.05:
                    break
                log(f"    {r.term:<12} {str(r.term_name)[:34]:<36} "
                    f"FE={r.fold_enrichment:<6} FDR={r.fdr_bh:.2e}")

            idx = df.set_index("term")
            for go, gname in STRESS_GO.items():
                if go not in idx.index:
                    continue
                r = idx.loc[go]
                fe_rows.append({
                    "subgenome": label, "trajectory": code,
                    "term": go, "term_name": gname,
                    "n_bg": int(r.n_bg), "n_fg": int(r.n_fg),
                    "fold_enrichment": r.fold_enrichment,
                    "p_hyper": r.p_hyper, "fdr_bh": r.fdr_bh,
                    "sig": bool(r.fdr_bh < 0.05),
                })

    pd.DataFrame(counts).to_csv(
        f"{OUT_DIR}/traj_enrich_counts.tsv", sep="\t", index=False)
    fe = pd.DataFrame(fe_rows)
    fe.to_csv(f"{OUT_DIR}/traj_enrich_stress_FE.tsv", sep="\t", index=False)
    log(f"\nsaved: {OUT_DIR}/traj_enrich_stress_FE.tsv ({len(fe)} 行)")

    with open(f"{OUT_DIR}/traj_enrich_summary.txt", "w", encoding="utf-8") as fh:
        w = fh.write
        w("四类 3D 轨迹的 GO 富集: sT / sA / ALL\n")
        w("=" * 74 + "\n\n")
        w("[目的]\n")
        w("  07c 只做 TP。要判断某功能是否 TP 特有, 必须同时看 PR / IR / GR,\n")
        w("  排除该功能在全部轨迹都富集的可能。本脚本补齐四类轨迹。\n\n")
        w("[方法]\n")
        w("  超几何上尾 + BH FDR, 与 07c 完全一致 (直接复用其 enrich 函数)。\n")
        w("  前景 = 该亚基因组该轨迹的基因 (有 >= 1 条 GO 注释者);\n")
        w("  背景 = 该亚基因组全部可分类基因 (有 >= 1 条 GO 注释者)。\n\n")
        w("[前景占背景比例与显著条目数]\n")
        for r in counts:
            w(f"  {r['subgenome']:<4} {r['trajectory']:<3} "
              f"fg={r['n_fg']:>6} ({r['fg_fraction']:.1%}) "
              f"FDR<0.05={r['n_sig']}\n")
        w("\n[读法警示]\n")
        w("  PR 前景占背景比例极高 (>60%), 前景近似等于背景, FE 必然逼近 1,\n")
        w("  无富集可测, 这是数学必然而非生物学结论。GR 前景极小, 统计力低,\n")
        w("  FE 虽可极端但多不显著。故跨轨迹比较应以 FE 的方向与量级为主,\n")
        w("  不可仅凭 '是否显著' 判断功能归属。\n\n")
        w("[胁迫相关 GO 的跨轨迹 FE 对照]\n")
        w("  见 traj_enrich_stress_FE.tsv。\n")
        if len(fe):
            for go, gname in STRESS_GO.items():
                sub = fe[(fe.term == go) & (fe.subgenome == "ALL")]
                if not len(sub):
                    continue
                w(f"\n  {go} {gname}\n")
                for _, r in sub.iterrows():
                    star = " *" if r.sig else ""
                    w(f"    {r['trajectory']:<3} FE={r['fold_enrichment']:<7} "
                      f"fg={r['n_fg']:>4}/{r['n_bg']:<5} "
                      f"FDR={r['fdr_bh']:.2e}{star}\n")
    log(f"saved: {OUT_DIR}/traj_enrich_summary.txt")
    log("=== done ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())