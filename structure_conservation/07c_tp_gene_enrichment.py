#!/usr/bin/env python3
"""
07c_tp_gene_enrichment.py — sT 与 sA 的 TP 基因分别做功能富集.

问题
----
07_trajectory_classification.py 把 sT 与 sA 的基因各自分为四类轨迹
(Parental_Retention / Immediate_Remodeling / Transient_Perturbation /
Gradual_Remodeling)。本脚本对两个亚基因组 **分别** 做功能富集:

    sT: 前景 = sT 被判为 Transient_Perturbation 的基因 (主分析 6041)
        背景 = sT 全部可分类基因 (27298)
    sA: 前景 = sA 被判为 Transient_Perturbation 的基因 (主分析 4002)
        背景 = sA 全部可分类基因 (30034)

两个亚基因组各自独立, 互不借用, 也不涉及同源对。

统计单元: 基因。背景 = 该亚基因组全部基因 (有 >= 1 条 GO 注释者)。

检验: 超几何分布上尾 (hypergeom.sf), Benjamini-Hochberg FDR 控制。
      N = 背景基因数, M = 前景基因数,
      n = 背景中注释某条目的基因数, K = 前景中注释该条目的基因数,
      p = P(X >= K) = hypergeom.sf(K-1, N, n, M)

GO / KEGG_Pathway / PFAM 三套注释体系分别检验, 各自独立做 BH 校正。

串联簇与位置审计
----------------
超几何检验假设各观测单元独立, 串联复制簇违反此假设 (如 chitinase 在 sT2
的 17.58-17.60 Mb 的一串基因会夸大独立观测数)。谨慎的串联簇折叠需亚基因组
内部平行同源数据, 仓库无现成结果, 本脚本不擅自重算 mmseqs。故对每个富集
条目做位置审计: 列出前景基因坐标, 并按 sT/sA 侧 100 kb 窗口统计。窗口数
小 (如 1-2) 说明富集由单个串联基因家族驱动, 应表述为该基因家族在 TP 中
富集, 而非广泛的通路级富集。

GO alt_id 处理
--------------
eggNOG 注释可能直接给出 GO 的 alt_id (如 GO:0017144, 规范条目为
GO:0006805)。脚本把 alt_id 映射回规范 id, 并剔除 obsolete 条目。

输入 (只读)
-----------
* trajectory_sT_Ath.tsv, trajectory_sA_Aar.tsv
    03.pairwise_subgenome_conservation/results/trajectory/
* Asu_emapper.emapper.annotations
    05.function_annotation/results/emapper_out/Asu/
* go-basic.obo
    05.function_annotation/resources/
* Asu.gff
    01.single/Asu/ref/   (仅用于取基因坐标做位置审计)

输出
----
results/trajectory/tp_gene_enrich/
  tp_enrich_{sT,sA}_GO.tsv / _KEGG.tsv / _PFAM.tsv   富集结果
  tp_genes_{sT,sA}.tsv                                该亚基因组 TP 基因明细 (坐标+注释)
  tp_enrich_summary.txt                               文字汇总
"""
from __future__ import annotations

import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

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
OUT_DIR = os.path.join(TRJ_DIR, "tp_gene_enrich")

TP = "Transient_Perturbation"
CLUSTER_WIN = 100_000
MIN_BG = 3
MAX_LIST = 60

NS_CODE = {"biological_process": "BP",
           "molecular_function": "MF",
           "cellular_component": "CC"}


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------- #
# 注释解析
# --------------------------------------------------------------------------- #
def parse_obo(path):
    """返回 (names, obsolete, alt2canon)."""
    names, obsolete, alt2canon = {}, set(), {}
    if not os.path.isfile(path):
        return names, obsolete, alt2canon
    with open(path, encoding="utf-8") as fh:
        in_term = False
        gid = name = ns = ""
        alts = []
        is_obs = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("[") and line.endswith("]"):
                if in_term:
                    if is_obs:
                        obsolete.add(gid)
                        obsolete.update(alts)
                    else:
                        names[gid] = (name, ns)
                        for a in alts:
                            alt2canon[a] = gid
                in_term = (line == "[Term]")
                gid = name = ns = ""
                alts = []
                is_obs = False
                continue
            if not in_term:
                continue
            if line.startswith("id:"):
                gid = line.split("id:", 1)[1].strip()
            elif line.startswith("name:"):
                name = line.split("name:", 1)[1].strip()
            elif line.startswith("namespace:"):
                ns = line.split("namespace:", 1)[1].strip()
            elif line.startswith("alt_id:"):
                alts.append(line.split("alt_id:", 1)[1].strip())
            elif line.startswith("is_obsolete:") and "true" in line.lower():
                is_obs = True
        if in_term:
            if is_obs:
                obsolete.add(gid)
                obsolete.update(alts)
            else:
                names[gid] = (name, ns)
                for a in alts:
                    alt2canon[a] = gid
    return names, obsolete, alt2canon


def canonicalize_go(gos, alt2canon, obsolete):
    out = set()
    for g in gos:
        if g in obsolete:
            continue
        c = alt2canon.get(g, g)
        if c in obsolete or c == "":
            continue
        out.add(c)
    return out


def parse_emapper(path, alt2canon, obsolete):
    """返回 (gene2go, gene2kegg, gene2pfam, gene2desc, gene2name)."""
    header = {}
    gene2go = defaultdict(set)
    gene2kegg = defaultdict(set)
    gene2pfam = defaultdict(set)
    gene2desc, gene2name = {}, {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                fields = line.lstrip("#").rstrip("\n").split("\t")
                header = {f: i for i, f in enumerate(fields)}
                continue
            if not header:
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < max(header.values()) + 1:
                continue
            gene = cols[header["query"]].split(".")[0]

            def multi(col):
                raw = cols[header[col]]
                if raw in ("-", ""):
                    return set()
                return {x.strip() for x in raw.split(",") if x.strip()}

            gene2go[gene] |= canonicalize_go(
                {g for g in multi("GOs") if g.startswith("GO:")},
                alt2canon, obsolete)
            gene2kegg[gene] |= {x for x in multi("KEGG_Pathway")
                                if x.startswith("ko")}
            gene2pfam[gene] |= multi("PFAMs")
            desc = cols[header["Description"]]
            if desc not in ("-", "") and gene not in gene2desc:
                gene2desc[gene] = desc
            pname = cols[header["Preferred_name"]]
            if pname not in ("-", "") and gene not in gene2name:
                gene2name[gene] = pname
    return dict(gene2go), dict(gene2kegg), dict(gene2pfam), gene2desc, gene2name


def parse_gff_genes(path):
    """返回 {gene_id: (chrom, start, end)}, 仅取 feature == gene."""
    coord = {}
    if not os.path.isfile(path):
        return coord
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            gid = None
            for kv in c[8].split(";"):
                if kv.startswith("ID="):
                    gid = kv[3:].split(".")[0]
                    break
            if gid:
                coord[gid] = (c[0], int(c[3]), int(c[4]))
    return coord


# --------------------------------------------------------------------------- #
# BH FDR
# --------------------------------------------------------------------------- #
def bh_fdr(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    running = 1.0
    for k, idx in enumerate(reversed(order), start=1):
        rank = n - k + 1
        running = min(running, pvals[idx] * n / rank)
        adj[idx] = running
    return [min(v, 1.0) for v in adj]


# --------------------------------------------------------------------------- #
# 富集
# --------------------------------------------------------------------------- #
def enrich(fg, universe, gene2term, coord, term_label=None,
           min_bg=MIN_BG):
    """超几何富集, 每个条目的前景基因附位置审计 (100 kb 簇标签)."""
    N = len(universe)
    M = len(fg)
    term2bg = defaultdict(set)
    for g in universe:
        for t in gene2term.get(g, ()):
            term2bg[t].add(g)

    rows = []
    for t, bg_genes in term2bg.items():
        n = len(bg_genes)
        if n < min_bg:
            continue
        hit = bg_genes & fg
        K = len(hit)
        if K == 0:
            continue
        expected = n * M / N
        p = float(hypergeom.sf(K - 1, N, n, M))
        name, ns = term_label.get(t, ("", "")) if term_label else ("", "")

        clusters = defaultdict(int)
        audit = []
        for g in sorted(hit):
            chrom, start = coord.get(g, ("-", 0))[:2]
            if chrom != "-":
                win = f"{chrom}:{int(start) // CLUSTER_WIN * CLUSTER_WIN // 1_000_000}Mb"
                clusters[win] += 1
                audit.append(f"{g}@{chrom}:{start // 1_000}kb")
            else:
                clusters["-"] += 1
                audit.append(f"{g}@-")
        top_win = dict(sorted(clusters.items(), key=lambda kv: -kv[1])[:5])

        rows.append({
            "term": t,
            "term_name": name,
            "namespace": NS_CODE.get(ns, ns),
            "n_bg": n,
            "n_fg": K,
            "expected": round(expected, 3),
            "fold_enrichment": round(K / expected, 3) if expected > 0 else np.nan,
            "p_hyper": p,
            "n_100kb_windows": len(clusters),
            "distribution": ";".join(f"{w}:{c}" for w, c in top_win.items()),
            "fg_genes": ",".join(audit[:MAX_LIST]) +
                        (",..." if len(audit) > MAX_LIST else ""),
        })

    cols = ["term", "term_name", "namespace", "n_bg", "n_fg", "expected",
            "fold_enrichment", "p_hyper", "fdr_bh", "n_100kb_windows",
            "distribution", "fg_genes"]
    if not rows:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(rows)
    df["fdr_bh"] = bh_fdr(df["p_hyper"].tolist())
    df = df.sort_values(["fdr_bh", "p_hyper"],
                        ascending=True).reset_index(drop=True)
    return df[cols]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    log("=== 07c_tp_gene_enrichment ===")

    st = pd.read_csv(f"{TRJ_DIR}/trajectory_sT_Ath.tsv", sep="\t")
    sa = pd.read_csv(f"{TRJ_DIR}/trajectory_sA_Aar.tsv", sep="\t")
    log(f"  sT genes: {len(st):,}  TP: {(st.trajectory == TP).sum():,}")
    log(f"  sA genes: {len(sa):,}  TP: {(sa.trajectory == TP).sum():,}")

    log("  载入注释与坐标 ...")
    go_names, go_obsolete, alt2canon = parse_obo(OBO)
    gene2go, gene2kegg, gene2pfam, gene2desc, gene2name = parse_emapper(
        EMAPPER, alt2canon, go_obsolete)
    coord = parse_gff_genes(GFF)
    log(f"  obo: {len(go_names):,}  obsolete: {len(go_obsolete):,}"
        f"  alt2canon: {len(alt2canon):,}")
    log(f"  emapper loci: {len(gene2go):,}  GFF coords: {len(coord):,}")

    # 每亚基因组独立富集 + 亚基因组合并 (ALL) 富集
    all_traj = pd.concat([st, sa], ignore_index=True)
    results = {}
    for label, dfg in [("sT", st), ("sA", sa), ("ALL", all_traj)]:
        fg_all = set(dfg[dfg.trajectory == TP].gene_id)
        bg_all = set(dfg.gene_id)
        log(f"\n  === {label} ===  foreground(TP)={len(fg_all):,}  "
            f"background={len(bg_all):,}")

        # 明细表
        tp = dfg[dfg.trajectory == TP].copy()
        tp["chrom"] = tp.gene_id.map(lambda g: coord.get(g, ("-", 0, 0))[0])
        tp["start"] = tp.gene_id.map(lambda g: coord.get(g, ("-", 0, 0))[1])
        tp["end"] = tp.gene_id.map(lambda g: coord.get(g, ("-", 0, 0))[2])
        tp["preferred_name"] = tp.gene_id.map(gene2name).fillna("-")
        tp["description"] = tp.gene_id.map(gene2desc).fillna("-")
        tp["GOs"] = tp.gene_id.map(
            lambda g: ",".join(sorted(gene2go.get(g, []))) or "-")
        tp["KEGG"] = tp.gene_id.map(
            lambda g: ",".join(sorted(gene2kegg.get(g, []))) or "-")
        tp["PFAMs"] = tp.gene_id.map(
            lambda g: ",".join(sorted(gene2pfam.get(g, []))) or "-")
        (tp[["gene_id", "nat_corr", "syn_corr", "chrom", "start", "end",
             "preferred_name", "description", "GOs", "KEGG", "PFAMs"]]
         .sort_values(["chrom", "start"])
         .to_csv(f"{OUT_DIR}/tp_genes_{label}.tsv", sep="\t", index=False))
        log(f"  写出 TP 基因明细: {len(tp):,}")

        for tag, mapping, labels in [("GO", gene2go, go_names),
                                     ("KEGG", gene2kegg, None),
                                     ("PFAM", gene2pfam, None)]:
            universe = {g for g in bg_all if mapping.get(g)}
            fg = fg_all & universe
            if not fg or not universe:
                log(f"  [{tag}] 无可检验, 跳过")
                continue
            df = enrich(fg, universe, mapping, coord, labels)
            df.to_csv(f"{OUT_DIR}/tp_enrich_{label}_{tag}.tsv",
                      sep="\t", index=False)
            n_sig = int((df.fdr_bh < 0.05).sum())
            log(f"  [{tag}] universe={len(universe):,} fg={len(fg):,} "
                f"terms={len(df):,} FDR<0.05={n_sig}")
            results[(label, tag)] = (df, len(universe), len(fg), n_sig)
            for _, r in df.head(14).iterrows():
                if r.fdr_bh >= 0.05:
                    break
                log(f"    {r.term:<12} {str(r.term_name)[:34]:<36} "
                    f"fg={r.n_fg:>3}/{r.n_bg:<5} FE={r.fold_enrichment:<6} "
                    f"FDR={r.fdr_bh:.2e} 簇={r.n_100kb_windows} "
                    f"[{r.distribution[:44]}]")

    # 文字汇总
    with open(f"{OUT_DIR}/tp_enrich_summary.txt", "w", encoding="utf-8") as fh:
        w = fh.write
        w("sT 与 sA 的 Transient_Perturbation 基因: 分别功能富集\n")
        w("=" * 78 + "\n\n")
        w("[定义]\n")
        w("  sT 与 sA 各自独立判为 TP 的基因, 分别做富集, 互不借用。\n")
        w("  前景 = 该亚基因组 TP 基因; 背景 = 该亚基因组全部可分类基因。\n\n")
        w("[计数]\n")
        w(f"  sT: 全部 {len(st):,}, TP {(st.trajectory == TP).sum():,}\n")
        w(f"  sA: 全部 {len(sa):,}, TP {(sa.trajectory == TP).sum():,}\n\n")
        w("[方法]\n")
        w("  超几何上尾检验 + Benjamini-Hochberg FDR。\n")
        w(f"  背景中注释基因数 < {MIN_BG} 的条目不检验。\n")
        w("  obsolete GO 已剔除; alt_id 映射回规范主条目。\n")
        w("  GO / KEGG / PFAM 各自独立做 BH 校正。\n\n")
        w("[串联簇位置审计]\n")
        w("  超几何假设单元独立, 串联复制簇违反之。谨慎折叠需亚基因组内部\n")
        w(f"  平行同源, 仓库无现成数据, 不重算。按 {CLUSTER_WIN // 1_000:,} kb\n")
        w("  窗口审计前景基因分布: 窗口数小说明富集由单个串联基因家族驱动,\n")
        w("  应表述为该家族在 TP 中富集, 而非通路级富集。\n\n")
        for label in ["sT", "sA"]:
            w(f"[{label}]\n")
            for tag in ["GO", "KEGG", "PFAM"]:
                if (label, tag) not in results:
                    w(f"  {tag}: 无结果\n")
                    continue
                df, nu, nf, n_sig = results[(label, tag)]
                w(f"  {tag}: universe={nu:,} fg={nf:,} terms={len(df):,} "
                  f"FDR<0.05={n_sig}\n")
                sig = df[df.fdr_bh < 0.05]
                if len(sig):
                    for _, r in sig.head(30).iterrows():
                        w(f"    {r.term:<12} [{r.namespace or '-':<2}] "
                          f"{str(r.term_name)[:44]:<46} "
                          f"fg={r.n_fg:>3}/{r.n_bg:<5} FE={r.fold_enrichment:<6} "
                          f"FDR={r.fdr_bh:.1e} 簇={r.n_100kb_windows}\n")
                        w(f"      {CLUSTER_WIN // 1_000}kb 窗口: {r.distribution}\n")
                        w(f"      前景基因: "
                          f"{r.fg_genes[:140]}{'...' if len(r.fg_genes) > 140 else ''}\n")
                else:
                    w(f"    无 FDR<0.05。p 最小 5 条:\n")
                    for _, r in df.head(5).iterrows():
                        w(f"      {r.term:<12} {str(r.term_name)[:44]:<46} "
                          f"{r.n_fg:>3}/{r.n_bg:<5} FE={r.fold_enrichment:<6} "
                          f"p={r.p_hyper:.1e} FDR={r.fdr_bh:.2f}\n")
            w("\n")

    log(f"\n  输出目录: {OUT_DIR}")
    log("=== done ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())