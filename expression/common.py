#!/usr/bin/env python3
"""
09.expression_3d_integration 公共变量与小工具.

提供:
  - BASE_DIR, 各子模块路径
  - load_tpm(species): 读 TPM 矩阵, 返回 (gene_id 列表, sample 列表, ndarray)
  - tpm_summary(): 每个基因的平均 / 中位 / log1p mean, 用于下游分析
  - normalize_gene_id(): 把 emapper query 的后缀剥掉, 对齐到 08 矩阵
"""
import os
import re
import numpy as np
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
EXPR_DIR  = os.path.join(BASE_DIR, "08.expression/04.count/matrix")
GENE3D_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/gene_3d")
GENES_BED_DIR = os.path.join(BASE_DIR, "04.cds_synteny_3d_conservation/results/cds")
PWCONS_DIR = os.path.join(BASE_DIR, "03.pairwise_subgenome_conservation/results")
EMAPPER_DIR = os.path.join(BASE_DIR, "05.function_annotation/results/emapper_out")
TE_DENSITY_DIR = os.path.join(BASE_DIR, "07.subgenome_TE_analysis/results/density")
COLL_DIR = os.path.join(BASE_DIR, "02.collinearity/04_collinear_pairs")
OUT_BASE = os.path.join(BASE_DIR, "09.expression_3d_integration/results")

# 染色体 -> 亚基因组
ST_CHROMS = {f"sT{i}" for i in range(1, 6)}
SA_CHROMS = {f"sA{i}" for i in range(6, 14)}


def load_tpm(species):
    """读 08 TPM 矩阵; 返回 DataFrame, index=gene_id, columns=sample."""
    path = os.path.join(EXPR_DIR, f"{species}.tpm.tsv")
    df = pd.read_csv(path, sep="\t").rename(columns={"Gene_Id": "gene_id"})
    df = df.set_index("gene_id")
    return df


TPM_MIN = 0.01


def tpm_summary(df):
    """每基因表达概要: mean, median, log_mean (log1p 后再求 mean), expressed (mean>=TPM_MIN)."""
    out = pd.DataFrame(index=df.index)
    out["mean_tpm"]   = df.mean(axis=1)
    out["median_tpm"] = df.median(axis=1)
    out["log_mean"]   = np.log1p(df.values).mean(axis=1)
    out["expressed"]  = (out["mean_tpm"] >= TPM_MIN).astype(int)
    return out


def normalize_gene_id(qid, species):
    """emapper 的 #query 字段还原成 08 矩阵里的 gene_id.
    Asu:  Asu01G000001.mRNA1  ->  Asu01G000001
    Ath:  AT1G01010.1          ->  AT1G01010
    Aar:  rna-Aar_LOCUS1       ->  Aar_LOCUS1
    """
    if species == "Asu":
        return re.sub(r"\.mRNA\d+$", "", qid)
    if species == "Ath":
        return re.sub(r"\.\d+$", "", qid)
    if species == "Aar":
        return qid.removeprefix("rna-") if qid.startswith("rna-") else qid
    return qid


def load_emapper_orthologs(species, want_species_prefix=None):
    """解析 emapper.orthologs (4 列: query orth_type species orthologs).
    返回 long-format DataFrame: query_id(已 normalize), orth_type, ortho_species, ortho_id

    want_species_prefix: 字符串或可迭代; 若给定, 只保留 ortho_species 以其开头的行,
        显著节省内存 (Asu emapper 完整展开 ~6M 行, 过滤后 <1M).
    """
    path = os.path.join(EMAPPER_DIR, species, f"{species}_emapper.emapper.orthologs")
    if want_species_prefix is None:
        prefixes = None
    elif isinstance(want_species_prefix, str):
        prefixes = (want_species_prefix,)
    else:
        prefixes = tuple(want_species_prefix)

    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            q, otype, osp, olist = parts[:4]
            if prefixes is not None and not osp.startswith(prefixes):
                continue
            q_norm = normalize_gene_id(q, species)
            # olist 形如 "*XP_010522873,*XP_010556656"
            for o in olist.split(","):
                rows.append((q_norm, otype, osp, o.lstrip("*")))
    return pd.DataFrame(rows, columns=["query_id", "orth_type", "ortho_species", "ortho_id"])


def load_genes_bed(species):
    """04/cds/{species}_genes.bed -> DataFrame chr/start/end/gene_id/strand"""
    path = os.path.join(GENES_BED_DIR, f"{species}_genes.bed")
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chr", "start", "end", "gene_id", "strand"])
    return df


def asu_subgenome(asu_chr):
    """sT1-sT5 -> 'sT'; sA6-sA13 -> 'sA'"""
    if asu_chr in ST_CHROMS: return "sT"
    if asu_chr in SA_CHROMS: return "sA"
    return None
