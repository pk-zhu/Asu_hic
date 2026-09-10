"""亚基因组→亲本染色体全局映射 (基于 PAF matches 总量, 不受分辨率影响)

从 sT_vs_Ath.paf 和 sA_vs_Aar.paf 的 matches 总量 greedy 分配:
  sT1->Chr1, sT2->Chr2, sT3->Chr3, sT4->Chr4, sT5->Chr5
  sA6->chr_1, sA7->chr_2, sA8->chr_3, sA9->chr_4,
  sA10->chr_5, sA11->chr_6, sA12->chr_7, sA13->chr_8
"""

# sT 亚基因组 -> Ath
ST_TO_ATH = {
    'sT1': 'Chr1',
    'sT2': 'Chr2',
    'sT3': 'Chr3',
    'sT4': 'Chr4',
    'sT5': 'Chr5',
}

# sA 亚基因组 -> Aar
SA_TO_AAR = {
    'sA6': 'chr_1',
    'sA7': 'chr_2',
    'sA8': 'chr_3',
    'sA9': 'chr_4',
    'sA10': 'chr_5',
    'sA11': 'chr_6',
    'sA12': 'chr_7',
    'sA13': 'chr_8',
}


def get_chr_mapping(filepath):
    """根据共线性文件路径返回对应的染色体映射表"""
    if 'Asu_Ath' in filepath:
        return ST_TO_ATH
    elif 'Asu_Aar' in filepath:
        return SA_TO_AAR
    else:
        raise ValueError(f"无法识别共线性文件类型: {filepath}")
