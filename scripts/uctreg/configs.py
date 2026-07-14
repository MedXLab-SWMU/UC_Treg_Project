from __future__ import annotations

from dataclasses import dataclass


DATASET_ORDER = [
    "1_GSE182270",
    "2_GSE214695",
    "3_GSE150115",
    "4_GSE116222",
    "5_GSE134649",
    "6_GSE125527_tissue",
    "6_GSE125527_pbmc",
    "7_GSE231993",
    "8_SCP259",
    "9_GSE114374",
    "10_GSE296969",
    "11_GSE235663",
]


DATASET_SPECS = {
    "1_GSE182270": {"loader": "10x", "raw": "1_GSE182270"},
    "2_GSE214695": {"loader": "10x", "raw": "2_GSE214695"},
    "3_GSE150115": {"loader": "matrix_t", "raw": "3_GSE150115"},
    "4_GSE116222": {"loader": "gse116222", "raw": "4_GSE116222"},
    "5_GSE134649": {"loader": "10x", "raw": "5_GSE134649"},
    "6_GSE125527_tissue": {
        "loader": "matrix",
        "raw": "6_GSE125527/Tissue",
    },
    "6_GSE125527_pbmc": {
        "loader": "matrix",
        "raw": "6_GSE125527/PBMC",
    },
    "7_GSE231993": {"loader": "10x", "raw": "7_GSE231993"},
    "8_SCP259": {"loader": "scp259", "raw": "8_SCP259/SCP259_combined.h5ad"},
    "9_GSE114374": {"loader": "gse114374", "raw": "9_GSE114374"},
    "10_GSE296969": {"loader": "10x", "raw": "10_GSE296969"},
    "11_GSE235663": {"loader": "gse235663", "raw": "11_GSE235663"},
}


QC_THRESHOLDS = {
    key: dict(min_genes=200, max_genes=6000, min_umis=500, max_umis=20000, max_mito=0.1)
    for key in DATASET_ORDER
}
QC_THRESHOLDS["5_GSE134649"]["max_mito"] = 0.2
QC_THRESHOLDS["9_GSE114374"]["max_genes"] = 5000
QC_THRESHOLDS["10_GSE296969"]["max_genes"] = 5000


GSE116222_ID_MAP = {
    "A1": "HC11", "A2": "UC14", "A3": "UC11",
    "B1": "HC12", "B2": "UC15", "B3": "UC12",
    "C1": "HC13", "C2": "UC16", "C3": "UC13",
}
GSE116222_GSM_MAP = {
    "UC11": "GSM3214201", "UC14": "GSM3214202", "HC11": "GSM3214203",
    "UC12": "GSM3214204", "UC15": "GSM3214205", "HC12": "GSM3214206",
    "UC13": "GSM3214207", "UC16": "GSM3214208", "HC13": "GSM3214209",
}


SCP259_SAMPLE_TO_ID = {
    "N10.EpiA":"HC47", "N10.LPA":"HC48", "N10.LPB":"HC49",
    "N11.LPA":"HC50", "N11.LPB":"HC51", "N13.LPA":"HC52",
    "N13.LPB":"HC53", "N15.EpiA":"HC54", "N15.LPA":"HC55",
    "N15.LPB":"HC56", "N16.LPA":"HC57", "N16.LPB":"HC58",
    "N17.EpiB":"HC59", "N17.LPA":"HC60", "N17.LPB":"HC61",
    "N18.EpiB":"HC62", "N18.LPA":"HC63", "N18.LPB":"HC64",
    "N20.LPA":"HC65", "N20.LPB":"HC66", "N21.LPA":"HC67",
    "N21.LPB":"HC68", "N8.LPA":"HC69", "N8.LPB":"HC70",
    "N10.EpiB":"HC71", "N11.EpiA":"HC72", "N11.EpiB":"HC73",
    "N13.EpiA":"HC74", "N13.EpiB":"HC75", "N15.EpiB":"HC76",
    "N16.EpiA":"HC77", "N16.EpiB":"HC78", "N17.EpiA":"HC79",
    "N18.EpiA":"HC80", "N20.EpiA":"HC81", "N20.EpiB":"HC82",
    "N21.EpiA":"HC83", "N21.EpiB":"HC84", "N8.EpiA":"HC85",
    "N8.EpiB":"HC86", "N12.LPB":"UC38", "N14.LPB":"UC39",
    "N19.LPB":"UC40", "N23.LPB":"UC41", "N26.LPB":"UC42",
    "N7.LPB":"UC43", "N9.EpiB":"UC44", "N9.LPB":"UC45",
    "N12.EpiB":"UC46", "N14.EpiB":"UC47", "N19.EpiB":"UC48",
    "N23.EpiB":"UC49", "N26.EpiB":"UC50", "N7.EpiB":"UC51",
    "N12.LPA":"UC52", "N14.LPA":"UC53", "N19.LPA":"UC54",
    "N23.LPA":"UC55", "N26.EpiA":"UC56", "N26.LPA":"UC57",
    "N7.LPA":"UC58", "N9.EpiA":"UC59", "N9.LPA":"UC60",
    "N12.EpiA":"UC61", "N14.EpiA":"UC62", "N19.EpiA":"UC63",
    "N23.EpiA":"UC64", "N7.EpiA":"UC65",
}


L1_CORE = {
    "Immune": ["PTPRC","CD3D","CD3E","CD2","CD79A","MS4A1","NKG7","GNLY","PRF1","LYZ"],
    "Epithelial": ["EPCAM","KRT8","KRT18","KRT19","KRT20","MUC2","TFF3"],
    "Stromal": ["COL1A1","COL3A1","DCN","LUM","PDGFRA","ACTA2","TAGLN","MYH11","PECAM1","VWF"],
    "Other": ["RBFOX3","TUBB3","SNAP25","MAP2","UCHL1","PLP1","SOX10","S100B","MBP"],
}
L1_EXT = {
    "Immune": ["PTPRC","CD3D","CD3E","CD2","CD247","IL7R","CCR7","CD8A","CD4","CD79A","CD79B","MS4A1","BANK1","CD19","CD74","NKG7","GNLY","PRF1","GZMB","GZMH","KLRD1","KLRK1","FGFBP2","LYZ","LST1","S100A8","S100A9","FCGR3A","TYROBP","CTSS","ITGAM","XBP1","JCHAIN","MZB1","IGKC","IGHG1"],
    "Epithelial": ["EPCAM","CDH1","CLDN3","CLDN4","OCLN","KRT8","KRT18","KRT19","KRT20","MUC2","MUC13","TFF3","ALPI","CA1","FABP1","FABP2","REG4","BEST4","AGR2","VIL1","SLC26A3","SPINK4","CEACAM5"],
    "Stromal": ["COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","DCN","LUM","FBLN1","PDGFRA","PDGFRB","ACTA2","TAGLN","MYH11","MYL9","TPM2","CNN1","VIM","DES","SPARC","FAP","PECAM1","VWF","CDH5","ENG","ESAM","KDR","FLT1","PROX1","LYVE1","CLDN5"],
    "Other": ["RBFOX3","TUBB3","SNAP25","SYT1","STMN2","UCHL1","MAP2","NEFL","NEFM","RBFOX1","ELAVL3","SLC17A6","GAD1","GAD2","PLP1","SOX10","MPZ","MBP","PMP22","ERBB3","FABP7","CDH19","GFAP","S100B","PHOX2B","HAND2","RET","ASCL1","CHAT","SLC18A3","SLC5A7","NOS1","VIP","CALB1","CALB2","TAC1","PENK"],
}


@dataclass(frozen=True)
class L2Config:
    l1_label: str
    prefix: str
    panels: dict[str, list[str]]
    resolutions: tuple[float, ...]
    main_resolution: float
    categories: tuple[str, ...]
    tissue_name: str
    marker_gap_quantile: float = 0.2
    score_gap_quantile: float = 0.2


L2_CONFIGS = {
    "immune": L2Config(
        "Immune", "Immune|",
        {
            "T_NK": ["PTPRC","CD3D","CD3E","CD2","TRAC","IL7R","CD4","CD8A","CD8B","NCR1","KLRD1","NKG7","GNLY","PRF1","GZMB","GZMH","KLRK1","FGFBP2","TYROBP"],
            "B": ["MS4A1","CD79A","CD79B","CD74","BANK1","CD19","CD37"],
            "Plasma": ["JCHAIN","MZB1","XBP1","IGHG1","IGKC","IGLC2","SDC1"],
            "Myeloid": ["LYZ","LST1","S100A8","S100A9","FCGR3A","VCAN","FCN1","CTSS","ITGAM","ITGAX"],
        },
        (0.2, 0.4, 0.6), 0.4, ("T_NK","B","Plasma","Myeloid","Other"), "immune cells",
    ),
    "stromal": L2Config(
        "Stromal", "Stromal|",
        {
            "Fibroblast": ["COL1A1","COL1A2","DCN","LUM","COL6A1","COL6A2","PDGFRA","THY1","FBLN1","TAGLN2"],
            "Endothelial": ["PECAM1","VWF","KDR","CDH5","CLDN5","TIE1","KLF2","EMCN","ESAM","ENG","ROBO4"],
            "Pericyte_SMC": ["RGS5","PDGFRB","NOTCH3","MCAM","ACTA2","TAGLN","MYH11","CNN1","MYLK","DES"],
        },
        (0.2, 0.4, 0.6), 0.4, ("Fibroblast","Endothelial","Pericyte_SMC","Other"), "stromal",
    ),
    "epithelial": L2Config(
        "Epithelial", "Epithelial|",
        {
            "Absorptive": ["CA1","CA2","SLC26A3","ALPI","BEST4","OTOP2","GUCA2A","AQP8","MS4A12","CEACAM5","CEACAM6"],
            "Goblet": ["MUC2","TFF3","FCGBP","SPINK4","AGR2"],
            "Enteroendocrine": ["CHGA","CHGB","NEUROD1","PAX6","PCSK1","TPH1","PYY","SST","GCG"],
            "Tuft": ["DCLK1","POU2F3","TRPM5","AVIL","GFI1B"],
            "Stem_TA_prolif": ["LGR5","ASCL2","OLFM4","MKI67","TOP2A","UBE2C","TK1","STMN1","PCNA"],
        },
        (0.3, 0.5, 0.7), 0.5,
        ("Absorptive","Goblet","Enteroendocrine","Tuft","Stem_TA_prolif","Other"),
        "Colonic epithelium",
    ),
    "other": L2Config(
        "Other", "ENS|",
        {
            "Glia": ["GPM6B","PMP22","MPZ","MAL","S100B","CRYAB","PSAP","NCAM1","CADM1","DAG1","SPARC","TIMP3","IGFBP7","FBLN2","COL18A1"],
            "Active_Glia": ["HSPA1A","HSPA1B","HSP90AA1","HSP90AB1","DNAJB1","DNAJA1","FOS","FOSB","JUN","EGR1","IER3","ATF3","HES1","AHR","CCL2","FGL2","LGALS1","IFITM3"],
            "Glia-like_Stromal": ["COL4A1","COL1A2","FBLN2","HSPG2","MCAM","MATN2","SPARC","IGFBP7"],
        },
        (0.1, 0.2, 0.3), 0.2,
        ("Glia","Active_Glia","Glia-like_Stromal","Other"), "Enteric Nervous System",
        marker_gap_quantile=0.0,
        score_gap_quantile=0.0,
    ),
}


# Final cluster decisions recorded in the reviewed L2 notebook outputs. These
# keep the default pipeline deterministic when the external GPT step is off.
L2_REVIEWED_CLUSTER_MAPS = {
    "immune": {
        "0":"T_NK", "1":"Plasma", "2":"T_NK", "3":"B",
        "4":"Plasma", "5":"Other", "6":"Myeloid", "7":"Other",
        "8":"Other", "9":"Plasma", "10":"Other", "11":"Plasma",
    },
    "stromal": {
        "0":"Fibroblast", "1":"Fibroblast", "2":"Fibroblast", "3":"Other",
        "4":"Other", "5":"Pericyte_SMC", "6":"Fibroblast", "7":"Endothelial",
        "8":"Other", "9":"Pericyte_SMC", "10":"Fibroblast", "11":"Fibroblast",
    },
    "epithelial": {
        "0":"Other", "1":"Other", "2":"Absorptive", "3":"Other",
        "4":"Stem_TA_prolif", "5":"Absorptive", "6":"Absorptive", "7":"Goblet",
        "8":"Goblet", "9":"Absorptive", "10":"Tuft", "11":"Other",
    },
    "other": {"0":"Active_Glia", "1":"Glia", "2":"Active_Glia", "3":"Active_Glia"},
}


L3_IMMUNE = {
    "tnk": {
        "subset": "T_NK", "hvg": 4000,
        "resolutions": (0.3, 0.5, 0.7, 0.9), "cluster_key": "leiden_0.9_TNK",
        "stem": "TNK", "mapping": {
            "6":"Treg", "3":"CD8+ Tem", "8":"CD8+ Trm", "5":"CD8+ Teff",
            "10":"NK", "4":"Naive T", "1":"CD4+ Memory", "7":"Tex",
            "2":"T_Stress", "12":"T_Stress", "0":"Low Quality",
            "15":"Low Quality", "16":"Low Quality", "9":"Contaminants",
            "11":"Contaminants", "13":"Contaminants", "14":"Contaminants",
        },
    },
    "myeloid": {
        "subset": "Myeloid", "hvg": 2000,
        "resolutions": (0.2, 0.4), "cluster_key": "leiden_0.4_Myeloid",
        "stem": "Myeloid", "mapping": {
            "1":"Mac Resident", "3":"Mac Inflammatory", "4":"Mac Activated",
            "9":"Mac Remodeling", "2":"Mono Classical", "6":"DC Mature",
            "10":"pDC", "5":"Prolif Myeloid", "0":"Low Quality",
            "7":"Contaminants", "8":"Contaminants",
        },
    },
    "b": {
        "subset": "B", "hvg": 2000,
        "resolutions": (0.2, 0.4), "cluster_key": "leiden_0.4_B",
        "stem": "B", "mapping": {
            "0":"Naive B", "1":"Plasma IgA", "2":"Low Quality",
            "3":"B_Stress", "4":"GC B", "5":"Contaminants",
        },
    },
    "plasma": {
        "subset": "Plasma", "hvg": 2000,
        "resolutions": (0.2, 0.4), "cluster_key": "leiden_0.4_Plasma",
        "stem": "Plasma", "mapping": {
            "2":"Plasma IgG", "1":"Plasma IgA", "6":"Plasma IgA",
            "0":"Plasma Stressed", "3":"Contaminants", "5":"Contaminants",
            "4":"Low Quality", "7":"Low Quality",
        },
    },
    "mast": {
        "subset": "Mast", "hvg": 2000,
        "resolutions": (0.1, 0.2), "cluster_key": "leiden_0.2_Mast",
        "stem": "Mast", "mapping": None,
    },
}


L3_STROMAL = {
    "fibroblast": {
        "subset":"Fibroblast", "hvg":2000, "resolutions":(0.2,0.4),
        "cluster_key":"leiden_0.2_fibroblast", "stem":"fibroblast",
        "mapping": {"0":"Inflammatory Fibro","1":"Villus-top Fibro","2":"Crypt-bottom Fibro","3":"Contaminants","4":"Contaminants","5":"Myofibroblasts","6":"Contaminants","7":"Cycling Fibroblasts"},
    },
    "endothelial": {
        "subset":"Endothelial", "hvg":1000, "resolutions":(0.2,0.1),
        "cluster_key":"leiden_0.2_Endo", "stem":"endo", "cluster_suffix":"Endo",
        "mapping": {"0":"Arterial ECs","1":"Contaminants","2":"Inflammatory Homing ECs","3":"Contaminants","4":"Capillary ECs"},
    },
    "pericyte_smc": {
        "subset":"Pericyte_SMC", "hvg":2000, "resolutions":(0.2,0.4),
        "cluster_key":"leiden_0.4_SMC", "stem":"SMC",
        "mapping": {"0":"Visceral SMCs","1":"Pericytes","2":"Activated Pericytes","3":"Inflammatory SMCs","4":"Vascular SMCs"},
    },
}
