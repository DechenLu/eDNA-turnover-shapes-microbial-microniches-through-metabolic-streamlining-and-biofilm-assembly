#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DNA 降解酶 / DNA degradation 相关基因详细统计脚本（final core 高置信主结果版）
================================================================================

用途
----
本脚本参照 HGT 汇总脚本的输入、递归读取、最佳注释表选择和输出方式，
用于从同一套 MAG/基因组功能注释表中筛选高置信 DNA 降解酶相关基因。

本版比 publicdb / strict / main_only 版本进一步收紧：
1. 不保留 REVIEW_OPTIONAL 和 REVIEW_WEAK_EVIDENCE 类条目；
2. 默认剔除 DNA repair / recombination-only nuclease、RNA/tRNA/rRNA 相关、
   polymerase exonuclease domain、广谱 hydrolase/phosphoesterase、methylase/specificity
   subunit-only、CRISPR/Cas/transposase-associated nuclease 和普通代谢酶；
3. 只保留有明确 DNA 降解/切割证据的条目，例如：
   DNase / deoxyribonuclease / exodeoxyribonuclease / endodeoxyribonuclease / 
   extracellular DNA degradation protein / RecJ ssDNA exonuclease / Exonuclease VII / 
   SbcB-ExoI / restriction endonuclease / type I-IV restriction enzyme R protein / Mrr。

输入
----
python dna_degradation_detail_summary_final_core.py <input_file_or_dir> <output_prefix>

支持 csv / tsv / txt / tsvm / xlsx / xls。
如果输入为目录，会递归读取子目录；同一子目录中优先选择文件名含
karen_rundbcan_prokka_egg_nog > karen_rundbcan > karen 的版本。

输出
----
- *_detailed.csv
- *_DNA_degradation_hits_only.csv
- *_sample_summary.csv
- *_group_summary.csv（如果存在 group/treatment/condition 列）
- *_pivot_main_class.csv
- *_pivot_main_category.csv
- *_pivot_expanded_class.csv
- *_pivot_expanded_abundance.csv
- *_sample_summary.json
"""

import os
import re
import json
import glob
import argparse
import warnings
from typing import Dict, List, Tuple, Iterable, Set

import pandas as pd

warnings.filterwarnings("ignore")


# =========================
# 1) 分类定义
# =========================
DNADEG_CLASS_TO_CATEGORY = {
    "extracellular_dnase": "dna_degradation_enzyme",
    "dnase_nuclease": "dna_degradation_enzyme",
    "endonuclease_exonuclease": "dna_degradation_enzyme",
    "restriction_endonuclease": "dna_degradation_enzyme",
}
CLASS_PRIORITY = [
    "extracellular_dnase",
    "dnase_nuclease",
    "restriction_endonuclease",
    "endonuclease_exonuclease",
]
CATEGORY_PRIORITY = ["dna_degradation_enzyme"]
CORE_DNADEG_CLASSES = set(DNADEG_CLASS_TO_CATEGORY.keys())
CORE_DNADEG_CATEGORIES = set(DNADEG_CLASS_TO_CATEGORY.values())


# =========================
# 2) 高置信白名单
# =========================
# 只保留直接 DNA 降解/切割证据。
CORE_KEYWORD_PATTERNS = {
    "extracellular_dnase": [
        r"\bextracellular\b.{0,80}\bdna\b.{0,80}\bdegradation\b",
        r"\bextracel+ular dna degradation protein\b",
        r"\bextracellular dna degradation protein\b",
        r"\bextracellular\b.{0,80}\bdeoxyribonuclease\b",
        r"\bextracellular\b.{0,80}\bdnase\b",
        r"\bedd[ab]\b",
    ],
    "dnase_nuclease": [
        r"\bdnase\b",
        r"\bdnase\s*[i1]\b",
        r"\bdnase\s*ii\b",
        r"\bdeoxyribonuclease\b",
        r"\bdeoxyribonuclease\s*[i1]\b",
        r"\bdeoxyribonuclease\s*ii\b",
        r"\bexodeoxyribonuclease\b",
        # generic endodeoxyribonuclease is not used alone because it captures RuvC/Holliday repair nucleases
        r"\bdna/rna non-specific endonuclease\b",
        r"\bdna-rna non-specific endonuclease\b",
        r"\bendonuclease g\b",
        r"\bstaphylococcal nuclease\b",
        r"\bstaphylococcal nuclease domain\b",
        r"\bsnase\b",
        r"\bthermonuclease\b",
        r"\bsite-specific deoxyribonuclease\b",
    ],
    "endonuclease_exonuclease": [
        r"\bexodeoxyribonuclease\s*(?:i|1|vii|7)\b",
        r"\bexodeoxyribonuclease\s*7\b",
        r"\bexonuclease\s*vii\b",
        r"\bexonuclease\s*7\b",
        r"\bexodeoxyribonuclease vii\b",
        r"\bexodeoxyribonuclease 7 large subunit\b",
        r"\bexodeoxyribonuclease 7 small subunit\b",
        r"\bexodeoxyribonuclease vii large subunit\b",
        r"\bexodeoxyribonuclease vii small subunit\b",
        r"\bsingle-stranded-dna-specific exonuclease recj\b",
        r"\brecj-like\b.{0,80}\bdna\b.{0,80}\bexonuclease\b",
        r"\bsbcb\b.{0,80}\bexodeoxyribonuclease\b",
        r"\bexoi\b.{0,80}\bexodeoxyribonuclease\b",
        r"\bexodeoxyribonuclease i\b",
    ],
    "restriction_endonuclease": [
        r"\brestriction endonuclease\b",
        r"\btype[- ]?i restriction endonuclease\b",
        r"\btype[- ]?ii restriction endonuclease\b",
        r"\btype[- ]?iii restriction endonuclease\b",
        r"\btype[- ]?iv restriction endonuclease\b",
        r"\btype[- ]?[ivx]+ restriction enzyme r protein\b",
        r"\btype[- ]?[ivx]+ restriction enzyme r subunit\b",
        r"\btype[- ]?[ivx]+ site-specific restriction-modification system r \(restriction\) subunit\b",
        r"\bpredicted type iv restriction endonuclease\b",
        r"\b5-methylcytosine-specific restriction enzyme\b",
        r"\bmethylcytosine-specific restriction enzyme\b",
        r"\bmrr restriction system protein\b",
        r"\bmrr restriction protein\b",
    ],
}

# 高置信 KO。只保留直接 DNA degradation / DNase / exodeoxyribonuclease / restriction endonuclease。
KO_TO_CLASS = {
    "K07004": "extracellular_dnase",       # eddB; extracellular DNA degradation protein EddB
    "K01113": "extracellular_dnase",       # eddA / DNase-like extracellular DNA degradation protein
    "K01173": "dnase_nuclease",            # endonuclease G / DNA-RNA non-specific endonuclease
    "K15979": "dnase_nuclease",            # staphylococcal nuclease domain-containing protein
    "K01141": "endonuclease_exonuclease",  # sbcB/exoI; exodeoxyribonuclease I
    "K01144": "endonuclease_exonuclease",  # exodeoxyribonuclease VII small subunit-like, retained only with direct text guard
    "K01153": "endonuclease_exonuclease",  # exonuclease-related, retained only with direct text guard
    "K07451": "restriction_endonuclease",  # 5-methylcytosine-specific restriction enzyme
    "K07504": "restriction_endonuclease",  # predicted type IV restriction endonuclease
}
KO_NOTES = {
    "K07004": "eddB; extracellular DNA degradation protein EddB",
    "K01113": "eddA; extracellular DNA degradation protein EddA / DNase-like",
    "K01173": "endonuclease G / DNA-RNA non-specific endonuclease",
    "K15979": "staphylococcal nuclease domain-containing protein",
    "K01141": "sbcB/exoI; exodeoxyribonuclease I",
    "K01144": "exodeoxyribonuclease-like KO; kept only with direct text evidence",
    "K01153": "exonuclease-like KO; kept only with direct text evidence",
    "K07451": "5-methylcytosine-specific restriction enzyme",
    "K07504": "predicted type IV restriction endonuclease",
}
# 这些 KO 只有在同一行有直接文字证据时才保留，避免 KO 版本差异带来的误收。
KO_REQUIRE_DIRECT_TEXT = {"K01144", "K01153"}

# 高置信 COG。删除 COG0507/0732/宽泛 COG，只保留直接 exonuclease / restriction enzyme。
COG_TO_CLASS = {
    "COG0608": "endonuclease_exonuclease",  # RecJ-like single-stranded DNA exonuclease
    "COG1570": "endonuclease_exonuclease",  # Exodeoxyribonuclease VII large subunit / XseA-like
    "COG2925": "endonuclease_exonuclease",  # Exodeoxyribonuclease I / SbcB-like
    "COG2810": "restriction_endonuclease",  # predicted type IV restriction endonuclease
    "COG1715": "restriction_endonuclease",  # Mrr restriction system protein
}

# 高置信 PFAM/domain。删除 HNH/GIY-YIG/PDDEXK/Exo_endo_phos/Metallophos/DHH/DHHA1/Pol exonuclease 等 domain-only。
PFAM_TO_CLASS = {
    "SNase": "dnase_nuclease",
    "SNase_OB": "dnase_nuclease",
    "DNase-RNase": "dnase_nuclease",
    "DNase_I": "dnase_nuclease",
    "DNase_II": "dnase_nuclease",
    "Exonuc_VII_L": "endonuclease_exonuclease",
    "Exonuc_VII_S": "endonuclease_exonuclease",
    "RE_BstXI": "restriction_endonuclease",
    "Mrr_N": "restriction_endonuclease",
    "Mrr_cat": "restriction_endonuclease",
}

# 高置信 gene/preferred。删除 HsdM/HsdS、RuvC、RecBCD、SbcCD、Nfo/Nth/Nfi/UvrC 等 optional/repair。
GENE_TO_CLASS = {
    "eddA": "extracellular_dnase",
    "eddB": "extracellular_dnase",
    "recJ": "endonuclease_exonuclease",
    "xseA": "endonuclease_exonuclease",
    "xseB": "endonuclease_exonuclease",
    "sbcB": "endonuclease_exonuclease",
    "exoI": "endonuclease_exonuclease",
    "ENDOG": "dnase_nuclease",
    "DNASE1": "dnase_nuclease",
    "DNASE1L1": "dnase_nuclease",
    "DNASE1L2": "dnase_nuclease",
    "DNASE1L3": "dnase_nuclease",
    "DNASE2": "dnase_nuclease",
    "DNASE2B": "dnase_nuclease",
    "mrr": "restriction_endonuclease",
    "hsdR": "restriction_endonuclease",
}
# gene-only 需要直接文字辅助的列表，避免单基因名误判。
GENE_REQUIRE_DIRECT_TEXT = {"ENDOG", "DNASE1", "DNASE1L1", "DNASE1L2", "DNASE1L3", "DNASE2", "DNASE2B", "hsdR"}


# =========================
# 3) 强假阳性 / 边缘项删除规则
# =========================
FALSE_POSITIVE_PATTERNS = [
    # RNA / rRNA / tRNA related
    r"\bribonuclease\b", r"\brnase\b", r"\boligoribonuclease\b", r"\brna degradation\b",
    r"\brna-binding\b", r"\brna helicase\b", r"\brna methyltransferase\b",
    r"\bribosomal rna\b", r"\brrna\b", r"\b23s rrna\b", r"\b16s rrna\b",
    r"\btrna\b", r"\baminoacyl-trna\b", r"\bd-aminoacyl-trna deacylase\b",
    r"\banticodon\b", r"\btrna_anti\b", r"\bpap phosphatase\b", r"\bnrna\b",

    # DNA repair / recombination / replication-only
    r"\bdna repair\b", r"\bnucleotide excision repair\b", r"\bbase excision repair\b",
    r"\bmismatch repair\b", r"\brecombinational repair\b", r"\brecombination protein\b",
    r"\bholliday junction\b", r"\bcrossover junction\b", r"\bruvc\b", r"\brecu\b",
    r"\bgen1\b", r"\byen1\b", r"\bzranb3\b", r"\breca\b", r"\brecb\b", r"\brecc\b", r"\brecd\b",
    r"\baddab\b", r"\badda\b", r"\baddb\b", r"\bsbcc\b", r"\bsbcd\b", r"\bmre11\b", r"\brad50\b",
    r"\buvra\b", r"\buvrb\b", r"\buvrc\b", r"\buvrd\b", r"\buvrabc\b",
    r"\bendonuclease iii\b", r"\bendonuclease iv\b", r"\bendonuclease v\b",
    r"\bap endonuclease\b", r"\bdna glycosylase\b", r"\buv dna damage endonuclease\b", r"\bnucs\b",
    r"\bdna polymerase\b", r"\bpolymerase exonuclease\b", r"\bdna_pol_[ab]_exo\b",

    # methylase / specificity subunit only, not direct cleavage
    r"\bmethylase\b", r"\bdna methyltransferase\b", r"\bmethyltransferase\b", r"\bn6[_ -]?mtase\b",
    r"\bmethylase_s\b", r"\bhsdm\b", r"\bhsds\b", r"\bspecificity subunit\b",

    # broad hydrolase / phosphoesterase / phosphodiesterase / non-specific domains
    r"\bmetallophos\b", r"\bmetallophosphoesterase\b", r"\bphosphoesterase\b", r"\bphosphodiesterase\b",
    r"\bmetal-dependent hydrolase\b", r"\bputative metal-dependent hydrolase\b",
    r"\bdhh\b", r"\bdhha1\b", r"\bduf2161\b", r"\bpho[dD]\b", r"\bpyrophosphatase\b",
    r"\bendonuclease/exonuclease/phosphatase family\b", r"\bexo_endo_phos\b",

    # CRISPR / transposase / mobile element nucleases not DNA degradation enzymes
    r"\bcrispr\b", r"\bcas[0-9a-z_-]*\b", r"\bcas_cas4\b", r"\btransposase\b", r"\bintegrase\b",
    r"\bmobile element\b", r"\bpd-\(?d/e\)?xk nuclease family transposase\b",

    # ordinary regulators, ATPases, transporters, metabolic enzymes
    r"\btranscriptional regulator\b", r"\btranscriptional repressor\b", r"\bdna-binding transcriptional\b",
    r"\bhelix-turn-helix\b", r"\bwinged helix\b", r"\batpase\b", r"\bsmc\b", r"\brara\b",
    r"\btransporter\b", r"\bpermease\b", r"\bbeta-hexosaminidase\b", r"\bprotease\b", r"\bpeptidase\b",
    r"\bglycosidase\b", r"\bphospholipase\b", r"\btransferase\b",
]



# 3.1) 最终强删除规则（不可被 direct_override 覆盖）
# 根据 ultra_strict 结果的二次质控加入，目标是只保留 KEEP_CORE 主结果。
# 这些类型即使含有 nuclease/deoxyribonuclease 字样，也默认不进入 DNA 降解酶主分析。
FINAL_HARD_REJECT_PATTERNS = [
    # RecBCD/AddAB/SbcCD 等 DNA repair/recombination nuclease，不作为 DNA degradation enzyme 主结果
    r"\brecb\b", r"\brecc\b", r"\brecd\b", r"\brecbcd\b",
    r"\batp-dependent helicase/deoxyribonuclease subunit [bc]\b",
    r"\bexodeoxyribonuclease v\b", r"\bexonuclease v\b",
    r"\baddab\b", r"\badda\b", r"\baddb\b",
    r"\bsbcc\b", r"\bsbcd\b", r"\bsbc[cd]\b",
    r"\bmre11\b", r"\brad50\b",
    r"\bholliday junction\b", r"\bcrossover junction\b", r"\bruvc\b", r"\brecu\b",
    r"\bgen1\b", r"\byen1\b", r"\bzranb3\b",
    r"\bnucleotide excision repair\b", r"\bbase excision repair\b", r"\bmismatch repair\b",
    r"\buvr[a-d]\b", r"\buvrabc\b", r"\bnucs\b",
    r"\bendonuclease iii\b", r"\bendonuclease iv\b", r"\bendonuclease v\b",
    r"\bap endonuclease\b", r"\bdna glycosylase\b", r"\buv dna damage endonuclease\b",

    # TatD/YcfH/YjjV/metal hydrolase ambiguous，上一轮质控建议整体剔除
    r"\btatd\b", r"\btatd-related\b",
    r"\bycfh\b", r"\byjjv\b", r"\bputative metal-dependent hydrolase\b",
    r"\bmetal-dependent hydrolase\b", r"\bmetallophos\b", r"\bmetallophosphoesterase\b",
    r"\bdhh\b", r"\bdhha1\b", r"\bduf2161\b", r"\bexo_endo_phos\b",
    r"\bendonuclease/exonuclease/phosphatase family\b", r"\bphosphoesterase\b", r"\bphosphodiesterase\b",

    # RNA/tRNA/rRNA related
    r"\bribonuclease\b", r"\brnase\b", r"\boligoribonuclease\b", r"\brna degradation\b",
    r"\bribosomal rna\b", r"\brrna\b", r"\b23s rrna\b", r"\b16s rrna\b",
    r"\btrna\b", r"\baminoacyl-trna\b", r"\bd-aminoacyl-trna deacylase\b",
    r"\btrna_anti\b", r"\banticodon\b", r"\bpap phosphatase\b",

    # Polymerase / replication / proofreading
    r"\bdna polymerase\b", r"\bpolymerase exonuclease\b", r"\bdna_pol_[ab]_exo\b",
    r"\bproofreading\b", r"\breplication-associated\b",

    # Restriction-modification non-cleavage components
    r"\bhsdm\b", r"\bhsds\b", r"\bspecificity subunit\b", r"\bs subunit\b",
    r"\bmethylase\b", r"\bmethylase_s\b", r"\bn6[_ -]?mtase\b",
    r"\bdna methyltransferase\b", r"\bmethyltransferase\b",

    # CRISPR/Cas, transposase/mobile element nuclease
    r"\bcrispr\b", r"\bcas[0-9a-z_-]*\b", r"\bcas_cas4\b",
    r"\btransposase\b", r"\bintegrase\b", r"\bmobile element\b",
    r"\bpd-\(?d/e\)?xk nuclease family transposase\b",

    # obvious non-DNA-degradation enzymes / cell proteins
    r"\bperiplasmic serine endoprotease\b", r"\bdegp\b", r"\bdegq\b",
    r"\blipopolysaccharide assembly protein\b", r"\bcell division coordinator\b", r"\bcpob\b",
    r"\bcardiolipin synthase\b", r"\bbeta-hexosaminidase\b", r"\bprotease\b", r"\bpeptidase\b",
    r"\bglycosidase\b", r"\bphospholipase\b", r"\btransferase\b", r"\btransporter\b",
    r"\btranscriptional regulator\b", r"\bdna-binding transcriptional\b", r"\batpase\b", r"\bsmc\b", r"\brara\b",
]

# 如果同行包含这些强直接证据，可以覆盖部分假阳性词。例如 DNase/RNase 双功能酶含 RNase 字样但仍应保留。
DIRECT_OVERRIDE_PATTERNS = [
    r"\bdnase\b", r"\bdeoxyribonuclease\b", r"\bexodeoxyribonuclease\b",
    r"\bdna/rna non-specific endonuclease\b", r"\bdna-rna non-specific endonuclease\b",
    r"\bextracellular dna degradation protein\b", r"\bextracel+ular dna degradation protein\b",
    r"\bexonuclease vii\b", r"\bexodeoxyribonuclease vii\b", r"\bexodeoxyribonuclease 7\b",
    r"\bsingle-stranded-dna-specific exonuclease recj\b", r"\bexodeoxyribonuclease i\b",
    r"\brestriction endonuclease\b", r"\btype[- ]?[ivx]+ restriction enzyme r protein\b",
    r"\btype[- ]?[ivx]+ restriction enzyme r subunit\b", r"\b5-methylcytosine-specific restriction enzyme\b",
    r"\bmrr restriction system protein\b", r"\bmrr restriction protein\b",
]


def first_nonempty(*vals):
    for v in vals:
        if v is not None:
            s = str(v).strip()
            if s and s.lower() != "nan" and s != "-":
                return s
    return ""


def normalize_text(x) -> str:
    s = str(x or "")
    s = s.replace("′", "'").replace("’", "'").replace("–", "-").replace("—", "-")
    s = re.sub(r"\s+", " ", s).strip()
    return s


_REGEX_ANY_CACHE = {}

def has_regex(patterns: Iterable[str], text: str) -> bool:
    pats = tuple([p for p in patterns if p])
    if not pats:
        return False
    rgx = _REGEX_ANY_CACHE.get(pats)
    if rgx is None:
        rgx = re.compile("|".join(f"(?:{p})" for p in pats), flags=re.I)
        _REGEX_ANY_CACHE[pats] = rgx
    return rgx.search(text) is not None


def extract_kos(text: str) -> List[str]:
    return sorted(set(re.findall(r"K\d{5}", str(text).upper())))


def extract_cogs(text: str) -> List[str]:
    return sorted(set(re.findall(r"COG\d{4}", str(text).upper())))


def token_pattern(term: str) -> str:
    return r"(?<![A-Za-z0-9_])" + re.escape(term) + r"(?![A-Za-z0-9_])"


def standardize_columns(df):
    mapping = {}
    for col in df.columns:
        low = str(col).lower().strip()
        if "filename" in low or low in {"sample", "sample_id", "bin"}:
            mapping["fileName"] = col
        elif low == "genome" or low == "genome_id" or "genome_id" in low:
            mapping.setdefault("fileName", col)
        elif low in {"group", "treatment", "condition", "group1"}:
            mapping["group"] = col
        elif low == "gene" or low.startswith("gene_") or "locus_tag" in low:
            mapping["gene"] = col
        elif "preferred_name" in low or "preferred" in low:
            mapping["preferred"] = col
        elif "description" in low:
            mapping["description"] = col
        elif low == "product" or "product" in low:
            mapping["product"] = col
        elif "uniprot" in low:
            mapping["uniprot"] = col
        elif low == "ec_number" or low == "ec" or "ec___" in low:
            mapping["ec"] = col
        elif "seed_ortholog" in low:
            mapping["seed_ortholog"] = col
        elif "eggnog_ogs" in low:
            mapping["eggnog_ogs"] = col
        elif "kegg_ko" in low or low == "ko" or "ko_id" in low:
            mapping["ko"] = col
        elif low == "cog" or ("cog" in low and "category" not in low):
            mapping["cog"] = col
        elif "pfam" in low:
            mapping["pfam"] = col
        elif "kegg_tc" in low:
            mapping["kegg_tc"] = col
        elif low == "go" or "gos" in low:
            mapping["go"] = col
        elif "run_dbcan" in low:
            mapping["dbcan"] = col
        elif "cazy" in low:
            mapping["cazy"] = col
        elif "karen" in low:
            mapping["karen"] = col
        elif "abundance" in low or "rpm" in low or "tpm" in low or "rpkm" in low or low == "count":
            mapping["abundance"] = col
    return mapping


def require_minimum_columns(col_map):
    useful = {"gene", "preferred", "description", "product", "ko", "cog", "pfam", "uniprot", "seed_ortholog"}
    found = useful.intersection(set(col_map.keys()))
    if len(found) < 2:
        raise ValueError(
            f"可识别注释列过少，仅识别到: {sorted(found)}。"
            "请至少提供 gene / description / product / preferred / KO / COG / PFAM / UniProtKB / seed_ortholog 中的若干列。"
        )


def get_abundance(row, col_map):
    if "abundance" not in col_map:
        return 1.0
    try:
        return float(row.get(col_map["abundance"], 1.0))
    except Exception:
        return 1.0


def collect_annotation_text(row, col_map):
    fields = [
        "gene", "preferred", "description", "product", "uniprot", "seed_ortholog",
        "eggnog_ogs", "karen", "ko", "cog", "pfam", "kegg_tc", "go", "ec", "dbcan", "cazy"
    ]
    parts = []
    for key in fields:
        col = col_map.get(key)
        if col is not None:
            parts.append(first_nonempty(row.get(col)))
    return normalize_text(" ".join([x for x in parts if x])).lower()


def pfam_hits(text: str) -> List[str]:
    hits = []
    for term in PFAM_TO_CLASS:
        if re.search(token_pattern(term), text, flags=re.I):
            hits.append(term)
    return sorted(set(hits))


def gene_hits(text: str) -> List[str]:
    hits = []
    for term in GENE_TO_CLASS:
        if re.search(r"\b" + re.escape(term) + r"\b", text, flags=re.I):
            hits.append(term)
    return sorted(set(hits))


def is_false_positive_context(text: str) -> bool:
    if not has_regex(FALSE_POSITIVE_PATTERNS, text):
        return False
    # 强直接证据覆盖部分假阳性词。但若是明确 RNA/tRNA/rRNA 或 polymerase/repair-only，上面强阳性一般不会同时出现。
    return not has_regex(DIRECT_OVERRIDE_PATTERNS, text)


def classify_dna_degradation(full_info: str, ko_text: str = ""):
    text = normalize_text(full_info).lower()
    hits: Set[str] = set()
    source_map: Dict[str, List[str]] = {}

    # final core: hard reject first, no override
    if has_regex(FINAL_HARD_REJECT_PATTERNS, text):
        return []

    if is_false_positive_context(text):
        return []

    has_direct_text = has_regex(DIRECT_OVERRIDE_PATTERNS, text)

    # 关键词命中：只使用直接高置信关键词。
    for cls, pats in CORE_KEYWORD_PATTERNS.items():
        for patt in pats:
            if re.search(patt, text, flags=re.I):
                hits.add(cls)
                source_map.setdefault(cls, []).append("direct_keyword/description/product")
                break

    # KO 命中。
    for ko in extract_kos((ko_text or "") + " " + text):
        if ko not in KO_TO_CLASS:
            continue
        if ko in KO_REQUIRE_DIRECT_TEXT and not has_direct_text:
            continue
        cls = KO_TO_CLASS[ko]
        hits.add(cls)
        source_map.setdefault(cls, []).append(f"KEGG_KO({ko})")
        if ko in KO_NOTES:
            source_map.setdefault(cls, []).append(KO_NOTES[ko])

    # COG 命中。
    for cog in extract_cogs(text):
        if cog in COG_TO_CLASS:
            cls = COG_TO_CLASS[cog]
            hits.add(cls)
            source_map.setdefault(cls, []).append(f"COG({cog})")

    # PFAM 命中：只保留强直接结构域；Mrr domain 如遇 methylase/specificity-only 会被 FP 规则挡掉。
    for pf in pfam_hits(text):
        cls = PFAM_TO_CLASS[pf]
        hits.add(cls)
        source_map.setdefault(cls, []).append(f"PFAM({pf})")

    # gene/preferred 命中。部分基因需要直接文字证据辅助。
    for g in gene_hits(text):
        if g in GENE_REQUIRE_DIRECT_TEXT and not has_direct_text:
            continue
        cls = GENE_TO_CLASS[g]
        hits.add(cls)
        source_map.setdefault(cls, []).append(f"gene_or_preferred({g})")

    # 二次强删除：即使命中，也删掉 repair/recombination-only 和 RNA-only。RecJ/XseA/SbcB 等被强直接词保护。
    repair_optional = [
        r"\bholliday junction\b", r"\bcrossover junction\b", r"\brecombination\b", r"\bruvc\b", r"\brecu\b",
        r"\brecb\b", r"\brecc\b", r"\brecd\b", r"\baddab\b", r"\buvr[a-d]\b", r"\bnucs\b",
        r"\bendonuclease iii\b", r"\bendonuclease iv\b", r"\bendonuclease v\b", r"\bdna polymerase\b",
    ]
    repair_keep = [
        r"\bdnase\b", r"\bdeoxyribonuclease\b", r"\bexodeoxyribonuclease\b",
        r"\bexonuclease vii\b", r"\bexodeoxyribonuclease vii\b", r"\bsingle-stranded-dna-specific exonuclease recj\b",
        r"\bsbcb\b.{0,80}\bexodeoxyribonuclease\b", r"\bexoi\b.{0,80}\bexodeoxyribonuclease\b",
        r"\brestriction endonuclease\b", r"\btype[- ]?[ivx]+ restriction enzyme r\b",
    ]
    if has_regex(repair_optional, text) and not has_regex(repair_keep, text):
        return []

    # Exodeoxyribonuclease III / Xth is primarily a DNA repair AP endonuclease/exonuclease;
    # it is excluded from the main DNA degradation enzyme result.
    if has_regex([r"\bxth\b", r"\bexodeoxyribonuclease iii\b", r"\bexodeoxyribonuclease 3\b"], text):
        return []

    # Remove restriction-modification associated non-cleavage components even if wording contains
    # "restriction" or an ambiguous KO. Keep R protein / R subunit / Mrr / true restriction endonuclease.
    rm_non_cleavage = [
        r"\bspecificity subunit\b", r"\bs subunit\b", r"\bs subunits\b", r"\bhsds\b",
        r"\bmethylase\b", r"\bmethylase_s\b", r"\bn6[_ -]?mtase\b", r"\bdna methyltransferase\b", r"\bhsdm\b",
    ]
    rm_keep = [
        r"\btype[- ]?[ivx]+ restriction enzyme r protein\b",
        r"\btype[- ]?[ivx]+ restriction enzyme r subunit\b", r"\bmrr restriction\b",
        r"\b5-methylcytosine-specific restriction enzyme\b",
    ]
    if has_regex(rm_non_cleavage, text) and not has_regex(rm_keep, text):
        return []

    if not hits:
        return []

    return [(cls, ";".join(sorted(set(source_map.get(cls, []))))) for cls in sorted(hits)]


def summarize_categories(classes):
    return sorted({DNADEG_CLASS_TO_CATEGORY[c] for c in classes if c in DNADEG_CLASS_TO_CATEGORY})


def score_dna_degradation_classes(classes):
    cls = set(classes)
    score = 0
    if "extracellular_dnase" in cls:
        score += 4
    if "dnase_nuclease" in cls:
        score += 3
    if "restriction_endonuclease" in cls:
        score += 3
    if "endonuclease_exonuclease" in cls:
        score += 3
    if len(cls) >= 2:
        score += 1
    return score


def classify_row(row, col_map):
    ko = first_nonempty(row.get(col_map.get("ko")))
    full_info = collect_annotation_text(row, col_map)
    results = classify_dna_degradation(full_info, ko)
    if not results:
        return pd.Series({
            "dna_deg_category": None,
            "dna_deg_class": None,
            "dna_deg_evidence": None,
            "dna_deg_all_categories": None,
            "dna_deg_all_classes": None,
            "has_DNA_degradation": 0,
            "dna_deg_risk_score": 0,
        })
    all_classes = [x[0] for x in results]
    all_categories = summarize_categories(all_classes)
    all_evidences = [f"{x[0]}:{x[1]}" for x in results]
    best_class = next((x for x in CLASS_PRIORITY if x in all_classes), sorted(all_classes)[0])
    return pd.Series({
        "dna_deg_category": DNADEG_CLASS_TO_CATEGORY.get(best_class),
        "dna_deg_class": best_class,
        "dna_deg_evidence": " | ".join(all_evidences),
        "dna_deg_all_categories": ";".join(all_categories),
        "dna_deg_all_classes": ";".join(sorted(all_classes)),
        "has_DNA_degradation": 1,
        "dna_deg_risk_score": score_dna_degradation_classes(all_classes),
    })


def read_any_table(path):
    lp = str(path).lower()
    if lp.endswith(".xlsx") or lp.endswith(".xls"):
        return pd.read_excel(path)
    if lp.endswith(".tsv") or lp.endswith(".tsvm") or lp.endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def file_priority_score(path: str) -> Tuple[int, int]:
    base = os.path.basename(path)
    low = base.lower()
    if "karen_rundbcan_prokka_egg_nog" in low:
        score = 4
    elif "karen_rundbcan" in low:
        score = 3
    elif "karen" in low:
        score = 2
    else:
        score = 1
    return (score, len(base))


def choose_best_files(root_dir: str, files: List[str]) -> List[str]:
    groups: Dict[str, List[str]] = {}
    for f in files:
        rel = os.path.relpath(f, root_dir)
        parent = os.path.dirname(rel)
        groups.setdefault(parent, []).append(f)
    chosen = []
    for _, flist in groups.items():
        chosen.append(sorted(flist, key=file_priority_score, reverse=True)[0])
    return sorted(chosen)


def collect_input_tables(input_path):
    all_dfs = []
    if os.path.isdir(input_path):
        pats = ["*.csv", "*.tsv", "*.tsvm", "*.txt", "*.xlsx", "*.xls"]
        files = []
        for p in pats:
            files.extend(glob.glob(os.path.join(input_path, "**", p), recursive=True))
        files = sorted(set(files))
        if not files:
            raise ValueError(f"目录中未找到可读取文件: {input_path}\n当前支持格式: {', '.join(pats)}")
        chosen_files = choose_best_files(input_path, files)
        print(f">>> 在目录及子目录中共找到 {len(files)} 个文件，实际选用 {len(chosen_files)} 个最佳版本文件")
        for f in chosen_files:
            try:
                df = read_any_table(f)
                if "fileName" not in df.columns:
                    df["fileName"] = os.path.relpath(f, input_path)
                all_dfs.append(df)
                print(f"    [OK] {os.path.relpath(f, input_path)}")
            except Exception as e:
                print(f"    [跳过] {os.path.relpath(f, input_path)} -> {e}")
        if not all_dfs:
            raise ValueError(f"目录中虽然找到了文件，但都读取失败: {input_path}")
    else:
        df = read_any_table(input_path)
        if "fileName" not in df.columns:
            df["fileName"] = os.path.basename(input_path)
        all_dfs.append(df)
    full_df = pd.concat(all_dfs, ignore_index=True)
    full_df = full_df[~full_df["fileName"].astype(str).str.contains("*", regex=False, na=False)].copy()
    return full_df


def build_sample_summary(df_out, col_map):
    sample_rows = []
    for sample, grp in df_out.groupby("fileName"):
        hit_grp = grp[grp["has_DNA_degradation"] == 1]
        class_set = set()
        category_set = set()
        for x in hit_grp["dna_deg_all_classes"].dropna():
            class_set.update([y for y in str(x).split(";") if y])
        for x in hit_grp["dna_deg_all_categories"].dropna():
            category_set.update([y for y in str(x).split(";") if y])
        class_counts = {cls: 0 for cls in CLASS_PRIORITY}
        if not hit_grp.empty:
            exp_cls = []
            for x in hit_grp["dna_deg_all_classes"].dropna():
                exp_cls.extend([y for y in str(x).split(";") if y])
            for c in exp_cls:
                if c in class_counts:
                    class_counts[c] += 1
        sample_rows.append({
            "Sample": sample,
            "Group": first_nonempty(grp[col_map["group"]].iloc[0]) if "group" in col_map else "NA",
            "Total_Features": len(grp),
            "DNA_Degradation_Hit_Count": int((grp["has_DNA_degradation"] == 1).sum()),
            "DNA_Degradation_Abundance": float(hit_grp["feature_abundance"].sum()),
            "Mean_DNA_Degradation_Risk_Score": round(float(hit_grp["dna_deg_risk_score"].mean()), 3) if len(hit_grp) > 0 else 0.0,
            "Max_DNA_Degradation_Risk_Score": int(hit_grp["dna_deg_risk_score"].max()) if len(hit_grp) > 0 else 0,
            "Unique_DNA_Degradation_Category_Count": len(category_set),
            "Unique_DNA_Degradation_Categories": ";".join(sorted(category_set)),
            "Category_Coverage_%": round(len(category_set.intersection(CORE_DNADEG_CATEGORIES)) / len(CORE_DNADEG_CATEGORIES) * 100, 1) if category_set else 0.0,
            "Unique_DNA_Degradation_Class_Count": len(class_set),
            "Unique_DNA_Degradation_Classes": ";".join(sorted(class_set)),
            "Core_DNA_Degradation_Coverage_%": round(len(class_set.intersection(CORE_DNADEG_CLASSES)) / len(CORE_DNADEG_CLASSES) * 100, 1) if class_set else 0.0,
            **{f"{cls}_Hit_Count": class_counts.get(cls, 0) for cls in CLASS_PRIORITY},
        })
    return pd.DataFrame(sample_rows)


def run_dna_degradation_detail_summary(input_path, output_prefix="dna_degradation_detail_summary_final_core"):
    full_df = collect_input_tables(input_path)
    col_map = standardize_columns(full_df)
    require_minimum_columns(col_map)
    print(f">>> 开始进行 DNA 降解酶 final core 高置信主结果统计，样本数: {full_df['fileName'].nunique()}")
    print(">>> 当前版本默认删除：RNA/tRNA/rRNA、DNA repair/recombination-only、TatD/YcfH/YjjV ambiguous、polymerase、broad hydrolase/phosphoesterase、methylase/specificity-only、CRISPR/transposase、普通非DNA酶等。")

    cls = full_df.apply(lambda r: classify_row(r, col_map), axis=1)
    df_out = pd.concat([full_df, cls], axis=1)
    df_out["feature_abundance"] = df_out.apply(lambda r: get_abundance(r, col_map), axis=1)
    df_out.to_csv(f"{output_prefix}_detailed.csv", index=False, encoding="utf-8-sig")

    df_hits = df_out[df_out["has_DNA_degradation"] == 1].copy()
    df_hits.to_csv(f"{output_prefix}_DNA_degradation_hits_only.csv", index=False, encoding="utf-8-sig")

    df_sample = build_sample_summary(df_out, col_map)
    df_sample.to_csv(f"{output_prefix}_sample_summary.csv", index=False, encoding="utf-8-sig")

    if "Group" in df_sample.columns:
        df_group = (
            df_sample.groupby("Group", dropna=False)
            .agg(
                Sample_N=("Sample", "nunique"),
                Mean_DNA_Degradation_Hit_Count=("DNA_Degradation_Hit_Count", "mean"),
                Mean_DNA_Degradation_Abundance=("DNA_Degradation_Abundance", "mean"),
                Mean_DNA_Degradation_Risk_Score=("Mean_DNA_Degradation_Risk_Score", "mean"),
                Mean_Category_Coverage=("Category_Coverage_%", "mean"),
                Mean_Core_DNA_Degradation_Coverage=("Core_DNA_Degradation_Coverage_%", "mean"),
            )
            .reset_index()
        )
        df_group.to_csv(f"{output_prefix}_group_summary.csv", index=False, encoding="utf-8-sig")

    if not df_hits.empty:
        pd.crosstab(df_hits["fileName"], df_hits["dna_deg_class"]).to_csv(f"{output_prefix}_pivot_main_class.csv", encoding="utf-8-sig")
        pd.crosstab(df_hits["fileName"], df_hits["dna_deg_category"]).to_csv(f"{output_prefix}_pivot_main_category.csv", encoding="utf-8-sig")

    expanded_rows = []
    for _, row in df_hits.iterrows():
        for cls_name in str(row["dna_deg_all_classes"]).split(";"):
            cls_name = cls_name.strip()
            if cls_name:
                expanded_rows.append({
                    "fileName": row["fileName"],
                    "dna_deg_class_expanded": cls_name,
                    "feature_abundance": row["feature_abundance"],
                })
    if expanded_rows:
        df_exp = pd.DataFrame(expanded_rows)
        pd.crosstab(df_exp["fileName"], df_exp["dna_deg_class_expanded"]).to_csv(f"{output_prefix}_pivot_expanded_class.csv", encoding="utf-8-sig")
        (
            df_exp.groupby(["fileName", "dna_deg_class_expanded"], as_index=False)["feature_abundance"]
            .sum()
            .pivot(index="fileName", columns="dna_deg_class_expanded", values="feature_abundance")
            .fillna(0)
            .to_csv(f"{output_prefix}_pivot_expanded_abundance.csv", encoding="utf-8-sig")
        )

    with open(f"{output_prefix}_sample_summary.json", "w", encoding="utf-8") as f:
        json.dump(df_sample.to_dict(orient="records"), f, indent=2, ensure_ascii=False)

    print(f"✅ 完成，输出前缀: {output_prefix}")
    print(f">>> 总基因/功能行数: {len(df_out)}")
    print(f">>> 高置信 DNA 降解酶命中数: {len(df_hits)}")


def main():
    parser = argparse.ArgumentParser(description="DNA 降解酶 final core 高置信主结果筛选汇总脚本。")
    parser.add_argument("input_path", help="输入注释文件或目录")
    parser.add_argument("output_prefix", nargs="?", default="dna_degradation_detail_summary_final_core", help="输出文件前缀")
    args = parser.parse_args()
    run_dna_degradation_detail_summary(args.input_path, args.output_prefix)


if __name__ == "__main__":
    main()
