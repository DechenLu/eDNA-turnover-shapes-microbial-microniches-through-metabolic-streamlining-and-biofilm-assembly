# Engineering eDNA Turnover with Core–Shell Hydrogel Microreactors to Mitigate Environmental Antimicrobial Resistance

This repository contains custom analysis scripts associated with the manuscript:

**Engineering eDNA Turnover with Core–Shell Hydrogel Microreactors to Mitigate Environmental Antimicrobial Resistance**

The code provided here supports the bioinformatic analyses described in the manuscript and its supplementary information.

## Repository contents

### DNA degradation-related gene analysis

`dna_degradation_detail_summary.py`

This custom Python script was used to identify and summarize high-confidence genes related to DNA degradation from MAG/genome functional annotation tables.

The analysis classifies DNA degradation-related genes into four major functional categories:

- extracellular DNase
- DNase/nuclease
- restriction endonuclease
- endonuclease/exonuclease

The script uses annotation information from multiple sources, including gene names, functional descriptions, KEGG orthologs, COGs, PFAM domains, and related annotation fields.

To reduce false-positive assignments, annotations primarily associated with functions not directly indicative of DNA degradation are excluded from the main analysis. These include, among others:

- DNA repair- or recombination-specific nucleases
- RNA-, rRNA-, or tRNA-associated nucleases
- DNA polymerase proofreading domains
- nonspecific hydrolases and phosphoesterases
- methylase or restriction-modification specificity subunits without direct cleavage evidence
- CRISPR/Cas-associated nucleases
- transposases and mobile-element-associated nucleases

The resulting high-confidence assignments were used to quantify DNA degradation-related gene categories in the metagenome-assembled genomes analyzed in the study.

## Requirements

The script requires Python 3 and the following Python package:

```bash
pip install pandas openpyxl
```

Most other dependencies are included in the Python standard library.

## Usage

The script can be run using:

```bash
python dna_degradation_detail_summary.py <input_file_or_directory> <output_prefix>
```

For example:

```bash
python dna_degradation_detail_summary.py annotation_results/ dna_degradation_results
```

The input can be either a single annotation file or a directory containing multiple annotation files.

Supported input formats include:

- `.csv`
- `.tsv`
- `.txt`
- `.tsvm`
- `.xlsx`
- `.xls`

When a directory is provided, the script recursively searches supported annotation files and selects the preferred annotation file version within each subdirectory.

## Main outputs

The script generates several output files, including:

- `*_detailed.csv`  
  Complete annotation table with DNA degradation-related classifications.

- `*_DNA_degradation_hits_only.csv`  
  High-confidence DNA degradation-related gene hits retained for downstream analysis.

- `*_sample_summary.csv`  
  Sample-level summary of DNA degradation-related gene categories.

- `*_group_summary.csv`  
  Group-level summary when treatment or group information is available.

- `*_pivot_main_class.csv`  
  Sample-by-gene-class count matrix.

- `*_pivot_main_category.csv`  
  Sample-by-functional-category count matrix.

- `*_pivot_expanded_class.csv`  
  Expanded gene-class count matrix.

- `*_pivot_expanded_abundance.csv`  
  Expanded abundance matrix.

- `*_sample_summary.json`  
  Sample-level summary in JSON format.

## DNA degradation gene categories

The main analysis retains four predefined high-confidence DNA degradation gene classes:

| Gene class | General interpretation |
| --- | --- |
| `extracellular_dnase` | Extracellular DNA degradation proteins and extracellular DNases |
| `dnase_nuclease` | DNases and nucleases with direct DNA cleavage-related annotation |
| `restriction_endonuclease` | Restriction endonucleases and related restriction enzymes |
| `endonuclease_exonuclease` | DNA-specific endonucleases and exonucleases retained under the high-confidence criteria |

The classification rules implemented in the script are intentionally conservative to minimize inclusion of proteins whose primary functions are DNA repair, RNA processing, replication, or unrelated hydrolysis.

## Reproducibility

The script included in this repository corresponds to the custom DNA degradation-related gene screening workflow described in the manuscript.

For reproducibility, the code should be used together with the input functional annotation tables described in the Methods and Supplementary Information of the paper.

## Citation

If you use this code, please cite the associated publication:

**Lu et al. Engineering eDNA Turnover with Core–Shell Hydrogel Microreactors to Mitigate Environmental Antimicrobial Resistance. Trends in Biotechnology.**

The final publication details and DOI will be added after publication.

## Code availability

The custom Python code used for DNA degradation-related gene identification and summarization is provided in this repository.

A DOI-archived version of the code will be deposited through Zenodo for permanent preservation and citation in the published article.

## License

Please contact the authors regarding reuse of the code until a formal software license is added to this repository.
