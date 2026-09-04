# eDNA turnover shapes microbial microniches through metabolic streamlining and biofilm assembly

This repository contains the bioinformatic and ecological-analysis code associated with the study **“eDNA turnover shapes microbial microniches through metabolic streamlining and biofilm assembly.”** It provides reproducible workflows for microbial-community analysis, metagenomic read processing, reconstruction and quality filtering of metagenome-assembled genomes (MAGs), species-level clustering, taxonomic classification and functional annotation.

## Overview

Extracellular DNA (eDNA) is both a nutrient-rich component of microbial necromass and a structural constituent of biofilms. This study investigates how localized eDNA turnover by the marine bacterium *Marinifilum nucleolyticum* D737T influences nutrient availability, metabolic interactions, community composition and biofilm organization within sediment-associated microbial microniches.

The study combines an *in situ* microencapsulation system with isolate experiments, synthetic-community assays, genome-scale metabolic modeling, amplicon profiling and metagenomic analysis. Local MAGs recovered from microcapsule surfaces and surrounding sediments were used to evaluate treatment-associated changes in predicted nucleotide-biosynthetic potential. A quality-controlled comparison of 3,436 MAGs from 15 selected habitat categories was used to place the local observations in a broader environmental context, rather than to represent a comprehensive global survey.

## Study framework

The principal analytical components are:

1. **Localized eDNA turnover:** experimental deployment of D737T within core–shell microcapsules in marine-sediment communities.
2. **Community responses:** alpha diversity, beta diversity, taxonomic shifts, niche breadth and microbial association networks.
3. **Metabolic interactions:** isolate and synthetic-community experiments combined with genome-scale metabolic modeling.
4. **Local metagenomic profiling:** reconstruction of MAGs from microcapsule surfaces and surrounding sediments.
5. **Nucleotide-biosynthetic potential:** evaluation of *de novo* dATP, dGTP, dCTP and dTTP pathway completeness in quality-filtered MAGs.
6. **Selected-habitat comparison:** environmental contextualization using 3,436 high-quality MAGs from 15 habitat categories.

## Repository structure

```text
.
├── README.md
├── GITHUB_ABOUT.md
├── CHANGELOG.md
├── CITATION.cff.template
├── LICENSE_GUIDANCE.md
├── config/
│   └── config.example.env
├── metadata/
│   ├── samples.example.tsv
│   └── custom_taxonomy.example.tsv
├── scripts/
│   ├── 00_check_setup.sh
│   ├── 01_quality_filter.sh
│   ├── 02_assemble.sh
│   ├── 03_bin.sh
│   ├── 04_refine_bins.sh
│   ├── 05_collect_refined_bins.sh
│   ├── 06_checkm.sh
│   ├── 07_filter_high_quality_mags.py
│   ├── 08_drep_ani95.sh
│   ├── 08_export_ani95_clusters.py
│   ├── 09_prokka.sh
│   ├── 10_eggnog_mapper.sh
│   ├── 11_gtdbtk_classify.sh
│   ├── 12_gtdbtk_tree.sh
│   └── 13_align_16s_mafft.sh
├── analysis/
│   └── README.md
├── docs/
│   ├── methods_mapping.md
│   ├── original_command_provenance.md
│   ├── release_checklist.md
│   └── code_availability_template.md
├── envs/
│   └── software_versions.tsv
├── tests/
│   └── test_filter_high_quality_mags.py
├── diversity.R
├── β-diversity.R
├── Levins' B C-T.R
├── SpiecEasi_network.R
└── plot_network.R
```

Large sequencing files, intermediate assemblies, external databases and study-specific source-data tables are not included in the repository.

## Community and ecological analyses

The root-level R scripts provide the following analyses:

| Script | Purpose |
|---|---|
| `diversity.R` | Alpha-diversity calculations and mixed-model comparisons |
| `β-diversity.R` | Bray–Curtis beta diversity, ordination and PERMANOVA |
| `Levins' B C-T.R` | Levins’ niche-breadth analysis |
| `SpiecEasi_network.R` | SPIEC-EASI microbial association-network inference |
| `plot_network.R` | Network visualization and export |

Representative R packages include `vegan`, `tidyverse`, `data.table`, `lme4`, `lmerTest`, `emmeans`, `glmmTMB`, `igraph` and `SpiecEasi`. Input filenames and group variables are documented at the beginning of each script and should be adjusted to match the corresponding source-data tables.

## MAG reconstruction and annotation workflow

The metagenomic workflow comprises:

1. Adapter removal and paired-end read-quality filtering with BBDuk/BBTools.
2. Per-sample assembly with MEGAHIT through metaWRAP.
3. Independent genome binning with MetaBAT2, MaxBin2 and CONCOCT.
4. Consolidation of bins using `metaWRAP bin_refinement`.
5. Genome-quality estimation with CheckM.
6. Retention of MAGs with completeness ≥90% and contamination ≤10%.
7. Genome dereplication and species-level clustering at 95% secondary ANI using dRep.
8. Export of the dRep `secondary_cluster` assignments as ANI95 cluster IDs and linkage to the final high-quality MAG subset.
9. Genome annotation with Prokka and eggNOG-mapper.
10. Taxonomic classification and phylogenomic reconstruction with GTDB-Tk.
11. Optional alignment of 16S rRNA sequences with MAFFT.

## Quick start

Clone the repository and prepare the configuration and sample-manifest files:

```bash
git clone https://github.com/DechenLu/eDNA-turnover-governs-microbial-niche-construction.git
cd eDNA-turnover-governs-microbial-niche-construction

cp config/config.example.env config/config.env
cp metadata/samples.example.tsv metadata/samples.tsv
```

Edit `config/config.env` to specify computational resources and database paths. Edit `metadata/samples.tsv` so that each row contains a sample identifier and its paired FASTQ files.

Validate the configuration:

```bash
bash scripts/00_check_setup.sh config/config.env
```

Run the complete core workflow:

```bash
bash run_pipeline.sh config/config.env
```

The pipeline can also be executed one stage at a time:

```bash
bash scripts/01_quality_filter.sh config/config.env
bash scripts/02_assemble.sh config/config.env
bash scripts/03_bin.sh config/config.env
bash scripts/04_refine_bins.sh config/config.env
bash scripts/05_collect_refined_bins.sh config/config.env
bash scripts/06_checkm.sh config/config.env

python scripts/07_filter_high_quality_mags.py \
  --checkm-table results/06_checkm/checkm_qa.tsv \
  --genome-dir results/05_refined_bins_combined \
  --output-dir results/07_high_quality_mags \
  --completeness 90 \
  --contamination 10 \
  --link

bash scripts/08_drep_ani95.sh config/config.env
python scripts/08_export_ani95_clusters.py \
  --cdb results/08_drep_ani95/data_tables/Cdb.csv \
  --hq-table results/07_high_quality_mags/high_quality_mags.tsv \
  --output results/08_drep_ani95/ani95_cluster_assignments.tsv
bash scripts/09_prokka.sh config/config.env
bash scripts/10_eggnog_mapper.sh config/config.env
bash scripts/11_gtdbtk_classify.sh config/config.env
```

Phylogenomic reconstruction requires a study-specific taxonomy file and is run separately:

```bash
cp metadata/custom_taxonomy.example.tsv metadata/custom_taxonomy.tsv
bash scripts/12_gtdbtk_tree.sh config/config.env
```

## Input data

### Metagenomic sample manifest

`metadata/samples.tsv` must be tab-separated and contain the following columns:

| Column | Description |
|---|---|
| `sample_id` | Unique filesystem-safe sample identifier |
| `read1` | Path to the forward FASTQ or FASTQ.GZ file |
| `read2` | Path to the reverse FASTQ or FASTQ.GZ file |

### Community-analysis inputs

The R scripts generally require an OTU/ASV abundance table, sample metadata and, where applicable, taxonomy and phylogenetic-tree files. Required filenames and columns are specified at the beginning of each script.

### External databases

GTDB-Tk and eggNOG-mapper databases are not distributed with this repository. Their local paths must be defined using `GTDBTK_DATA_PATH` and `EGGNOG_DATA_DIR` in `config/config.env`.

## Key parameters

The supplied workflow uses the following principal settings:

- BBDuk: `ktrim=r`, `k=28`, `mink=12`, `hdist=1`, `tbo=t`, `tpe=t`, `qtrim=rl`, `trimq=20`, `minlen=100`.
- Assembly: MEGAHIT through `metaWRAP assembly --megahit`.
- Binning: MetaBAT2, MaxBin2 and CONCOCT; minimum contig length 500 bp.
- Bin refinement: minimum completeness 50% and maximum contamination 10%.
- Final MAG filtering: CheckM completeness ≥90% and contamination ≤10%.
- ANI dereplication and species-level clustering: dRep `-comp 50 -con 10 -sa 0.95` on the refined MAG collection.
- Final analytical MAG subset: CheckM completeness ≥90% and contamination ≤10%, linked to the dRep cluster assignments by genome identifier.
- eggNOG-mapper: DIAMOND mode with a seed-ortholog E-value threshold of 1 × 10⁻⁵.
- GTDB-Tk: classification against GTDB release 214, subject to final version verification.

The 50% completeness threshold used during bin refinement and dRep clustering is distinct from the final ≥90% completeness criterion used to retain MAGs for the principal downstream analyses. The 95%-ANI cluster identifier corresponds to the dRep `secondary_cluster` field in `data_tables/Cdb.csv`.

## Software

The principal command-line tools are BBTools/BBMap, MEGAHIT, metaWRAP, MetaBAT2, MaxBin2, CONCOCT, CheckM, dRep, Prokka, eggNOG-mapper, GTDB-Tk, FastTree, MAFFT and GNU Parallel. Reported and command-record versions are summarized in `envs/software_versions.tsv`.

## Validation

Shell syntax and the high-quality-MAG filtering utility can be checked without downloading the external bioinformatic databases:

```bash
find scripts -name '*.sh' -print0 | xargs -0 -n1 bash -n
python -m unittest discover -s tests -v
```

## Reproducibility notes

- Machine-specific paths and usernames have been replaced with configuration variables.
- MEGAHIT has been confirmed as the assembler used for the final analysis.
- The final MAG filter uses inclusive thresholds: completeness ≥90% and contamination ≤10%.
- dRep dereplication and clustering use the supplied settings `-comp 50 -con 10 -sa 0.95`; ANI95 cluster assignments are exported from `Cdb.csv`.
- The exact GTDB-Tk version should be confirmed from the final analysis environment before the repository is archived.
- Custom pathway-completeness, statistical and figure-generation scripts should be placed under `analysis/`. These files were not included in the supplied command archive and should be added before stating that all figure-generation code is publicly available.

See `docs/release_checklist.md` for the remaining checks.

## Data availability

Raw sequencing data are not stored in this GitHub repository. Replace the following placeholders before release:

- Sequence Read Archive/BioProject accession: `[TO BE ADDED]`
- Processed source-data repository: `[TO BE ADDED]`
- Archived code release DOI: `[TO BE ADDED]`

## Citation

Associated manuscript:

> Yu, W.-X., Lu, D.-C. *et al.* **eDNA turnover shapes microbial microniches through metabolic streamlining and biofilm assembly.**

After publication, add the journal, year, DOI and complete citation. Complete `CITATION.cff.template`, rename it to `CITATION.cff`, and archive the final tagged release with Zenodo or an equivalent repository.

## Code availability

The scripts used for microbial-community analysis, metagenomic processing, MAG reconstruction, genome-quality filtering, species-level clustering, taxonomic classification and functional annotation are provided in this repository. A manuscript-ready statement is available in `docs/code_availability_template.md`.

## License

No software license has been assigned automatically because license selection may require author or institutional approval. See `LICENSE_GUIDANCE.md` before public release.
