#!/usr/bin/env bash
set -Eeuo pipefail

CONFIG_FILE="${1:-config/config.env}"

bash scripts/00_check_setup.sh "$CONFIG_FILE"
bash scripts/01_quality_filter.sh "$CONFIG_FILE"
bash scripts/02_assemble.sh "$CONFIG_FILE"
bash scripts/03_bin.sh "$CONFIG_FILE"
bash scripts/04_refine_bins.sh "$CONFIG_FILE"
bash scripts/05_collect_refined_bins.sh "$CONFIG_FILE"
bash scripts/06_checkm.sh "$CONFIG_FILE"

# shellcheck disable=SC1090
source "$CONFIG_FILE"
python scripts/07_filter_high_quality_mags.py \
  --checkm-table "${RESULTS_DIR}/06_checkm/checkm_qa.tsv" \
  --genome-dir "${RESULTS_DIR}/05_refined_bins_combined" \
  --output-dir "${RESULTS_DIR}/07_high_quality_mags" \
  --completeness "$HQ_COMPLETENESS" \
  --contamination "$HQ_CONTAMINATION" \
  --link

bash scripts/08_drep_ani95.sh "$CONFIG_FILE"
python scripts/08_export_ani95_clusters.py \
  --cdb "${RESULTS_DIR}/08_drep_ani95/data_tables/Cdb.csv" \
  --hq-table "${RESULTS_DIR}/07_high_quality_mags/high_quality_mags.tsv" \
  --output "${RESULTS_DIR}/08_drep_ani95/ani95_cluster_assignments.tsv"
bash scripts/09_prokka.sh "$CONFIG_FILE"
bash scripts/10_eggnog_mapper.sh "$CONFIG_FILE"
bash scripts/11_gtdbtk_classify.sh "$CONFIG_FILE"
