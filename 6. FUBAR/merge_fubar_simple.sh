#!/bin/bash
#SBATCH -J merge_fubar
#SBATCH -t 01:00:00
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH -o logs/merge_fubar_%j.out
#SBATCH -e logs/merge_fubar_%j.err

###############################################################################
# Merge per-OG FUBAR TSVs into ONE master TSV
#
# INPUT (per OG):
#   OGxxxxxx/FUBAR_per_codon_withGaps_realcoords.tsv
#
# OUTPUT:
#   FUBAR_MASTER_12Feb2026.tsv   (written in BASE)
#
# OUTPUT FORMAT:
#   Adds OG as first column:
#     OG  full_pos  codon  alpha  beta  omega  P_pos  P_neg  BF_pos  beta_minus_alpha
#
# HOW TO RUN:
#   cd /.../Results_Oct16_6/OG_paml
#   mkdir -p logs
#   sbatch merge_fubar_master_simple_12Feb2026.sbatch
###############################################################################

BASE="/n/holylabs/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix/OrthoFinder/Results_Oct16_6/OG_paml"
OUT="${BASE}/FUBAR_MASTER_12Feb2026.tsv"

mkdir -p "${BASE}/logs"

echo "[INFO] Start: $(date)"
echo "[INFO] BASE:  ${BASE}"
echo "[INFO] OUT:   ${OUT}"
echo

# fresh output every time
echo -e "OG\tfull_pos\tcodon\talpha\tbeta\tomega\tP_pos\tP_neg\tBF_pos\tbeta_minus_alpha" > "${OUT}"

n_scanned=0
n_added=0
n_missing=0

for OGDIR in "${BASE}"/OG*; do
  [[ -d "${OGDIR}" ]] || continue
  OG="$(basename "${OGDIR}")"
  TSV="${OGDIR}/FUBAR_per_codon_withGaps_realcoords.tsv"

  ((n_scanned++)) || true

  if [[ ! -s "${TSV}" ]]; then
    ((n_missing++)) || true
    continue
  fi

  # append all data rows, prefix with OG
  # (skip header with tail -n +2)
  tail -n +2 "${TSV}" | sed "s/^/${OG}\t/" >> "${OUT}"
  ((n_added++)) || true
done

echo
echo "[INFO] Done: $(date)"
echo "[INFO] Summary:"
echo "[INFO]   OG dirs scanned : ${n_scanned}"
echo "[INFO]   OGs merged      : ${n_added}"
echo "[INFO]   Missing/empty   : ${n_missing}"
echo "[INFO]   Master lines    : $(wc -l < "${OUT}")  (includes header)"
echo "[INFO] Output: ${OUT}"
