#!/bin/bash
#SBATCH -J merge_fel
#SBATCH -t 01:00:00
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH -o logs/merge_fel_%j.out
#SBATCH -e logs/merge_fel_%j.err

###############################################################################
# Merge per-OG FEL TSVs into ONE master TSV
#
# INPUT (per OG):
#   OGxxxxxx/OGxxxxxx.FEL.tsv
#     (produced by parse_fel.py)
#
# OUTPUT:
#   FEL_MASTER_27Feb2026.tsv   (written in BASE)
#
# OUTPUT FORMAT:
#   Adds OG as first column:
#     OG  codon  partition  alpha  beta  LRT  p_value  dnds_mle  dnds_lb  dnds_ub  selection  dnds_ci
#
# HOW TO RUN:
#   cd /.../Results_Oct16_6/OG_paml
#   mkdir -p logs
#   sbatch merge_fel_master_simple_27Feb2026.sbatch
###############################################################################

BASE="/n/holylabs/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix/OrthoFinder/Results_Oct16_6/OG_paml"
OUT="${BASE}/FEL_MASTER_27Feb2026.tsv"

mkdir -p "${BASE}/logs"

echo "[INFO] Start: $(date)"
echo "[INFO] BASE:  ${BASE}"
echo "[INFO] OUT:   ${OUT}"
echo

# fresh output every time (header matches your parse_fel.py output + OG prefix)
echo -e "OG\tcodon\tpartition\talpha\tbeta\tLRT\tp_value\tdnds_mle\tdnds_lb\tdnds_ub\tselection\tdnds_ci" > "${OUT}"

n_scanned=0
n_added=0
n_missing=0

for OGDIR in "${BASE}"/OG*; do
  [[ -d "${OGDIR}" ]] || continue
  OG="$(basename "${OGDIR}")"
  TSV="${OGDIR}/${OG}.FEL.tsv"

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
