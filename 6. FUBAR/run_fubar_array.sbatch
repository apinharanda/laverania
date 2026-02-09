#!/usr/bin/env bash
###############################################################################
# run_fubar_array.sbatch
#
# Purpose
# -------
# Run HyPhy FUBAR independently for each orthogroup (OG) using a SLURM array.
# Each array task processes exactly ONE orthogroup.
#
# Inputs (per orthogroup directory: OG_paml/OGxxxx/)
# --------------------------------------------------
#   - OGxxxx.codon.fa     Codon alignment (FASTA, in-frame, gaps allowed)
#   - OGxxxx_tree.txt    Newick gene tree with branch lengths
#
# Outputs (written to the SAME OG directory)
# ------------------------------------------
#   - OGxxxx.FUBAR.json  HyPhy FUBAR result (JSON)
#
# Logs
# ----
#   logs/OGxxxx.FUBAR.log        stdout/stderr from HyPhy
#   logs/fubar_status.tsv       One-line summary per OG (appended)
#
# Status file columns
# -------------------
# timestamp  
# OG          Orthogroup ID
# status      OK | FAIL
# reason      short failure reason or "ran"
# hyphy_exit  HyPhy exit code (NA if not run)
#
# Dependencies
# ------------
# - HyPhy >= 2.5 available as: `hyphy`
# - SLURM workload manager
#
# Assumptions
# -----------
# - OG list exists at: OG_paml/og_dirs.txt (one OG per line)
# - Tree and alignment names match OG directory name
# - No attempt is made to "fix" bad inputs
# - If inputs are missing or HyPhy fails, OG is marked FAIL
#
###############################################################################

#SBATCH --job-name=fubar_OG
#SBATCH --time=36:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=36G
#SBATCH --array=1-4987


# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
BASE_DIR="OG_paml"
OG_LIST="${BASE_DIR}/og_dirs.txt"
LOG_DIR="logs"
STATUS_FILE="${LOG_DIR}/fubar_status.tsv"

mkdir -p "${LOG_DIR}"

# ------------------------------------------------------------------
# Get orthogroup for this array task
# ------------------------------------------------------------------
OG="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OG_LIST}")"
[[ -z "${OG}" ]] && exit 0

OG_DIR="${BASE_DIR}/${OG}"
ALN="${OG_DIR}/${OG}.codon.fa"
TREE="${OG_DIR}/${OG}_tree.txt"
JSON="${OG_DIR}/${OG}.FUBAR.json"
LOG="${LOG_DIR}/${OG}.FUBAR.log"

# ------------------------------------------------------------------
# Initialize status file header (safe if repeated)
# ------------------------------------------------------------------
if [[ ! -s "${STATUS_FILE}" ]]; then
  printf "timestamp\tOG\tstatus\treason\thyphy_exit\n" > "${STATUS_FILE}"
fi

timestamp="$(date -Is)"

# ------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------
if [[ ! -s "${ALN}" ]]; then
  printf "%s\t%s\tFAIL\tmissing_alignment\tNA\n" \
    "${timestamp}" "${OG}" >> "${STATUS_FILE}"
  exit 0
fi

if [[ ! -s "${TREE}" ]]; then
  printf "%s\t%s\tFAIL\tmissing_tree\tNA\n" \
    "${timestamp}" "${OG}" >> "${STATUS_FILE}"
  exit 0
fi

# ------------------------------------------------------------------
# Run HyPhy FUBAR
# ------------------------------------------------------------------
set +e
hyphy fubar \
  --alignment "${ALN}" \
  --tree "${TREE}" \
  --output "${JSON}" \
  > "${LOG}" 2>&1
hyphy_exit=$?
set -e

# ------------------------------------------------------------------
# Check results
# ------------------------------------------------------------------
if [[ "${hyphy_exit}" -ne 0 ]]; then
  printf "%s\t%s\tFAIL\thyphy_exit\t%s\n" \
    "${timestamp}" "${OG}" "${hyphy_exit}" >> "${STATUS_FILE}"
  exit 0
fi

if [[ ! -s "${JSON}" ]]; then
  printf "%s\t%s\tFAIL\tjson_missing\t%s\n" \
    "${timestamp}" "${OG}" "${hyphy_exit}" >> "${STATUS_FILE}"
  exit 0
fi

# ------------------------------------------------------------------
# Success
# ------------------------------------------------------------------
printf "%s\t%s\tOK\tran\t%s\n" \
  "${timestamp}" "${OG}" "${hyphy_exit}" >> "${STATUS_FILE}"
