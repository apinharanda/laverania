#!/usr/bin/env bash
###############################################################################
# parse_all_fel_json_to_tsv.sh
#
# Purpose
# -------
# Batch-parse HyPhy FEL JSON outputs across OG_paml/OG*/ directories using
# parse_fel.py, producing one TSV per OG plus a master summary log.
#
# Assumptions (per OG directory)
# ------------------------------
#   OG_paml/OGXXXX/OGXXXX.FEL.json      input (must exist and be non-empty)
# Output (per OG directory)
# -------------------------
#   OG_paml/OGXXXX/OGXXXX.FEL.tsv       output (created/overwritten)
# Logs
# ----
#   logs/fel_parse_status.tsv           one row per OG with status and counts
#
# Usage
# -----
#   cd /n/holylabs/.../Results_Oct16_6
#   bash parse_all_fel_json_to_tsv.sh
#
# Notes
# -----
# - Codon count is estimated from the TSV rows (excluding header).
# - Failure reasons are kept short and practical.
###############################################################################

# ----------------------------
# paths
# ----------------------------
BASE_DIR="OG_paml"
PARSER="./parse_fel.py"

LOG_DIR="logs"
STATUS_FILE="${LOG_DIR}/fel_parse_status.tsv"

mkdir -p "${LOG_DIR}"

# ----------------------------
# Check parser exists
# ----------------------------
if [[ ! -x "${PARSER}" && ! -f "${PARSER}" ]]; then
  echo "[ERROR] Cannot find parser at: ${PARSER}"
  echo "Run this script from Results_Oct16_6 (so ./parse_fel.py exists), or edit PARSER."
  exit 1
fi

# ----------------------------
# Initialize status file (safe if repeated)
# ----------------------------
if [[ ! -s "${STATUS_FILE}" ]]; then
  printf "timestamp\tOG\tstatus\treason\tfel_json\ttsv_out\tn_codon_rows\tparser_exit\n" > "${STATUS_FILE}"
fi

timestamp="$(date -Is)"

# ----------------------------
# Iterate OG directories
# ----------------------------
shopt -s nullglob
OG_DIRS=( "${BASE_DIR}"/OG* )

if [[ ${#OG_DIRS[@]} -eq 0 ]]; then
  echo "[ERROR] No OG directories found under: ${BASE_DIR}/OG*"
  exit 1
fi

for d in "${OG_DIRS[@]}"; do
  OG="$(basename "$d")"
  JSON="${d}/${OG}.FEL.json"
  OUT="${d}/${OG}.FEL.tsv"

  # ---- Input checks
  if [[ ! -s "${JSON}" ]]; then
    printf "%s\t%s\tFAIL\tmissing_or_empty_json\t%s\t%s\tNA\tNA\n" \
      "${timestamp}" "${OG}" "${JSON}" "${OUT}" >> "${STATUS_FILE}"
    continue
  fi

  # ---- Run parser
  set +e
  python3 "${PARSER}" "${JSON}" "${OUT}" > "${LOG_DIR}/${OG}.parse_fel.log" 2>&1
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    # keep reason short; detailed error is in logs/OGxxxx.parse_fel.log
    printf "%s\t%s\tFAIL\tparser_error\t%s\t%s\tNA\t%s\n" \
      "${timestamp}" "${OG}" "${JSON}" "${OUT}" "${rc}" >> "${STATUS_FILE}"
    continue
  fi

  if [[ ! -s "${OUT}" ]]; then
    printf "%s\t%s\tFAIL\ttsv_missing_or_empty\t%s\t%s\tNA\t%s\n" \
      "${timestamp}" "${OG}" "${JSON}" "${OUT}" "${rc}" >> "${STATUS_FILE}"
    continue
  fi

  # ---- Count codon rows (TSV lines minus header)
  n_lines=$(wc -l < "${OUT}" | tr -d ' ')
  if [[ "${n_lines}" -ge 1 ]]; then
    n_codons=$(( n_lines - 1 ))
  else
    n_codons=0
  fi

  printf "%s\t%s\tOK\tran\t%s\t%s\t%s\t%s\n" \
    "${timestamp}" "${OG}" "${JSON}" "${OUT}" "${n_codons}" "${rc}" >> "${STATUS_FILE}"
done

echo "[OK] Done. Status log: ${STATUS_FILE}"
echo "[OK] Per-OG parser logs: ${LOG_DIR}/OG*.parse_fel.log"
