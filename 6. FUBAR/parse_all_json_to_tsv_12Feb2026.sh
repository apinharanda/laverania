#!/bin/bash
#SBATCH -J parse_fubar
#SBATCH -t 48:00:00
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH -o logs/parse_fubar_%j.out
#SBATCH -e logs/parse_fubar_%j.err

###############################################################################
# SCRIPT: parse_all_json_to_tsv_12Feb2026.sh
#
# WHAT THIS DOES
#   Batch-parses HyPhy FUBAR JSON results (one JSON per orthogroup directory)
#   into a tab-separated table (TSV) with one row per codon.
#
# WHY THIS EXISTS / WHERE IT FITS IN THE PIPELINE
#   Typical workflow for each orthogroup OGxxxx:
#     (1) Build a codon alignment (e.g. OG0003384.codon.fa)
#     (2) Run HyPhy FUBAR -> produces OG0003384.FUBAR.json
#     (3) THIS SCRIPT: convert JSON -> FUBAR_per_codon_withGaps_realcoords.tsv
#     (4) Downstream plotting/summary scripts read the TSV (faster than JSON).
#
# IMPORTANT DEFINITIONS (COORDINATE SYSTEMS)
#   This parser writes:
#     - codon    : 1..N codon index in the alignment
#     - full_pos : set equal to codon (same index)
#
#   "full_pos == codon" means:
#     we are *not* mapping to Pf3D7 genomic coordinates here.
#     we are also *not* doing any gap->real coordinate conversion here.
#
#
# EXPECTED DIRECTORY STRUCTURE (Feb 2026 run)
#   BASE = /.../Results_Oct16_6/OG_paml/
#     OG0003384/
#       OG0003384.FUBAR.json
#       (output will be written here)
#     OG0004504/
#       OG0004504.FUBAR.json
#       ...
#
# OUTPUT (per OG dir)
#   FUBAR_per_codon_withGaps_realcoords.tsv
#
# SKIP LOGIC (IMPORTANT)
#   - If OUT TSV already exists and is non-empty, we skip that OG.
#   - If JSON missing or empty, we skip that OG.
#
# WHAT TO CHECK IF "IT RUNS BUT DOES NOTHING"
#   1) Look at the log: logs/parse_fubar_<jobid>.out
#      If everything is [SKIP], it’s usually because JSON name pattern is wrong.
#   2) Confirm JSON files are named EXACTLY: OGxxxxx.FUBAR.json
#      If not, change JSON=... below.
#
# HOW TO RUN
#   cd /.../Results_Oct16_6/OG_paml
#   mkdir -p logs
#   sbatch parse_all_json_to_tsv_12Feb2026.sh
###############################################################################

# --- REQUIRED PATHS ---
# BASE: folder containing OG directories (OG0000001, OG0000002, ...)
BASE="/n/holylabs/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix/OrthoFinder/Results_Oct16_6/OG_paml"

# PARSER: python script that converts one JSON into one TSV
PARSER="${BASE}/parse_fubar_json_to_tsv_withGaps_fullposIsCodon.py"
# ----------------------

# Ensure a logs folder exists at BASE for SLURM stdout/stderr (SBATCH uses logs/...)
mkdir -p "${BASE}/logs"

echo "[INFO] Starting parse job at: $(date)"
echo "[INFO] BASE:   ${BASE}"
echo "[INFO] PARSER: ${PARSER}"
echo "[INFO] Host:   $(hostname)"
echo

# --- FAIL EARLY if the basic setup is wrong ---
if [[ ! -d "${BASE}" ]]; then
  echo "[ERROR] BASE directory does not exist: ${BASE}" >&2
  exit 1
fi

if [[ ! -s "${PARSER}" ]]; then
  echo "[ERROR] Parser script not found or empty: ${PARSER}" >&2
  exit 1
fi

# --- MAIN LOOP over OG directories ---
for OGDIR in "${BASE}"/OG*; do
  [[ -d "${OGDIR}" ]] || continue

  OG="$(basename "${OGDIR}")"

  # IMPORTANT: This assumes JSON file is exactly named like:
  #   OG0003384/OG0003384.FUBAR.json
  JSON="${OGDIR}/${OG}.FUBAR.json"

  # Output TSV name (alignment coords unless parser maps)
  OUT="${OGDIR}/FUBAR_per_codon_withGaps_realcoords.tsv"

  # If TSV exists, skip (prevents overwrite and saves time)
  if [[ -s "${OUT}" ]]; then
    echo "[SKIP] ${OG} -> TSV already exists (${OUT})"
    continue
  fi

  # If JSON missing or empty, skip
  if [[ ! -s "${JSON}" ]]; then
    echo "[SKIP] ${OG} -> missing/empty JSON (${JSON})"
    continue
  fi

  echo "[RUN ] ${OG}"
  echo "       JSON: ${JSON}"
  echo "       OUT : ${OUT}"

  # Run parser. If it fails, warn and continue to next OG.
  if ! python3 "${PARSER}" "${JSON}" "${OUT}"; then
    echo "[WARN] ${OG} -> parse failed; removing partial OUT (if any)."
    rm -f "${OUT}"
    continue
  fi

  # Sanity check output is non-empty
  if [[ ! -s "${OUT}" ]]; then
    echo "[WARN] ${OG} -> parser finished but OUT is empty; removing."
    rm -f "${OUT}"
    continue
  fi

  echo "[DONE] ${OG}"
  echo
done

echo "[INFO] Finished parse job at: $(date)"
