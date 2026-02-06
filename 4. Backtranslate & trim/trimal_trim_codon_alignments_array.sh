#!/usr/bin/env bash
###############################################################################
# Script name : trimal_trim_codon_alignments_array.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 06 February 2026
#
# Purpose
# -------
# Trim codon alignments produced by PAL2NAL (Codons_raw/*.codon.fa) using trimAl,
# trimming in CODON units (i.e., whole codons). This avoids breaking reading
# frames and is the recommended order: PAL2NAL first, then trim codons.
#
# Inputs
# ------
# Codon alignments (FASTA) from PAL2NAL:
#   Codons_raw/<OG>.codon.fa
#
# Outputs
# -------
# Trimmed codon alignments (FASTA):
#   Codons_trim/<OG>.codon.trim.fa
#
# Logs
# ----
# logs/trimal_codon_<OG>.log
# logs/trimal_codon_<OG>.status.tsv
#
# How to run (test)
# -----------------
# sbatch --array=1-20 03b_trimal_trim_codon_alignments_array.sh
#
# How to run (full; update array max to your file count)
# ------------------------------------------------------
# sbatch --array=1-5079 trimal_trim_codon_alignments_array.sh
#
# Requirements
# ------------
# trimal
#
###############################################################################

#SBATCH --job-name=trimal_codon
#SBATCH --time=06:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --array=1-5079

IN_DIR="Codons_raw"
OUT_DIR="Codons_trim"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Pick the Nth codon alignment for this array task
codon=$(ls -1 "$IN_DIR"/*.codon.fa 2>/dev/null | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [[ -z "$codon" ]]; then
  exit 0
fi

OG="$(basename "$codon" .codon.fa)"
out="${OUT_DIR}/${OG}.codon.trim.fa"
tmp="${OUT_DIR}/${OG}.codon.trim.tmp.fa"
log="${LOG_DIR}/trimal_codon_${OG}.log"
status="${LOG_DIR}/trimal_codon_${OG}.status.tsv"

ts="$(date -Is)"
echo "START  ${ts}  OG=${OG}  input=${codon}" > "$log"

# trimAl options used:
#   -automated1 : trimAl's automated heuristic trimming
#   -codon      : trim in codon units (whole codons)
trimal -in "$codon" -out "$tmp" -automated1 -codon 2>>"$log"
tri_exit=$?

out_headers=$(grep -c '^>' "$tmp" 2>/dev/null || echo 0)

codon_ok=1
if [[ "$out_headers" -eq 0 ]]; then
  codon_ok=0
else
  modset=$(awk '
    /^>/ {if(seq!=""){gsub(/-/,"",seq); print length(seq)%3}; seq=""; next}
    {gsub(/[ \t\r\n]/,""); seq=seq$0}
    END{if(seq!=""){gsub(/-/,"",seq); print length(seq)%3}}
  ' "$tmp" | sort -u | tr '\n' ' ')
  echo "CODON_LEN_MOD3_SET=${modset}" >> "$log"
  if echo "$modset" | grep -qvE '^(0[[:space:]]*)$'; then
    codon_ok=0
  fi
fi

ts_end="$(date -Is)"

if [[ "$tri_exit" -eq 0 && "$out_headers" -gt 0 && "$codon_ok" -eq 1 ]]; then
  mv -f "$tmp" "$out"
  echo "END    ${ts_end}  OG=${OG}  STATUS=OK  out=${out}" >> "$log"
  printf "%s\t%s\tOK\t%s\t%s\n" "$ts_end" "$OG" "$out_headers" "$out" > "$status"
else
  rm -f "$tmp"
  echo "END    ${ts_end}  OG=${OG}  STATUS=FAIL  trimal_exit=${tri_exit}  out_headers=${out_headers}  codon_ok=${codon_ok}" >> "$log"
  printf "%s\t%s\tFAIL\t%s\t%s\n" "$ts_end" "$OG" "$out_headers" "$out" > "$status"
fi
