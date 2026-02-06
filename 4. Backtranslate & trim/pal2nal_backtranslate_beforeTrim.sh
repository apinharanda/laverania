#!/usr/bin/env bash
###############################################################################
# Script name : pal2nal_backtranslate_beforeTrim.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 06 February 2026
#
# Purpose
# -------
# Back-translate OrthoFinder protein alignments to codon alignments using PAL2NAL.
#
# IMPORTANT
# ---------
# This script runs PAL2NAL on the UNTRIMMED protein alignments from OrthoFinder
# (MultipleSequenceAlignments/*.fa). Trimming happens later on the codon
# alignments (and is codon-aware)
#
# Inputs
# ------
# Protein MSAs (FASTA):
#   MultipleSequenceAlignments/*.fa
#
# Combined standardized CDS FASTA (headers must match protein IDs):
#   all_species.cds.std.fa   (set via CDS_ALL below)
#
# Outputs
# -------
# Codon alignments (FASTA):
#   Codons_raw/<OG>.codon.fa
#
# Per-OG logs:
#   logs/pal2nal_<OG>.log
#   logs/pal2nal_<OG>.status.tsv
#
# SLURM
# -----
# Run as an array job where each task processes one orthogroup alignment.
#
# How to run (test)
# -----------------
# sbatch --array=1-20 03a_pal2nal_backtranslate_array.sh
#
# How to run (full; update array max to your file count)
# ------------------------------------------------------
# sbatch --array=1-5079 03a_pal2nal_backtranslate_array.sh
#
# Requirements
# ------------
# pal2nal.pl
# seqkit
#
###############################################################################

#SBATCH --job-name=pal2nal
#SBATCH --time=06:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --array=1-5079

PROT_DIR="MultipleSequenceAlignments"
OUT_DIR="Codons_raw"
LOG_DIR="logs"

# Path to your combined CDS FASTA (edit if needed)
CDS_ALL="/n/holylabs/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix/OrthoFinder/Results_Sep23_1/all_species.cds.std.fa"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Pick the Nth alignment for this array task
aln=$(ls -1 "$PROT_DIR"/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [[ -z "$aln" ]]; then
  exit 0
fi

OG="$(basename "$aln" .fa)"

ids="${OG}.ids"
cds="${OG}.cds.fa"
out="${OUT_DIR}/${OG}.codon.fa"
log="${LOG_DIR}/pal2nal_${OG}.log"
status="${LOG_DIR}/pal2nal_${OG}.status.tsv"

# Make a reproducible start line in the log
ts="$(date -Is)"
echo "START  ${ts}  OG=${OG}  input=${aln}" > "$log"
echo "CDS_ALL=${CDS_ALL}" >> "$log"

# 1) IDs from the protein alignment (match headers exactly, without the '>')
grep '^>' "$aln" | sed 's/^>//' > "$ids"
ids_n=$(wc -l < "$ids")
echo "IDS_N=${ids_n}  ids_file=${ids}" >> "$log"

# 2) Extract matching CDS records
seqkit grep -f "$ids" "$CDS_ALL" > "$cds"
cds_n=$(grep -c '^>' "$cds" 2>/dev/null || echo 0)
echo "CDS_HEADERS_N=${cds_n}  cds_file=${cds}" >> "$log"

# 3) PAL2NAL (FASTA output; keep gaps; standard genetic code)
#    -output fasta      : FASTA codon alignment output
#    -codontable 1      : standard genetic code
pal2nal.pl "$aln" "$cds" -output fasta -codontable 1 > "$out" 2>> "$log"
pal_exit=$?

# 4) Basic success criteria: file exists, has at least 1 header, and codon-length sanity
out_headers=$(grep -c '^>' "$out" 2>/dev/null || echo 0)
out_bytes=$(wc -c < "$out" 2>/dev/null || echo 0)

# Check that every sequence length (excluding gaps) is a multiple of 3
# If output is empty, codon_ok becomes 0
codon_ok=1
if [[ "$out_headers" -eq 0 ]]; then
  codon_ok=0
else
  modset=$(awk '
    /^>/ {if(seq!=""){gsub(/-/,"",seq); print length(seq)%3}; seq=""; next}
    {gsub(/[ \t\r\n]/,""); seq=seq$0}
    END{if(seq!=""){gsub(/-/,"",seq); print length(seq)%3}}
  ' "$out" | sort -u | tr '\n' ' ')
  # Expect modset to be "0 " (or "0")
  echo "CODON_LEN_MOD3_SET=${modset}" >> "$log"
  if echo "$modset" | grep -qvE '^(0[[:space:]]*)$'; then
    codon_ok=0
  fi
fi

ts_end="$(date -Is)"

if [[ "$pal_exit" -eq 0 && "$out_headers" -gt 0 && "$codon_ok" -eq 1 ]]; then
  echo "END    ${ts_end}  OG=${OG}  STATUS=OK" >> "$log"
  printf "%s\t%s\tOK\t%s\t%s\t%s\t%s\n" "$ts_end" "$OG" "$ids_n" "$cds_n" "$out_headers" "$out" > "$status"
else
  echo "END    ${ts_end}  OG=${OG}  STATUS=FAIL  pal_exit=${pal_exit}  out_headers=${out_headers}  codon_ok=${codon_ok}  bytes=${out_bytes}" >> "$log"
  # Keep a placeholder output path for consistency
  printf "%s\t%s\tFAIL\t%s\t%s\t%s\t%s\n" "$ts_end" "$OG" "$ids_n" "$cds_n" "$out_headers" "$out" > "$status"
fi
