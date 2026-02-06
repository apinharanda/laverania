#!/usr/bin/env bash
###############################################################################
# Script name : 01_trim_protein_alignments_trimal_sbatch.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 06 February 2026
#
# Purpose
# -------
# Trim OrthoFinder protein multiple sequence alignments using trimAl.
# Each alignment is processed independently using a SLURM job array
# Trimming is performed at the protein level only (column removal).
# Resulting trimmed protein alignments are used downstream for codon
# back-translation
#
# Inputs
# ------
# Protein multiple sequence alignments (FASTA):
#   MultipleSequenceAlignments/*.fa #uploaded in Zenodo
#
# Outputs
# -------
# Trimmed protein alignments:
#   trimAl/<OG>.trimmed.fa
#
# Log files:
#   trimAl/trimAl_warnings_array.log
#     - trimAl warnings (e.g. all-gap sequences removed)
#
#   trimAl/trimAl_empty_alignments_array.txt
#     - Orthogroups that became empty after trimming
#
# SLURM
# -----
# One orthogroup alignment per array task:
#   --array=1-5079
#
# Requirements
# ------------
# trimAl
#
# How to run
# ----------
# sbatch 01_trim_protein_alignments_trimal_sbatch.sh
#
###############################################################################

#SBATCH --job-name=trimAl
#SBATCH --time=06:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH --array=1-5079

PROT_DIR="MultipleSequenceAlignments"
TRIM_DIR="trimAl"
mkdir -p "$TRIM_DIR"

WARN_LOG="$TRIM_DIR/trimAl_warnings_array.log"
EMPTY_LIST="$TRIM_DIR/trimAl_empty_alignments_array.txt"

# one alignment per array task
aln=$(ls -1 "$PROT_DIR"/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")

# if task id exceeds number of files, exit quietly
if [[ -z "$aln" ]]; then
  exit 0
fi

OG="$(basename "$aln" .fa)"
out="$TRIM_DIR/${OG}.trimmed.fa"
tmp="$TRIM_DIR/${OG}.trimmed.tmp.fa"

trimal -in "$aln" -out "$tmp" -automated1 2>>"$WARN_LOG"

if [[ ! -s "$tmp" || $(grep -c '^>' "$tmp") -eq 0 ]]; then
  echo "$OG" >> "$EMPTY_LIST"
  rm -f "$tmp"
  exit 0
fi

mv -f "$tmp" "$out"
