#!/usr/bin/env bash
###############################################################################
# Script name : add_species_prefix.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 23 September 2025
#
# Purpose
# -------
# Add a species-specific prefix to FASTA headers prior to orthology inference.
# Each sequence header is prefixed with the species name inferred from the
# input filename, allowing faster identity assignment when plotting/analysis
#
# Expected input filenames:
#   <Species>_longest.fa
#
# Example:
#   Padleri_longest.fa  →  prefixed_Padleri.fa
#   >XP_001234567       →  >Padleri|XP_001234567
#
# Usage
# -----
# bash add_species_prefix.sh *.fa
#
# Requirements
# ------------
# bash, awk
#
###############################################################################

for f in "$@"; do
    species=$(basename "$f" _longest.fa)
    out="prefixed_${species}.fa"

    echo "[INFO] Processing $f -> $out"

    awk -v sp="$species" '
        /^>/ {
            sub(/^>/, ">" sp "|", $0)
        }
        { print }
    ' "$f" > "$out"
done
