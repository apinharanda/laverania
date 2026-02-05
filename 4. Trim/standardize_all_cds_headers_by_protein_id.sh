#!/usr/bin/env bash
###############################################################################
# Script name : standardize_all_cds_headers_by_protein_id.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 30 September 2025
#
# Purpose
# -------
# Standardize CDS FASTA headers for all Laverania species so that each CDS record
# is keyed by the corresponding protein_id. This is used to ensure PAL2NAL can
# match CDS sequences to protein alignment IDs (e.g. OrthoFinder inputs with
# headers like: >prefixed_Species|PROTEIN_ID).
#
# Output header format
# --------------------
#   >prefixed_<Species>|<protein_id>
#
# Input
# -----
# NCBI CDS FASTA files (e.g. *_cds_from_genomic.fna), where headers contain
# protein_id attributes and/or a "_cds_" token. Uploaded to Zenodo
#
# Output
# ------
# One standardized CDS FASTA per species, written to CDS_std/:
#   CDS_std/prefixed_<Species>.cds.std.fa also in Zenodo
#
# Notes
# -----
# Extraction priority for protein_id:
#   1) protein_id=TOKEN
#   2) protein_id="TOKEN"
#   3) fallback: token after "_cds_"
# Records without a detectable protein_id are skipped.
#
# How to run
# ----------
# Edit the CDS= paths below if needed, then run:
#   bash standardize_all_cds_headers_by_protein_id.sh
#
###############################################################################

mkdir -p CDS_std

# ---------------------------------------------------------------------------
# Inputs (edit paths here if they move)
# ---------------------------------------------------------------------------
CDS_Pfalciparum="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Pfalciparum/GCF_000002765.6_GCA_000002765_cds_from_genomic.fna"
CDS_Ppraefalciparum="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Ppraefalciparum/GCA_900095595.1_PPRFG01_cds_from_genomic.fna"
CDS_Preichenowi="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Preichenowi/GCF_001601855.1_ASM160185v1_cds_from_genomic.fna"
CDS_Pgaboni="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Pgaboni/GCF_001602025.1_ASM160202v1_cds_from_genomic.fna"
CDS_Padleri="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Padleri/GCF_900097015.1_PADLG01_cds_from_genomic.fna"
CDS_Pblacklocki="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Pblacklocki/GCA_900097035.1_PBLACG01_cds_from_genomic.fna"
CDS_Pbillcollinsi="/n/holylabs/LABS/neafsey_lab/Lab/shared_resources/Pbillcollinsi/GCA_900257145.2_Plasmodium_billcollinsi_cds_from_genomic.fna"

# ---------------------------------------------------------------------------
# Body (repeat block per species; intentionally "stupid simple", no functions)
# ---------------------------------------------------------------------------

# Pfalciparum
SPEC="prefixed_Pfalciparum"
CDS="$CDS_Pfalciparum"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"

# Ppraefalciparum
SPEC="prefixed_Ppraefalciparum"
CDS="$CDS_Ppraefalciparum"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"

# Preichenowi
SPEC="prefixed_Preichenowi"
CDS="$CDS_Preichenowi"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"

# Pgaboni
SPEC="prefixed_Pgaboni"
CDS="$CDS_Pgaboni"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"

# Padleri
SPEC="prefixed_Padleri"
CDS="$CDS_Padleri"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"

# Pblacklocki
SPEC="prefixed_Pblacklocki"
CDS="$CDS_Pblacklocki"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"

# Pbillcollinsi
SPEC="prefixed_Pbillcollinsi"
CDS="$CDS_Pbillcollinsi"
OUT="CDS_std/${SPEC}.cds.std.fa"
awk -v spec="$SPEC" '
  /^>/ {
    pid=""
    if (match($0, /protein_id=([^] ]+)/, m)) pid=m[1]
    else if (match($0, /protein_id="([^"]+)"/, m)) pid=m[1]
    if (pid=="") { if (match($0, /_cds_([^_ ]+)/, m2)) pid=m2[1] }
    if (pid=="") { skip=1; next }
    skip=0
    print ">" spec "|" pid
    next
  }
  skip!=1 { print }
' "$CDS" > "$OUT"
echo "$SPEC  input=$(grep -c '^>' "$CDS")  output=$(grep -c '^>' "$OUT")  file=$OUT"
