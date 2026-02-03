#!/usr/bin/env bash
###############################################################################
# Script name : pick_longest_manual_run.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 23 September 2025
#
# Purpose
# -------
# Manual-run variant used for assemblies where GPFF protein accessions are not
# captured by strict underscore-based patterns (e.g. SCQ12544.1).
#
#   1) Build protein→gene map from GPFF (GeneID/locus_tag/gene), emitting
#      accession keys with and without version suffix.
#   2) Compute protein lengths and FASTA order.
#   3) Join and select the longest protein per gene (ties: earliest in FASTA).
#   4) Extract selected sequences to OUT.
#
# Requirements
# ------------
# seqkit, awk, sort, join
#
###############################################################################

usage() {
cat <<'EOF'
Usage:
  pick_longest_manual_run.sh -f <proteins.faa> -g <proteins.gpff> -o <output.fa>

Notes:
  - Requires: seqkit, awk, sort, join
  - Set KEEP_TMP=1 to keep the temporary directory for debugging
EOF
}

FAA=""
GPFF=""
OUT=""

if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--faa)   FAA="$2"; shift 2;;
    -g|--gpff)  GPFF="$2"; shift 2;;
    -o|--out)   OUT="$2"; shift 2;;
    -h|--help)  usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit not found in PATH"; exit 1; }
[[ -r "$FAA"  ]] || { echo "ERROR: Cannot read FAA: $FAA"; exit 1; }
[[ -r "$GPFF" ]] || { echo "ERROR: Cannot read GPFF: $GPFF"; exit 1; }
[[ -n "$OUT"  ]] || { echo "ERROR: Missing -o/--out"; exit 1; }

TMP="$(mktemp -d)"
if [[ "${KEEP_TMP:-0}" != "1" ]]; then
  trap 'rm -rf "$TMP"' EXIT
fi

echo "Input FAA  : $(realpath "$FAA")"
echo "Input GPFF : $(realpath "$GPFF")"
echo "Output FAA : $(realpath "$OUT")"
echo "TMPDIR     : $TMP"

# 1) Build protein → gene map from GPFF; emit keys with and without version
awk -v RS="//" '
function emit(acc,key, base){
  if(!acc||!key) next
  print acc "\t" key
  base=acc; sub(/\.[0-9]+$/,"",base)
  print base "\t" key
}
{
  acc="";
  if (match($0, /\nVERSION[ \t]+([A-Za-z0-9_]+(\.[0-9]+)?)/, m)) acc=m[1];
  else if (match($0, /\nACCESSION[ \t]+([A-Za-z0-9_]+)/, a))    acc=a[1];

  if      (match($0, /db_xref="GeneID:([0-9]+)"/, g)) key="GeneID:" g[1];
  else if (match($0, /locus_tag="([^"]+)"/,       l)) key="locus:"  l[1];
  else if (match($0, /\ngene="([^"]+)"/,          s)) key="gene:"   s[1];

  emit(acc, key)
}' "$GPFF" | sort -u > "$TMP/prot2gene.tsv"

echo "prot2gene rows: $(wc -l < "$TMP/prot2gene.tsv")"

# 2) FASTA order and lengths (ID = header up to first space)
seqkit fx2tab -n -l "$FAA" \
  | awk -F'\t' '{id=$1; sub(/ .*/,"",id); print NR"\t"id"\t"$2}' \
  > "$TMP/prot_len.tsv"

# 3) Join and pick longest per gene; ties → earliest in FASTA
join -t $'\t' -j 1 \
  <(sort -k2,2 "$TMP/prot_len.tsv" | awk -F'\t' '{print $2"\t"$1"\t"$3}') \
  <(sort -k1,1 "$TMP/prot2gene.tsv") \
| awk -F'\t' '{print $1,$3,$2,$4}' OFS='\t' \
| sort -k4,4 -k2,2nr -k3,3n \
| awk -F'\t' '!seen[$4]++{print $1}' > "$TMP/keep.txt"

echo "to keep: $(wc -l < "$TMP/keep.txt") proteins"

# 4) Extract sequences (match header name up to first space)
awk -v keep="$TMP/keep.txt" '
BEGIN{
  while((getline k < keep)>0){ sub(/\r$/,"",k); ids[k]=1 }
  close(keep)
}
function hid(s){ sub(/^>/,"",s); sub(/ .*/,"",s); return s }
BEGINFILE{ hit=0 }
/^>/ { hit = ids[hid($0)] }
/^>/ || hit { if(hit) print }
' "$FAA" > "$OUT"

ls -lh "$OUT"
grep -c '^>' "$OUT"
