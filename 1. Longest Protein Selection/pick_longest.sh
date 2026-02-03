#!/usr/bin/env bash
###############################################################################
# Script name : pick_longest.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 23 September 2025
#
# Purpose
# -------
# Select the longest protein isoform per gene using:
#   - a protein FASTA (.faa) for protein sequences and lengths
#   - a protein GenPept flat file (.gpff) to map protein accessions to gene keys
#
# Tie-break rule
# -------------
# If multiple proteins map to the same gene key and have identical lengths,
# keep the protein that appears earliest in the FASTA file.
#
# Requirements
# ------------
# seqkit, awk, sort, join
#
###############################################################################

# pick_longest.sh — longest protein per gene using FAA+GPFF
# Tie-break: if equal length, keep earliest appearance in FASTA.

usage() {
cat <<'EOF'
Usage:
  pick_longest.sh PREFIX
    # Uses: PREFIX_protein.faa and PREFIX_protein.gpff in CWD
  pick_longest.sh -f proteins.faa -g proteins.gpff [-o output.faa]

Notes:
  - Requires: seqkit, awk, sort, join
  - Tie-breaker: earliest appearance in FASTA wins
  - Set KEEP_TMP=1 to keep temp files for debugging
EOF
}

FAA=""; GPFF=""; OUT=""
if [[ $# -eq 0 ]]; then usage; exit 1; fi

if [[ $# -eq 1 && "$1" != "-"* ]]; then
  PREFIX="$1"
  FAA="${PREFIX}_protein.faa"
  GPFF="${PREFIX}_protein.gpff"
  OUT="${PREFIX}.longest_per_gene.faa"
else
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -f|--faa)   FAA="$2"; shift 2;;
      -g|--gpff)  GPFF="$2"; shift 2;;
      -o|--out)   OUT="$2"; shift 2;;
      -h|--help)  usage; exit 0;;
      *) echo "Unknown arg: $1"; usage; exit 1;;
    esac
  done
  [[ -z "${OUT}" ]] && OUT="$(basename "${FAA%.faa}").longest_per_gene.faa"
fi

command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit not found in PATH"; exit 1; }
[[ -r "$FAA"  ]] || { echo "ERROR: Cannot read FAA: $FAA"; exit 1; }
[[ -r "$GPFF" ]] || { echo "ERROR: Cannot read GPFF: $GPFF"; exit 1; }

echo "Input FAA  : $(realpath "$FAA")"
echo "Input GPFF : $(realpath "$GPFF")"
echo "Output FAA : $(realpath "$OUT")"
echo

tmpdir="$(mktemp -d)"
if [[ "${KEEP_TMP:-0}" != "1" ]]; then
  trap 'rm -rf "$tmpdir"' EXIT
fi
echo "TMPDIR     : $tmpdir"

prot2gene="$tmpdir/prot2gene.tsv"
protlen="$tmpdir/prot_len.tsv"
joined="$tmpdir/joined.tsv"
keep="$tmpdir/keep_proteins.txt"
keep_pub="$(dirname "$OUT")/$(basename "$OUT").ids.txt"

# 1) Build prot->gene map (robust; mawk/gawk-safe); emit versioned + versionless keys
awk -v RS="//" '
function emit(acc,key,   base){
  if (acc=="" || key=="") return
  print acc "\t" key
  base=acc; sub(/\.[0-9]+$/,"",base)
  print base "\t" key
}
{
  acc=""; key=""
  if (match($0, /\nVERSION[ \t]+([A-Z]{2}_[0-9]+(\.[0-9]+)?)/, mV)) acc=mV[1]
  if (acc=="" && match($0, /\nACCESSION[ \t]+([A-Z]{2}_[0-9]+)/, mA)) acc=mA[1]
  if (match($0, /db_xref="GeneID:([0-9]+)"/, mG))      key="GeneID:" mG[1]
  else if (match($0, /locus_tag="([^"]+)"/, mL))       key="locus:"  mL[1]
  else if (match($0, /\ngene="([^"]+)"/, mS))          key="gene:"   mS[1]
  emit(acc, key)
}' "$GPFF" | sort -u > "$prot2gene"

echo "prot2gene rows : $(wc -l < "$prot2gene")"
[[ -s "$prot2gene" ]] || { echo "ERROR: No protein→gene mappings found."; exit 1; }

# 2) Protein order + lengths (id is before first space), add line# for tie-break
seqkit fx2tab -n -l "$FAA" \
  | awk -F'\t' '{id=$1; sub(/ .*/,"",id); print NR, id, $2}' OFS='\t' > "$protlen"

# 3) Join prot_id -> gene_key
# joined.raw: prot_id  line#  length  gene_key
join -t $'\t' -j 1 \
  <(sort -k2,2 "$protlen" | awk -F'\t' '{print $2, $1, $3}' OFS='\t') \
  <(sort -k1,1 "$prot2gene") > "$joined".raw

n_join=$(wc -l < "$joined".raw)
echo "joined rows   : $n_join"
[[ $n_join -gt 0 ]] || { echo "ERROR: No FASTA IDs matched GPFF IDs."; exit 1; }

# 4) Pick longest per gene; ties → smallest line# (first in FASTA)
awk -F'\t' '{print $1,$3,$2,$4}' OFS='\t' "$joined".raw \
| sort -k4,4 -k2,2nr -k3,3n \
| awk -F'\t' '!seen[$4]++{print $1}' > "$keep"

cp -f "$keep" "$keep_pub"  # keep list beside OUT for visibility
echo "keep list     : $keep_pub (n=$(wc -l < "$keep_pub"))"

# 5) Extract sequences with AWK (match header name up to first space)
awk -v keep="$keep" -v FAA="$FAA" '
BEGIN{
  while ((getline line < keep) > 0) { gsub(/\r$/,"",line); ids[line]=1 }
  close(keep)
}
function header_id(s,   t){ t=s; sub(/^>/,"",t); sub(/ .*/,"",t); return t }
BEGINFILE{ printing=0 }
FNR==1 && ARGC>0 { }   # placeholder
/^>/ { printing = ids[ header_id($0) ] ? 1 : 0 }
printing { print }
' "$FAA" > "$OUT"

echo "Done."
echo "Total input proteins : $(grep -c '^>' "$FAA")"
echo "Unique genes kept    : $(wc -l < "$keep")"
echo "Output written to    : $(realpath "$OUT")"
ls -lh "$OUT" || true
grep -n '^>' "$OUT" | head || true
