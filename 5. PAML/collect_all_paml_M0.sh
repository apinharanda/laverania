#!/usr/bin/env bash
###############################################################################
# Script name : collect_all_paml_M0.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 09 February 2026
#
# Purpose
# -------
# Summarize codeml (PAML) M0 results across OGxxxx directories into one TSV.
# IMPORTANT: values are parsed ONLY from codeml output (*.paml.out).
#
# Expected directory layout
# -------------------------
# OG_paml/
#   OG0000000/OG0000000.paml.out
#   OG0000001/OG0000001.paml.out
#   ...
#
# Output
# ------
# paml_M0_summary_allOG.tsv
# paml_M0_summary_allOG_sorted.tsv
#
# Columns
# -------
# OG, omega_M0, kappa, lnL, S, N,
# tree_length, tree_length_dN, tree_length_dS,
# nseq, aln_len_codons, status
#
###############################################################################

OUT="paml_M0_summary_allOG.tsv"
OUT_SORTED="paml_M0_summary_allOG_sorted.tsv"
NUM='-?[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?'

printf "OG\tomega_M0\tkappa\tlnL\tS\tN\ttree_length\ttree_length_dN\ttree_length_dS\tnseq\taln_len_codons\tstatus\n" > "$OUT"

for d in OG*; do
  [[ -d "$d" ]] || continue
  og="$d"
  out="${d}/${og}.paml.out"

  omega="NA"; kappa="NA"; lnL="NA"
  S="NA"; N="NA"
  tree_length="NA"; tree_length_dN="NA"; tree_length_dS="NA"
  nseq="NA"; aln_len_codons="NA"
  status="FAIL"

  if [[ -s "$out" ]]; then
    # omega (M0)
    omega=$(awk -F'=' '
      /^[[:space:]]*omega[[:space:]]*\(dN\/dS\)[[:space:]]*=/ {
        gsub(/[ \t]/,"",$2); print $2; exit
      }' "$out" 2>/dev/null || true)
    [[ -n "${omega:-}" ]] || omega="NA"

    # kappa
    kappa=$(awk -F'=' '
      /^[[:space:]]*kappa[[:space:]]*\(ts\/tv\)[[:space:]]*=/ {
        gsub(/[ \t]/,"",$2); print $2; exit
      }' "$out" 2>/dev/null || true)
    [[ -n "${kappa:-}" ]] || kappa="NA"

    # lnL (the numeric value after "):")
    lnL=$(awk '
      /^[[:space:]]*lnL\(ntime:/ {
        if (match($0, /\):[[:space:]]*('"$NUM"')/, m)) { print m[1]; exit }
      }' "$out" 2>/dev/null || true)
    [[ -n "${lnL:-}" ]] || lnL="NA"

    # codeml-reported tree lengths (re-estimated under M0)
    tree_length=$(awk '
      /^[[:space:]]*tree length[[:space:]]*=/ {
        if (match($0, /= *('"$NUM"')/, m)) { print m[1]; exit }
      }' "$out" 2>/dev/null || true)
    [[ -n "${tree_length:-}" ]] || tree_length="NA"

    tree_length_dN=$(awk '
      /^[[:space:]]*tree length for dN[[:space:]]*:/ {
        if (match($0, /: *('"$NUM"')/, m)) { print m[1]; exit }
      }' "$out" 2>/dev/null || true)
    [[ -n "${tree_length_dN:-}" ]] || tree_length_dN="NA"

    tree_length_dS=$(awk '
      /^[[:space:]]*tree length for dS[[:space:]]*:/ {
        if (match($0, /: *('"$NUM"')/, m)) { print m[1]; exit }
      }' "$out" 2>/dev/null || true)
    [[ -n "${tree_length_dS:-}" ]] || tree_length_dS="NA"

    # nseq + alignment length (codons) from the header line "<ns> <ls>"
    read -r ns_tmp ls_tmp < <(
      awk '
        /^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]*$/ { print $1, $2; exit }
      ' "$out" 2>/dev/null || true
    )
    if [[ -n "${ns_tmp:-}" && -n "${ls_tmp:-}" ]]; then
      nseq="$ns_tmp"
      aln_len_codons="$ls_tmp"
    fi

    # N and S from the first branch row in the "dN & dS for each branch" table
    read -r N_tmp S_tmp < <(
      awk '
        /dN[[:space:]]*&[[:space:]]*dS[[:space:]]*for/ { inblk=1; next }
        inblk && $1 ~ /^[0-9]+\.\.[0-9]+$/ && NF>=7 { print $3, $4; exit }
      ' "$out" 2>/dev/null || true
    )
    if [[ -n "${N_tmp:-}" && -n "${S_tmp:-}" ]]; then
      N="$N_tmp"
      S="$S_tmp"
    fi

    # success criterion: omega line present
    [[ "$omega" != "NA" ]] && status="OK"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$og" "$omega" "$kappa" "$lnL" "$S" "$N" \
    "$tree_length" "$tree_length_dN" "$tree_length_dS" \
    "$nseq" "$aln_len_codons" "$status" >> "$OUT"
done

# Sort (OG IDs are zero-padded; plain sort works)
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -k1,1) > "$OUT_SORTED"

echo "[OK] Wrote: $OUT_SORTED"
