#!/usr/bin/env bash
###############################################################################
# Script name : 04_make_paml_M0_dirs.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 06 February 2026
#
# Purpose
# -------
# Create a PAML (codeml) M0 directory per orthogroup (OG000xxxx), containing:
#   1) a symlink to the codon alignment
#   2) a symlink to the corresponding gene tree (_tree.txt)
#   3) a codeml control file (codeml.<OG>.ctl) whose seqfile/treefile/outfile
#      names match the files in that directory.
#
# Inputs
# ------
# Codon alignments (PAML PHYLIP):
#   /n/holylabs/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix/OrthoFinder/Results_Oct16_6/Codons_raw
#
# Gene trees (Newick):
#   Gene_Trees/OG000xxxx_tree.txt
#
# Output
# ------
# OG_paml/OG000xxxx/ directory per OG, containing:
#   OG000xxxx.codon.fa -> (symlink)
#   OG000xxxx_tree.txt     -> (symlink)
#   codeml.OG000xxxx.ctl   (new file)
#
# Notes on gaps (important)
# -------------------------
# - Do NOT set cleandata = 1 for gappy alignments unless to explicitly
#   drop every column containing any gap (very aggressive).
# - Here we set: cleandata = 0 (keep columns; PAML treats gaps as missing data).
#
# How to run
# ----------
# From: /n/holylabs/.../OrthoFinder/Results_Sep23_1
#   bash make_paml_M0_dirs.sh
#
###############################################################################


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CODON_DIR="Codons_raw"
TREE_DIR="Gene_Trees"
OUT_BASE="OG_paml"

mkdir -p "$OUT_BASE"

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
[[ -d "$CODON_DIR" ]] || { echo "ERROR: Missing codon dir: $CODON_DIR"; exit 1; }
[[ -d "$TREE_DIR"  ]] || { echo "ERROR: Missing tree dir:  $TREE_DIR"; exit 1; }

# ---------------------------------------------------------------------------
# Loop over all codon alignments and build per-OG folders
# ---------------------------------------------------------------------------
shopt -s nullglob

n_total=0
n_made=0
n_missing_tree=0

for phy in "$CODON_DIR"/OG*.codon.fa; do
  ((n_total++))

  og="$(basename "$phy" .codon.fa)"          # OG000xxxx
  tree="$TREE_DIR/${og}_tree.txt"

  # skip if tree missing
  if [[ ! -s "$tree" ]]; then
    echo "[WARN] Missing tree for ${og}: ${tree}"
    ((n_missing_tree++))
    continue
  fi

  og_dir="${OUT_BASE}/${og}"
  mkdir -p "$og_dir"

  # symlinks (force update)
  ln -sf "$(realpath "$phy")"  "${og_dir}/${og}.codon.fa"
  ln -sf "$(realpath "$tree")" "${og_dir}/${og}_tree.txt"

  # write ctl (overwrite to keep deterministic)
  cat > "${og_dir}/codeml.${og}.ctl" <<EOF
      seqfile = ${og}.codon.fa     * alignment file
      treefile = ${og}_tree.txt        * tree file (Newick)
      outfile = ${og}.paml.out         * main result file

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1          * 1:codons
    CodonFreq = 2          * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

        clock = 0
       aaDist = 0
   model = 0               * 0:one-ratio (M0)
  NSsites = 0

      icode = 0            * 0:universal genetic code
  fix_kappa = 0
      kappa = 2
  fix_omega = 0
      omega = 0.2

  fix_alpha = 1
      alpha = 0
     Malpha = 0
      ncatG = 8

      getSE = 0
RateAncestor = 0

    cleandata = 0          * 0:keep columns w/ gaps as missing; 1:remove any-gap columns

  fix_blength = 0
       method = 0
EOF

  ((n_made++))
done

echo "[DONE] Found codon alignments: ${n_total}"
echo "[DONE] Built OG directories:  ${n_made}"
echo "[DONE] Missing trees:         ${n_missing_tree}"
echo "[DONE] Output base:           ${OUT_BASE}"
