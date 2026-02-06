#!/bin/bash
#SBATCH --job-name=raxml_species_tree
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=raxml_species_tree.%j.out
#SBATCH --error=raxml_species_tree.%j.err

###############################################################################
# Script name : run_raxml_species_tree.sbatch
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 21 Oct 2025
#
# Purpose
# -------
# Infer a species tree using RAxML-NG from a concatenated alignment of
# single-copy orthologues, with bootstrap support.
#
# Input
# -----
# species_concat.fa
#
# Output
# ------
# ../Species_Tree/raxml_species_bs.*
#
# Requirements
# ------------
# raxml-ng available in PATH
#
###############################################################################

# Number of threads (must match SLURM allocation)
THREADS=16

# Run RAxML-NG
raxml-ng --all \
  --msa SpeciesTreeAlignment.fa \
  --model LG+G \
  --bs-trees 200 \
  --bs-metric fbp \
  --seed 42 \
  --threads ${THREADS} \
  --prefix ../Species_Tree/raxml_species_bs
