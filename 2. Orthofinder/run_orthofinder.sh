#!/usr/bin/env bash
###############################################################################
# Script name : run_orthofinder.sh
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 16 October 2025
#
# Purpose
# -------
# Run OrthoFinder on Laverania protein datasets after species-specific
# prefixing of FASTA headers. This step infers orthogroups, multiple sequence
# alignments, and gene trees for downstream comparative analyses.
#
# The input directory is expected to contain one protein FASTA file per species,
# with headers already prefixed by species name (see add_species_prefix.sh).
#
# Execution environment
# ---------------------
# SLURM-managed HPC cluster
#
# Requirements
# ------------
# orthofinder
# mafft
# raxml-ng
#
###############################################################################

#SBATCH --job-name=laverania
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000

###############################################################################
# OrthoFinder parameters
#
# -f : input directory containing protein FASTA files (one per species)
# -t : total number of threads available to OrthoFinder
# -a : number of threads per alignment / tree-building job
# -M : run in MSA mode (build multiple sequence alignments)
# -A : aligner used for MSAs (mafft)
# -T : tree inference program (raxml-ng)
# -I : inflation parameter for MCL clustering
#
###############################################################################

orthofinder \
  -f /n/holylabs/LABS/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix \
  -t 16 \
  -a 8 \
  -M msa \
  -A mafft \
  -T raxml-ng \
  -I 1.5
