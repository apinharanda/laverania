#!/usr/bin/env Rscript
###############################################################################
# Script name : orthofinder_diagnostic_heatmaps_and_missingness.R
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 04 February 2026
#
# Purpose
# -------
# This script performs simple, transparent diagnostics on OrthoFinder output
# to understand discrepancies between:
#   (i) the inferred species tree (e.g. RAxML)
#   (ii) clustering based on orthogroup sharing
#
# Specifically, it:
#   1) Computes orthogroup missingness per species from Orthogroups.tsv
#   2) Plots the number of missing orthogroups per species, annotating
#      the percentage of orthogroups missing
#   3) Computes pairwise Jaccard similarity between species based on
#      orthogroup presence/absence
#   4) Plots Jaccard similarity heatmaps using:
#        - all orthogroups
#        - orthogroups present in >= (N-1) species
#
# The code is intentionally explicit and verbose to aid readability and
# reproducibility.
#
# Input
# -----
# Orthogroups/Orthogroups.tsv
#
# Outputs (date-stamped)
# ----------------------
# orthogroup_missingness_<date>.tsv
# orthogroup_missingness_<date>.pdf
# orthogroup_sharing_ALL_<date>.pdf
# orthogroup_sharing_SHARED_Nminus1_<date>.pdf
#
###############################################################################

# Load required library (only pheatmap is used)
suppressPackageStartupMessages(library(pheatmap))

# -----------------------------
# Define input and date string
# -----------------------------

# Path to OrthoFinder orthogroup assignment table
og_file <- "Orthogroups/Orthogroups.tsv"

# Date string used for all output filenames
today <- format(Sys.Date(), "%Y-%m-%d")

# Fail early if input is missing
if (!file.exists(og_file)) {
  stop("Input file not found: ", og_file)
}

# -----------------------------
# Read Orthogroups.tsv
# -----------------------------

# Read the orthogroup table exactly as written by OrthoFinder.
# First column = orthogroup ID
# Remaining columns = species, with protein IDs or empty cells
og <- read.delim(
  og_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Extract species names from column headers (excluding orthogroup ID)
species <- colnames(og)[-1]

# Remove the "prefixed_" string added prior to OrthoFinder
species_clean <- sub("^prefixed_", "", species)

# Total number of species
n_species <- length(species_clean)

# -----------------------------
# Build presence/absence matrix
# -----------------------------

# Convert orthogroup assignments into a binary matrix:
#   rows    = orthogroups
#   columns = species
#   values  = 1 (present) or 0 (absent)
pa <- as.matrix(
  sapply(
    og[, -1, drop = FALSE],
    function(x) as.integer(!(is.na(x) | x == ""))
  )
)

# Assign clean species names to columns
colnames(pa) <- species_clean

# Assign orthogroup IDs to rows
rownames(pa) <- og[[1]]

# -----------------------------
# Quantify orthogroup missingness
# -----------------------------

# Total number of orthogroups in the analysis
orthogroups_total <- nrow(pa)

# Count how many orthogroups are absent in each species
orthogroups_absent <- colSums(pa == 0)

# Count how many orthogroups are present in each species
orthogroups_present <- orthogroups_total - orthogroups_absent

# Fraction of orthogroups missing per species
frac_missing <- orthogroups_absent / orthogroups_total

# Assemble missingness summary table
missingness <- data.frame(
  species = names(orthogroups_absent),
  orthogroups_total = orthogroups_total,
  orthogroups_present = as.integer(orthogroups_present),
  orthogroups_absent = as.integer(orthogroups_absent),
  frac_missing = as.numeric(frac_missing),
  stringsAsFactors = FALSE
)

# Sort species by decreasing number of missing orthogroups
missingness <- missingness[order(-missingness$orthogroups_absent), ]

# Write missingness table to disk
miss_tsv <- paste0("orthogroup_missingness_", today, ".tsv")
write.table(
  missingness,
  miss_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Plot missingness bar plot
# -----------------------------

miss_pdf <- paste0("orthogroup_missingness_", today, ".pdf")
pdf(miss_pdf, width = 7, height = 4)

# Increase bottom margin for long species names
par(mar = c(8, 4, 2, 1) + 0.1)

# Add extra headroom so percentage labels are not clipped
ymax <- max(missingness$orthogroups_absent)
ylim_top <- ymax * 1.12

# Draw bar plot: number of absent orthogroups per species
bp <- barplot(
  missingness$orthogroups_absent,
  names.arg = missingness$species,
  las = 2,
  ylab = "Absent orthogroups (count)",
  main = "Orthogroup missingness per species",
  ylim = c(0, ylim_top)
)

# Add percentage of orthogroups missing above each bar
text(
  x = bp,
  y = missingness$orthogroups_absent,
  labels = sprintf("%.1f%%", 100 * missingness$frac_missing),
  pos = 3,
  offset = 0.3,
  cex = 0.8
)

dev.off()

# Print missingness table to stdout for logging
cat("\nMissingness summary (sorted by absent orthogroups):\n")
print(missingness, row.names = FALSE)

# -----------------------------
# Jaccard similarity: ALL orthogroups
# -----------------------------

# Initialise empty square matrix for Jaccard similarities
sp <- colnames(pa)
jac_all <- matrix(
  NA_real_,
  nrow = length(sp),
  ncol = length(sp),
  dimnames = list(sp, sp)
)

# Compute pairwise Jaccard similarity between species
# Jaccard = |A ∩ B| / |A ∪ B|
for (i in sp) {
  for (j in sp) {
    a <- pa[, i] == 1
    b <- pa[, j] == 1
    jac_all[i, j] <- sum(a & b) / sum(a | b)
  }
}

# Check whether matrix is informative (not constant)
rng_all <- range(jac_all, na.rm = TRUE)

if (is.finite(rng_all[1]) && is.finite(rng_all[2]) && diff(rng_all) > 0) {

  out_all <- paste0("orthogroup_sharing_ALL_", today, ".pdf")

  # Convert similarity to distance for clustering
  d_all <- as.dist(1 - jac_all)

  pdf(out_all, width = 6, height = 6)
  pheatmap(
    jac_all,
    clustering_distance_rows = d_all,
    clustering_distance_cols = d_all,
    clustering_method = "average",
    border_color = NA,
    display_numbers = FALSE,
    main = "Orthogroup sharing (all orthogroups)"
  )
  dev.off()

} else {
  cat("\nSkipping ALL orthogroup heatmap (constant matrix)\n")
}

# -----------------------------
# Jaccard similarity: shared by >= (N-1) species
# -----------------------------

# Count how many species each orthogroup is present in
present_per_og <- rowSums(pa)

# Retain orthogroups present in all but at most one species
min_shared <- max(2, n_species - 1)
pa_shared <- pa[present_per_og >= min_shared, , drop = FALSE]

# Initialise Jaccard matrix for shared orthogroups
jac_shared <- matrix(
  NA_real_,
  nrow = length(sp),
  ncol = length(sp),
  dimnames = list(sp, sp)
)

# Compute Jaccard similarities using the filtered matrix
for (i in sp) {
  for (j in sp) {
    a <- pa_shared[, i] == 1
    b <- pa_shared[, j] == 1
    jac_shared[i, j] <- sum(a & b) / sum(a | b)
  }
}

# Check for non-constant matrix
rng_shared <- range(jac_shared, na.rm = TRUE)

if (is.finite(rng_shared[1]) && is.finite(rng_shared[2]) && diff(rng_shared) > 0) {

  out_shared <- paste0("orthogroup_sharing_SHARED_Nminus1_", today, ".pdf")
  d_shared <- as.dist(1 - jac_shared)

  pdf(out_shared, width = 6, height = 6)
  pheatmap(
    jac_shared,
    clustering_distance_rows = d_shared,
    clustering_distance_cols = d_shared,
    clustering_method = "average",
    border_color = NA,
    display_numbers = FALSE,
    main = paste0(
      "Orthogroup sharing (shared ≥ ",
      min_shared,
      " species)"
    )
  )
  dev.off()

} else {
  cat("\nSkipping SHARED orthogroup heatmap (constant matrix)\n")
}

# -----------------------------
# End of script
# -----------------------------
cat("\nAnalysis complete.\n")
