#!/usr/bin/env Rscript
###############################################################################
# Script name : orthofinder_core_heatmap_and_missingness.R
#
# Author      : Ana Pinharanda
# Email       : app@hsph.harvard.edu
#
# Date        : 04 February 2026
#
# Purpose
# -------
# Use OrthoFinder's Comparative_Genomics_Statistics outputs to generate:
#   1) Orthogroup sharing heatmap (Jaccard) from Orthogroups_SpeciesOverlaps.tsv
#   2) Orthogroup missingness per species from:
#        - Statistics_Overall.tsv (total orthogroups)
#        - Statistics_PerSpecies.tsv ("Number of orthogroups containing species")
#
# Inputs (fixed paths)
# --------------------
# Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv
# Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv
# Comparative_Genomics_Statistics/Statistics_Overall.tsv
#
# Outputs (date-stamped)
# ----------------------
# orthogroup_sharing_heatmap_from_SpeciesOverlaps_<date>.pdf
# orthogroup_missingness_from_stats_<date>.tsv
# orthogroup_missingness_from_stats_<date>.pdf
#
###############################################################################

suppressPackageStartupMessages(library(pheatmap))

today <- format(Sys.Date(), "%Y-%m-%d")

stats_dir <- "Comparative_Genomics_Statistics"
overlaps_file <- file.path(stats_dir, "Orthogroups_SpeciesOverlaps.tsv")
per_species_file <- file.path(stats_dir, "Statistics_PerSpecies.tsv")
overall_file <- file.path(stats_dir, "Statistics_Overall.tsv")

if (!file.exists(overlaps_file)) stop("Missing required file: ", overlaps_file)
if (!file.exists(per_species_file)) stop("Missing required file: ", per_species_file)
if (!file.exists(overall_file)) stop("Missing required file: ", overall_file)

# =============================================================================
# PART A: Jaccard heatmap from Orthogroups_SpeciesOverlaps.tsv
# =============================================================================

ov <- read.delim(overlaps_file, sep = "\t", header = TRUE, check.names = FALSE)

# First column = rownames (species); first header cell is blank
rownames(ov) <- ov[[1]]
ov <- ov[, -1, drop = FALSE]

# Clean species names for plotting
rownames(ov) <- sub("^prefixed_", "", rownames(ov))
colnames(ov) <- sub("^prefixed_", "", colnames(ov))

ov_mat <- as.matrix(ov)
storage.mode(ov_mat) <- "numeric"

present <- diag(ov_mat)
names(present) <- rownames(ov_mat)

sp <- rownames(ov_mat)

jac <- matrix(NA_real_, nrow = length(sp), ncol = length(sp),
              dimnames = list(sp, sp))

for (i in sp) {
  for (j in sp) {
    shared <- ov_mat[i, j]
    denom <- present[i] + present[j] - shared
    jac[i, j] <- shared / denom
  }
}

heat_pdf <- paste0("orthogroup_sharing_heatmap_from_SpeciesOverlaps_", today, ".pdf")
pdf(heat_pdf, width = 6, height = 6)
d <- as.dist(1 - jac)
pheatmap(
  jac,
  clustering_distance_rows = d,
  clustering_distance_cols = d,
  clustering_method = "average",
  border_color = NA,
  display_numbers = FALSE,
  main = "Orthogroup sharing (Jaccard; OrthoFinder SpeciesOverlaps)"
)
dev.off()

# =============================================================================
# PART B: Missingness from Statistics_Overall.tsv + Statistics_PerSpecies.tsv
# =============================================================================

# ---- total orthogroups from Statistics_Overall.tsv ----
overall <- read.delim(overall_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(overall) <- c("metric", "value")
overall$metric <- trimws(overall$metric)

total_og <- as.integer(overall$value[overall$metric == "Number of orthogroups"])
if (length(total_og) != 1 || is.na(total_og)) {
  stop("Could not read 'Number of orthogroups' from: ", overall_file)
}

# ---- per-species table from Statistics_PerSpecies.tsv ----
ps <- read.delim(per_species_file, sep = "\t", header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE)

# First column = metric names
metric_names <- trimws(ps[[1]])

# Remaining columns = species
ps <- ps[, -1, drop = FALSE]
colnames(ps) <- sub("^prefixed_", "", colnames(ps))

# Locate the exact row for orthogroups present in each species
idx <- which(metric_names == "Number of orthogroups containing species")
if (length(idx) != 1) {
  stop("Could not find row 'Number of orthogroups containing species' in: ", per_species_file)
}

# Extract counts and KEEP species names
orthogroups_present <- as.integer(ps[idx, ])
names(orthogroups_present) <- colnames(ps)

if (length(orthogroups_present) == 0) stop("orthogroups_present is empty (unexpected).")
if (any(is.na(orthogroups_present))) stop("NA values found in orthogroups_present (unexpected).")

missingness <- data.frame(
  species = names(orthogroups_present),
  orthogroups_total = rep(total_og, length(orthogroups_present)),
  orthogroups_present = as.integer(orthogroups_present),
  stringsAsFactors = FALSE
)

missingness$orthogroups_absent <- missingness$orthogroups_total - missingness$orthogroups_present
missingness$frac_missing <- missingness$orthogroups_absent / missingness$orthogroups_total

missingness <- missingness[order(-missingness$orthogroups_absent), ]

miss_tsv <- paste0("orthogroup_missingness_from_stats_", today, ".tsv")
write.table(missingness, miss_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nMissingness summary (sorted by absent orthogroups):\n")
print(missingness, row.names = FALSE)

miss_pdf <- paste0("orthogroup_missingness_from_stats_", today, ".pdf")
pdf(miss_pdf, width = 7, height = 4)
par(mar = c(8, 4, 2, 1) + 0.1)

ymax <- max(missingness$orthogroups_absent)
ylim_top <- ymax * 1.20

bp <- barplot(
  missingness$orthogroups_absent,
  names.arg = missingness$species,
  las = 2,
  ylab = "Absent orthogroups (count)",
  main = "Orthogroup missingness per species",
  ylim = c(0, ylim_top)
)

text(
  x = bp,
  y = missingness$orthogroups_absent,
  labels = sprintf("%.1f%%", 100 * missingness$frac_missing),
  pos = 3,
  offset = 0.3,
  cex = 0.85
)

dev.off()

cat("\nOutputs written:\n")
cat("  ", heat_pdf, "\n", sep = "")
cat("  ", miss_tsv, "\n", sep = "")
cat("  ", miss_pdf, "\n", sep = "")
