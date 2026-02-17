#!/usr/bin/env Rscript
###############################################################################
# Script: make_piNpiS_mean_across_pops.R
#
# Purpose
# -------
# Compute per-gene mean piN and mean piS across populations,
# then compute piN/piS from those means and plot distribution.
#
# Optional:
#   Apply piS cutoff AFTER averaging across populations.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------
INFILE <- "/n/holylabs/neafsey_lab/Lab/apinharanda/laverania/proteins_to_run_orthofinder/prefix/OrthoFinder/Results_Sep23_1/edit_sarah_table.csv"
BINS   <- 30

# Optional cutoff (set to NULL to disable)
PIS_CUTOFF <- NULL       # e.g. 1e-4 

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
df <- read.csv(INFILE, stringsAsFactors = FALSE)

df_clean <- df %>%
  mutate(
    pi_ns = as.numeric(pi_ns),
    pi_s  = as.numeric(pi_s)
  ) %>%
  filter(is.finite(pi_ns), is.finite(pi_s))

# ------------------------------------------------------------------
# Compute per-gene MEAN across populations
# ------------------------------------------------------------------
gene_mean <- df_clean %>%
  group_by(GENE) %>%
  summarise(
    mean_piN = mean(pi_ns),
    mean_piS = mean(pi_s),
    .groups = "drop"
  ) %>%
  mutate(
    pnps_mean = mean_piN / mean_piS
  ) %>%
  filter(is.finite(pnps_mean))

# ------------------------------------------------------------------
# Optional piS cutoff
# ------------------------------------------------------------------
if (!is.null(PIS_CUTOFF)) {
  gene_mean <- gene_mean %>%
    filter(mean_piS >= PIS_CUTOFF)
  cat("Applied piS cutoff:", PIS_CUTOFF, "\n")
}

stopifnot(nrow(gene_mean) > 0)

cat("\nTotal genes:", nrow(gene_mean), "\n")

# ------------------------------------------------------------------
# Plot distribution
# ------------------------------------------------------------------
p <- ggplot(gene_mean, aes(x = pnps_mean)) +
  geom_histogram(
    bins = BINS,
    fill = "grey80",
    color = "grey30"
  ) +
  labs(
    x = expression("Mean " * pi[N] / pi[S] * " across populations"),
    y = "Number of genes"
  ) +
  theme_classic(base_size = 16)

pdf(paste0("piN_piS_meanAcrossPops_", Sys.Date(), ".pdf"),
    width = 6, height = 4)
print(p)
dev.off()

# ------------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------------
summary_tbl <- tibble(
  n_genes     = nrow(gene_mean),
  mean_pnps   = mean(gene_mean$pnps_mean),
  median_pnps = median(gene_mean$pnps_mean),
  sd_pnps     = sd(gene_mean$pnps_mean),
  q05         = quantile(gene_mean$pnps_mean, 0.05),
  q25         = quantile(gene_mean$pnps_mean, 0.25),
  q75         = quantile(gene_mean$pnps_mean, 0.75),
  q95         = quantile(gene_mean$pnps_mean, 0.95),
  min_pnps    = min(gene_mean$pnps_mean),
  max_pnps    = max(gene_mean$pnps_mean)
)

write.table(
  summary_tbl,
  file = paste0("piN_piS_meanAcrossPops_summary_", Sys.Date(), ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSummary statistics:\n")
print(summary_tbl)

cat("\nDone.\n")
