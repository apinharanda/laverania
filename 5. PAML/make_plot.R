#!/usr/bin/env Rscript
###############################################################################
# Script: make_plot.R
#
# Purpose
# -------
# Plot the distribution of genome-wide dN/dS (omega) estimates obtained under
# the PAML codeml M0 (one-ratio) model across orthogroups, and compute basic
# summary statistics.
#
# Input
# -----
# paml_M0_summary_allOG_sorted.tsv
#   Tab-delimited table with one row per orthogroup, including:
#     - omega_M0 : M0 dN/dS estimate
#     - status   : "OK" if codeml converged and omega was estimated
#
# Output
# ------
# 1) Histogram of M0 omega values (30 bins):
#    M0_omega_distribution_<DATE>.pdf
#
# 2) Summary statistics table:
#    M0_omega_summary_stats_<DATE>.tsv
#
# Notes
# -----
# - Only orthogroups with status == "OK" are included.
# - No upper bound or truncation is applied to omega values.
# - Histogram binning is fixed at 30 bins for interpretability.
# - The x-axis label explicitly indicates omega = dN/dS under M0.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
df <- read.delim("paml_M0_summary_allOG_sorted.tsv",
                 stringsAsFactors = FALSE)

# Keep only successful M0 runs and valid omega estimates
df_ok <- df %>%
  filter(status == "OK") %>%
  mutate(omega_M0 = as.numeric(omega_M0)) %>%
  filter(is.finite(omega_M0))

stopifnot(nrow(df_ok) > 0)

# ------------------------------------------------------------------
# Plot: distribution of M0 omega values
# ------------------------------------------------------------------
# x-axis: omega (dN/dS) estimated under the M0 model
# y-axis: number of orthogroups (gene families)
p <- ggplot(df_ok, aes(x = omega_M0)) +
  geom_histogram(
    bins = 30,               # fixed number of bins
    fill = "grey80",
    color = "grey30"
  ) +
  labs(
    x = expression(omega~"(dN/dS, M0)"),
    y = "Number of orthogroups"
  ) +
  theme_classic()

pdf(paste0("M0_omega_distribution_", Sys.Date(), ".pdf"),
    width = 6, height = 4)
print(p)
dev.off()

# ------------------------------------------------------------------
# Summary statistics for M0 omega
# ------------------------------------------------------------------
summary_tbl <- tibble(
  n_genes      = nrow(df_ok),
  mean_omega   = mean(df_ok$omega_M0),
  median_omega = median(df_ok$omega_M0),
  sd_omega     = sd(df_ok$omega_M0),
  q05          = quantile(df_ok$omega_M0, 0.05),
  q10          = quantile(df_ok$omega_M0, 0.10),
  q25          = quantile(df_ok$omega_M0, 0.25),
  q75          = quantile(df_ok$omega_M0, 0.75),
  q90          = quantile(df_ok$omega_M0, 0.90),
  q95          = quantile(df_ok$omega_M0, 0.95),
  min_omega    = min(df_ok$omega_M0),
  max_omega    = max(df_ok$omega_M0)
)

write.table(
  summary_tbl,
  file = paste0("M0_omega_summary_stats_", Sys.Date(), ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Print summary to stdout for logs
cat("\nM0 omega summary (successful orthogroups only):\n")
print(summary_tbl)
