#!/usr/bin/env Rscript
###############################################################################
# Script: make_piN_plot.R
#
# Purpose
# -------
# Plot the distribution of site-level nonsynonymous nucleotide diversity (piN)
# values from the Pf3D7 master per-codon file (drc_pi_ns_by_codon.csv), and
# compute basic summary statistics.
#
# Input
# -----
# /n/holylabs/LABS/neafsey_lab/Lab/perkins/data/drc_pi_ns_by_codon.csv
#   CSV with columns:
#     - gene            (e.g., PF3D7_0100100)
#     - genome_position (bp, optional for this script)
#     - codon_position  (0-based codon index)
#     - pi_ns           (piN value per codon)
#
# Output
# ------
# 1) Histogram (30 bins):
#    piN_distribution_<DATE>.pdf
#
# 2) Summary statistics table:
#    piN_summary_stats_<DATE>.tsv
#
# Notes
# -----
# - This is a genome-wide site-level distribution (all codons across all genes).
# - Only finite pi_ns values are included.
# - No truncation/winsorization is applied.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------
PIN_FILE <- "/n/holylabs/LABS/neafsey_lab/Lab/perkins/data/drc_pi_ns_by_codon.csv"
BINS <- 30

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
stopifnot(file.exists(PIN_FILE))

df <- read.csv(PIN_FILE, stringsAsFactors = FALSE)

need_cols <- c("gene", "codon_position", "pi_ns")
if (!all(need_cols %in% names(df))) {
  stop("Input file missing required columns: ",
       paste(setdiff(need_cols, names(df)), collapse = ", "),
       "\nFound: ", paste(names(df), collapse = ", "))
}

df_ok <- df %>%
  mutate(
    gene = as.character(gene),
    codon_position = suppressWarnings(as.integer(codon_position)),
    piN = suppressWarnings(as.numeric(pi_ns))
  ) %>%
  filter(is.finite(piN))

stopifnot(nrow(df_ok) > 0)

# ------------------------------------------------------------------
# Plot: distribution of piN values (all codons)
# ------------------------------------------------------------------
p <- ggplot(df_ok, aes(x = piN)) +
  geom_histogram(
    bins = BINS,
    fill = "grey80",
    color = "grey30"
  ) +
  labs(
    x = "piN (per codon)",
    y = "Number of codons"
  ) +
  theme_classic()

pdf(paste0("piN_distribution_", Sys.Date(), ".pdf"),
    width = 6, height = 4)
print(p)
dev.off()

# ------------------------------------------------------------------
# Summary statistics for piN
# ------------------------------------------------------------------
summary_tbl <- tibble(
  n_codons    = nrow(df_ok),
  n_genes     = dplyr::n_distinct(df_ok$gene),
  mean_piN    = mean(df_ok$piN),
  median_piN  = median(df_ok$piN),
  sd_piN      = sd(df_ok$piN),
  q05         = as.numeric(quantile(df_ok$piN, 0.05)),
  q10         = as.numeric(quantile(df_ok$piN, 0.10)),
  q25         = as.numeric(quantile(df_ok$piN, 0.25)),
  q75         = as.numeric(quantile(df_ok$piN, 0.75)),
  q90         = as.numeric(quantile(df_ok$piN, 0.90)),
  q95         = as.numeric(quantile(df_ok$piN, 0.95)),
  min_piN     = min(df_ok$piN),
  max_piN     = max(df_ok$piN),
  frac_zero   = mean(df_ok$piN == 0, na.rm = TRUE)
)

write.table(
  summary_tbl,
  file = paste0("piN_summary_stats_", Sys.Date(), ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Print summary to stdout for logs
cat("\npiN summary (all finite codons across all genes):\n")
print(summary_tbl)
