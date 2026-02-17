#!/usr/bin/env Rscript
###############################################################################
# Script: make_piN_per_gene_plot.R
#
# Purpose
# -------
# Summarise site-level piN (pi_ns) per gene, then plot the distribution across
# genes and write summary statistics.
#
# Input
# -----
# /n/holylabs/LABS/neafsey_lab/Lab/perkins/data/drc_pi_ns_by_codon.csv
#   Columns: gene, codon_position (0-based), pi_ns
#
# Output
# ------
# 1) Histogram of per-gene piN (30 bins):
#    piN_per_gene_distribution_<DATE>.pdf
#
# 2) Per-gene summary table (one row per gene):
#    piN_per_gene_table_<DATE>.tsv
#
# 3) Global summary stats across genes:
#    piN_per_gene_summary_stats_<DATE>.tsv
#
# Notes
# -----
# - By default uses the MEAN piN per gene (mean across codons).
# - Also writes median piN, gene length in codons, and fraction of zero codons.
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
# Choose which per-gene value to plot: "mean" or "median"
PLOT_STAT <- "mean"

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
  filter(!is.na(gene) & nzchar(gene) & is.finite(piN))

stopifnot(nrow(df_ok) > 0)

# ------------------------------------------------------------------
# Per-gene summaries
# ------------------------------------------------------------------
gene_tbl <- df_ok %>%
  group_by(gene) %>%
  summarise(
    n_codons = dplyr::n(),
    mean_piN = mean(piN),
    median_piN = median(piN),
    sd_piN = sd(piN),
    frac_zero = mean(piN == 0),
    .groups = "drop"
  )

stopifnot(nrow(gene_tbl) > 0)

# Pick plotting value
if (PLOT_STAT == "median") {
  gene_tbl <- gene_tbl %>% mutate(piN_gene = median_piN)
  xlab_txt <- "Per-gene piN (median across codons)"
} else {
  gene_tbl <- gene_tbl %>% mutate(piN_gene = mean_piN)
  xlab_txt <- "Per-gene piN (mean across codons)"
}

# ------------------------------------------------------------------
# Plot distribution across genes
# ------------------------------------------------------------------
p <- ggplot(gene_tbl, aes(x = piN_gene)) +
  geom_histogram(
    bins = BINS,
    fill = "grey80",
    color = "grey30"
  ) +
  labs(
    x = xlab_txt,
    y = "Number of genes"
  ) +
  theme_classic()

pdf(paste0("piN_per_gene_distribution_", Sys.Date(), ".pdf"),
    width = 6, height = 4)
print(p)
dev.off()

# ------------------------------------------------------------------
# Write per-gene table
# ------------------------------------------------------------------
write.table(
  gene_tbl %>% arrange(desc(piN_gene)),
  file = paste0("piN_per_gene_table_", Sys.Date(), ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------------
# Summary stats across genes
# ------------------------------------------------------------------
summary_tbl <- tibble(
  n_genes      = nrow(gene_tbl),
  mean_piN     = mean(gene_tbl$piN_gene),
  median_piN   = median(gene_tbl$piN_gene),
  sd_piN       = sd(gene_tbl$piN_gene),
  q05          = as.numeric(quantile(gene_tbl$piN_gene, 0.05)),
  q10          = as.numeric(quantile(gene_tbl$piN_gene, 0.10)),
  q25          = as.numeric(quantile(gene_tbl$piN_gene, 0.25)),
  q75          = as.numeric(quantile(gene_tbl$piN_gene, 0.75)),
  q90          = as.numeric(quantile(gene_tbl$piN_gene, 0.90)),
  q95          = as.numeric(quantile(gene_tbl$piN_gene, 0.95)),
  min_piN      = min(gene_tbl$piN_gene),
  max_piN      = max(gene_tbl$piN_gene)
)

write.table(
  summary_tbl,
  file = paste0("piN_per_gene_summary_stats_", Sys.Date(), ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nPer-gene piN summary across genes:\n")
print(summary_tbl)
