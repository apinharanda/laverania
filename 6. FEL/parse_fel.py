#!/usr/bin/env python3
###############################################################################
# parse_fel.py
#
# Author: Ana Pinharanda (with ChatGPT helper)
# Date:   2026-02-27
#
# Goal
# ----
# Take ONE HyPhy FEL JSON output (produced by `hyphy fel --ci Yes`) and extract
# the per-codon table into a TSV that matches what you see in the HyPhy log:
#
#   Codon | Partition | alpha | beta | LRT | p-value | dN/dS MLE | LB | UB | selection
#
# The TSV is designed to be:
#   - easy to import into R (separate numeric columns)
#   - also includes a human-readable dN/dS CI string like: 0.000(0.00-0.46)
#
# Assumptions (keep it simple)
# ----------------------------
# This script assumes the FEL JSON has exactly this structure
#
#   json["MLE"]["headers"]  : list of [name, description]
#   json["MLE"]["content"]  : dict keyed by partition index as strings ("0","1",...)
#                             each value is a 2D list (rows=coding sites, cols=fields)
#
# If JSON breaks this schema, the script should fail loudly rather than
# trying lots of fallbacks.
#
# Why "partition 0" in TSV but "partition 1" in HyPhy log?
# --------------------------------------------------------
# The JSON uses 0-based partition indexing ("0"). Some HyPhy console tables show
# partitions as 1-based. The biology is the same; it is only a label convention.
# 
#
# Usage
# -----
#   python3 parse_fel_json_to_tsv.py OG0004504.FEL.json OG0004504.FEL.tsv
#
###############################################################################

import json
import sys


# ----------------------------
# Columns we will extract
# ----------------------------
# These names MUST match the column names found in json["MLE"]["headers"].
COL_ALPHA = "alpha"        # synonymous rate estimate at the site
COL_BETA  = "beta"         # nonsynonymous rate estimate at the site
COL_LRT   = "LRT"          # likelihood ratio statistic for beta != alpha
COL_PVAL  = "p-value"      # asymptotic p-value for beta != alpha
COL_DNDS  = "dN/dS MLE"    # point estimate of dN/dS at the site
COL_LB    = "dN/dS LB"     # lower CI bound for dN/dS (requires --ci Yes)
COL_UB    = "dN/dS UB"     # upper CI bound for dN/dS (requires --ci Yes)

# Output TSV column order (kept stable for downstream scripts)
OUT_COLS = [
    "codon",
    "partition",
    "alpha",
    "beta",
    "LRT",
    "p_value",
    "dnds_mle",
    "dnds_lb",
    "dnds_ub",
    "selection",
    "dnds_ci",
]


def load_json(path):
    """
    Read a JSON file and return the parsed Python object (dict/list).
    """
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def get_header_names(mle_obj):
    """
    Extract column names from FEL:
      mle_obj["headers"] is a list of [name, description]

    We return a simple list of names, like:
      ["alpha", "beta", "alpha=beta", "LRT", "p-value", ..., "dN/dS UB"]
    """
    headers = mle_obj["headers"]
    names = []
    for h in headers:
        # h is like ["alpha", "Synonymous substitution rate at a site"]
        names.append(h[0])
    return names


def make_column_index(header_names):
    """
    Build a dict mapping:
      column name -> integer index in each data row

    Example:
      idx["alpha"] = 0
      idx["beta"]  = 1
      idx["LRT"]   = 3
      ...
    """
    idx = {}
    for i, name in enumerate(header_names):
        idx[name] = i

    # Fail early if FEL output does not contain the columns we rely on.
    required = [COL_ALPHA, COL_BETA, COL_LRT, COL_PVAL, COL_DNDS, COL_LB, COL_UB]
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(
            "FEL headers missing required columns: " + ", ".join(missing) +
            "\nFound headers: " + str(header_names)
        )

    return idx


def as_float(x):
    """
    Convert a value from the JSON table into a float.

    HyPhy writes numbers as JSON numbers, so in practice x is already int/float.
    We still cast with float(x) for safety and to standardize types.
    """
    if x is None:
        return None
    return float(x)


def selection_label(alpha, beta, p_value, p_cutoff=0.1):
    """
    Recreate the "Selection detected?" idea from the log, but in a simple label.

    The FEL log prints things like:
      "Neg. p = 0.0749"
      "Pos. p = 0.0423"

    We implement:
      - if p_value <= p_cutoff and beta > alpha  => "Pos"
      - if p_value <= p_cutoff and beta < alpha  => "Neg"
      - else                                     => "None"

    Notes:
      - This is *not* used for the numeric outputs.
      - Change p_cutoff if you want 0.05 or something else.
    """
    if p_value is None or alpha is None or beta is None:
        return "None"
    if p_value > p_cutoff:
        return "None"
    if beta > alpha:
        return "Pos"
    if beta < alpha:
        return "Neg"
    return "None"


def format_dnds_ci(dnds_mle, dnds_lb, dnds_ub):
    """
    Produce a human-readable string matching the log style:
      0.000(0.00-0.46)

    Formatting choices:
      - MLE: 3 decimals (like log table typically shows)
      - LB/UB: 2 decimals
    """
    if dnds_mle is None or dnds_lb is None or dnds_ub is None:
        return "NA"
    return f"{dnds_mle:.3f}({dnds_lb:.2f}-{dnds_ub:.2f})"


def parse_fel_json(in_json_path):
    """
    Main extraction logic:
      1) Load JSON
      2) Read MLE headers -> find column indices
      3) Iterate over partitions in MLE.content
      4) For each row (codon), pull alpha/beta/LRT/p/etc and write an output row

    Returns:
      list of dicts (one dict per codon per partition)
    """
    j = load_json(in_json_path)

    # FEL site table lives under the "MLE" object (per HyPhy JSON schema).
    mle = j["MLE"]

    # (1) Convert headers into column names
    header_names = get_header_names(mle)

    # (2) Map each needed name -> column index in each row of content
    col = make_column_index(header_names)

    # (3) FEL tables are stored by partition under mle["content"]
    #     Example: mle["content"]["0"] is a list of rows (codons) for partition 0.
    content_by_partition = mle["content"]

    out_rows = []

    #  iterate partitions in numeric order: "0", "1", ...
    partition_keys = sorted(content_by_partition.keys(), key=lambda k: int(k))

    for part_key in partition_keys:
        partition_id = int(part_key)          # partition index as stored in JSON (0-based)
        rows = content_by_partition[part_key] # 2D array: rows=codons, cols=fields

        # If you want the TSV to show partitions as 1-based (like some logs), change here:
        partition_out = partition_id          # or: partition_id + 1

        # codon is simply the row number (1-based) in this table
        for i, row in enumerate(rows, start=1):
            # Pull each metric by its column index.
            # This mirrors: row[col["alpha"]], row[col["beta"]], etc.
            alpha = as_float(row[col[COL_ALPHA]])
            beta  = as_float(row[col[COL_BETA]])
            lrt   = as_float(row[col[COL_LRT]])
            pval  = as_float(row[col[COL_PVAL]])
            dnds  = as_float(row[col[COL_DNDS]])
            lb    = as_float(row[col[COL_LB]])
            ub    = as_float(row[col[COL_UB]])

            out_rows.append({
                "codon": i,
                "partition": partition_out,
                "alpha": alpha,
                "beta": beta,
                "LRT": lrt,
                "p_value": pval,
                "dnds_mle": dnds,
                "dnds_lb": lb,
                "dnds_ub": ub,
                "selection": selection_label(alpha, beta, pval, p_cutoff=0.1),
                "dnds_ci": format_dnds_ci(dnds, lb, ub),
            })

    return out_rows


def write_tsv(rows, out_path):
    """
    Write the extracted rows to a TSV with a fixed column order.

    Numeric formatting:
      - We keep numeric columns as plain numbers
    """
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\t".join(OUT_COLS) + "\n")

        for r in rows:
            fields = []
            for c in OUT_COLS:
                v = r[c]
                if v is None:
                    fields.append("NA")
                elif isinstance(v, float):
                    fields.append(f"{v:.10g}")
                else:
                    fields.append(str(v))
            f.write("\t".join(fields) + "\n")


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python3 parse_fel_json_to_tsv.py <in.FEL.json> <out.tsv>\n")
        sys.exit(2)

    in_json = sys.argv[1]
    out_tsv = sys.argv[2]

    rows = parse_fel_json(in_json)
    write_tsv(rows, out_tsv)

    sys.stderr.write(f"[OK] Parsed {len(rows)} rows -> {out_tsv}\n")


if __name__ == "__main__":
    main()
