#!/usr/bin/env python3
###############################################################################
# SCRIPT: parse_fubar_json_to_tsv_withGaps_fullposIsCodon.py
# Author: Ana Pinharanda (with ChatGPT helper) 
# WHAT THIS DOES
#   Convert HyPhy FUBAR JSON output into a TSV table (one row per codon/site).
#
# WHY
#   HyPhy JSON is annoying
# INPUT
#   A FUBAR JSON file produced by HyPhy, typically named:
#     OG0003384.FUBAR.json
#
# OUTPUT
#   A TSV file with columns:
#     full_pos    : alignment codon index (1..N)
#     codon       : alignment codon index (1..N)  [same as full_pos here]
#     alpha       : synonymous rate proxy (HyPhy)
#     beta        : nonsynonymous rate proxy (HyPhy)
#     omega       : beta/alpha (NaN if alpha==0)
#     P_pos       : Prob[alpha<beta] (evidence for positive selection at site)
#     P_neg       : Prob[alpha>beta] (evidence for negative selection at site)
#     BF_pos      : BayesFactor[alpha<beta]
#     beta_minus_alpha : (beta-alpha)
#
# COORDINATE WARNING (IMPORTANT)
#   This script DOES NOT map to genomic coordinates.
#   "full_pos" here is *alignment codon index*.
#   It is called "full_pos" because downstream plotting scripts expects
#   that column name.
#
#  Pf3D7 coordinates / "real coords" is a separate step.
#
# USAGE
#   python3 parse_fubar_json_to_tsv_withGaps_fullposIsCodon.py <fubar.json> <out.tsv>
#
# EXAMPLE
#   python3 parse_fubar_json_to_tsv_withGaps_fullposIsCodon.py \
#       OG0003384.FUBAR.json FUBAR_per_codon_withGaps_realcoords.tsv
#
# COMMON FAILURE MODES
#   1) JSON is empty/corrupt -> JSONDecodeError
#   2) HyPhy JSON schema differs -> missing MLE/headers/content keys
#   3) Column names differ -> script exits with "Missing column"
#
# NOTE ON "WITH GAPS"
#   "withGaps" in the filename refers to the MSA. parser itself does not “add gaps”.
###############################################################################

import json
import sys

# Columns we need from the HyPhy table.
# These names must match the JSON headers exactly.
NEED = [
    "alpha",
    "beta",
    "Prob[alpha<beta]",
    "Prob[alpha>beta]",
    "BayesFactor[alpha<beta]",
    "beta-alpha",
]


def die(msg: str):
    """Exit immediately with an informative message."""
    raise SystemExit(f"[ERROR] {msg}")


def load_rows(j):
    """
    Pull out:
      - the header list (MLE.headers)
      - a lookup dict from header -> column index
      - the rows (MLE.content[<key>]) where each row is a list of values

    HyPhy sometimes stores headers as either:
      ["alpha","beta",...]
    or:
      [["alpha",...],["beta",...],...]
    so we normalize them to strings.
    """
    if "MLE" not in j:
        die("Missing top-level key 'MLE' in JSON.")
    if "headers" not in j["MLE"]:
        die("Missing key 'MLE.headers' in JSON.")
    if "content" not in j["MLE"]:
        die("Missing key 'MLE.content' in JSON.")

    hdr_raw = j["MLE"]["headers"]
    hdr = [h[0] if isinstance(h, list) else h for h in hdr_raw]  # normalize

    idx = {name: i for i, name in enumerate(hdr)}

    # Ensure required columns exist
    for k in NEED:
        if k not in idx:
            die(f"Missing column in JSON: {k}. Have: {hdr}")

    # HyPhy content is a dict of (usually) {"0": rows}
    content = j["MLE"]["content"]
    if not isinstance(content, dict) or len(content) == 0:
        die("MLE.content is empty or not a dict.")

    key = "0" if "0" in content else next(iter(content.keys()))
    rows = content[key]

    if not isinstance(rows, list) or len(rows) == 0:
        die(f"MLE.content['{key}'] is empty or not a list of rows.")

    return rows, idx


def main():
    # --------------------
    # Parse arguments
    # --------------------
    if len(sys.argv) != 3:
        die("Usage: parse_fubar_json_to_tsv_withGaps_fullposIsCodon.py <fubar.json> <out.tsv>")

    json_path, out_tsv = sys.argv[1], sys.argv[2]

    # --------------------
    # Load JSON
    # --------------------
    try:
        with open(json_path, "r") as f:
            j = json.load(f)
    except FileNotFoundError:
        die(f"JSON not found: {json_path}")
    except json.JSONDecodeError as e:
        die(f"JSON decode error in {json_path}: {e}")

    # --------------------
    # Extract table rows
    # --------------------
    rows, idx = load_rows(j)

    # --------------------
    # Write TSV
    # --------------------
    with open(out_tsv, "w") as out:
        out.write("full_pos\tcodon\talpha\tbeta\tomega\tP_pos\tP_neg\tBF_pos\tbeta_minus_alpha\n")

        #  codons 1..N (alignment codon index)
        for i, r in enumerate(rows, start=1):
            a = float(r[idx["alpha"]])
            b = float(r[idx["beta"]])
            ppos = float(r[idx["Prob[alpha<beta]"]])
            pneg = float(r[idx["Prob[alpha>beta]"]])
            bf = float(r[idx["BayesFactor[alpha<beta]"]])
            bma = float(r[idx["beta-alpha"]])

            # omega here is simply beta/alpha (HyPhy-style)
            omega = (b / a) if a > 0 else float("nan")

            # "full_pos" is defined as alignment codon index (same as codon)
            out.write(f"{i}\t{i}\t{a}\t{b}\t{omega}\t{ppos}\t{pneg}\t{bf}\t{bma}\n")

    print(f"[OK] Wrote {out_tsv} with {len(rows)} codons (full_pos==codon)")


if __name__ == "__main__":
    main()
