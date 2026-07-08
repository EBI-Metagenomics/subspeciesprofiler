#!/usr/bin/env python3
"""
SPP Eligibility from a per-species genome QC table
==================================================
Parses one species' per-genome completeness/contamination CSV (the `qc_csv`
column of the pipeline samplesheet), classifies each genome as HQ / MQ /
DISCARDED, and assigns an SPP eligibility label to the species group.

Writes a single-row `spp_eligibility_report.tsv` with columns:
    species_name  spp_label  n_hq  n_mq  n_eff  hq_ratio
    n_total_passing_qc  n_discarded

Usage:
    spp_eligibility_from_qc.py \
        --species_name "Escherichia coli" \
        --qc_csv       ecoli_qc.csv \
        --output       spp_eligibility_report.tsv

QC CSV format (comma-separated, one row per genome, header required):
    genome,completeness,contamination
    genome1,99.5,0.5
    genome2,98.2,1.1

    - The completeness and contamination columns are matched case-insensitively
      by name (any column whose header contains 'completeness'/'contamination'),
      so CheckM/CheckM2-style headers work directly.

Classification thresholds:
    - HQ: completeness >= 90 and contamination <= 1%
    - MQ: completeness >= 80 and contamination <= 5%
    - otherwise DISCARDED (a >=90%-complete genome with 1-5% contamination is
      demoted to MQ rather than dropped).

The classification thresholds and labelling rules are kept in sync with
`bin/spp_species_eligibility.py`.
"""

import argparse
import csv
import sys

# ---------------------------------------------------------------------------
# SPP thresholds (kept in sync with bin/spp_species_eligibility.py)
# ---------------------------------------------------------------------------
MIN_COMP = 80.0
HQ_COMP  = 90.0
HQ_CONT  = 1.0   # max contamination allowed for HQ; MAX_CONT (5.0) is the MQ ceiling
MAX_CONT = 5.0

STRONG_N_HQ  = 75
STRONG_N_EFF = 100
ACCEPT_N_HQ  = 50
ACCEPT_N_EFF = 80

BOUNDARY_HQ_MARGIN  = 5
BOUNDARY_EFF_MARGIN = 8
NEAR_THRESH_N_HQ    = 40
NEAR_THRESH_N_EFF   = 64


# ---------------------------------------------------------------------------
# QC and inclusion logic (mirrors spp_species_eligibility.py)
# ---------------------------------------------------------------------------

def classify_genome(comp, cont, min_comp, hq_comp, hq_cont, max_cont):
    try:
        comp = float(comp)
        cont = float(cont)
    except (ValueError, TypeError):
        return "DISCARDED"
    if comp >= hq_comp and cont <= hq_cont:      # HQ: >=90 completeness, <=1% contamination
        return "HQ"
    elif comp >= min_comp and cont <= max_cont:  # MQ: >=80 completeness, <=5% contamination
        return "MQ"
    return "DISCARDED"


def compute_label(n_hq, n_mq):
    n_eff    = n_hq + 0.5 * n_mq
    hq_ratio = n_hq / (n_hq + n_mq) if (n_hq + n_mq) > 0 else 0.0

    if n_hq >= STRONG_N_HQ and n_eff >= STRONG_N_EFF:
        base = "STRONG"
    elif n_hq >= ACCEPT_N_HQ and n_eff >= ACCEPT_N_EFF:
        base = "ACCEPTABLE"
    elif n_hq >= ACCEPT_N_HQ:
        base = "NOT_ELIGIBLE"
    else:
        base = "DISCARDED"

    suffix = ""
    if base in ("STRONG", "ACCEPTABLE"):
        thr_hq  = STRONG_N_HQ  if base == "STRONG" else ACCEPT_N_HQ
        thr_eff = STRONG_N_EFF if base == "STRONG" else ACCEPT_N_EFF
        if (n_hq - thr_hq) < BOUNDARY_HQ_MARGIN or (n_eff - thr_eff) < BOUNDARY_EFF_MARGIN:
            suffix = "_BOUNDARY"
    elif base == "DISCARDED":
        if n_hq >= NEAR_THRESH_N_HQ and n_eff >= NEAR_THRESH_N_EFF:
            suffix = "_NEAR_THRESHOLD"

    return base + suffix, round(n_eff, 1), round(hq_ratio, 3)


# ---------------------------------------------------------------------------
# QC table parsing
# ---------------------------------------------------------------------------

def resolve_columns(fieldnames):
    """Find the completeness and contamination columns case-insensitively."""
    comp_col = cont_col = None
    for name in fieldnames or []:
        low = name.strip().lower()
        if comp_col is None and "completeness" in low:
            comp_col = name
        if cont_col is None and "contamination" in low:
            cont_col = name
    missing = [n for n, c in (("completeness", comp_col), ("contamination", cont_col)) if c is None]
    if missing:
        sys.exit(
            "ERROR: QC CSV must contain completeness and contamination columns; "
            f"missing {missing}. Found columns: {list(fieldnames or [])}"
        )
    return comp_col, cont_col


def count_qc(qc_csv, min_comp, hq_comp, hq_cont, max_cont):
    n_hq = n_mq = n_disc = 0
    with open(qc_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        comp_col, cont_col = resolve_columns(reader.fieldnames)
        for row in reader:
            label = classify_genome(row.get(comp_col), row.get(cont_col), min_comp, hq_comp, hq_cont, max_cont)
            if label == "HQ":
                n_hq += 1
            elif label == "MQ":
                n_mq += 1
            else:
                n_disc += 1
    return n_hq, n_mq, n_disc


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Assess SPP eligibility for one species from its genome QC table."
    )
    parser.add_argument("--species_name", required=True,
                        help="Species name (the samplesheet 'species' value).")
    parser.add_argument("--qc_csv", required=True,
                        help="Per-genome completeness/contamination CSV for this species.")
    parser.add_argument("--output", default="spp_eligibility_report.tsv",
                        help="Output report TSV (default: spp_eligibility_report.tsv).")
    parser.add_argument("--min_completeness",  type=float, default=MIN_COMP)
    parser.add_argument("--hq_completeness",   type=float, default=HQ_COMP)
    parser.add_argument("--hq_contamination",  type=float, default=HQ_CONT)
    parser.add_argument("--max_contamination", type=float, default=MAX_CONT)
    args = parser.parse_args()

    n_hq, n_mq, n_disc = count_qc(
        args.qc_csv, args.min_completeness, args.hq_completeness,
        args.hq_contamination, args.max_contamination
    )
    spp_label, n_eff, hq_ratio = compute_label(n_hq, n_mq)
    n_passing = n_hq + n_mq

    columns = [
        "species_name", "spp_label", "n_hq", "n_mq", "n_eff",
        "hq_ratio", "n_total_passing_qc", "n_discarded",
    ]
    values = [
        args.species_name, spp_label, n_hq, n_mq, n_eff,
        hq_ratio, n_passing, n_disc,
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(columns)
        writer.writerow(values)

    print(
        f"{args.species_name}: {n_hq} HQ, {n_mq} MQ, {n_disc} discarded "
        f"(n_eff={n_eff}, hq_ratio={hq_ratio}) -> {spp_label}",
        file=sys.stderr,
    )
    print(f"Report written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
