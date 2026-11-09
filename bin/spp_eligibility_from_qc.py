#!/usr/bin/env python3
"""
SPP Eligibility from a per-species genome QC table
==================================================
Parses one species' per-genome completeness/contamination CSV (the `qc_csv`
column of the pipeline samplesheet), classifies each genome as HQ / MQ /
DISCARDED, drops the HQ/MQ genomes that are too fragmented (see Contiguity),
and assigns an SPP eligibility label to the species group from the genomes
that pass both criteria. Every downstream branch uses this one set of genomes.

Writes a single-row `spp_eligibility_report.tsv` with columns:
    species_name  spp_label  n_hq  n_mq  n_eff  hq_ratio
    n_total_passing_qc  n_discarded  n_discarded_contiguity

Usage:
    spp_eligibility_from_qc.py \
        --species_name "Escherichia coli" \
        --qc_csv       ecoli_qc.csv \
        --stats        ecoli_seqkit_stats.tsv \
        --output       spp_eligibility_report.tsv

QC CSV format (comma-separated, one row per genome, header required):
    genome,completeness,contamination
    genome1.fna.gz,99.5,0.5
    genome2.fna.gz,98.2,1.1

    - The completeness and contamination columns are matched case-insensitively
      by name (any column whose header contains 'completeness'/'contamination'),
      so CheckM/CheckM2-style headers work directly.
    - The genome column holds the assembly *filename* (with extension); the
      PopPUNK sample name is that filename with the FASTA extension stripped.

Classification thresholds:
    - HQ: completeness >= 90 and contamination <= 1%
    - MQ: completeness >= 80 and contamination <= 5%
    - otherwise DISCARDED (a >=90%-complete genome with 1-5% contamination is
      demoted to MQ rather than dropped).

Contiguity (when --stats is given; `seqkit stats --all -N 90 --tabular` of the
species' genome files):
    - an HQ/MQ genome is DISCARDED if its N50 (default) or N90
      (--contiguity_metric N90) is below --min_contiguity bp (default 5000, the
      5-kb SynTracker window: shorter contigs give it no regions). A genome at
      the threshold passes. Without --stats there is no contiguity filter.
    - every genome of the QC CSV must have a row in the stats table.

In addition to the species eligibility report, this writes two per-genome files
covering the genomes that pass (DISCARDED excluded):
    - a PopPUNK r-file (`--r-files` input): no header, TAB-separated,
      `sample_name<TAB><genomes_prefix>/<filename>` (default prefix `./genomes`,
      matching the `stageAs: 'genomes'` used by the poppunk/createdb module).
    - a labels CSV (`genome,label`) consumed by the ANI evaluator.
and, with --genome_qc_output, an audit table of every genome of the QC CSV:
    genome  file  completeness  contamination  N50  N90  quality  label  reason
where `quality` is the class from completeness/contamination alone, `label` the
final class, and `reason` is `pass`, `quality` or `contiguity`.

The classification thresholds and labelling rules are kept in sync with
`bin/spp_species_eligibility.py`.
"""

import argparse
import csv
import os
import sys

# FASTA suffixes stripped to derive a sample name from an assembly filename.
# Kept in sync with normalise_genome_id() in bin/evaluate_poppunk_fastani.py.
FASTA_SUFFIXES = [
    ".fasta.gz", ".fa.gz", ".fna.gz",
    ".fasta", ".fa", ".fna",
    ".fas", ".contigs",
]


def strip_fasta_ext(filename):
    """Return the assembly filename with a single trailing FASTA extension removed."""
    name = os.path.basename(str(filename))
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name

# ---------------------------------------------------------------------------
# SPP thresholds (kept in sync with bin/spp_species_eligibility.py)
# ---------------------------------------------------------------------------
MIN_COMP = 80.0
HQ_COMP  = 90.0
HQ_CONT  = 1.0   # max contamination allowed for HQ; MAX_CONT (5.0) is the MQ ceiling
MAX_CONT = 5.0

CONTIGUITY_METRICS = ("N50", "N90")
MIN_CONTIGUITY = 5000  # bp; the SynTracker window length

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
    """Find the genome, completeness and contamination columns case-insensitively."""
    genome_col = comp_col = cont_col = None
    for name in fieldnames or []:
        low = name.strip().lower()
        if genome_col is None and low in ("genome", "genome_id", "sample", "id", "name"):
            genome_col = name
        if comp_col is None and "completeness" in low:
            comp_col = name
        if cont_col is None and "contamination" in low:
            cont_col = name
    # Fall back to the first column for the genome identifier.
    if genome_col is None and fieldnames:
        genome_col = fieldnames[0]
    missing = [n for n, c in (("completeness", comp_col), ("contamination", cont_col)) if c is None]
    if missing:
        sys.exit(
            "ERROR: QC CSV must contain completeness and contamination columns; "
            f"missing {missing}. Found columns: {list(fieldnames or [])}"
        )
    return genome_col, comp_col, cont_col


def read_stats(path):
    """seqkit stats --all -N 90 --tabular -> {file basename: {"N50": int, "N90": int}}."""
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            missing = [c for c in ("file", *CONTIGUITY_METRICS) if c not in (reader.fieldnames or [])]
            if missing:
                sys.exit(f"ERROR: seqkit stats '{path}' lacks columns {missing} (run seqkit stats with --all -N 90 --tabular)")
            stats = {}
            for row in reader:
                values = {}
                for metric in CONTIGUITY_METRICS:
                    try:
                        values[metric] = int(float(row[metric]))
                    except (TypeError, ValueError):
                        values[metric] = 0
                stats[os.path.basename(row["file"])] = values
    except FileNotFoundError as err:
        sys.exit(f"ERROR: cannot read seqkit stats '{path}': {err}")
    return stats


def count_qc(qc_csv, min_comp, hq_comp, hq_cont, max_cont, stats=None, metric="N50", min_contiguity=MIN_CONTIGUITY):
    """Classify every genome; return counts, the passing rows and the per-genome audit rows.

    Passing rows are tuples of (sample_name, filename, label) for the genomes that
    pass both quality and contiguity (DISCARDED excluded). With stats=None there
    is no contiguity filter.
    """
    n_hq = n_mq = n_disc = n_disc_contig = 0
    passing = []
    audit = []
    missing = []
    with open(qc_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        genome_col, comp_col, cont_col = resolve_columns(reader.fieldnames)
        for row in reader:
            filename = row.get(genome_col)
            comp, cont = row.get(comp_col), row.get(cont_col)
            quality = classify_genome(comp, cont, min_comp, hq_comp, hq_cont, max_cont)
            nx = {}
            if stats is not None:
                nx = stats.get(os.path.basename(str(filename)))
                if nx is None:
                    missing.append(filename)
                    nx = {}
            label, reason = quality, "pass"
            if quality == "DISCARDED":
                reason = "quality"
            elif stats is not None and nx.get(metric, 0) < min_contiguity:
                label, reason = "DISCARDED", "contiguity"
            audit.append([strip_fasta_ext(filename), filename, comp, cont,
                          nx.get("N50", ""), nx.get("N90", ""), quality, label, reason])
            if label == "HQ":
                n_hq += 1
            elif label == "MQ":
                n_mq += 1
            else:
                n_disc += 1
                n_disc_contig += reason == "contiguity"
                continue
            passing.append((strip_fasta_ext(filename), filename, label))
    if missing:
        sys.exit(f"ERROR: {len(missing)} genomes of the QC CSV have no row in the seqkit stats "
                 f"(missing from genomes_dir?): {sorted(missing)[:5]}")
    return n_hq, n_mq, n_disc, n_disc_contig, passing, audit


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
    parser.add_argument("--rfile_output",
                        help="Output PopPUNK r-file (TAB-separated sample_name<TAB>path, "
                             "HQ/MQ genomes only). If omitted, no r-file is written.")
    parser.add_argument("--labels_output",
                        help="Output labels CSV (genome,label for HQ/MQ genomes) for the "
                             "ANI evaluator. If omitted, no labels file is written.")
    parser.add_argument("--genomes_prefix", default="./genomes",
                        help="Path prefix written before each filename in the r-file; must match "
                             "the poppunk/createdb `stageAs` directory (default: ./genomes).")
    parser.add_argument("--genome_qc_output",
                        help="Output per-genome audit TSV (every genome of the QC CSV with its "
                             "classes and the reason it passed or was discarded).")
    parser.add_argument("--stats",
                        help="seqkit stats --all -N 90 --tabular of the species' genome files. "
                             "If omitted, there is no contiguity filter.")
    parser.add_argument("--contiguity_metric", choices=CONTIGUITY_METRICS, default="N50",
                        help="Contiguity statistic an HQ/MQ genome must reach (default: N50).")
    parser.add_argument("--min_contiguity", type=int, default=MIN_CONTIGUITY,
                        help=f"Minimum --contiguity_metric, in bp (default: {MIN_CONTIGUITY}).")
    parser.add_argument("--min_completeness",  type=float, default=MIN_COMP)
    parser.add_argument("--hq_completeness",   type=float, default=HQ_COMP)
    parser.add_argument("--hq_contamination",  type=float, default=HQ_CONT)
    parser.add_argument("--max_contamination", type=float, default=MAX_CONT)
    args = parser.parse_args()

    stats = read_stats(args.stats) if args.stats else None
    n_hq, n_mq, n_disc, n_disc_contig, passing, audit = count_qc(
        args.qc_csv, args.min_completeness, args.hq_completeness,
        args.hq_contamination, args.max_contamination,
        stats, args.contiguity_metric, args.min_contiguity,
    )
    spp_label, n_eff, hq_ratio = compute_label(n_hq, n_mq)
    n_passing = n_hq + n_mq

    prefix = args.genomes_prefix.rstrip("/")
    if args.rfile_output:
        # PopPUNK r-file: no header, TAB-separated. Use "\n" line endings (not the
        # csv default "\r\n") so a trailing CR never corrupts the genome paths.
        with open(args.rfile_output, "w", newline="") as rf:
            writer = csv.writer(rf, delimiter="\t", lineterminator="\n")
            for sample_name, filename, _label in passing:
                writer.writerow([sample_name, f"{prefix}/{filename}"])

    if args.labels_output:
        with open(args.labels_output, "w", newline="") as lf:
            writer = csv.writer(lf, lineterminator="\n")
            writer.writerow(["genome", "label"])
            for sample_name, _filename, label in passing:
                writer.writerow([sample_name, label])

    if args.genome_qc_output:
        with open(args.genome_qc_output, "w", newline="") as gf:
            writer = csv.writer(gf, delimiter="\t", lineterminator="\n")
            writer.writerow(["genome", "file", "completeness", "contamination",
                             "N50", "N90", "quality", "label", "reason"])
            writer.writerows(audit)

    columns = [
        "species_name", "spp_label", "n_hq", "n_mq", "n_eff",
        "hq_ratio", "n_total_passing_qc", "n_discarded", "n_discarded_contiguity",
    ]
    values = [
        args.species_name, spp_label, n_hq, n_mq, n_eff,
        hq_ratio, n_passing, n_disc, n_disc_contig,
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(values)

    print(
        f"{args.species_name}: {n_hq} HQ, {n_mq} MQ, {n_disc} discarded "
        f"({n_disc - n_disc_contig} quality, {n_disc_contig} contiguity"
        + (f": {args.contiguity_metric} < {args.min_contiguity} bp" if stats is not None else ", no stats given")
        + f") (n_eff={n_eff}, hq_ratio={hq_ratio}) -> {spp_label}",
        file=sys.stderr,
    )
    print(f"Report written to: {args.output}", file=sys.stderr)
    if args.rfile_output:
        print(f"PopPUNK r-file written to: {args.rfile_output} ({n_passing} passing HQ/MQ genomes)", file=sys.stderr)
    if args.labels_output:
        print(f"Labels written to: {args.labels_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
