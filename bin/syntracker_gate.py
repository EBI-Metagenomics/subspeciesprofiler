#!/usr/bin/env python3
"""
Contiguity gate for SynTracker targets
======================================
SynTracker compares genomes in 5-kb windows around the reference's regions, so contigs shorter
than a window give it nothing: a fragmented genome joins the run but keeps few regions, and its
APSS is biased. This gate keeps the HQ and MQ genomes of one species whose N50 (default) or N90
(--metric N90, stricter: 90% of the assembly in contigs of at least that length) is at least
--min-length bp (default 5000, the SynTracker window).

Outputs (`<p>` = --out-prefix), genomes named by their file names:
  <p>_syntracker_gate.tsv         genome, file, label, N50, N90, metric, min_length, status (pass/fail)
  <p>_syntracker_gate_stats.tsv   the seqkit table of the genomes that pass, plus a `label` column
                                  (the input of syntracker_select.py)

Usage:
    syntracker_gate.py --stats seqkit_stats.tsv --labels labels.csv \
        --metric N50 --min-length 5000 --out-prefix B_longum
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from evaluate_poppunk_fastani import normalise_genome_id, read_labels

METRICS = ("N50", "N90")


def read_stats(path: str) -> pd.DataFrame:
    """seqkit stats --all -N 90 --tabular: one row per FASTA file, with N50 and N90."""
    try:
        df = pd.read_csv(path, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError) as err:
        sys.exit(f"ERROR: cannot read seqkit stats '{path}': {err}")
    missing = {"file", *METRICS} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: seqkit stats '{path}' lacks columns {sorted(missing)} (run with --all -N 90 --tabular)")
    df["file"] = df["file"].map(lambda f: Path(str(f)).name)
    df["genome"] = df["file"].map(normalise_genome_id)
    for metric in METRICS:
        df[metric] = pd.to_numeric(df[metric], errors="coerce").fillna(0).astype(int)
    if df["genome"].duplicated().any():
        sys.exit(f"ERROR: two files map to the same genome name in '{path}': "
                 f"{sorted(df.loc[df['genome'].duplicated(), 'genome'])[:5]}")
    return df


def gate(stats: pd.DataFrame, labels: pd.DataFrame, metric: str, min_length: int) -> pd.DataFrame:
    """Add `label` and `status` (pass if `metric` >= min_length) to the seqkit table."""
    label_of = dict(zip(labels["genome_id"], labels["quality_status"]))
    out = stats.copy()
    out["label"] = out["genome"].map(label_of)
    unlabelled = out.loc[out["label"].isna(), "genome"].tolist()
    if unlabelled:
        sys.exit(f"ERROR: {len(unlabelled)} genomes in the seqkit stats have no label: {sorted(unlabelled)[:5]}")
    out["status"] = (out[metric] >= min_length).map({True: "pass", False: "fail"})
    return out


def main():
    parser = argparse.ArgumentParser(description="Keep the genomes contiguous enough for SynTracker.")
    parser.add_argument("--stats", required=True, help="seqkit stats --all -N 90 --tabular of the HQ/MQ genomes.")
    parser.add_argument("--labels", required=True, help="Per-genome labels CSV (genome,label).")
    parser.add_argument("--metric", choices=METRICS, default="N50", help="Contiguity statistic to gate on (default N50).")
    parser.add_argument("--min-length", type=int, default=5000, help="Minimum value of --metric, in bp (default 5000).")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    args = parser.parse_args()
    if args.min_length < 0:
        sys.exit("ERROR: --min-length must be >= 0")

    table = gate(read_stats(args.stats), read_labels(args.labels), args.metric, args.min_length)

    report = table[["genome", "file", "label", "N50", "N90"]].copy()
    report["metric"] = args.metric
    report["min_length"] = args.min_length
    report["status"] = table["status"]
    report.sort_values(["status", "label", args.metric, "genome"], ascending=[False, True, False, True]) \
        .to_csv(f"{args.out_prefix}_syntracker_gate.tsv", sep="\t", index=False)
    passed = table[table["status"] == "pass"].drop(columns=["genome", "status"])
    passed.to_csv(f"{args.out_prefix}_syntracker_gate_stats.tsv", sep="\t", index=False)

    counts = table.groupby(["label", "status"]).size()
    summary = ", ".join(f"{label} {counts.get((label, 'pass'), 0)}/{counts[label].sum()}"
                        for label in sorted(table["label"].unique()))
    print(f"{args.metric} >= {args.min_length} bp: {len(passed)} of {len(table)} genomes pass ({summary})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
