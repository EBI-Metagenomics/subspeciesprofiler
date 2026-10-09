#!/usr/bin/env python3
"""
SynTracker reference and target selection for one species
==========================================================
Picks the SynTracker reference genome and the target list from the species' HQ
genomes, using the all-vs-all FastANI the pipeline already computes and per-genome
N50 (seqkit stats).

Reference: the HQ genome with the highest species-wide centrality (mean ANI to every
other HQ genome) among those with N50 >= --min-n50. Centrality matters most: SynTracker
only compares regions of the reference, and a target loses regions the further it is
from the reference (BLAST keeps hits at >= 97% identity), so a central reference spreads
that loss evenly across lineages instead of penalising every lineage but its own. The
N50 floor keeps a fragmented genome from being chosen. If no genome passes the floor,
the genome with the highest N50 is used, with a warning.

Targets: every HQ genome, the reference included. Above --max-targets (SynTracker's cost
grows with the square of the targets), the genomes with the highest N50 are kept: contig
breaks remove regions and bias APSS, and the reference is always kept.

Outputs (`<p>` = --out-prefix), genomes named by their file names:
  <p>_syntracker_reference.txt   the reference file name
  <p>_syntracker_targets.txt     the target file names, one per line
  <p>_syntracker_selection.tsv   genome, file, N50, centrality, role (reference/target/capped)

Usage:
    syntracker_select.py --ani fastani.txt --stats seqkit_stats.tsv \
        --min-n50 100000 --max-targets 150 --out-prefix B_longum
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_poppunk_fastani import normalise_genome_id, read_fastani


def read_stats(path: str) -> pd.DataFrame:
    """seqkit stats --all --tabular: one row per FASTA file, with N50."""
    try:
        df = pd.read_csv(path, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError) as err:
        sys.exit(f"ERROR: cannot read seqkit stats '{path}': {err}")
    missing = {"file", "N50"} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: seqkit stats '{path}' lacks columns {sorted(missing)} (run with --all --tabular)")
    out = pd.DataFrame({"file": df["file"].map(lambda f: Path(str(f)).name)})
    out["genome"] = out["file"].map(normalise_genome_id)
    out["N50"] = pd.to_numeric(df["N50"], errors="coerce").fillna(0).astype(int)
    if out["genome"].duplicated().any():
        sys.exit(f"ERROR: two files map to the same genome name in '{path}': "
                 f"{sorted(out.loc[out['genome'].duplicated(), 'genome'])[:5]}")
    return out


def centrality(fastani: pd.DataFrame, genomes: list) -> pd.Series:
    """Mean ANI (0-1) of each genome to every other genome in `genomes`; pairs in both directions are averaged."""
    keep = fastani["query"].isin(genomes) & fastani["reference"].isin(genomes) & (fastani["query"] != fastani["reference"])
    pairs = fastani.loc[keep, ["query", "reference", "ani"]]
    g1 = np.where(pairs["query"] < pairs["reference"], pairs["query"], pairs["reference"])
    g2 = np.where(pairs["query"] < pairs["reference"], pairs["reference"], pairs["query"])
    pairs = pd.DataFrame({"g1": g1, "g2": g2, "ani": pairs["ani"].to_numpy()}).groupby(["g1", "g2"], as_index=False)["ani"].mean()
    both = pd.concat([pairs.rename(columns={"g1": "genome"})[["genome", "ani"]],
                      pairs.rename(columns={"g2": "genome"})[["genome", "ani"]]])
    return both.groupby("genome")["ani"].mean().reindex(genomes)


def select(stats: pd.DataFrame, cent: pd.Series, min_n50: int, max_targets: int):
    """Return (reference genome, target genomes in order, table with roles)."""
    table = stats.copy()
    table["centrality"] = table["genome"].map(cent)

    eligible = table[(table["N50"] >= min_n50) & table["centrality"].notna()]
    if eligible.empty:
        print(f"WARNING: no HQ genome has N50 >= {min_n50} and FastANI values; "
              "the reference is the genome with the highest N50.", file=sys.stderr)
        reference = table.sort_values(["N50", "genome"], ascending=[False, True]).iloc[0]["genome"]
    else:
        reference = eligible.sort_values(["centrality", "N50", "genome"], ascending=[False, False, True]).iloc[0]["genome"]

    by_n50 = table.sort_values(["N50", "genome"], ascending=[False, True])["genome"].tolist()
    ordered = [reference] + [g for g in by_n50 if g != reference]
    targets = ordered[:max_targets]
    if len(ordered) > max_targets:
        print(f"WARNING: {len(ordered)} HQ genomes > --max-targets ({max_targets}); SynTracker runs on the "
              f"{max_targets} with the highest N50, leaving {len(ordered) - max_targets} out.", file=sys.stderr)

    kept = set(targets)
    table["role"] = np.where(table["genome"] == reference, "reference",
                             np.where(table["genome"].isin(kept), "target", "capped"))
    table["order"] = table["role"].map({"reference": 0, "target": 1, "capped": 2})
    return reference, targets, table.sort_values(["order", "N50", "genome"], ascending=[True, False, True])


def main():
    parser = argparse.ArgumentParser(description="Choose the SynTracker reference and targets for one species.")
    parser.add_argument("--ani", required=True, help="All-vs-all FastANI output (percent ANI).")
    parser.add_argument("--stats", required=True, help="seqkit stats --all --tabular of the HQ genomes.")
    parser.add_argument("--min-n50", type=int, default=100000, help="Minimum N50 of the reference (default 100000).")
    parser.add_argument("--max-targets", type=int, default=150, help="Maximum number of targets (default 150).")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    args = parser.parse_args()
    if args.max_targets < 2:
        sys.exit("ERROR: --max-targets must be at least 2")

    stats = read_stats(args.stats)
    if len(stats) < 2:
        sys.exit(f"ERROR: SynTracker needs at least 2 HQ genomes, got {len(stats)}")
    cent = centrality(read_fastani(args.ani), stats["genome"].tolist())
    reference, targets, table = select(stats, cent, args.min_n50, args.max_targets)

    file_of = dict(zip(stats["genome"], stats["file"]))
    Path(f"{args.out_prefix}_syntracker_reference.txt").write_text(file_of[reference] + "\n")
    Path(f"{args.out_prefix}_syntracker_targets.txt").write_text("".join(file_of[g] + "\n" for g in targets))
    table[["genome", "file", "N50", "centrality", "role"]].to_csv(
        f"{args.out_prefix}_syntracker_selection.tsv", sep="\t", index=False, float_format="%.5f")
    row = table.set_index("genome").loc[reference]
    print(f"reference: {reference} (centrality {row['centrality']:.4f}, N50 {row['N50']}); "
          f"{len(targets)} targets", file=sys.stderr)


if __name__ == "__main__":
    main()
