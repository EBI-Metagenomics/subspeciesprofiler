#!/usr/bin/env python3
"""
SynTracker APSS -> de novo subspecies clusters (Leiden on an APSS graph)
=======================================================================
Turns the pairwise APSS (average pairwise synteny score) that SynTracker computes
against one reference genome into `Taxon,Cluster` tables, one per
(min_apss, resolution) point of a sweep, in the format of PopPUNK's
`*_clusters.csv` so poppunk/evaluate scores them like any PopPUNK fit.

Input is the all-regions APSS table (`avg_synteny_scores_all_regions.csv`): each pair's
mean synteny score over every region both genomes share. It is deterministic and the
most precise: SynTracker's subsampled tables (n = 40 ... 200 regions per pair) add noise
of about 0.23/sqrt(n) to every pair, which on E. lenta was as large as the real spread
between pairs at n = 40.

Steps:
  1. Read the pairs (one reference genome asserted), keep one row per unordered pair and
     drop pairs compared on fewer than --min-regions regions.
  2. Genome coverage: each target's median Compared_regions over its pairs. Targets
     below --min-genome-coverage x the species median are excluded (`low_coverage`):
     a genome keeps fewer regions when it is fragmented or far from the reference, and
     the regions it keeps are not a random sample, which biases its APSS.
  3. Graph: one node per remaining target, one edge per pair, weight = APSS; edges with
     APSS < min_apss are pruned. Leiden (modularity, weights = APSS, `resolution`, fixed
     seed) gives the clusters; a target left without edges is a singleton.
  4. Noise: the per-region synteny-score SD, estimated from the difference between a
     subsampled table (--apss-subsampled, n regions per pair) and the all-regions table,
     whose variance is SD^2 (1/n - 1/R); the APSS standard error of a typical pair is
     SD / sqrt(median regions per pair). The evaluator uses it as the gap margin.

Outputs (`<p>` = --out-prefix):
  <p>_syntracker_leiden_apss<t>_r<res>_clusters.csv   Taxon,Cluster (cluster 1 = largest)
  <p>_syntracker_genome_coverage.tsv                   genome, median_regions, relative_coverage, status
  <p>_syntracker_noise.tsv                             region_sd, median_regions, apss_se, ...

APSS input columns (SynTracker 1.4.0): Ref_genome, Sample1, Sample2, APSS,
Compared_regions. Names are reconciled with `normalise_genome_id` from
`evaluate_poppunk_fastani.py`, the single place that defines the convention.

Usage:
    syntracker_apss_clusters.py \
        --apss avg_synteny_scores_all_regions.csv --apss-subsampled avg_synteny_scores_40_regions.csv \
        --targets targets.txt --min-apss 0.70,0.75,0.80 --resolutions 0.5,1.0,2.0 --out-prefix B_longum
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_poppunk_fastani import normalise_genome_id

DEFAULT_REGION_SD = 0.226  # per-region synteny-score SD measured on E. lenta (0.226) and B. longum (0.227)


def read_pairs(path: str) -> pd.DataFrame:
    """Pairs (g1 < g2, apss, regions) from one APSS table; one reference genome asserted."""
    try:
        df = pd.read_csv(path)
    except (FileNotFoundError, pd.errors.EmptyDataError) as err:
        sys.exit(f"ERROR: cannot read APSS file '{path}': {err}")
    missing = {"Sample1", "Sample2", "APSS", "Compared_regions"} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: APSS file '{path}' lacks columns {sorted(missing)}; found {list(df.columns)}")
    if df.empty:
        sys.exit(f"ERROR: APSS file '{path}' has no pairs (did SynTracker's synteny step fail?)")
    if "Ref_genome" in df.columns and df["Ref_genome"].nunique() > 1:
        sys.exit(f"ERROR: APSS file '{path}' mixes {df['Ref_genome'].nunique()} reference genomes; "
                 "clusters are only meaningful within one reference.")

    s1 = df["Sample1"].map(normalise_genome_id)
    s2 = df["Sample2"].map(normalise_genome_id)
    pairs = pd.DataFrame({
        "g1": s1.where(s1 <= s2, s2),
        "g2": s2.where(s1 <= s2, s1),
        "apss": pd.to_numeric(df["APSS"], errors="coerce"),
        "regions": pd.to_numeric(df["Compared_regions"], errors="coerce"),
    })
    pairs = pairs[pairs["g1"] != pairs["g2"]].dropna(subset=["apss", "regions"])
    return pairs.drop_duplicates(subset=["g1", "g2"], keep="first").reset_index(drop=True)


def genome_coverage(pairs: pd.DataFrame, targets: list, min_relative: float) -> pd.DataFrame:
    """Each target's median Compared_regions over its pairs, relative to the species median."""
    both = pd.concat([pairs[["g1", "regions"]].rename(columns={"g1": "genome"}),
                      pairs[["g2", "regions"]].rename(columns={"g2": "genome"})])
    median = both.groupby("genome")["regions"].median().reindex(targets)
    species = float(median.median()) if median.notna().any() else np.nan
    cov = pd.DataFrame({"genome": targets, "median_regions": median.to_numpy()})
    cov["relative_coverage"] = cov["median_regions"] / species
    low = cov["relative_coverage"].isna() | (cov["relative_coverage"] < min_relative)
    cov["status"] = np.where(low, "low_coverage", "ok")
    return cov


def estimate_noise(pairs: pd.DataFrame, subsampled_path) -> dict:
    """Per-region SD from a subsampled table vs the all-regions table, and the typical pair's APSS SE."""
    median_regions = float(pairs["regions"].median())
    region_sd, source, n = DEFAULT_REGION_SD, "default", np.nan
    if subsampled_path:
        sub = read_pairs(subsampled_path)
        n = float(sub["regions"].median())
        both = sub.merge(pairs, on=["g1", "g2"], suffixes=("_n", "_all"))
        both = both[both["regions_all"] > both["regions_n"]]
        scale = (1.0 / both["regions_n"] - 1.0 / both["regions_all"]).mean()
        estimate = float(np.sqrt((both["apss_n"] - both["apss_all"]).var() / scale)) if len(both) >= 10 and scale > 0 else np.nan
        if np.isfinite(estimate) and estimate > 0:
            region_sd, source = estimate, Path(subsampled_path).name
        else:
            print(f"WARNING: too few pairs to estimate the per-region SD from '{subsampled_path}'; "
                  f"using the default {DEFAULT_REGION_SD}.", file=sys.stderr)
    return {
        "region_sd": region_sd,
        "median_regions": median_regions,
        "apss_se": region_sd / np.sqrt(median_regions),
        "sd_source": source,
        "subsampled_n": n,
    }


def leiden(pairs: pd.DataFrame, targets: list, min_apss: float, resolution: float, seed: int) -> list:
    """Leiden communities (one label per target, in `targets` order) on the APSS graph pruned at min_apss."""
    import igraph as ig

    index = {genome: i for i, genome in enumerate(targets)}
    kept = pairs[pairs["apss"] >= min_apss]
    graph = ig.Graph(n=len(targets), edges=list(zip(kept["g1"].map(index), kept["g2"].map(index))))
    graph.es["weight"] = kept["apss"].astype(float).tolist()
    random.seed(seed)  # igraph draws from Python's `random`
    partition = graph.community_leiden(
        objective_function="modularity", weights="weight", resolution=resolution, n_iterations=-1
    )
    return list(partition.membership)


def renumber(genomes, labels) -> dict:
    """Map each genome to a cluster id: 1 = largest cluster, ties broken by first genome name."""
    groups = {}
    for genome, label in zip(genomes, labels):
        groups.setdefault(label, []).append(genome)
    ordered = sorted(groups.values(), key=lambda members: (-len(members), min(members)))
    return {genome: cid for cid, members in enumerate(ordered, start=1) for genome in members}


def parse_floats(text: str, name: str) -> list:
    tokens = [t.strip() for t in text.split(",") if t.strip()]
    try:
        values = [(t, float(t)) for t in tokens]
    except ValueError:
        sys.exit(f"ERROR: {name} must be comma-separated numbers, got '{text}'")
    if not values:
        sys.exit(f"ERROR: {name} is empty")
    return values


def main():
    parser = argparse.ArgumentParser(description="Cluster SynTracker targets from all-regions APSS with Leiden over a parameter sweep.")
    parser.add_argument("--apss", required=True, help="SynTracker avg_synteny_scores_all_regions.csv.")
    parser.add_argument("--apss-subsampled", default=None,
                        help="A subsampled table (avg_synteny_scores_<n>_regions.csv) from the same run, "
                             f"used only to estimate the per-region score SD (default {DEFAULT_REGION_SD} without it).")
    parser.add_argument("--targets", default=None,
                        help="SynTracker target names, one per line; every target is a graph node. "
                             "Default: the genomes in the APSS table.")
    parser.add_argument("--min-regions", type=int, default=100,
                        help="Drop pairs compared on fewer regions (default 100).")
    parser.add_argument("--min-genome-coverage", type=float, default=0.5,
                        help="Exclude targets whose median regions per pair is below this fraction of "
                             "the species median (default 0.5).")
    parser.add_argument("--min-apss", default="0.75", help="Comma-separated edge-pruning floors (default 0.75).")
    parser.add_argument("--resolutions", default="1.0", help="Comma-separated Leiden resolutions (default 1.0).")
    parser.add_argument("--seed", type=int, default=42, help="Leiden random seed.")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix.")
    args = parser.parse_args()

    floors = parse_floats(args.min_apss, "--min-apss")
    resolutions = parse_floats(args.resolutions, "--resolutions")
    if not 0 <= args.min_genome_coverage <= 1:
        sys.exit(f"ERROR: --min-genome-coverage must be in [0, 1], got {args.min_genome_coverage}")

    pairs = read_pairs(args.apss)
    if args.targets:
        targets = sorted({normalise_genome_id(line.strip()) for line in open(args.targets) if line.strip()})
    else:
        targets = sorted(set(pairs["g1"]) | set(pairs["g2"]))
    outside = ~(pairs["g1"].isin(targets) & pairs["g2"].isin(targets))
    if outside.any():
        print(f"WARNING: {int(outside.sum())} pair(s) involve non-target genomes; ignored.", file=sys.stderr)
        pairs = pairs[~outside]

    noise = estimate_noise(pairs, args.apss_subsampled)
    short = pairs["regions"] < args.min_regions
    print(f"pairs: {len(pairs)}, {int(short.sum())} below {args.min_regions} regions dropped", file=sys.stderr)
    pairs = pairs[~short]
    if pairs.empty:
        sys.exit(f"ERROR: no APSS pairs with at least {args.min_regions} regions")

    cov = genome_coverage(pairs, targets, args.min_genome_coverage)
    cov.to_csv(f"{args.out_prefix}_syntracker_genome_coverage.tsv", sep="\t", index=False, float_format="%.4f")
    kept = cov.loc[cov["status"].eq("ok"), "genome"].tolist()
    if len(kept) < len(targets):
        print(f"WARNING: {len(targets) - len(kept)} low-coverage target(s) excluded: "
              f"{', '.join(cov.loc[cov['status'].ne('ok'), 'genome'][:10])}", file=sys.stderr)
    if len(kept) < 2:
        sys.exit("ERROR: fewer than 2 targets with enough coverage to cluster")
    pairs = pairs[pairs["g1"].isin(kept) & pairs["g2"].isin(kept)]

    noise["n_pairs"] = len(pairs)
    noise["n_targets"] = len(targets)
    noise["n_low_coverage"] = len(targets) - len(kept)
    pd.DataFrame([noise]).to_csv(f"{args.out_prefix}_syntracker_noise.tsv", sep="\t", index=False, float_format="%.5f")
    print(f"noise: per-region SD {noise['region_sd']:.3f} ({noise['sd_source']}), "
          f"APSS SE {noise['apss_se']:.4f} at {noise['median_regions']:.0f} regions", file=sys.stderr)

    for t_token, t_value in floors:
        for r_token, r_value in resolutions:
            assignment = renumber(kept, leiden(pairs, kept, t_value, r_value, args.seed))
            model = f"syntracker_leiden_apss{t_token}_r{r_token}"
            out = pd.DataFrame(sorted(assignment.items()), columns=["Taxon", "Cluster"])
            out.to_csv(f"{args.out_prefix}_{model}_clusters.csv", index=False)
            sizes = out["Cluster"].value_counts()
            print(f"{model}: {len(sizes)} clusters, {int((sizes == 1).sum())} singletons, "
                  f"largest {list(sizes.head(4))}, {len(out)} genomes", file=sys.stderr)


if __name__ == "__main__":
    main()
