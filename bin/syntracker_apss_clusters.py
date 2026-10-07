#!/usr/bin/env python3
"""
SynTracker APSS -> clustering tables
====================================
Turns SynTracker's average pairwise synteny scores (APSS) into one
`Taxon,Cluster` table per APSS threshold, in the format PopPUNK writes, so each
table can be scored by `evaluate_poppunk_fastani.py` like any PopPUNK fit.

Methods (`--method`):
  average    Average linkage (UPGMA) on APSS: repeatedly merge the two clusters
             with the highest mean pairwise APSS while that mean is >= t.
             The default: on real data (B. longum, 588 genomes) it recovers the
             subspecies structure across a wide band of t, because a few
             high-APSS "bridging" pairs between groups cannot merge them.
  connected  Single linkage: genomes linked by any pair with APSS >= t, clusters
             are the connected components. Kept for comparison; it chains
             through bridging pairs and only works in a narrow band of t.

A pair missing from the APSS table (too few comparable regions) is left out of
the averages rather than counted as dissimilar. A genome with no qualifying pair
is a singleton.

Optional dRep propagation (`--drep-cdb`): when SynTracker was run on dRep
representatives only, every genome in a representative's dRep secondary cluster
inherits that representative's cluster. Genomes whose dRep cluster has no
SynTracker target are left out of the tables.

APSS input columns (SynTracker 1.4.0): Ref_genome, Sample1, Sample2, APSS,
Compared_regions. Sample names are SynTracker's target file basenames without
extension; all names are reconciled with `normalise_genome_id` from
`evaluate_poppunk_fastani.py`, the single place that defines the convention.

Output: `<out-prefix>_syntracker_<avg|cc>_apss<t>_clusters.csv` per threshold,
with cluster 1 the largest.

Usage:
    syntracker_apss_clusters.py \
        --apss        avg_synteny_scores_all_regions.csv \
        --thresholds  0.70,0.72,0.74 \
        --out-prefix  B_longum \
        [--method average|connected] [--labels labels.csv] [--drep-cdb Cdb.csv]
"""

import argparse
import sys

import numpy as np
import pandas as pd

from evaluate_poppunk_fastani import normalise_genome_id

METHOD_TAGS = {"average": "avg", "connected": "cc"}


def read_apss(path: str):
    """Return (genomes, symmetric APSS matrix with NaN for missing pairs and 1.0 on the diagonal)."""
    try:
        df = pd.read_csv(path)
    except (FileNotFoundError, pd.errors.EmptyDataError) as err:
        sys.exit(f"ERROR: cannot read APSS file '{path}': {err}")
    missing = {"Sample1", "Sample2", "APSS"} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: APSS file '{path}' lacks columns {sorted(missing)}; found {list(df.columns)}")
    if df.empty:
        sys.exit(f"ERROR: APSS file '{path}' has no pairs (did SynTracker's synteny step fail?)")

    s1 = df["Sample1"].map(normalise_genome_id)
    s2 = df["Sample2"].map(normalise_genome_id)
    pairs = pd.DataFrame({
        "g1": np.where(s1 <= s2, s1, s2),
        "g2": np.where(s1 <= s2, s2, s1),
        "apss": pd.to_numeric(df["APSS"], errors="coerce"),
    })
    pairs = pairs[pairs["g1"] != pairs["g2"]].dropna(subset=["apss"])
    # Several reference genomes would give the same pair more than once: average them.
    pairs = pairs.groupby(["g1", "g2"], as_index=False)["apss"].mean()

    genomes = sorted(set(pairs["g1"]) | set(pairs["g2"]))
    index = pd.Index(genomes)
    i, j = index.get_indexer(pairs["g1"]), index.get_indexer(pairs["g2"])
    sim = np.full((len(genomes), len(genomes)), np.nan)
    sim[i, j] = sim[j, i] = pairs["apss"].to_numpy(dtype=float)
    np.fill_diagonal(sim, 1.0)
    return genomes, sim


def average_linkage_merges(sim: np.ndarray):
    """
    Greedy average-linkage agglomeration on a similarity matrix (NaN = missing pair).

    Returns the merge sequence [(a, b, mean_similarity), ...] where cluster b is
    absorbed into a. The average over two clusters uses only the pairs present,
    via running sums and counts. Cutting at threshold t replays the merges up to
    the first one whose similarity is below t, which is exactly what a separate
    greedy run stopping at t would do (even if the heights are not monotone).
    """
    n = len(sim)
    present = ~np.isnan(sim)
    total = np.where(present, sim, 0.0)
    count = present.astype(float)
    np.fill_diagonal(total, 0.0)
    np.fill_diagonal(count, 0.0)
    active = np.ones(n, dtype=bool)
    merges = []
    for _ in range(n - 1):
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = total / count
        mean[count == 0] = -np.inf
        mean[~active, :] = -np.inf
        mean[:, ~active] = -np.inf
        np.fill_diagonal(mean, -np.inf)
        k = int(np.argmax(mean))  # first maximum -> deterministic tie-break
        a, b = divmod(k, n)
        best = mean[a, b]
        if not np.isfinite(best):
            break  # remaining clusters share no APSS pairs
        if a > b:
            a, b = b, a
        merges.append((a, b, float(best)))
        total[a, :] += total[b, :]
        total[:, a] = total[a, :]
        count[a, :] += count[b, :]
        count[:, a] = count[a, :]
        total[a, a] = count[a, a] = 0.0
        active[b] = False
    return merges


def cut_average(n: int, merges, threshold: float) -> np.ndarray:
    labels = np.arange(n)
    for a, b, similarity in merges:
        if similarity < threshold:
            break
        labels[labels == labels[b]] = labels[a]
    return labels


def cut_connected(sim: np.ndarray, threshold: float) -> np.ndarray:
    n = len(sim)
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    linked = np.triu(np.nan_to_num(sim, nan=-np.inf) >= threshold, k=1)
    for a, b in zip(*np.nonzero(linked)):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    return np.array([find(x) for x in range(n)])


def renumber(genomes, labels) -> dict:
    """Map each genome to a cluster id: 1 = largest cluster, ties broken by first genome name."""
    groups = {}
    for genome, label in zip(genomes, labels):
        groups.setdefault(label, []).append(genome)
    ordered = sorted(groups.values(), key=lambda members: (-len(members), min(members)))
    return {genome: cid for cid, members in enumerate(ordered, start=1) for genome in members}


def read_drep_cdb(path: str) -> dict:
    """genome -> dRep secondary cluster, from dRep's Cdb.csv."""
    cdb = pd.read_csv(path)
    missing = {"genome", "secondary_cluster"} - set(cdb.columns)
    if missing:
        sys.exit(f"ERROR: dRep Cdb file '{path}' lacks columns {sorted(missing)}; found {list(cdb.columns)}")
    return dict(zip(cdb["genome"].map(normalise_genome_id), cdb["secondary_cluster"].astype(str)))


def read_label_genomes(path: str) -> list:
    labels = pd.read_csv(path, sep=None, engine="python")
    column = next((c for c in labels.columns if c.lower() in ("genome", "genome_id", "sample", "id", "name")), labels.columns[0])
    return labels[column].map(normalise_genome_id).tolist()


def propagate(assignment: dict, genome_to_cluster: dict) -> dict:
    """Give every genome in a target's dRep cluster the target's cluster id."""
    cluster_to_target = {}
    for genome, drep_cluster in genome_to_cluster.items():
        if genome in assignment:
            cluster_to_target.setdefault(drep_cluster, genome)
    out = dict(assignment)
    for genome, drep_cluster in genome_to_cluster.items():
        target = cluster_to_target.get(drep_cluster)
        if target is not None and genome not in out:
            out[genome] = assignment[target]
    return out


def main():
    parser = argparse.ArgumentParser(description="Cluster genomes from SynTracker APSS at a sweep of thresholds.")
    parser.add_argument("--apss", required=True, help="SynTracker avg_synteny_scores_*.csv")
    parser.add_argument("--thresholds", required=True, help="Comma-separated APSS thresholds, e.g. 0.70,0.72,0.74")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix.")
    parser.add_argument("--method", choices=sorted(METHOD_TAGS), default="average",
                        help="Clustering method (default: average linkage).")
    parser.add_argument("--labels", default=None,
                        help="Per-genome labels CSV (genome,label). Without --drep-cdb, labelled genomes "
                             "absent from the APSS table are reported as singletons.")
    parser.add_argument("--drep-cdb", default=None,
                        help="dRep Cdb.csv: propagate each target's cluster to its dRep secondary cluster.")
    args = parser.parse_args()

    tokens = [t.strip() for t in args.thresholds.split(",") if t.strip()]
    try:
        thresholds = [(t, float(t)) for t in tokens]
    except ValueError:
        sys.exit(f"ERROR: --thresholds must be comma-separated numbers, got '{args.thresholds}'")
    if not thresholds:
        sys.exit("ERROR: --thresholds is empty")

    genomes, sim = read_apss(args.apss)
    merges = average_linkage_merges(sim) if args.method == "average" else None
    genome_to_cluster = read_drep_cdb(args.drep_cdb) if args.drep_cdb else None
    extra_singletons = []
    if args.labels and genome_to_cluster is None:
        extra_singletons = sorted(set(read_label_genomes(args.labels)) - set(genomes))

    tag = METHOD_TAGS[args.method]
    for token, value in thresholds:
        labels = cut_average(len(genomes), merges, value) if merges is not None else cut_connected(sim, value)
        assignment = renumber(genomes, labels)
        if genome_to_cluster is not None:
            assignment = propagate(assignment, genome_to_cluster)
        next_id = max(assignment.values(), default=0) + 1
        for offset, genome in enumerate(extra_singletons):
            assignment[genome] = next_id + offset

        out = pd.DataFrame(sorted(assignment.items()), columns=["Taxon", "Cluster"])
        out.to_csv(f"{args.out_prefix}_syntracker_{tag}_apss{token}_clusters.csv", index=False)
        sizes = out["Cluster"].value_counts()
        print(
            f"apss>={token} ({args.method}): {len(sizes)} clusters, {int((sizes == 1).sum())} singletons, "
            f"largest {list(sizes.head(4))}, {len(out)} genomes",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
