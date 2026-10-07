#!/usr/bin/env python3
"""
SynTracker APSS -> de novo subspecies clusters (Leiden on an APSS graph)
=======================================================================
Turns the pairwise APSS (average pairwise synteny score) that SynTracker computes
against one reference genome into `Taxon,Cluster` tables, one per
(min_apss, resolution) point of a sweep, in the format of PopPUNK's
`*_clusters.csv` so poppunk/evaluate scores them like any PopPUNK fit.

Recipe (SynTracker paper's network analysis):
  0. One APSS table at one subsampling depth `n` (`avg_synteny_scores_<n>_regions.csv`);
     depths are never mixed, and one table holds exactly one reference genome.
     With `--depth auto` the depth is chosen from the data: for each `n`, count the
     targets with at least one pair and the pairs; take the highest `n` reached before
     either count drops below `--min-retention` x its value at the lowest `n` (the
     cliff). The retention table and plot record the choice.
  1. Clean the pair list: drop missing APSS, keep one row per unordered pair. In a
     per-`n` table every pair must have `Compared_regions >= n` (SynTracker only reports
     pairs with at least `n` regions); that is asserted, not filtered.
  2. Build a graph: one node per SynTracker target, one edge per pair, weight = APSS.
     A pair SynTracker could not compare is a missing edge, not a fake distance.
  3. Prune edges with APSS < min_apss.
  4. Leiden community detection (modularity, weights = APSS, `resolution`, fixed seed).
     Each community is one cluster.
  5. A target left without edges is its own cluster (singleton). Genomes that were not
     SynTracker targets are absent from the tables.
  6. Cluster 1 is the largest. A per-cluster QC table gives size, mean intra-cluster
     APSS (cohesion), max inter-cluster APSS (separation) and `low_confidence` (size <= 2),
     both computed on the cleaned, unpruned pairs.

Known limitation (follow-up): after pruning every weight lies in [min_apss, 1], so
modularity is driven mostly by which edges exist. Rescaling the weights to
(APSS - min_apss) / (1 - min_apss) would restore their dynamic range; not done here.

APSS input columns (SynTracker 1.4.0): Ref_genome, Sample1, Sample2, APSS,
Compared_regions. Sample names are SynTracker's target file basenames without
extension; all names are reconciled with `normalise_genome_id` from
`evaluate_poppunk_fastani.py`, the single place that defines the convention.

Outputs (per species, `<p>` = --out-prefix):
  <p>_syntracker_leiden_n<n>_apss<t>_r<res>_clusters.csv    Taxon,Cluster
  <p>_syntracker_leiden_n<n>_apss<t>_r<res>_cluster_qc.tsv  per-cluster QC
  <p>_syntracker_depth_retention.tsv / .png                 depth choice (always written)

Usage:
    syntracker_apss_clusters.py \
        --apss        avg_synteny_scores_40_regions.csv avg_synteny_scores_60_regions.csv ... \
        --targets     targets.txt \
        --depth       auto \
        --min-apss    0.70,0.75,0.80 \
        --resolutions 0.5,1.0,2.0 \
        --out-prefix  B_longum
"""

import argparse
import random
import re
import sys
from pathlib import Path

import pandas as pd

from evaluate_poppunk_fastani import normalise_genome_id

DEPTH_RE = re.compile(r"avg_synteny_scores_(\d+|all)_regions\.csv$")


def table_depth(path: str):
    """Subsampling depth from a SynTracker APSS file name: an int, or 'all'."""
    m = DEPTH_RE.search(Path(path).name)
    if not m:
        sys.exit(f"ERROR: '{path}' is not a SynTracker avg_synteny_scores_<n|all>_regions.csv table")
    return m.group(1) if m.group(1) == "all" else int(m.group(1))


def read_pairs(path: str, depth, allow_empty: bool = False) -> pd.DataFrame:
    """
    Cleaned pair list (g1 < g2, apss) from one APSS table.

    Asserts one reference genome and, for a per-`n` table, Compared_regions >= n.
    """
    try:
        df = pd.read_csv(path)
    except (FileNotFoundError, pd.errors.EmptyDataError) as err:
        sys.exit(f"ERROR: cannot read APSS file '{path}': {err}")
    missing = {"Sample1", "Sample2", "APSS"} - set(df.columns)
    if depth != "all":
        missing |= {"Compared_regions"} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: APSS file '{path}' lacks columns {sorted(missing)}; found {list(df.columns)}")
    if df.empty and not allow_empty:
        sys.exit(f"ERROR: APSS file '{path}' has no pairs (did SynTracker's synteny step fail?)")

    if "Ref_genome" in df.columns and df["Ref_genome"].nunique() > 1:
        sys.exit(
            f"ERROR: APSS file '{path}' mixes {df['Ref_genome'].nunique()} reference genomes; "
            "clusters are only meaningful within one reference, so cluster each one separately."
        )
    if depth != "all":
        regions = pd.to_numeric(df["Compared_regions"], errors="coerce")
        short = df[~(regions >= depth)]
        if not short.empty:
            sys.exit(
                f"ERROR: APSS file '{path}' is the {depth}-region table but {len(short)} pair(s) have "
                f"Compared_regions < {depth} (e.g. {short.iloc[0]['Sample1']} vs {short.iloc[0]['Sample2']}); "
                "this table is not what SynTracker writes for that depth."
            )

    s1 = df["Sample1"].map(normalise_genome_id)
    s2 = df["Sample2"].map(normalise_genome_id)
    pairs = pd.DataFrame({
        "g1": s1.where(s1 <= s2, s2),
        "g2": s2.where(s1 <= s2, s1),
        "apss": pd.to_numeric(df["APSS"], errors="coerce"),
    })
    pairs = pairs[pairs["g1"] != pairs["g2"]].dropna(subset=["apss"])
    return pairs.drop_duplicates(subset=["g1", "g2"], keep="first").reset_index(drop=True)


def restrict_to_targets(pairs: pd.DataFrame, targets: list, path: str) -> pd.DataFrame:
    keep = pairs["g1"].isin(targets) & pairs["g2"].isin(targets)
    if not keep.all():
        print(f"WARNING: {int((~keep).sum())} pair(s) in '{path}' involve non-target genomes; ignored.",
              file=sys.stderr)
    return pairs[keep].reset_index(drop=True)


def retention(pairs: pd.DataFrame) -> tuple:
    """(targets with at least one pair, pairs)."""
    return len(set(pairs["g1"]) | set(pairs["g2"])), len(pairs)


def select_depth(tables: dict, n_targets: int, min_retention: float, out_prefix: str) -> int:
    """
    Highest per-`n` depth before retention falls off a cliff; writes the retention TSV and plot.

    tables: {n: cleaned pairs}. Depths are walked upwards from the lowest; the walk stops at
    the first depth where retained targets or pairs drop below min_retention x the lowest
    depth's value.
    """
    depths = sorted(tables)
    counts = {n: retention(tables[n]) for n in depths}
    base_samples, base_pairs = counts[depths[0]]
    if base_pairs == 0:
        sys.exit(f"ERROR: no APSS pairs at any subsampling depth ({', '.join(map(str, depths))})")

    chosen = depths[0]
    for n in depths[1:]:
        samples, pairs = counts[n]
        if samples < min_retention * base_samples or pairs < min_retention * base_pairs:
            break
        chosen = n

    rows = [
        {"n": n, "retained_samples": counts[n][0], "retained_pairs": counts[n][1],
         "total_targets": n_targets, "selected": n == chosen}
        for n in depths
    ]
    pd.DataFrame(rows).to_csv(f"{out_prefix}_syntracker_depth_retention.tsv", sep="\t", index=False)
    plot_retention(rows, n_targets, chosen, min_retention, f"{out_prefix}_syntracker_depth_retention.png")
    print(f"depth: n={chosen} chosen from {depths} (min retention {min_retention} of n={depths[0]})",
          file=sys.stderr)
    return chosen


def plot_retention(rows: list, n_targets: int, chosen: int, min_retention: float, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    possible = max(n_targets * (n_targets - 1) / 2, 1)
    xs = [r["n"] for r in rows]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(xs, [r["retained_samples"] / max(n_targets, 1) for r in rows], marker="o", label="targets with >= 1 pair")
    ax.plot(xs, [r["retained_pairs"] / possible for r in rows], marker="s", label="pairs (of all target pairs)")
    ax.axvline(chosen, color="grey", linestyle="--", label=f"chosen n = {chosen}")
    ax.set_xlabel("regions subsampled per pair (n)")
    ax.set_ylabel("fraction retained")
    ax.set_ylim(0, 1.05)
    ax.set_xticks(xs)
    ax.set_title(f"SynTracker depth retention (cliff: < {min_retention:g} of lowest n)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


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


def cluster_qc(assignment: dict, pairs: pd.DataFrame, reference: str, depth) -> pd.DataFrame:
    """Per cluster: size, mean intra-cluster APSS, max inter-cluster APSS, low_confidence."""
    c1 = pairs["g1"].map(assignment)
    c2 = pairs["g2"].map(assignment)
    intra = pairs[c1 == c2].groupby(c1[c1 == c2])["apss"].mean()
    inter = pd.concat([
        pd.Series(pairs["apss"][c1 != c2].values, index=c1[c1 != c2].values),
        pd.Series(pairs["apss"][c1 != c2].values, index=c2[c1 != c2].values),
    ])
    inter = inter.groupby(level=0).max() if not inter.empty else pd.Series(dtype=float)
    sizes = pd.Series(assignment).value_counts().sort_index()
    return pd.DataFrame({
        "cluster": sizes.index,
        "size": sizes.values,
        "mean_intra_apss": [intra.get(c, float("nan")) for c in sizes.index],
        "max_inter_apss": [inter.get(c, float("nan")) for c in sizes.index],
        "low_confidence": sizes.values <= 2,
        "reference_genome": reference,
        "n": depth,
    })


def parse_floats(text: str, name: str) -> list:
    tokens = [t.strip() for t in text.split(",") if t.strip()]
    try:
        values = [(t, float(t)) for t in tokens]
    except ValueError:
        sys.exit(f"ERROR: {name} must be comma-separated numbers, got '{text}'")
    if not values:
        sys.exit(f"ERROR: {name} is empty")
    return values


def read_reference(path: str) -> str:
    df = pd.read_csv(path, usecols=lambda c: c == "Ref_genome")
    return str(df["Ref_genome"].iloc[0]) if "Ref_genome" in df.columns and not df.empty else ""


def main():
    parser = argparse.ArgumentParser(description="Cluster SynTracker targets from APSS with Leiden over a parameter sweep.")
    parser.add_argument("--apss", required=True, nargs="+",
                        help="SynTracker avg_synteny_scores_<n|all>_regions.csv table(s) from one run.")
    parser.add_argument("--targets", default=None,
                        help="SynTracker target names, one per line; every target becomes a node. "
                             "Default: the genomes in the APSS table.")
    parser.add_argument("--depth", default="auto",
                        help="'auto' (choose from the per-n tables), a per-n depth (40, 60, ...) or 'all'.")
    parser.add_argument("--min-retention", type=float, default=0.9,
                        help="Auto depth: retained targets and pairs must stay >= this fraction of the lowest n.")
    parser.add_argument("--min-apss", default="0.75", help="Comma-separated edge-pruning floors (default 0.75).")
    parser.add_argument("--resolutions", default="1.0", help="Comma-separated Leiden resolutions (default 1.0).")
    parser.add_argument("--seed", type=int, default=42, help="Leiden random seed.")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix.")
    args = parser.parse_args()

    floors = parse_floats(args.min_apss, "--min-apss")
    resolutions = parse_floats(args.resolutions, "--resolutions")
    if not 0 < args.min_retention <= 1:
        sys.exit(f"ERROR: --min-retention must be in (0, 1], got {args.min_retention}")

    by_depth = {}
    for path in args.apss:
        depth = table_depth(path)
        if depth in by_depth:
            sys.exit(f"ERROR: two APSS tables for depth {depth}: '{by_depth[depth]}' and '{path}'")
        by_depth[depth] = path

    targets = None
    if args.targets:
        targets = sorted({normalise_genome_id(line.strip()) for line in open(args.targets) if line.strip()})

    if args.depth == "auto":
        per_n = {n: p for n, p in by_depth.items() if n != "all"}
        if not per_n:
            sys.exit("ERROR: --depth auto needs the per-n tables (avg_synteny_scores_<n>_regions.csv)")
        tables = {n: read_pairs(p, n, allow_empty=True) for n, p in per_n.items()}
        if targets is None:
            targets = sorted(set().union(*[set(t["g1"]) | set(t["g2"]) for t in tables.values()]))
        tables = {n: restrict_to_targets(t, targets, per_n[n]) for n, t in tables.items()}
        depth = select_depth(tables, len(targets), args.min_retention, args.out_prefix)
        pairs = tables[depth]
    else:
        depth = args.depth if args.depth == "all" else int(args.depth)
        if depth not in by_depth:
            sys.exit(f"ERROR: --depth {args.depth} but no avg_synteny_scores_{depth}_regions.csv among --apss")
        pairs = read_pairs(by_depth[depth], depth)
        if targets is None:
            targets = sorted(set(pairs["g1"]) | set(pairs["g2"]))
        pairs = restrict_to_targets(pairs, targets, by_depth[depth])
    if pairs.empty:
        sys.exit(f"ERROR: no APSS pairs between targets at depth {depth}")
    reference = read_reference(by_depth[depth])

    for t_token, t_value in floors:
        for r_token, r_value in resolutions:
            labels = leiden(pairs, targets, t_value, r_value, args.seed)
            assignment = renumber(targets, labels)
            model = f"syntracker_leiden_n{depth}_apss{t_token}_r{r_token}"

            out = pd.DataFrame(sorted(assignment.items()), columns=["Taxon", "Cluster"])
            out.to_csv(f"{args.out_prefix}_{model}_clusters.csv", index=False)
            qc = cluster_qc(assignment, pairs, reference, depth)
            qc.to_csv(f"{args.out_prefix}_{model}_cluster_qc.tsv", sep="\t", index=False)

            sizes = out["Cluster"].value_counts()
            print(
                f"{model}: {len(sizes)} clusters, {int((sizes == 1).sum())} singletons, "
                f"largest {list(sizes.head(4))}, {len(out)} genomes",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
