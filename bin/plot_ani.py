#!/usr/bin/env python3
"""
Per-species ANI overview: pairwise ANI histogram + clustered ANI heatmap
========================================================================
One figure per species from the all-vs-all FastANI the pipeline already computes, with two
panels side by side:

  left   histogram of pairwise ANI (%) -- one count per unordered genome pair -- with a black
         kernel density estimate scaled to counts, to show at a glance how the pairwise
         identities are distributed (one mode = a continuum; several = discrete groups).
  right  heatmap of ANI between every pair of genomes, rows and columns ordered by UPGMA
         hierarchical clustering on 100 - ANI, with the dendrograms on top and left and a
         strip marking each genome HQ or MQ (from the SPECIESQC labels).

Pairs: self-hits are dropped and the two FastANI directions of a pair are averaged (as for the
SynTracker reference centrality). FastANI reports nothing for pairs below ~80% ANI; those are
"missing": grey in the heatmap, counted in the title and the summary, and treated as the most
distant observed pair for the clustering.

Names and the ANI scale come from `read_fastani` / `normalise_genome_id` in
`evaluate_poppunk_fastani.py`, the single place that defines them.

Outputs (`<p>` = --out-prefix):
  <p>_ani_overview.png   the figure (200 dpi)
  <p>_ani_overview.pdf   the same figure, vector
  <p>_ani_overview.tsv   n_genomes, n_pairs, n_missing_pairs, ani_min, ani_median, ani_max,
                         colour_min (all ANI in %)

Usage:
    plot_ani.py --ani B_longum.ani.txt --labels B_longum_labels.csv --title B_longum --out-prefix B_longum
"""

import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from evaluate_poppunk_fastani import normalise_genome_id, read_fastani  # noqa: E402

MAX_LABELLED = 50  # genome names on the heatmap axes up to this many genomes
LABEL_COLOURS = {"HQ": "#1b5e20", "MQ": "#a5d6a7", "other": "#bdbdbd"}


def pairwise(fastani: pd.DataFrame) -> pd.DataFrame:
    """One row per unordered pair (g1 < g2) with the mean ANI (%) of the two directions."""
    df = fastani[fastani["query"] != fastani["reference"]]
    g1 = np.where(df["query"] < df["reference"], df["query"], df["reference"])
    g2 = np.where(df["query"] < df["reference"], df["reference"], df["query"])
    pairs = pd.DataFrame({"g1": g1, "g2": g2, "ani": df["ani"].to_numpy() * 100.0})
    return pairs.groupby(["g1", "g2"], as_index=False)["ani"].mean()


def ani_matrix(pairs: pd.DataFrame, genomes: list) -> pd.DataFrame:
    """Symmetric ANI (%) matrix over `genomes`, diagonal 100, NaN where FastANI reported no pair."""
    m = pd.DataFrame(np.nan, index=genomes, columns=genomes)
    for g1, g2, ani in pairs[["g1", "g2", "ani"]].itertuples(index=False):
        m.at[g1, g2] = ani
        m.at[g2, g1] = ani
    np.fill_diagonal(m.values, 100.0)
    return m


def cluster_order(matrix: pd.DataFrame):
    """UPGMA on 100 - ANI (missing pairs = the most distant observed pair): (linkage, leaf order)."""
    observed = matrix.values[~np.isnan(matrix.values)]
    dist = 100.0 - matrix.fillna(observed.min()).values
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    z = linkage(squareform(np.clip(dist, 0.0, None), checks=False), method="average",
                optimal_ordering=len(matrix) <= 1000)
    order = dendrogram(z, no_plot=True)["leaves"]
    return z, order


def read_labels(path) -> dict:
    """genome -> HQ/MQ from the SPECIESQC labels CSV (`genome,label`); {} without a file."""
    if not path:
        return {}
    labels = pd.read_csv(path)
    if not {"genome", "label"} <= set(labels.columns):
        sys.exit(f"ERROR: labels file '{path}' needs columns genome,label; found {list(labels.columns)}")
    return dict(zip(labels["genome"].map(normalise_genome_id), labels["label"].astype(str)))


def histogram(ax, values: np.ndarray) -> None:
    lo, hi = float(values.min()), float(values.max())
    n_bins = int(np.clip(np.ceil((hi - lo) / 0.1), 10, 200))
    counts, edges, _ = ax.hist(values, bins=n_bins, range=(lo, hi) if hi > lo else (lo - 0.5, hi + 0.5),
                               color="#90a4ae", edgecolor="white", linewidth=0.3)
    if len(np.unique(values)) > 1:
        width = edges[1] - edges[0]
        grid = np.linspace(edges[0], edges[-1], 512)
        ax.plot(grid, gaussian_kde(values)(grid) * len(values) * width, color="black", linewidth=1.5)
    ax.set_xlabel("ANI (%)")
    ax.set_ylabel("Pairwise comparisons")
    ax.set_title("Pairwise ANI distribution", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)


def clustered_heatmap(fig, spec, matrix: pd.DataFrame, labels: dict, colour_min: float) -> None:
    z, order = cluster_order(matrix)
    ordered = matrix.iloc[order, order]
    n = len(ordered)

    labelled = n <= MAX_LABELLED
    # columns: left dendrogram, HQ/MQ strip, heatmap, gap (room for row names), colour bar
    grid = GridSpecFromSubplotSpec(2, 5, subplot_spec=spec,
                                   width_ratios=[0.18, 0.03, 1.0, 0.16 if labelled else 0.03, 0.04],
                                   height_ratios=[0.18, 1.0], wspace=0.02, hspace=0.02)
    ax_corner = fig.add_subplot(grid[0, 0:2])
    ax_top = fig.add_subplot(grid[0, 2])
    ax_left = fig.add_subplot(grid[1, 0])
    ax_strip = fig.add_subplot(grid[1, 1])
    ax_heat = fig.add_subplot(grid[1, 2])
    ax_bar = fig.add_subplot(grid[1, 4])

    # scipy draws leaf i at 10 i + 5: pin the dendrograms to the heatmap cells
    black = {"link_color_func": lambda _: "black", "no_labels": True}
    # thin lines for large species, or hundreds of leaves merge into a solid block
    with plt.rc_context({"lines.linewidth": 1.0 if n <= 100 else 0.25}):
        dendrogram(z, ax=ax_top, orientation="top", **black)
        dendrogram(z, ax=ax_left, orientation="left", **black)
    ax_top.set_xlim(0, 10 * n)
    ax_left.set_ylim(10 * n, 0)  # leaf 0 at the top, like the heatmap rows
    for ax in (ax_corner, ax_top, ax_left):
        ax.axis("off")
    ax_top.set_title(f"ANI between genomes (UPGMA, {n} genomes)", fontsize=11)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e0e0e0")
    image = ax_heat.imshow(np.ma.masked_invalid(ordered.values), cmap=cmap, vmin=colour_min, vmax=100.0,
                           aspect="auto", interpolation="nearest")
    if labelled:
        size = max(4, min(8, 400 // n))
        ax_heat.set_xticks(range(n), ordered.columns, rotation=90, fontsize=size)
        ax_heat.set_yticks(range(n), ordered.index, fontsize=size)
        ax_heat.yaxis.tick_right()
    else:
        ax_heat.set_xticks([])
        ax_heat.set_yticks([])

    classes = [labels.get(g, "other") for g in ordered.index]
    classes = [c if c in LABEL_COLOURS else "other" for c in classes]
    present = [c for c in LABEL_COLOURS if c in classes]
    strip = np.array([[present.index(c)] for c in classes])
    ax_strip.imshow(strip, cmap=ListedColormap([LABEL_COLOURS[c] for c in present]), aspect="auto",
                    interpolation="nearest", vmin=-0.5, vmax=len(present) - 0.5)
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    ax_corner.legend(handles=[Patch(color=LABEL_COLOURS[c], label=c) for c in present], title="Genome",
                     loc="center", fontsize=8, title_fontsize=8, frameon=False)

    fig.colorbar(image, cax=ax_bar, label="ANI (%)")


def main():
    parser = argparse.ArgumentParser(description="Per-species ANI histogram and clustered ANI heatmap.")
    parser.add_argument("--ani", required=True, help="All-vs-all FastANI output (percent ANI).")
    parser.add_argument("--labels", default=None, help="SPECIESQC labels CSV (genome,label) for the HQ/MQ strip.")
    parser.add_argument("--title", default="", help="Figure title (the species name).")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    args = parser.parse_args()

    fastani = read_fastani(args.ani)
    genomes = sorted(set(fastani["query"]) | set(fastani["reference"]))
    pairs = pairwise(fastani)
    n_possible = len(genomes) * (len(genomes) - 1) // 2
    n_missing = n_possible - len(pairs)
    values = pairs["ani"].to_numpy()

    summary = {"n_genomes": len(genomes), "n_pairs": len(pairs), "n_missing_pairs": n_missing,
               "ani_min": np.nan, "ani_median": np.nan, "ani_max": np.nan, "colour_min": np.nan}
    fig = plt.figure(figsize=(16, 7))
    if len(genomes) < 2 or pairs.empty:
        fig.text(0.5, 0.5, f"{args.title}: not enough genomes or ANI pairs to plot "
                           f"({len(genomes)} genomes, {len(pairs)} pairs)", ha="center", va="center")
        print(f"WARNING: {len(genomes)} genomes and {len(pairs)} pairs; nothing to plot.", file=sys.stderr)
    else:
        colour_min = float(np.floor(values.min() * 2) / 2)
        summary.update(ani_min=values.min(), ani_median=float(np.median(values)), ani_max=values.max(),
                       colour_min=colour_min)
        layout = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.35], wspace=0.18)
        histogram(fig.add_subplot(layout[0, 0]), values)
        clustered_heatmap(fig, layout[0, 1], ani_matrix(pairs, genomes), read_labels(args.labels), colour_min)
        missing = f", {n_missing} pair{'s' if n_missing > 1 else ''} without a FastANI hit (grey)" if n_missing else ""
        fig.suptitle(f"{args.title} — {len(genomes)} genomes, {len(pairs)} pairwise comparisons{missing}",
                     fontsize=13)

    for ext, kwargs in (("png", {"dpi": 200}), ("pdf", {})):
        fig.savefig(f"{args.out_prefix}_ani_overview.{ext}", bbox_inches="tight", **kwargs)
    pd.DataFrame([summary]).to_csv(f"{args.out_prefix}_ani_overview.tsv", sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
