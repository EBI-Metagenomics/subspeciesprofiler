#!/usr/bin/env python3
"""
Compare two clusterings of the same genomes
===========================================
Compares two ``Taxon,Cluster`` tables, for example a PopPUNK fit and a SynTracker
clustering of one species, and reports how well they agree. You run it after the
pipeline; the pipeline does not run it.

Only the genomes in both tables are compared. Genome names go through
``normalise_genome_id`` (bin/evaluate_poppunk_fastani.py), so PopPUNK and SynTracker
names match.

AGREEMENT MEASURES
------------------
fraction_2_in_1        Each method-2 cluster is matched to the method-1 cluster that
                       shares most of its genomes (tie: the larger method-1 cluster,
                       then the first name). A genome agrees when its method-2 cluster
                       matches its own method-1 cluster. Computed over the genomes in
                       method-2 clusters of at least 2 genomes: a singleton always
                       matches, so singletons are counted separately.
fraction_1_in_2        The same rule in the other direction. A split keeps
                       fraction_2_in_1 high and lowers fraction_1_in_2; a merge does
                       the opposite.
pair_agreement_2_in_1  Of the genome pairs that method 2 puts together, the fraction
                       that method 1 also puts together (and _1_in_2 the reverse).
ARI                    Adjusted Rand index: 1 = identical partitions, ~0 = random.

OUTPUTS (written to --outdir, prefixed with --prefix)
-------
  <prefix>_sankey.png / .pdf          Sankey plot, method 1 left and method 2 right
  <prefix>_agreement_summary.tsv      one row with all the measures and counts
  <prefix>_cluster_matches.tsv        the best match of each cluster, both directions
  <prefix>_contingency.tsv            genome count per (cluster_1, cluster_2) pair
  <prefix>_genome_placement.tsv       genome, cluster_1, cluster_2, agreement
                                      (agreed / not_agreed / singleton)

DEPENDENCIES
------------
python3, numpy, pandas, matplotlib (all in the pipeline's PopPUNK image)

EXAMPLES
--------
  # PopPUNK fit vs SynTracker clustering
  compare_clusterings.py \
      --clusters1 poppunk/fitmodel/dbscan_D3_mcp0.001/sp_fitmodel/sp_fitmodel_clusters.csv \
      --clusters2 syntracker/clusters/sp_syntracker_leiden_apss0.971_r1.0_clusters.csv \
      --label1 "PopPUNK dbscan_D3_mcp0.001" --label2 "SynTracker apss0.971_r1.0" \
      --outdir comparison --prefix poppunk_vs_syntracker

  # two PopPUNK fits
  compare_clusterings.py --clusters1 dbscan_clusters.csv --clusters2 refine_clusters.csv \
      --label1 dbscan --label2 refine --outdir comparison
"""

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import PathPatch, Rectangle  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_poppunk_fastani import normalise_genome_id  # noqa: E402

logger = logging.getLogger(__name__)

SINGLETONS = "singletons"
SMALL = "small clusters"
GREY = "#bdbdbd"


def read_clusters(path: Path) -> pd.Series:
    """Read a clustering table.

    :param path: CSV with a header; column 1 = genome, column 2 = cluster.
    :return: The cluster of each genome (index = normalised genome name), as strings.
    :raises SystemExit: If the table has fewer than 2 columns or lists a genome twice.
    """
    table = pd.read_csv(path, dtype=str)
    if table.shape[1] < 2:
        sys.exit(f"ERROR: clustering table '{path}' needs 2 columns (genome, cluster).")
    genomes = table.iloc[:, 0].astype(str).str.strip().map(normalise_genome_id)
    clusters = pd.Series(table.iloc[:, 1].astype(str).str.strip().to_numpy(), index=genomes)
    duplicated = sorted(set(clusters.index[clusters.index.duplicated()]))
    if duplicated:
        sys.exit(f"ERROR: clustering table '{path}' lists genomes more than once: {duplicated[:10]}")
    return clusters


def best_matches(contingency: pd.DataFrame) -> dict[str, tuple[str, int]]:
    """Match each column cluster to the row cluster that shares most of its genomes.

    :param contingency: Genome counts; rows = clusters to match to, columns = clusters to match.
    :return: For each column cluster, its best row cluster and the shared genome count.
        A tie goes to the larger row cluster, then to the first name.
    """
    row_sizes = contingency.sum(axis=1)
    matches = {}
    for column in contingency.columns:
        shared = contingency[column]
        best = min(shared.index, key=lambda row: (-shared[row], -row_sizes[row], str(row)))
        matches[column] = (best, int(shared[best]))
    return matches


def natural_key(label: str) -> tuple[int, float, str]:
    """Sort key that puts numeric cluster names in number order (2 before 10)."""
    try:
        return (0, float(label), "")
    except ValueError:
        return (1, 0.0, str(label))


def comb2(values: np.ndarray) -> float:
    """Return the total number of pairs, the sum of n(n-1)/2 over the values."""
    values = np.asarray(values, dtype=float)
    return float(np.sum(values * (values - 1) / 2))


def adjusted_rand_index(contingency: pd.DataFrame) -> float:
    """Compute the adjusted Rand index from a contingency table.

    :param contingency: Genome counts per (cluster_1, cluster_2) pair.
    :return: The ARI; 1.0 for two identical trivial partitions, NaN if undefined.
    """
    n = float(contingency.to_numpy().sum())
    index = comb2(contingency.to_numpy().ravel())
    rows = comb2(contingency.sum(axis=1).to_numpy())
    columns = comb2(contingency.sum(axis=0).to_numpy())
    total = n * (n - 1) / 2
    expected = rows * columns / total if total else 0.0
    maximum = (rows + columns) / 2
    if maximum == expected:
        return 1.0 if index == maximum else float("nan")
    return (index - expected) / (maximum - expected)


def agreement(clusters1: pd.Series, clusters2: pd.Series) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Compare two clusterings of the same genomes.

    :param clusters1: Method-1 cluster of each shared genome.
    :param clusters2: Method-2 cluster of each shared genome (same index).
    :return: The contingency table, the summary measures, and the genome placement table.
    """
    contingency = pd.crosstab(clusters1, clusters2)
    sizes1, sizes2 = clusters1.value_counts(), clusters2.value_counts()
    match_2_to_1 = best_matches(contingency)
    match_1_to_2 = best_matches(contingency.T)

    def fraction(own, other, sizes, matches):
        kept = own[own.map(sizes) >= 2]
        if kept.empty:
            return float("nan")
        return float(np.mean([matches[c][0] == other[g] for g, c in kept.items()]))

    placement = pd.DataFrame({"genome": clusters1.index, "cluster_1": clusters1.to_numpy(),
                              "cluster_2": clusters2.to_numpy()})
    placement["agreement"] = [
        "singleton" if sizes2[c2] == 1 else "agreed" if match_2_to_1[c2][0] == c1 else "not_agreed"
        for c1, c2 in zip(placement["cluster_1"], placement["cluster_2"])
    ]
    order = sorted(range(len(placement)), key=lambda i: (
        natural_key(placement.at[i, "cluster_1"]), natural_key(placement.at[i, "cluster_2"]),
        placement.at[i, "genome"]))
    placement = placement.loc[order].reset_index(drop=True)

    together = comb2(contingency.to_numpy().ravel())
    pairs2, pairs1 = comb2(sizes2.to_numpy()), comb2(sizes1.to_numpy())
    summary = {
        "n_clusters_1": len(sizes1),
        "n_clusters_2": len(sizes2),
        "n_singletons_1": int((sizes1 == 1).sum()),
        "n_singletons_2": int((sizes2 == 1).sum()),
        "fraction_2_in_1": fraction(clusters2, clusters1, sizes2, match_2_to_1),
        "fraction_1_in_2": fraction(clusters1, clusters2, sizes1, match_1_to_2),
        "pair_agreement_2_in_1": together / pairs2 if pairs2 else float("nan"),
        "pair_agreement_1_in_2": together / pairs1 if pairs1 else float("nan"),
        "ARI": adjusted_rand_index(contingency),
    }
    return contingency, summary, placement


def cluster_match_table(contingency: pd.DataFrame) -> pd.DataFrame:
    """List the best match of every cluster in both directions.

    :param contingency: Genome counts per (cluster_1, cluster_2) pair.
    :return: Columns direction, cluster, size, best_match, n_shared, fraction_shared.
    """
    rows = []
    for direction, table in (("2_in_1", contingency), ("1_in_2", contingency.T)):
        sizes = table.sum(axis=0)
        for cluster, (best, shared) in best_matches(table).items():
            rows.append({"direction": direction, "cluster": cluster, "size": int(sizes[cluster]),
                         "best_match": best, "n_shared": shared,
                         "fraction_shared": round(shared / sizes[cluster], 4)})
    return pd.DataFrame(rows).sort_values(["direction", "size", "cluster"], ascending=[False, False, True])


def collapse_small(clusters: pd.Series, min_size: int) -> pd.Series:
    """Merge singletons and clusters smaller than `min_size` into two grey blocks of the plot.

    :param clusters: The cluster of each genome.
    :param min_size: The smallest cluster drawn on its own.
    :return: The cluster label to draw for each genome.
    """
    sizes = clusters.map(clusters.value_counts())
    return clusters.where(sizes >= min_size, np.where(sizes == 1, SINGLETONS, SMALL))


def node_order(sizes: pd.Series) -> list[str]:
    """Order the blocks of one column: largest cluster first, the grey blocks last."""
    real = sorted((c for c in sizes.index if c not in (SMALL, SINGLETONS)), key=lambda c: (-sizes[c], natural_key(c)))
    return real + [c for c in (SMALL, SINGLETONS) if c in sizes.index]


def plot_sankey(clusters1: pd.Series, clusters2: pd.Series, labels: tuple[str, str], title: str,
                path_base: Path) -> None:
    """Draw a two-column Sankey plot: method-1 blocks left, method-2 blocks right.

    Block and band heights are genome counts. Bands take the colour of their method-1 block.
    The method-2 blocks follow the order of their best-matching method-1 block, so that
    fewer bands cross.

    :param clusters1: Method-1 cluster to draw for each genome (after `collapse_small`).
    :param clusters2: Method-2 cluster to draw for each genome (same index).
    :param labels: The names of method 1 and method 2.
    :param title: The plot title.
    :param path_base: Output path without the extension (``.png`` and ``.pdf`` are written).
    """
    flows = pd.crosstab(clusters1, clusters2)
    left = node_order(flows.sum(axis=1))
    position = {c: i for i, c in enumerate(left)}
    right_sizes = flows.sum(axis=0)
    grey_rank = {SMALL: 1, SINGLETONS: 2}
    right = sorted(right_sizes.index, key=lambda c: (
        grey_rank.get(c, 0), position[flows[c].idxmax()], -right_sizes[c], str(c)))

    n = float(len(clusters1))
    gap = max(n * 0.01, 0.5)
    cmap = plt.get_cmap("tab20")
    real_left = [c for c in left if c not in grey_rank]
    colour = {c: cmap(i % cmap.N) for i, c in enumerate(real_left)}
    colour.update({c: GREY for c in grey_rank})

    def layout(nodes, sizes):
        tops, y = {}, 0.0
        for node in nodes:
            tops[node] = y
            y += sizes[node] + gap
        return tops, y - gap

    left_top, left_height = layout(left, flows.sum(axis=1))
    right_top, right_height = layout(right, right_sizes)
    height = max(left_height, right_height)
    shift_left, shift_right = (height - left_height) / 2, (height - right_height) / 2

    fig_height = min(max(5.0, 0.2 * max(len(left), len(right))), 40.0)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    x0, x1, width = 0.0, 1.0, 0.04
    font = max(4.0, min(9.0, 300.0 / max(len(left), len(right))))

    left_offset = {c: left_top[c] + shift_left for c in left}
    right_offset = {c: right_top[c] + shift_right for c in right}
    for a in left:
        for b in right:
            count = flows.at[a, b]
            if count == 0:
                continue
            ya, yb = left_offset[a], right_offset[b]
            xm = (x0 + width + x1) / 2
            vertices = [(x0 + width, ya), (xm, ya), (xm, yb), (x1, yb), (x1, yb + count),
                        (xm, yb + count), (xm, ya + count), (x0 + width, ya + count), (x0 + width, ya)]
            codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.LINETO,
                     MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=colour[a], edgecolor="none", alpha=0.45))
            left_offset[a] += count
            right_offset[b] += count

    # Left blocks take their band colour; right blocks are dark grey; the merged blocks are light grey.
    for nodes, tops, shift, x, align, sizes in (
            (left, left_top, shift_left, x0, "right", flows.sum(axis=1)),
            (right, right_top, shift_right, x1, "left", right_sizes)):
        for node in nodes:
            y = tops[node] + shift
            face = GREY if node in grey_rank else colour[node] if align == "right" else "#4d4d4d"
            ax.add_patch(Rectangle((x, y), width, sizes[node], facecolor=face, edgecolor="black", lw=0.4))
            text_x = x - 0.01 if align == "right" else x + width + 0.01
            ax.text(text_x, y + sizes[node] / 2, f"{node} (n={sizes[node]})", ha=align, va="center",
                    fontsize=font)

    ax.text(x0 + width / 2, -gap * 2, labels[0], ha="center", va="bottom", fontsize=10, weight="bold")
    ax.text(x1 + width / 2, -gap * 2, labels[1], ha="center", va="bottom", fontsize=10, weight="bold")
    ax.set_xlim(-0.45, 1.5)
    ax.set_ylim(height + gap, -gap * 6)
    ax.axis("off")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{path_base}.png", dpi=300)
    fig.savefig(f"{path_base}.pdf")
    plt.close(fig)


def main() -> None:
    """Compare two clustering tables and write the plot and the tables."""
    parser = argparse.ArgumentParser(
        description="Compare two Taxon,Cluster tables: Sankey plot and agreement measures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--clusters1", required=True, type=Path, help="Method-1 clustering (genome, cluster).")
    parser.add_argument("--clusters2", required=True, type=Path, help="Method-2 clustering (genome, cluster).")
    parser.add_argument("--label1", default="method 1", help="Name of method 1 in the plot and the summary.")
    parser.add_argument("--label2", default="method 2", help="Name of method 2 in the plot and the summary.")
    parser.add_argument("--min-size", type=int, default=3,
                        help="Clusters smaller than this are merged into one grey block in the plot.")
    parser.add_argument("--outdir", default=".", type=Path)
    parser.add_argument("--prefix", default="comparison")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    all1, all2 = read_clusters(args.clusters1), read_clusters(args.clusters2)
    shared = sorted(set(all1.index) & set(all2.index))
    only1 = sorted(set(all1.index) - set(all2.index))
    only2 = sorted(set(all2.index) - set(all1.index))
    for name, genomes in ((args.label1, only1), (args.label2, only2)):
        if genomes:
            logger.warning("%d genome(s) are only in %s and are left out: %s%s", len(genomes), name,
                           genomes[:10], " ..." if len(genomes) > 10 else "")
    if len(shared) < 2:
        sys.exit("ERROR: fewer than 2 genomes are in both tables; check the genome names.")

    clusters1, clusters2 = all1.loc[shared], all2.loc[shared]
    contingency, measures, placement = agreement(clusters1, clusters2)
    summary = {"label_1": args.label1, "label_2": args.label2, "n_genomes_1": len(all1),
               "n_genomes_2": len(all2), "n_shared": len(shared), "n_only_1": len(only1),
               "n_only_2": len(only2), **measures}

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / args.prefix
    pd.DataFrame([summary]).to_csv(f"{out}_agreement_summary.tsv", sep="\t", index=False,
                                   float_format="%.4f")
    cluster_match_table(contingency).to_csv(f"{out}_cluster_matches.tsv", sep="\t", index=False)
    long = contingency.stack().rename("n_genomes").reset_index()
    long.columns = ["cluster_1", "cluster_2", "n_genomes"]
    long = long[long["n_genomes"] > 0]
    long = long.loc[sorted(long.index, key=lambda i: (natural_key(long.at[i, "cluster_1"]), -long.at[i, "n_genomes"]))]
    long.to_csv(f"{out}_contingency.tsv", sep="\t", index=False)
    placement.to_csv(f"{out}_genome_placement.tsv", sep="\t", index=False)

    title = (f"{len(shared)} shared genomes; fraction_2_in_1 = {measures['fraction_2_in_1']:.2f}, "
             f"fraction_1_in_2 = {measures['fraction_1_in_2']:.2f}, ARI = {measures['ARI']:.2f}")
    plot_sankey(collapse_small(clusters1, args.min_size), collapse_small(clusters2, args.min_size),
                (args.label1, args.label2), title, Path(f"{out}_sankey"))

    for key, value in summary.items():
        print(f"{key:<24}{value:.4f}" if isinstance(value, float) else f"{key:<24}{value}")


if __name__ == "__main__":
    main()
