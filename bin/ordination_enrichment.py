#!/usr/bin/env python3
"""
ordination_enrichment.py

Post-pipeline analysis of a gene presence/absence pangenome matrix.

Given a gene_presence_absence.Rtab (from ggCaller or Panaroo), a clustering
table assigning each genome to a cluster, and a per-gene annotation table,
this script produces:

  1) An ordination plot (PCoA on Jaccard distance of the presence/absence
     matrix) with genomes coloured by their cluster. Singleton clusters
     (size 1) are drawn as hollow markers with a dark-grey rim. Each non-singleton cluster with >= 3 genomes
     is drawn with a 90% confidence ellipse.
  2) A pairwise ANOSIM table (Analysis of Similarities) between groups,
     with the R statistic, permutation p-value and BH-adjusted FDR.
  3) A gene-enrichment table based on Fisher's exact tests, with odds ratio,
     p-value, BH-adjusted FDR and direction (enriched / depleted). See the
     enrichment logic below.
  4) An iTOL TREE_COLORS annotation file: one coloured range per genome, in the
     cluster colours of the ordination plot (singletons left uncoloured), to drop onto the
     Newick tree (e.g. the pipeline's <species>_ani_nj_tree.nwk) in iTOL.
  5) With --metadata, for each metadata column: one more ordination plot with the
     genomes coloured by category (cluster ellipses as grey outlines), and one iTOL
     DATASET_COLORSTRIP file, which iTOL draws as an outer ring of the tree.

METADATA
--------
The metadata table has a header row. Column 1 is the genome ID; each other column is
one category variable (e.g. disease_status, genome_type, genome_source). IDs must
equal the Rtab genome names. A genome that the table does not list, or an empty
cell, gets the category 'NA'. 'NA' and 'unknown' are light grey. Each column has its
own colours (tab10, or tab20 above 10 categories), the same in its plot and its ring.

ITOL NODE NAMES
---------------
iTOL applies an annotation only to a leaf whose label equals the node ID exactly;
anything else is silently ignored. The pipeline's tree (bin/fastani_to_tree.py)
names leaves after the genome with every character outside [A-Za-z0-9_.-]
replaced by '_', so the same replacement is applied to the cluster names. With
--tree, every node ID is checked against the tree's leaf labels: unmatched
genomes and uncoloured leaves are reported, only matched genomes are written,
and the script fails if none match. Without --tree the names are not checked.

ENRICHMENT LOGIC
----------------
* If the user supplies a focal cluster (--enrich-cluster):
    - if that cluster has >= --min-cluster-size genomes, run a one-vs-rest
      comparison (focal cluster vs all other genomes);
    - otherwise write 'no statistical power in group <id>' and run nothing.
* If no focal cluster is supplied:
    - compare the two largest (non-singleton) clusters against each other.

INPUTS
------
--rtab          gene_presence_absence.Rtab  (tab-separated; first column = gene
                name, remaining columns = genome IDs, values 0/1)
--clusters      CSV with columns: Taxon,Cluster  (genome id, cluster label)
--annotations   Annotation table. First column = ID, second column = annotation
                string. Delimiter auto-detected; header expected (use
                --anno-no-header if there is none). A key can have several rows
                (e.g. one VF hit and one ARG hit); their annotations are joined.
                Without --gene-families, the ID is an Rtab gene family name.
--gene-families (optional) Panaroo's gene_presence_absence.csv. With it, the
                --annotations ID is a protein ID as in the ggCaller .faa files,
                and each gene family gets the distinct annotations of its
                proteins with their counts, e.g. 'vfdb:fliC (48/51)' (48 of the
                family's 51 proteins carry that annotation).
--tree          (optional) Newick tree to check the iTOL node names against.
--metadata      (optional) Genome metadata table (see METADATA). Delimiter
                auto-detected; header expected.

OUTPUTS (written to --outdir, prefixed with --prefix)
-------
  <prefix>_ordination.png / .pdf   ordination plot with 90% ellipses
  <prefix>_ordination_coords.csv   per-genome PCoA coordinates + cluster (+ metadata)
  <prefix>_anosim.csv              pairwise ANOSIM table
  <prefix>_enrichment.csv          gene enrichment table (or a note if skipped)
  <prefix>_itol_colors.txt         iTOL TREE_COLORS file (coloured range per genome)
  With --metadata, per metadata column <col>:
  <prefix>_ordination_<col>.png / .pdf   ordination plot coloured by category
  <prefix>_itol_<col>.txt                iTOL DATASET_COLORSTRIP file (outer ring)

DEPENDENCIES
------------
python3, numpy, pandas, scipy, matplotlib

EXAMPLES
--------
  # default: ordination + pairwise ANOSIM + two-largest-group enrichment
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations anno.tsv --outdir results --prefix subspecies

  # also check the iTOL colour file against the tree it will be dropped on
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations anno.tsv --tree B_longum_ani_nj_tree.nwk --outdir results

  # one-vs-rest enrichment for a chosen cluster
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations anno.tsv --enrich-cluster 2 --outdir results

  # metadata plots and iTOL rings (disease_status, genome_type, ...)
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations anno.tsv --tree B_longum_ani_nj_tree.nwk \
      --metadata metadata.tsv --outdir results

  # protein-level annotations (e.g. VF/ARG hits keyed by .faa protein ID)
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations vf_arg_hits.tsv --gene-families gene_presence_absence.csv \
      --outdir results
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, chi2, rankdata

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.patches import Ellipse


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def eprint(*args):
    print(*args, file=sys.stderr)


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR, aligned to input order."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    fdr_sorted = ranked * n / (np.arange(n) + 1)
    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    fdr_sorted = np.clip(fdr_sorted, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = fdr_sorted
    return out


def pcoa(dist_square):
    """Principal Coordinates Analysis. Returns (coords n x k, prop per axis)."""
    d = np.asarray(dist_square, dtype=float)
    n = d.shape[0]
    a = -0.5 * (d ** 2)
    j = np.eye(n) - np.ones((n, n)) / n
    b = j @ a @ j
    b = (b + b.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(b)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos = eigvals > 1e-9
    eigvals_pos = eigvals[pos]
    coords = eigvecs[:, pos] * np.sqrt(eigvals_pos)
    prop = eigvals_pos / eigvals_pos.sum()
    return coords, prop


def confidence_ellipse(x, y, ax, conf=0.90, **kwargs):
    """Draw a `conf` confidence ellipse for 2D points (bivariate-normal)."""
    if len(x) < 3:
        return None
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    vals = np.clip(vals, 0, None)
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    scale = np.sqrt(chi2.ppf(conf, df=2))
    width, height = 2 * scale * np.sqrt(vals)
    ell = Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height,
                  angle=theta, **kwargs)
    ax.add_patch(ell)
    return ell


def anosim(dist_square, labels, n_perm=999, rng=None):
    """Two-or-more-group ANOSIM. Returns (R, p_value)."""
    if rng is None:
        rng = np.random.default_rng()
    labels = np.asarray(labels)
    n = len(labels)
    iu = np.triu_indices(n, k=1)
    ranks = rankdata(dist_square[iu])
    denom = n * (n - 1) / 4.0

    def r_stat(lab):
        within = lab[iu[0]] == lab[iu[1]]
        rw = ranks[within].mean()
        rb = ranks[~within].mean()
        return (rb - rw) / denom

    r_obs = r_stat(labels)
    count = 0
    for _ in range(n_perm):
        if r_stat(rng.permutation(labels)) >= r_obs:
            count += 1
    p = (count + 1) / (n_perm + 1)
    return r_obs, p


def itol_node_id(name):
    """Leaf name as bin/fastani_to_tree.py writes it: [^A-Za-z0-9_.-] -> '_'."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def newick_leaf_names(path):
    """Leaf labels of a Newick file (names right after '(' or ',').

    Internal node labels (after ')') and branch lengths are skipped; quoted
    labels are unquoted.
    """
    with open(path, encoding="utf-8") as fh:
        text = re.sub(r"\[[^\]]*\]", "", fh.read())  # drop [comments]
    names = re.findall(r"[(,]\s*('[^']*'|[^\s(),:;]+)", text)
    return [n[1:-1].replace("''", "'") if n.startswith("'") else n for n in names]


def match_tree_leaves(genomes, leaves=None):
    """Map each genome to its iTOL node ID and check the IDs against the tree.

    Fails if two genomes get the same node ID. With `leaves` (the tree's leaf
    labels), both kinds of mismatch are reported, and no match at all is an
    error. Returns (node ID per genome, unmatched genomes, leaves without a
    cluster); both lists are empty without `leaves`.
    """
    node_of = {g: itol_node_id(g) for g in genomes}
    clashes = pd.Series(node_of).loc[lambda s: s.duplicated(keep=False)]
    if not clashes.empty:
        eprint(f"ERROR: genomes collide on the same iTOL node ID: "
               f"{sorted(clashes.index)[:10]}")
        sys.exit(1)

    unmatched, uncoloured = [], []
    if leaves is not None:
        leaf_set = set(leaves)
        unmatched = sorted(g for g in genomes if node_of[g] not in leaf_set)
        uncoloured = sorted(leaf_set - set(node_of.values()))
        if len(unmatched) == len(genomes):
            eprint(f"ERROR: no genome of the clustering is a leaf of the tree "
                   f"(e.g. genome '{genomes[0]}' -> node '{node_of[genomes[0]]}'; "
                   f"tree leaves e.g. {sorted(leaf_set)[:3]}). iTOL would colour "
                   f"nothing.")
            sys.exit(1)
        if unmatched:
            eprint(f"WARNING: {len(unmatched)} genome(s) are not leaves of the tree "
                   f"and are left out of the iTOL file: {unmatched[:10]}"
                   + (" ..." if len(unmatched) > 10 else ""))
        if uncoloured:
            eprint(f"WARNING: {len(uncoloured)} tree leaves have no cluster and "
                   f"stay uncoloured: {uncoloured[:10]}"
                   + (" ..." if len(uncoloured) > 10 else ""))
    return node_of, unmatched, uncoloured


def write_itol_colors(path, genomes, cluster_of, sizes, singletons, color_of,
                      node_of, skip):
    """Write an iTOL TREE_COLORS file: one coloured range per genome.

    Colours are those of the ordination plot; singletons get no line, so iTOL
    leaves them uncoloured. Genomes in `skip` (not leaves of the tree) are left
    out. Returns the number of coloured genomes.
    """
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("TREE_COLORS\n"
                 "#Cluster colours of the ordination plot, one coloured range per "
                 "genome (singletons left uncoloured).\n"
                 "#Drop this file onto the tree in iTOL.\n"
                 "SEPARATOR TAB\n"
                 "DATA\n"
                 "#NODE_ID\tTYPE\tCOLOR\tLABEL\n")
        n = 0
        for g, c in zip(genomes, cluster_of):
            if g in skip or c in singletons:
                continue
            fh.write(f"{node_of[g]}\trange\t{to_hex(color_of[c])}\t{c} (n={sizes[c]})\n")
            n += 1
    return n


MISSING_CATEGORY = "NA"
GREY_CATEGORIES = {"na", "unknown"}
GREY = "#d9d9d9"


def read_metadata(path):
    """Read the genome metadata table.

    :param path: Table with a header row; column 1 = genome ID, the other columns
        are category variables. The delimiter is auto-detected.
    :return: One row per genome (index = genome ID), every value a string; empty
        cells are ``NA``.
    """
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    if df.shape[1] < 2:
        eprint(f"ERROR: metadata table '{path}' needs a genome column and at least "
               f"one category column.")
        sys.exit(1)
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str).str.strip()
    duplicated = sorted(set(df.index[df.index.duplicated()]))
    if duplicated:
        eprint(f"ERROR: metadata table '{path}' lists genomes more than once: {duplicated[:10]}")
        sys.exit(1)
    df = df.apply(lambda col: col.str.strip()).replace("", np.nan)
    return df.fillna(MISSING_CATEGORY)


def category_colours(values):
    """Give each category a colour; ``NA`` and ``unknown`` are light grey.

    :param values: The category of each genome.
    :return: Hex colour per category, in sorted category order.
    """
    categories = sorted(set(values))
    coloured = [c for c in categories if c.lower() not in GREY_CATEGORIES]
    if len(coloured) > 20:
        eprint(f"WARNING: {len(coloured)} categories; colours repeat after 20.")
    cmap = plt.get_cmap("tab10" if len(coloured) <= 10 else "tab20")
    colours = {c: to_hex(cmap(i % cmap.N)) for i, c in enumerate(coloured)}
    colours.update({c: GREY for c in categories if c.lower() in GREY_CATEGORIES})
    return {c: colours[c] for c in categories}


def plot_metadata(path_base, column, ax1, ax2, axis_labels, categories, colour_of,
                  cluster_of, sizes, singletons, conf):
    """Draw the ordination with the genomes coloured by one metadata column.

    Non-singleton genomes are filled with their category colour; singleton genomes
    are hollow with a rim in that colour. Each cluster of at least 3 genomes gets its
    confidence ellipse as a grey outline, with the cluster name at its centre.

    :param path_base: Output path without the extension (``.png`` and ``.pdf`` are written).
    :param column: The metadata column name, for the title and the legend.
    :param axis_labels: The x and y axis labels of the cluster plot.
    """
    fig, axp = plt.subplots(figsize=(8.5, 6.5))
    sing_mask = np.array([c in singletons for c in cluster_of])
    for category, colour in colour_of.items():
        m = categories == category
        n = int(m.sum())
        axp.scatter(ax1[m & ~sing_mask], ax2[m & ~sing_mask], color=colour, s=55,
                    edgecolor="black", linewidth=0.4, alpha=0.9,
                    label=f"{category} (n={n})")
        axp.scatter(ax1[m & sing_mask], ax2[m & sing_mask], s=55, facecolors="none",
                    edgecolors=colour, linewidth=1.2, alpha=0.9)
    if sing_mask.any():
        axp.scatter([], [], s=55, facecolors="none", edgecolors="#404040",
                    linewidth=1.2, label="hollow = singleton cluster")
    for cluster in sorted(c for c in sizes if c not in singletons):
        m = cluster_of == cluster
        if confidence_ellipse(ax1[m], ax2[m], axp, conf=conf, facecolor="none",
                              edgecolor="#707070", lw=1.0, ls="--") is not None:
            axp.text(np.mean(ax1[m]), np.mean(ax2[m]), cluster, ha="center",
                     va="center", fontsize=8, color="#404040")
    axp.set_xlabel(axis_labels[0])
    axp.set_ylabel(axis_labels[1])
    axp.set_title(f"PCoA of gene presence/absence, coloured by {column}")
    axp.axhline(0, color="grey", lw=0.5, ls="--")
    axp.axvline(0, color="grey", lw=0.5, ls="--")
    axp.legend(title=column, bbox_to_anchor=(1.02, 1), loc="upper left",
               fontsize=8, frameon=False)
    fig.text(0.01, 0.01, f"Dashed ellipses: clusters, {conf * 100:.0f}% confidence",
             fontsize=7, color="grey")
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(f"{path_base}.png", dpi=300)
    fig.savefig(f"{path_base}.pdf")
    plt.close(fig)


def write_itol_strip(path, column, genomes, categories, colour_of, node_of, skip):
    """Write an iTOL DATASET_COLORSTRIP file: one outer ring for one metadata column.

    :param column: The metadata column name (dataset label and legend title).
    :param categories: The category of each genome, in the order of `genomes`.
    :param colour_of: Hex colour per category.
    :param node_of: iTOL node ID per genome.
    :param skip: Genomes to leave out (not leaves of the tree).
    :return: The number of genomes written.
    """
    rows = [(node_of[g], c) for g, c in zip(genomes, categories) if g not in skip]
    legend = sorted({c for _, c in rows})
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("DATASET_COLORSTRIP\n"
                 "SEPARATOR TAB\n"
                 f"DATASET_LABEL\t{column}\n"
                 "COLOR\t#000000\n"
                 f"LEGEND_TITLE\t{column}\n"
                 "LEGEND_SHAPES\t" + "\t".join("1" for _ in legend) + "\n"
                 "LEGEND_COLORS\t" + "\t".join(colour_of[c] for c in legend) + "\n"
                 "LEGEND_LABELS\t" + "\t".join(legend) + "\n"
                 "STRIP_WIDTH\t25\n"
                 "MARGIN\t5\n"
                 "BORDER_WIDTH\t0\n"
                 "#Drop this file onto the tree in iTOL: one outer ring.\n"
                 "DATA\n"
                 "#NODE_ID\tCOLOR\tLABEL\n")
        for node, category in rows:
            fh.write(f"{node}\t{colour_of[category]}\t{category}\n")
    return len(rows)


def read_annotations(path, has_header):
    """Read the annotation table.

    :param path: Table with the ID in column 1 and the annotation in column 2.
    :param has_header: Whether the first row is a header.
    :return: The distinct non-empty annotations of each ID, in file order. Empty
        annotations (e.g. Panaroo genes without a product) are left out.
    """
    header = 0 if has_header else None
    df = pd.read_csv(path, sep=None, engine="python", header=header, dtype=str)
    if df.shape[1] < 2:
        eprint(f"ERROR: annotation table '{path}' needs >= 2 columns.")
        sys.exit(1)
    annotations = {}
    for key, value in zip(df.iloc[:, 0].astype(str).str.strip(), df.iloc[:, 1].fillna("").astype(str).str.strip()):
        if value and value not in annotations.setdefault(key, []):
            annotations[key].append(value)
    return {key: values for key, values in annotations.items() if values}


def family_annotations(gene_families, protein_annotations):
    """Give each gene family the annotations of its proteins.

    :param gene_families: Panaroo's gene_presence_absence.csv: column 1 = gene family
        (the Rtab gene name), columns 4 onwards = one column per genome, with the
        family's protein IDs in that genome (several separated by ';').
    :param protein_annotations: The annotations of each protein ID.
    :return: The annotation string of each annotated family, with its distinct
        annotations as 'annotation (proteins with it/proteins in the family)', most
        frequent first; and the number of annotated proteins found in the table.
    """
    family_map, found = {}, 0
    with open(gene_families, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader, [])
        if not header or header[0] != "Gene":
            eprint(f"ERROR: '{gene_families}' is not a Panaroo gene_presence_absence.csv "
                   f"(its first column must be 'Gene').")
            sys.exit(1)
        for row in reader:
            proteins = [p.strip() for cell in row[3:] if cell for p in cell.split(";") if p.strip()]
            counts = Counter()
            for protein in proteins:
                hits = protein_annotations.get(protein)
                if hits:
                    found += 1
                    counts.update(hits)
            if counts:
                family_map[row[0]] = "; ".join(
                    f"{label} ({n}/{len(proteins)})"
                    for label, n in sorted(counts.items(), key=lambda item: (-item[1], item[0])))
    return family_map, found


def fisher_enrichment(mat, anno_map, g_idx, r_idx, group_label, ref_label,
                      comparison_label):
    """Per-gene Fisher's exact test, group vs reference (index positions)."""
    M = mat.values
    genes = mat.index.to_numpy()
    n_g, n_r = len(g_idx), len(r_idx)
    pres_g = M[:, g_idx].sum(axis=1)
    pres_r = M[:, r_idx].sum(axis=1)
    total = pres_g + pres_r
    variable = (total > 0) & (total < (n_g + n_r))  # skip invariant genes

    records = []
    for gi in np.where(variable)[0]:
        a = int(pres_g[gi])
        b = n_g - a
        c = int(pres_r[gi])
        d = n_r - c
        orr, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        prop_g = a / n_g if n_g else 0.0
        prop_r = c / n_r if n_r else 0.0
        direction = ("enriched" if prop_g > prop_r
                     else "depleted" if prop_g < prop_r else "none")
        records.append({
            "Gene": genes[gi],
            "Annotation": anno_map.get(str(genes[gi]), ""),
            "Comparison": comparison_label,
            "Group": group_label,
            "Reference": ref_label,
            "group_size": n_g,
            "reference_size": n_r,
            "present_in_group": a,
            "prop_group": round(prop_g, 4),
            "present_in_reference": c,
            "prop_reference": round(prop_r, 4),
            "odds_ratio": orr,
            "p_value": p,
            "direction": direction,
        })
    enr = pd.DataFrame.from_records(records)
    if not enr.empty:
        enr["FDR"] = bh_fdr(enr["p_value"].values)
        enr = enr.sort_values(["FDR", "p_value"]).reset_index(drop=True)
        cols = ["Gene", "Annotation", "Comparison", "Group", "Reference",
                "group_size", "reference_size", "present_in_group",
                "prop_group", "present_in_reference", "prop_reference",
                "odds_ratio", "p_value", "FDR", "direction"]
        enr = enr[cols]
    return enr


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description="PCoA ordination (+90% ellipses), pairwise ANOSIM and "
                    "gene enrichment from a pangenome Rtab.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--annotations", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--prefix", default="pangenome")
    ap.add_argument("--metric", default="jaccard",
                    choices=["jaccard", "braycurtis"])
    ap.add_argument("--enrich-cluster", default=None,
                    help="Focal cluster for one-vs-rest enrichment. If omitted, "
                         "the two largest clusters are compared.")
    ap.add_argument("--min-cluster-size", type=int, default=5,
                    help="Minimum size for a group to have statistical power.")
    ap.add_argument("--ellipse-conf", type=float, default=0.90,
                    help="Confidence level for the ordination ellipses.")
    ap.add_argument("--anosim-min-size", type=int, default=3,
                    help="Minimum group size to include in pairwise ANOSIM.")
    ap.add_argument("--anosim-perm", type=int, default=999,
                    help="Number of ANOSIM permutations.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--anno-no-header", action="store_true")
    ap.add_argument("--gene-families", default=None,
                    help="Panaroo gene_presence_absence.csv. With it, the --annotations "
                         "IDs are protein IDs (as in the .faa files), mapped to their "
                         "gene families.")
    ap.add_argument("--tree", default=None,
                    help="Newick tree the iTOL colour file is for; its leaf "
                         "labels are checked against the genome names.")
    ap.add_argument("--metadata", default=None,
                    help="Genome metadata table: genome ID, then one column per "
                         "category variable. Adds one ordination plot and one "
                         "iTOL ring file per column.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---- Load presence/absence matrix (genes x genomes) ------------------- #
    mat = pd.read_csv(args.rtab, sep="\t", index_col=0)
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0)
    mat = (mat > 0).astype(int)

    # ---- Load clustering table -------------------------------------------- #
    clust = pd.read_csv(args.clusters)
    if clust.shape[1] < 2:
        eprint("ERROR: clustering table needs >= 2 columns (Taxon, Cluster).")
        sys.exit(1)
    tax_col, clu_col = clust.columns[:2]
    clust = clust[[tax_col, clu_col]].copy()
    clust.columns = ["Taxon", "Cluster"]
    clust["Taxon"] = clust["Taxon"].astype(str)
    clust["Cluster"] = clust["Cluster"].astype(str)

    # ---- Align genomes ---------------------------------------------------- #
    genomes_clu = set(clust["Taxon"])
    common = [g for g in mat.columns if g in genomes_clu]
    dropped_mat = sorted(set(mat.columns) - genomes_clu)
    dropped_clu = sorted(genomes_clu - set(mat.columns))
    if dropped_mat:
        eprint(f"WARNING: {len(dropped_mat)} matrix genome(s) lack a cluster and "
               f"are dropped: {dropped_mat[:10]}"
               + (" ..." if len(dropped_mat) > 10 else ""))
    if dropped_clu:
        eprint(f"WARNING: {len(dropped_clu)} clustered genome(s) not in matrix "
               f"are dropped: {dropped_clu[:10]}"
               + (" ..." if len(dropped_clu) > 10 else ""))
    if len(common) < 3:
        eprint("ERROR: fewer than 3 shared genomes; cannot ordinate.")
        sys.exit(1)

    mat = mat[common]
    clu_map = dict(zip(clust["Taxon"], clust["Cluster"]))
    cluster_of = np.array([clu_map[g] for g in common])
    sizes = pd.Series(cluster_of).value_counts().to_dict()
    singletons = {c for c, n in sizes.items() if n == 1}

    annotations = read_annotations(args.annotations, not args.anno_no_header)
    if not annotations:
        eprint(f"NOTE: the annotation table '{args.annotations}' has no non-empty "
               f"annotation; the enrichment table will have empty annotations.")
    if args.gene_families:
        anno_map, n_found = family_annotations(args.gene_families, annotations)
        if annotations and n_found == 0:
            eprint(f"WARNING: no protein ID of the annotation table is in "
                   f"'{args.gene_families}'. Use the protein IDs of the ggCaller .faa files.")
        elif annotations:
            eprint(f"Annotations: {n_found} of {len(annotations)} annotated proteins found, "
                   f"in {len(anno_map)} gene families.")
    else:
        anno_map = {key: "; ".join(values) for key, values in annotations.items()}
        if annotations and not any(gene in anno_map for gene in mat.index):
            eprint("WARNING: no ID of the annotation table is an Rtab gene family. If the "
                   "table uses protein IDs from the .faa files, add --gene-families "
                   "gene_presence_absence.csv.")

    # ---- Metadata (optional) ---------------------------------------------- #
    metadata = None
    if args.metadata:
        metadata = read_metadata(args.metadata)
        unlisted = [g for g in common if g not in metadata.index]
        n_extra = len(set(metadata.index) - set(common))
        metadata = metadata.reindex(common).fillna(MISSING_CATEGORY)
        if unlisted:
            eprint(f"WARNING: {len(unlisted)} genome(s) are not in the metadata table "
                   f"and get '{MISSING_CATEGORY}': {unlisted[:10]}"
                   + (" ..." if len(unlisted) > 10 else ""))
        if n_extra:
            eprint(f"WARNING: {n_extra} metadata row(s) are not analysed genomes and "
                   f"are ignored.")
        safe = {col: itol_node_id(col) for col in metadata.columns}
        if len(set(safe.values())) < len(safe):
            eprint(f"ERROR: metadata column names collide in file names: {sorted(safe)}")
            sys.exit(1)

    # ---- Ordination (PCoA) ------------------------------------------------ #
    X = mat.T.values.astype(int)  # genomes x genes
    dist = squareform(pdist(X, metric=args.metric))
    coords, prop = pcoa(dist)
    if coords.shape[1] < 2:
        eprint("ERROR: ordination produced < 2 axes (too little variation).")
        sys.exit(1)
    ax1, ax2 = coords[:, 0], coords[:, 1]
    pc1_pct, pc2_pct = prop[0] * 100, prop[1] * 100

    coords_df = pd.DataFrame({
        "Taxon": common, "PCo1": ax1, "PCo2": ax2, "Cluster": cluster_of,
        "Cluster_size": [sizes[c] for c in cluster_of],
        "Is_singleton": [c in singletons for c in cluster_of],
    })
    if metadata is not None:
        for col in metadata.columns:
            coords_df[col] = metadata[col].to_numpy()
    coords_path = os.path.join(args.outdir, f"{args.prefix}_ordination_coords.csv")
    coords_df.to_csv(coords_path, index=False)

    # ---- Plot ------------------------------------------------------------- #
    non_single = sorted([c for c in sizes if c not in singletons],
                        key=lambda x: (len(x), x))
    cmap = plt.get_cmap("tab20" if len(non_single) > 10 else "tab10")
    color_of = {c: cmap(i % cmap.N) for i, c in enumerate(non_single)}

    fig, axp = plt.subplots(figsize=(8.5, 6.5))
    sing_mask = np.array([c in singletons for c in cluster_of])
    if sing_mask.any():
        axp.scatter(ax1[sing_mask], ax2[sing_mask], s=55, facecolors="none",
                    edgecolors="#404040", linewidth=1.0, alpha=0.9,
                    label="Singleton")
    for c in non_single:
        m = cluster_of == c
        axp.scatter(ax1[m], ax2[m], color=color_of[c], s=55,
                    edgecolor="black", linewidth=0.4, alpha=0.9,
                    label=f"{c} (n={sizes[c]})")
        confidence_ellipse(ax1[m], ax2[m], axp, conf=args.ellipse_conf,
                           facecolor=color_of[c], edgecolor=color_of[c],
                           alpha=0.12, lw=1.2)
    axp.set_xlabel(f"PCo1 ({pc1_pct:.1f}%)")
    axp.set_ylabel(f"PCo2 ({pc2_pct:.1f}%)")
    axp.set_title(f"PCoA of gene presence/absence ({args.metric})")
    axp.axhline(0, color="grey", lw=0.5, ls="--")
    axp.axvline(0, color="grey", lw=0.5, ls="--")
    axp.legend(title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left",
               fontsize=8, frameon=False)
    fig.text(0.01, 0.01, f"Ellipses: {args.ellipse_conf*100:.0f}% confidence",
             fontsize=7, color="grey")
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    png_path = os.path.join(args.outdir, f"{args.prefix}_ordination.png")
    pdf_path = os.path.join(args.outdir, f"{args.prefix}_ordination.pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    # ---- iTOL colour ranges (same colours as the plot) -------------------- #
    itol_path = os.path.join(args.outdir, f"{args.prefix}_itol_colors.txt")
    leaves = newick_leaf_names(args.tree) if args.tree else None
    node_of, itol_unmatched, _ = match_tree_leaves(common, leaves)
    n_itol = write_itol_colors(itol_path, common, cluster_of, sizes, singletons,
                               color_of, node_of, set(itol_unmatched))

    # ---- Metadata plots and iTOL rings ------------------------------------ #
    metadata_files = []
    if metadata is not None:
        axis_labels = (f"PCo1 ({pc1_pct:.1f}%)", f"PCo2 ({pc2_pct:.1f}%)")
        for col in metadata.columns:
            categories = metadata[col].to_numpy()
            colour_of = category_colours(categories)
            base = os.path.join(args.outdir, f"{args.prefix}_ordination_{safe[col]}")
            plot_metadata(base, col, ax1, ax2, axis_labels, categories, colour_of,
                          cluster_of, sizes, singletons, args.ellipse_conf)
            strip_path = os.path.join(args.outdir, f"{args.prefix}_itol_{safe[col]}.txt")
            write_itol_strip(strip_path, col, common, categories, colour_of, node_of,
                             set(itol_unmatched))
            metadata_files.append((col, f"{base}.png", strip_path))

    # ---- Pairwise ANOSIM -------------------------------------------------- #
    anosim_path = os.path.join(args.outdir, f"{args.prefix}_anosim.csv")
    anosim_groups = sorted([c for c, n in sizes.items()
                            if n >= args.anosim_min_size and c not in singletons],
                           key=lambda x: (len(x), x))
    idx_of = {g: i for i, g in enumerate(common)}
    group_idx = {c: [idx_of[g] for g in common if cluster_of[idx_of[g]] == c]
                 for c in anosim_groups}
    anosim_rows = []
    for i in range(len(anosim_groups)):
        for k in range(i + 1, len(anosim_groups)):
            g1, g2 = anosim_groups[i], anosim_groups[k]
            sub = group_idx[g1] + group_idx[g2]
            labels = np.array([g1] * len(group_idx[g1]) + [g2] * len(group_idx[g2]))
            d_sub = dist[np.ix_(sub, sub)]
            r, p = anosim(d_sub, labels, n_perm=args.anosim_perm, rng=rng)
            anosim_rows.append({"Group1": g1, "Group2": g2,
                                "n1": len(group_idx[g1]), "n2": len(group_idx[g2]),
                                "ANOSIM_R": round(r, 4), "p_value": p})
    if anosim_rows:
        an = pd.DataFrame(anosim_rows)
        an["FDR"] = bh_fdr(an["p_value"].values)
        an = an.sort_values(["p_value", "ANOSIM_R"],
                            ascending=[True, False]).reset_index(drop=True)
        an.to_csv(anosim_path, index=False)
        n_anosim = len(an)
    else:
        with open(anosim_path, "w") as fh:
            fh.write(f"# Pairwise ANOSIM skipped: fewer than 2 groups with "
                     f">= {args.anosim_min_size} genomes.\n")
        n_anosim = 0

    # ---- Enrichment ------------------------------------------------------- #
    enrich_path = os.path.join(args.outdir, f"{args.prefix}_enrichment.csv")
    enr = None
    enrich_note = None

    def idx_for(cluster_id):
        return [i for i in range(len(common)) if cluster_of[i] == cluster_id]

    if args.enrich_cluster is not None:
        focal = str(args.enrich_cluster)
        if focal not in sizes:
            enrich_note = (f"Enrichment skipped: cluster '{focal}' not found. "
                           f"Available: {sorted(sizes)}")
        elif sizes[focal] < args.min_cluster_size:
            enrich_note = f"no statistical power in group {focal}"
        else:
            g_idx = idx_for(focal)
            r_idx = [i for i in range(len(common)) if cluster_of[i] != focal]
            enr = fisher_enrichment(mat, anno_map, g_idx, r_idx,
                                    group_label=focal, ref_label="rest",
                                    comparison_label=f"{focal}_vs_rest")
    else:
        # default: compare the two largest (non-singleton) clusters
        ranked = sorted([(n, c) for c, n in sizes.items() if c not in singletons],
                        key=lambda t: (-t[0], t[1]))
        if len(ranked) < 2:
            enrich_note = ("Enrichment skipped: fewer than 2 non-singleton "
                           "clusters to compare.")
        else:
            g1 = ranked[0][1]
            g2 = ranked[1][1]
            g_idx = idx_for(g1)
            r_idx = idx_for(g2)
            enr = fisher_enrichment(mat, anno_map, g_idx, r_idx,
                                    group_label=g1, ref_label=g2,
                                    comparison_label=f"{g1}_vs_{g2}")
            if sizes[g1] < args.min_cluster_size or sizes[g2] < args.min_cluster_size:
                eprint(f"WARNING: a compared group is smaller than "
                       f"{args.min_cluster_size}; low statistical power.")

    if enr is not None:
        enr.to_csv(enrich_path, index=False)
        n_tests = len(enr)
    else:
        eprint(enrich_note)
        with open(enrich_path, "w") as fh:
            fh.write("# " + enrich_note + "\n")
        n_tests = 0

    # ---- Summary ---------------------------------------------------------- #
    print("=" * 62)
    print(f"Genomes analysed       : {len(common)}")
    print(f"Clusters               : {len(sizes)} "
          f"({len(singletons)} singleton[s])")
    print(f"Ordination plot        : {png_path}")
    print(f"Coordinates table      : {coords_path}")
    print(f"ANOSIM groups (>= {args.anosim_min_size:>2}) : {anosim_groups}")
    print(f"Pairwise ANOSIM table  : {anosim_path} ({n_anosim} pair[s])")
    print(f"iTOL colour file       : {itol_path} ({n_itol} genome[s] coloured"
          + (f", {len(itol_unmatched)} not in the tree)" if args.tree
             else ", names not checked: no --tree)"))
    if enr is not None:
        print(f"Enrichment table       : {enrich_path} "
              f"({enr['Comparison'].iloc[0]}, {n_tests} tests)")
    else:
        print(f"Enrichment             : {enrich_note}")
    for col, png, strip in metadata_files:
        print(f"Metadata '{col}'".ljust(23) + f": {png}, {strip}")
    print("=" * 62)


if __name__ == "__main__":
    main()
