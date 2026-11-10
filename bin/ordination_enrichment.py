#!/usr/bin/env python3
"""
ordination_enrichment.py

Post-pipeline analysis of a gene presence/absence pangenome matrix.

Given a gene_presence_absence.Rtab (from ggCaller or Panaroo), a clustering
table assigning each genome to a cluster, and a per-gene annotation table,
this script produces:

  1) An ordination plot (PCoA on Jaccard distance of the presence/absence
     matrix) with genomes coloured by their cluster. Singleton clusters
     (size 1) are coloured grey. Each non-singleton cluster with >= 3 genomes
     is drawn with a 90% confidence ellipse.
  2) A pairwise ANOSIM table (Analysis of Similarities) between groups,
     with the R statistic, permutation p-value and BH-adjusted FDR.
  3) A gene-enrichment table based on Fisher's exact tests, with odds ratio,
     p-value, BH-adjusted FDR and direction (enriched / depleted). See the
     enrichment logic below.

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
--annotations   Table with one annotation string per gene. First column = gene
                id (matching Rtab gene names), second column = annotation string.
                Delimiter auto-detected; header expected (use --anno-no-header
                if there is none).

OUTPUTS (written to --outdir, prefixed with --prefix)
-------
  <prefix>_ordination.png / .pdf   ordination plot with 90% ellipses
  <prefix>_ordination_coords.csv   per-genome PCoA coordinates + cluster
  <prefix>_anosim.csv              pairwise ANOSIM table
  <prefix>_enrichment.csv          gene enrichment table (or a note if skipped)

DEPENDENCIES
------------
python3, numpy, pandas, scipy, matplotlib

EXAMPLES
--------
  # default: ordination + pairwise ANOSIM + two-largest-group enrichment
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations anno.tsv --outdir results --prefix subspecies

  # one-vs-rest enrichment for a chosen cluster
  python ordination_enrichment.py --rtab gpa.Rtab --clusters clusters.csv \
      --annotations anno.tsv --enrich-cluster 2 --outdir results
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, chi2, rankdata

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def read_annotations(path, has_header):
    header = 0 if has_header else None
    df = pd.read_csv(path, sep=None, engine="python", header=header, dtype=str)
    if df.shape[1] < 2:
        eprint(f"ERROR: annotation table '{path}' needs >= 2 columns.")
        sys.exit(1)
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))


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
        a = int(pres_g[gi]); b = n_g - a
        c = int(pres_r[gi]); d = n_r - c
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

    anno_map = read_annotations(args.annotations, not args.anno_no_header)

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
        axp.scatter(ax1[sing_mask], ax2[sing_mask], c="grey", s=55,
                    edgecolor="black", linewidth=0.4, alpha=0.8,
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
    if enr is not None:
        print(f"Enrichment table       : {enrich_path} "
              f"({enr['Comparison'].iloc[0]}, {n_tests} tests)")
    else:
        print(f"Enrichment             : {enrich_note}")
    print("=" * 62)


if __name__ == "__main__":
    main()
