#!/usr/bin/env python3

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def normalise_genome_id(x: str) -> str:
    """
    Normalise genome identifiers so that FastANI paths and PopPUNK sample names can match.

    This removes directory paths and common FASTA extensions.
    Adjust this function if your PopPUNK sample IDs intentionally include extensions.
    """
    x = str(x)
    name = Path(x).name

    for suffix in [
        ".fasta.gz", ".fa.gz", ".fna.gz",
        ".fasta", ".fa", ".fna",
        ".fas", ".contigs"
    ]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return name


def read_fastani(path: str) -> pd.DataFrame:
    """
    Read all-vs-all FastANI output.

    Expected columns:
      query, reference, ani, fragments_mapped, total_fragments

    FastANI ANI may be reported as 99.3 or 0.993.
    This function converts ANI to 0-1 scale.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["query", "reference", "ani", "fragments_mapped", "total_fragments"],
    )

    df["query"] = df["query"].map(normalise_genome_id)
    df["reference"] = df["reference"].map(normalise_genome_id)

    df["ani"] = pd.to_numeric(df["ani"], errors="coerce")
    df["fragments_mapped"] = pd.to_numeric(df["fragments_mapped"], errors="coerce")
    df["total_fragments"] = pd.to_numeric(df["total_fragments"], errors="coerce")

    # Convert 95-100 scale to 0-1 scale if needed.
    if df["ani"].dropna().median() > 1:
        df["ani"] = df["ani"] / 100.0

    df["alignment_fraction"] = df["fragments_mapped"] / df["total_fragments"]

    return df


def read_clusters(path: str) -> pd.DataFrame:
    """
    Read PopPUNK cluster assignment file.

    Expected common PopPUNK columns:
      Taxon, Cluster

    Also accepts:
      sample, cluster
      id, cluster
      genome_id, cluster
    """
    df = pd.read_csv(path)

    lower_to_original = {c.lower(): c for c in df.columns}

    sample_col = None
    for candidate in ["taxon", "sample", "id", "genome_id", "name"]:
        if candidate in lower_to_original:
            sample_col = lower_to_original[candidate]
            break

    cluster_col = None
    for candidate in ["cluster", "vlkc", "lineage"]:
        if candidate in lower_to_original:
            cluster_col = lower_to_original[candidate]
            break

    if sample_col is None or cluster_col is None:
        raise ValueError(
            f"Could not identify sample and cluster columns in {path}. "
            f"Found columns: {list(df.columns)}"
        )

    out = df[[sample_col, cluster_col]].copy()
    out.columns = ["genome_id", "cluster_id"]
    out["genome_id"] = out["genome_id"].map(normalise_genome_id)
    out["cluster_id"] = out["cluster_id"].astype(str)

    return out


def read_labels(path: str) -> pd.DataFrame:
    """
    Read the per-genome HQ/MQ labels file produced by
    `bin/spp_eligibility_from_qc.py` (`genome,label`).

    The genome column is matched among genome/genome_id/sample/id/name and the
    label column among label/quality_status/quality/qc_status. DISCARDED genomes
    are already excluded upstream, so only HQ/MQ genomes are expected here.
    """
    df = pd.read_csv(path, sep=None, engine="python")

    lower_to_original = {c.lower(): c for c in df.columns}

    genome_col = None
    for candidate in ["genome", "genome_id", "sample", "id", "name"]:
        if candidate in lower_to_original:
            genome_col = lower_to_original[candidate]
            break

    if genome_col is None:
        raise ValueError(
            f"Labels file must contain a genome/genome_id column. Found columns: {list(df.columns)}"
        )

    quality_col = None
    for candidate in ["label", "quality_status", "quality", "qc_status"]:
        if candidate in lower_to_original:
            quality_col = lower_to_original[candidate]
            break

    if quality_col is None:
        raise ValueError(
            f"Labels file must contain label, quality_status, quality, or qc_status. "
            f"Found columns: {list(df.columns)}"
        )

    out = df[[genome_col, quality_col]].copy()
    out.columns = ["genome_id", "quality_status"]
    out["genome_id"] = out["genome_id"].map(normalise_genome_id)
    out["quality_status"] = out["quality_status"].astype(str).str.upper()

    return out


def make_symmetric_ani_matrix(fastani: pd.DataFrame, genomes: list[str]) -> pd.DataFrame:
    """
    Construct a symmetric ANI matrix.

    FastANI all-vs-all may contain both A->B and B->A. If both exist, use their mean.
    Diagonal is set to 1.0.
    Missing values remain NaN.
    """
    pair_df = fastani.copy()

    pair_df["g1"] = pair_df[["query", "reference"]].min(axis=1)
    pair_df["g2"] = pair_df[["query", "reference"]].max(axis=1)

    pair_mean = (
        pair_df
        .groupby(["g1", "g2"], as_index=False)["ani"]
        .mean()
    )

    matrix = pd.DataFrame(np.nan, index=genomes, columns=genomes, dtype=float)

    for row in pair_mean.itertuples(index=False):
        g1 = row.g1
        g2 = row.g2
        ani = row.ani

        if g1 in matrix.index and g2 in matrix.columns:
            matrix.loc[g1, g2] = ani
            matrix.loc[g2, g1] = ani

    np.fill_diagonal(matrix.values, 1.0)

    return matrix


def safe_percentile(values, q):
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if len(values) == 0:
        return np.nan

    return float(np.percentile(values, q))


def score_gap(tail_gap, median_gap):
    if np.isnan(tail_gap) or np.isnan(median_gap):
        return np.nan
    if tail_gap >= 0.002:
        return 1.0
    if tail_gap >= 0 and median_gap >= 0.002:
        return 0.7
    if tail_gap < 0 and median_gap > 0:
        return 0.4
    return 0.0


def score_overlap(overlap_fraction):
    if np.isnan(overlap_fraction):
        return np.nan
    if overlap_fraction <= 0.05:
        return 1.0
    if overlap_fraction <= 0.10:
        return 0.7
    if overlap_fraction <= 0.25:
        return 0.4
    return 0.0


def score_cohesion(p5_intra_ani):
    if np.isnan(p5_intra_ani):
        return np.nan
    if p5_intra_ani >= 0.980:
        return 1.0
    if p5_intra_ani >= 0.975:
        return 0.7
    if p5_intra_ani >= 0.970:
        return 0.4
    return 0.0


def score_support(cluster_size, n_hq):
    if cluster_size >= 5 and n_hq >= 3:
        return 1.0
    if cluster_size >= 3 and n_hq >= 2:
        return 0.6
    if cluster_size >= 2:
        return 0.3
    return np.nan


def validity_status(score, cluster_size):
    if cluster_size == 1:
        return "Singleton"
    if np.isnan(score):
        return "Unresolved"
    if score >= 0.80:
        return "Accepted"
    if score >= 0.55:
        return "Borderline"
    return "Defective"


def confidence_status(cluster_size, n_hq, hq_ratio):
    if n_hq >= 3 and hq_ratio >= 0.40 and cluster_size >= 5:
        return "High_confidence"
    if n_hq >= 2 and cluster_size >= 3:
        return "Medium_confidence"
    if cluster_size >= 2 and n_hq >= 2:
        return "Low_confidence"
    if n_hq < 2:
        return "Unresolved_low_HQ"
    return "Low_confidence"


def evaluate_clusters(clusters: pd.DataFrame, metadata: pd.DataFrame, ani_matrix: pd.DataFrame) -> pd.DataFrame:
    merged = clusters.merge(metadata, on="genome_id", how="left")
    merged["quality_status"] = merged["quality_status"].fillna("UNKNOWN")

    cluster_to_genomes = (
        merged
        .groupby("cluster_id")["genome_id"]
        .apply(list)
        .to_dict()
    )

    rows = []

    for cluster_id, members in cluster_to_genomes.items():
        cluster_size = len(members)

        n_hq = int(
            merged.loc[
                merged["cluster_id"].eq(cluster_id)
                & merged["quality_status"].eq("HQ"),
                "genome_id"
            ].nunique()
        )

        hq_ratio = n_hq / cluster_size if cluster_size > 0 else np.nan

        if cluster_size == 1:
            rows.append({
                "cluster_id": cluster_id,
                "cluster_size": cluster_size,
                "N_HQ": n_hq,
                "HQ_ratio_cluster": hq_ratio,
                "nearest_external_cluster": np.nan,
                "p5_intra_ANI": np.nan,
                "median_intra_ANI": np.nan,
                "p95_nearest_inter_ANI": np.nan,
                "median_nearest_inter_ANI": np.nan,
                "tail_gap": np.nan,
                "median_gap": np.nan,
                "overlap_fraction": np.nan,
                "gap_score": np.nan,
                "overlap_score": np.nan,
                "cohesion_score": np.nan,
                "support_score": np.nan,
                "cluster_structure_score": np.nan,
                "validity_status": "Singleton",
                "confidence_status": confidence_status(cluster_size, n_hq, hq_ratio),
            })
            continue

        # Intra-cluster ANI values.
        intra_values = []
        for g1, g2 in itertools.combinations(members, 2):
            if g1 in ani_matrix.index and g2 in ani_matrix.columns:
                intra_values.append(ani_matrix.loc[g1, g2])

        p5_intra = safe_percentile(intra_values, 5)
        median_intra = safe_percentile(intra_values, 50)

        # Find nearest external cluster by highest median inter-cluster ANI.
        best_external_cluster = None
        best_external_median = -np.inf
        best_external_values = []

        for other_cluster_id, other_members in cluster_to_genomes.items():
            if other_cluster_id == cluster_id:
                continue

            inter_values = []
            for g1 in members:
                for g2 in other_members:
                    if g1 in ani_matrix.index and g2 in ani_matrix.columns:
                        inter_values.append(ani_matrix.loc[g1, g2])

            inter_values = np.asarray(inter_values, dtype=float)
            inter_values = inter_values[~np.isnan(inter_values)]

            if len(inter_values) == 0:
                continue

            inter_median = float(np.median(inter_values))

            if inter_median > best_external_median:
                best_external_median = inter_median
                best_external_cluster = other_cluster_id
                best_external_values = inter_values

        p95_inter = safe_percentile(best_external_values, 95)
        median_inter = safe_percentile(best_external_values, 50)

        tail_gap = p5_intra - p95_inter if not np.isnan(p5_intra) and not np.isnan(p95_inter) else np.nan
        median_gap = median_intra - median_inter if not np.isnan(median_intra) and not np.isnan(median_inter) else np.nan

        if not np.isnan(p5_intra) and len(best_external_values) > 0:
            overlap_fraction = float(np.mean(best_external_values > p5_intra))
        else:
            overlap_fraction = np.nan

        gap_s = score_gap(tail_gap, median_gap)
        overlap_s = score_overlap(overlap_fraction)
        cohesion_s = score_cohesion(p5_intra)
        support_s = score_support(cluster_size, n_hq)

        components = np.array([gap_s, overlap_s, cohesion_s, support_s], dtype=float)
        weights = np.array([0.40, 0.30, 0.20, 0.10], dtype=float)

        if np.any(np.isnan(components)):
            cluster_score = np.nan
        else:
            cluster_score = float(np.sum(weights * components))

        rows.append({
            "cluster_id": cluster_id,
            "cluster_size": cluster_size,
            "N_HQ": n_hq,
            "HQ_ratio_cluster": hq_ratio,
            "nearest_external_cluster": best_external_cluster,
            "p5_intra_ANI": p5_intra,
            "median_intra_ANI": median_intra,
            "p95_nearest_inter_ANI": p95_inter,
            "median_nearest_inter_ANI": median_inter,
            "tail_gap": tail_gap,
            "median_gap": median_gap,
            "overlap_fraction": overlap_fraction,
            "gap_score": gap_s,
            "overlap_score": overlap_s,
            "cohesion_score": cohesion_s,
            "support_score": support_s,
            "cluster_structure_score": cluster_score,
            "validity_status": validity_status(cluster_score, cluster_size),
            "confidence_status": confidence_status(cluster_size, n_hq, hq_ratio),
        })

    return pd.DataFrame(rows)


def calculate_silhouette_for_genome(genome, own_cluster, clusters_by_id, ani_matrix):
    """
    Calculate ANI-based silhouette for one genome.

    Distance = 1 - ANI.

    Returns NaN for singleton clusters or if distances are missing.
    """
    own_members = [g for g in clusters_by_id[own_cluster] if g != genome]

    if len(own_members) == 0:
        return np.nan, np.nan

    own_distances = []
    for g in own_members:
        if genome in ani_matrix.index and g in ani_matrix.columns:
            ani = ani_matrix.loc[genome, g]
            if not np.isnan(ani):
                own_distances.append(1.0 - ani)

    if len(own_distances) == 0:
        return np.nan, np.nan

    a = float(np.mean(own_distances))

    nearest_cluster = None
    b = np.inf

    for cluster_id, members in clusters_by_id.items():
        if cluster_id == own_cluster:
            continue

        distances = []
        for g in members:
            if genome in ani_matrix.index and g in ani_matrix.columns:
                ani = ani_matrix.loc[genome, g]
                if not np.isnan(ani):
                    distances.append(1.0 - ani)

        if len(distances) == 0:
            continue

        mean_dist = float(np.mean(distances))
        if mean_dist < b:
            b = mean_dist
            nearest_cluster = cluster_id

    if nearest_cluster is None or np.isinf(b):
        return np.nan, np.nan

    denom = max(a, b)
    if denom == 0:
        silhouette = 0.0
    else:
        silhouette = (b - a) / denom

    return silhouette, nearest_cluster


def evaluate_genomes(clusters: pd.DataFrame, metadata: pd.DataFrame, ani_matrix: pd.DataFrame) -> pd.DataFrame:
    merged = clusters.merge(metadata, on="genome_id", how="left")
    merged["quality_status"] = merged["quality_status"].fillna("UNKNOWN")

    clusters_by_id = (
        merged
        .groupby("cluster_id")["genome_id"]
        .apply(list)
        .to_dict()
    )

    rows = []

    for row in merged.itertuples(index=False):
        silhouette, nearest_cluster = calculate_silhouette_for_genome(
            row.genome_id,
            row.cluster_id,
            clusters_by_id,
            ani_matrix,
        )

        rows.append({
            "genome_id": row.genome_id,
            "cluster_id": row.cluster_id,
            "quality_status": row.quality_status,
            "nearest_external_cluster": nearest_cluster,
            "silhouette_ANI": silhouette,
            "is_negative_silhouette": (
                bool(silhouette < 0) if not np.isnan(silhouette) else np.nan
            ),
        })

    return pd.DataFrame(rows)


def weighted_tool_score(clusters: pd.DataFrame, cluster_metrics: pd.DataFrame, genome_subset=None):
    """
    Genome-weighted mean cluster score.

    Singleton genomes are excluded because their cluster score is undefined.
    """
    tmp = clusters.merge(
        cluster_metrics[["cluster_id", "cluster_structure_score", "validity_status"]],
        on="cluster_id",
        how="left",
    )

    if genome_subset is not None:
        tmp = tmp[tmp["genome_id"].isin(genome_subset)]

    tmp = tmp[tmp["validity_status"] != "Singleton"]
    tmp = tmp.dropna(subset=["cluster_structure_score"])

    if len(tmp) == 0:
        return np.nan

    return float(tmp["cluster_structure_score"].mean())


def fraction_status(clusters_with_metrics, status, genome_subset=None):
    tmp = clusters_with_metrics.copy()

    if genome_subset is not None:
        tmp = tmp[tmp["genome_id"].isin(genome_subset)]

    if len(tmp) == 0:
        return np.nan

    return float(np.mean(tmp["validity_status"].eq(status)))


def tool_metrics(clusters: pd.DataFrame, metadata: pd.DataFrame, cluster_metrics: pd.DataFrame, genome_metrics: pd.DataFrame) -> pd.DataFrame:
    merged = clusters.merge(metadata, on="genome_id", how="left")
    merged["quality_status"] = merged["quality_status"].fillna("UNKNOWN")

    merged = merged.merge(
        cluster_metrics[["cluster_id", "cluster_size", "validity_status", "cluster_structure_score"]],
        on="cluster_id",
        how="left",
    )

    hq_genomes = set(merged.loc[merged["quality_status"].eq("HQ"), "genome_id"])

    def rate_singleton(df):
        if len(df) == 0:
            return np.nan
        return float(np.mean(df["cluster_size"].eq(1)))

    def rate_tiny(df):
        if len(df) == 0:
            return np.nan
        return float(np.mean(df["cluster_size"] < 5))

    hq_df = merged[merged["quality_status"].eq("HQ")]
    total_df = merged

    gm_hq = genome_metrics[genome_metrics["quality_status"].eq("HQ")]
    gm_total = genome_metrics

    def silhouette_summary(df, column):
        values = pd.to_numeric(df["silhouette_ANI"], errors="coerce").dropna()
        if len(values) == 0:
            return np.nan
        if column == "median":
            return float(np.median(values))
        if column == "p10":
            return float(np.percentile(values, 10))
        if column == "negative_fraction":
            return float(np.mean(values < 0))
        raise ValueError(column)

    tool_structure_score_hq = weighted_tool_score(clusters, cluster_metrics, genome_subset=hq_genomes)
    tool_structure_score_total = weighted_tool_score(clusters, cluster_metrics, genome_subset=None)

    defective_hq_fraction = fraction_status(merged, "Defective", genome_subset=hq_genomes)
    singleton_rate_hq = rate_singleton(hq_df)
    tiny_cluster_rate_hq = rate_tiny(hq_df)
    neg_sil_hq = silhouette_summary(gm_hq, "negative_fraction")

    if (
        tool_structure_score_hq >= 0.80
        and tool_structure_score_total >= 0.70
        and defective_hq_fraction <= 0.05
        and singleton_rate_hq <= 0.10
        and tiny_cluster_rate_hq <= 0.15
        and neg_sil_hq <= 0.05
    ):
        status = "Strong"
    elif (
        tool_structure_score_hq >= 0.65
        and tool_structure_score_total >= 0.55
        and defective_hq_fraction <= 0.20
        and singleton_rate_hq <= 0.20
        and tiny_cluster_rate_hq <= 0.30
        and neg_sil_hq <= 0.15
    ):
        status = "Moderate"
    elif (
        tool_structure_score_hq < 0.65
        or tool_structure_score_total < 0.55
        or defective_hq_fraction > 0.20
        or singleton_rate_hq > 0.30
        or tiny_cluster_rate_hq > 0.40
        or neg_sil_hq > 0.25
    ):
        status = "Weak"
    else:
        status = "Mixed"

    if status in ["Strong", "Moderate"]:
        action = "ACCEPT"
    elif singleton_rate_hq > 0.20 and defective_hq_fraction <= 0.20:
        action = "TRY_MORE_PERMISSIVE_BOUNDARY"
    elif defective_hq_fraction > 0.20 or neg_sil_hq > 0.25:
        action = "TRY_MORE_STRICT_BOUNDARY"
    elif status == "Mixed":
        action = "REVIEW"
    else:
        action = "REJECT"

    row = {
        "tool_structure_score_HQ": tool_structure_score_hq,
        "tool_structure_score_total": tool_structure_score_total,

        "accepted_HQ_fraction": fraction_status(merged, "Accepted", genome_subset=hq_genomes),
        "accepted_total_fraction": fraction_status(merged, "Accepted"),

        "borderline_HQ_fraction": fraction_status(merged, "Borderline", genome_subset=hq_genomes),
        "borderline_total_fraction": fraction_status(merged, "Borderline"),

        "defective_HQ_fraction": defective_hq_fraction,
        "defective_total_fraction": fraction_status(merged, "Defective"),

        "singleton_rate_HQ": singleton_rate_hq,
        "singleton_rate_total": rate_singleton(total_df),

        "tiny_cluster_rate_HQ": tiny_cluster_rate_hq,
        "tiny_cluster_rate_total": rate_tiny(total_df),

        "median_silhouette_HQ": silhouette_summary(gm_hq, "median"),
        "p10_silhouette_HQ": silhouette_summary(gm_hq, "p10"),
        "negative_silhouette_HQ_fraction": neg_sil_hq,

        "median_silhouette_total": silhouette_summary(gm_total, "median"),
        "p10_silhouette_total": silhouette_summary(gm_total, "p10"),
        "negative_silhouette_total_fraction": silhouette_summary(gm_total, "negative_fraction"),

        "tool_status": status,
        "recommended_action": action,
    }

    return pd.DataFrame([row])


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate PopPUNK cluster assignments using all-vs-all FastANI results."
    )
    parser.add_argument("--fastani", required=True, help="All-vs-all FastANI output TSV.")
    parser.add_argument("--clusters", required=True, help="PopPUNK *_clusters.csv file.")
    parser.add_argument("--labels", required=True, help="Per-genome HQ/MQ labels CSV (genome,label) from spp_eligibility_from_qc.py.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")

    args = parser.parse_args()

    fastani = read_fastani(args.fastani)
    clusters = read_clusters(args.clusters)
    metadata = read_labels(args.labels)

    clusters = clusters.merge(metadata[["genome_id"]], on="genome_id", how="inner")

    if clusters.empty:
        raise ValueError(
            "No overlapping genome IDs between PopPUNK clusters and metadata after normalisation."
        )

    genomes = sorted(clusters["genome_id"].unique())

    fastani_genomes = set(fastani["query"]).union(set(fastani["reference"]))
    missing_fastani = sorted(set(genomes) - fastani_genomes)

    if missing_fastani:
        print(
            f"WARNING: {len(missing_fastani)} genomes from clusters/metadata are missing from FastANI results.",
            file=sys.stderr,
        )
        print(
            "First missing genomes: " + ", ".join(missing_fastani[:10]),
            file=sys.stderr,
        )

    ani_matrix = make_symmetric_ani_matrix(fastani, genomes)

    cluster_metrics = evaluate_clusters(clusters, metadata, ani_matrix)
    genome_metrics = evaluate_genomes(clusters, metadata, ani_matrix)
    summary_metrics = tool_metrics(clusters, metadata, cluster_metrics, genome_metrics)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    cluster_metrics.to_csv(f"{args.out_prefix}.cluster_metrics.tsv", sep="\t", index=False)
    genome_metrics.to_csv(f"{args.out_prefix}.genome_metrics.tsv", sep="\t", index=False)
    summary_metrics.to_csv(f"{args.out_prefix}.tool_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
