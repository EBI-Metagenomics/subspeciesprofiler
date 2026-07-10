#!/usr/bin/env python3

import argparse
import csv
import itertools
import sys
from collections import Counter
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


def _read_csv_or_exit(path: str, description: str, **kwargs) -> pd.DataFrame:
    """Read a CSV/TSV, failing with a clear message instead of a pandas traceback."""
    try:
        return pd.read_csv(path, **kwargs)
    except FileNotFoundError:
        sys.exit(f"ERROR: {description} file not found: '{path}'.")
    except (pd.errors.EmptyDataError, csv.Error):
        sys.exit(f"ERROR: {description} file is empty or unreadable: '{path}'.")


def read_fastani(path: str) -> pd.DataFrame:
    """
    Read all-vs-all FastANI output.

    Expected columns (tab-separated, no header):
      query, reference, ani, fragments_mapped, total_fragments

    FastANI reports ANI as a percentage (0-100). FastANI is run internally by the
    pipeline, so we rely on that contract and convert to the 0-1 scale used
    throughout. The values are validated defensively: rather than guessing the
    scale from the data, we fail loudly if the FastANI output ever stops looking
    like a percentage, instead of silently mis-scaling every threshold downstream.
    """
    df = _read_csv_or_exit(
        path,
        "FastANI",
        sep="\t",
        header=None,
        names=["query", "reference", "ani", "fragments_mapped", "total_fragments"],
    )

    df["query"] = df["query"].map(normalise_genome_id)
    df["reference"] = df["reference"].map(normalise_genome_id)

    df["ani"] = pd.to_numeric(df["ani"], errors="coerce")

    ani = df["ani"].dropna()
    if ani.empty:
        sys.exit(f"ERROR: no numeric ANI values found in FastANI output '{path}'.")
    hi = float(ani.max())
    if hi <= 1.5:
        sys.exit(
            f"ERROR: FastANI ANI looks like a 0-1 fraction (max={hi:.4f}), but percentage "
            f"(0-100) output is required from the FastANI step ('{path}')."
        )
    if hi > 100.5:
        sys.exit(
            f"ERROR: FastANI ANI exceeds 100 (max={hi:.4f}); not valid percentage output ('{path}')."
        )

    df["ani"] = df["ani"] / 100.0

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
    df = _read_csv_or_exit(path, "clusters")

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
        sys.exit(
            f"ERROR: could not identify sample and cluster columns in clusters file '{path}'. "
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
    df = _read_csv_or_exit(path, "labels", sep=None, engine="python")

    lower_to_original = {c.lower(): c for c in df.columns}

    genome_col = None
    for candidate in ["genome", "genome_id", "sample", "id", "name"]:
        if candidate in lower_to_original:
            genome_col = lower_to_original[candidate]
            break

    if genome_col is None:
        sys.exit(
            f"ERROR: labels file '{path}' must contain a genome/genome_id column. "
            f"Found columns: {list(df.columns)}"
        )

    quality_col = None
    for candidate in ["label", "quality_status", "quality", "qc_status"]:
        if candidate in lower_to_original:
            quality_col = lower_to_original[candidate]
            break

    if quality_col is None:
        sys.exit(
            f"ERROR: labels file '{path}' must contain label, quality_status, quality, or qc_status. "
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
    elif n_hq >= 2 and cluster_size >= 3:
        return "Medium_confidence"
    elif n_hq >= 2:
        return "Low_confidence"
    else:
        return "Unresolved_low_HQ"


def evaluate_clusters(clusters: pd.DataFrame, metadata: pd.DataFrame, ani_matrix: pd.DataFrame, min_comparison_cluster_size: int) -> pd.DataFrame:
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
        # Only clusters of at least `min_comparison_cluster_size` count as a
        # competing cluster -- singletons are not treated as separate groupings
        # for separation metrics (over-splitting is penalised via singleton_rate /
        # tiny_cluster_rate instead), so a lone outlier can't distort the score.
        best_external_cluster = None
        best_external_median = -np.inf
        best_external_values = []

        for other_cluster_id, other_members in cluster_to_genomes.items():
            if other_cluster_id == cluster_id:
                continue
            if len(other_members) < min_comparison_cluster_size:
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


def calculate_silhouette_for_genome(genome, own_cluster, clusters_by_id, ani_matrix, min_comparison_cluster_size):
    """
    Calculate ANI-based silhouette for one genome.

    Distance = 1 - ANI.

    Returns NaN for singleton clusters or if distances are missing. The nearest
    other cluster (b) is chosen only among clusters of at least
    `min_comparison_cluster_size`, so singletons are not treated as competing
    clusters (a lone outlier near a real cluster should not tank its silhouette;
    over-splitting is penalised separately via singleton_rate).
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
        if len(members) < min_comparison_cluster_size:
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


def evaluate_genomes(clusters: pd.DataFrame, metadata: pd.DataFrame, ani_matrix: pd.DataFrame, min_comparison_cluster_size: int) -> pd.DataFrame:
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
            min_comparison_cluster_size,
        )

        rows.append({
            "genome_id": row.genome_id,
            "cluster_id": row.cluster_id,
            "quality_status": row.quality_status,
            "nearest_external_cluster": nearest_cluster,
            "silhouette_ANI": silhouette,
        })

    result = pd.DataFrame(rows)
    # Derive the negative-silhouette flag vectorised, as a nullable boolean, so the
    # column has one clean dtype (True/False/<NA>) rather than mixing Python bool
    # with float NaN (object dtype). NaN silhouettes occur for singleton genomes
    # (no same-cluster neighbour) and single-cluster inputs (no external cluster).
    sil = result["silhouette_ANI"]
    result["is_negative_silhouette"] = (sil < 0).astype("boolean").mask(sil.isna())
    return result


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


def tool_metrics(clusters: pd.DataFrame, metadata: pd.DataFrame, cluster_metrics: pd.DataFrame, genome_metrics: pd.DataFrame, model_name: str, accept_status: set) -> pd.DataFrame:
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

    # The decision hinges on the HQ genomes -- the QC-trustworthy signal, expected
    # to be a solid fraction of the group. `tool_structure_score_hq` is NaN exactly
    # when no HQ genome sits in a non-singleton cluster, i.e. the HQ genomes show no
    # cohesive structure. This covers a single cluster, all singletons, AND the case
    # where only the MQ genomes cluster: MQ-only structure is suspect because contig/
    # gene fragmentation in lower-quality assemblies distorts the composition signal,
    # so a real biological signal should be recoverable by another model that also
    # resolves the HQ genomes. Any of these means this model is unsuitable -- label
    # it "Weak" immediately (rather than letting the NaN thresholds fall through to
    # "Mixed") so the model-selection loop moves on to the next model.
    n_clusters = int(len(cluster_metrics))
    n_nonsingleton = int((cluster_metrics["cluster_size"] > 1).sum()) if n_clusters else 0
    no_hq_structure = bool(np.isnan(tool_structure_score_hq))

    reason = None
    if no_hq_structure:
        status = "Weak"
        if n_clusters <= 1:
            reason = "no evaluable HQ structure: single cluster (no partitioning)"
        elif n_nonsingleton == 0:
            reason = "no evaluable HQ structure: all singletons (no cohesive clusters)"
        else:
            reason = (
                "no evaluable HQ structure: structure only among MQ genomes, HQ unresolved "
                "(likely a fragmentation artifact, not true subspecies signal)"
            )
    elif (
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

    # Human-readable reason for the profiler-history report (the no-structure
    # cases already set one above).
    if reason is None:
        if status in ("Strong", "Moderate"):
            reason = f"{status.lower()} clustering (HQ structure score {tool_structure_score_hq:.2f})"
        elif status == "Weak":
            fails = []
            if tool_structure_score_hq < 0.65:
                fails.append(f"HQ structure score {tool_structure_score_hq:.2f} < 0.65")
            if tool_structure_score_total < 0.55:
                fails.append(f"total structure score {tool_structure_score_total:.2f} < 0.55")
            if defective_hq_fraction > 0.20:
                fails.append(f"defective HQ fraction {defective_hq_fraction:.2f} > 0.20")
            if singleton_rate_hq > 0.30:
                fails.append(f"HQ singleton rate {singleton_rate_hq:.2f} > 0.30")
            if tiny_cluster_rate_hq > 0.40:
                fails.append(f"tiny-cluster HQ rate {tiny_cluster_rate_hq:.2f} > 0.40")
            if neg_sil_hq > 0.25:
                fails.append(f"negative-silhouette HQ fraction {neg_sil_hq:.2f} > 0.25")
            reason = "weak clustering: " + ("; ".join(fails) if fails else "below Moderate thresholds")
        else:  # Mixed: structure/separation reach Moderate quality, but a
            # fragmentation/misplacement metric sits in the band between the
            # Moderate and Weak thresholds -- genuinely mixed signals (good on some
            # axes, borderline on others), not a lack of information.
            borderline = []
            if singleton_rate_hq > 0.20:
                borderline.append(f"HQ singleton rate {singleton_rate_hq:.2f} (> Moderate's 0.20)")
            if tiny_cluster_rate_hq > 0.30:
                borderline.append(f"tiny-cluster HQ rate {tiny_cluster_rate_hq:.2f} (> Moderate's 0.30)")
            if neg_sil_hq > 0.15:
                borderline.append(f"negative-silhouette HQ fraction {neg_sil_hq:.2f} (> Moderate's 0.15)")
            reason = "mixed signal: acceptable structure but " + (
                "; ".join(borderline) if borderline else "borderline over-splitting/misplacement"
            )

    # Binomial control signal for the model-selection loop: ACCEPT (stop) when the
    # quality label is in the acceptance set, else TRY_NEXT_MODEL (keep searching).
    decision = "ACCEPT" if status in accept_status else "TRY_NEXT_MODEL"

    print(
        f"[{model_name or 'model'}] tool_status={status} decision={decision} :: {reason}",
        file=sys.stderr,
    )

    row = {
        "model": model_name,

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
        "decision": decision,
        "reason": reason,
    }

    return pd.DataFrame([row])


def check_species_ani(ani_matrix, genomes, min_ani):
    """
    Verify every pair of genomes in the species group is within-species close.

    Genomes of the same species are expected at ANI >= ~0.95 (the species
    boundary). A pair below `min_ani` -- or absent from the all-vs-all FastANI
    entirely, which means FastANI could not align them (ANI below its ~0.80
    detection limit) -- indicates a genome too distant to belong to this species
    group.

    Returns a list of offending (g1, g2, ani_or_None) pairs; `None` marks a pair
    missing from FastANI.
    """
    problems = []
    n = len(genomes)
    for i in range(n):
        for j in range(i + 1, n):
            g1, g2 = genomes[i], genomes[j]
            ani = ani_matrix.loc[g1, g2]
            if np.isnan(ani):
                problems.append((g1, g2, None))
            elif ani < min_ani:
                problems.append((g1, g2, float(ani)))
    return problems


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate PopPUNK cluster assignments using all-vs-all FastANI results."
    )
    parser.add_argument("--fastani", required=True, help="All-vs-all FastANI output TSV.")
    parser.add_argument("--clusters", required=True, help="PopPUNK *_clusters.csv file.")
    parser.add_argument("--labels", required=True, help="Per-genome HQ/MQ labels CSV (genome,label) from spp_eligibility_from_qc.py.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    parser.add_argument(
        "--min-ani", type=float, default=0.90,
        help="Minimum within-species ANI (0-1 scale). Fail if any genome pair in the "
             "species group is below this, or missing from FastANI (below its ~0.80 "
             "detection limit). Default 0.90: flags genomes likely misplaced during data "
             "preparation while tolerating borderline within-species cases (species "
             "boundary ~0.95).",
    )
    parser.add_argument(
        "--model-name", default="",
        help="Name of the PopPUNK model that produced these clusters (e.g. bgmm, dbscan, "
             "lineage), recorded in the output for the profiler-history report.",
    )
    parser.add_argument(
        "--min-comparison-cluster-size", type=int, default=2,
        help="Minimum cluster size for a cluster to count as a competing/nearest cluster "
             "in the separation metrics (silhouette b and nearest-external ANI). Default 2 "
             "excludes singletons -- they are not treated as genuine groupings for "
             "separation (over-splitting is penalised via singleton_rate instead).",
    )
    parser.add_argument(
        "--accept-status", default="Strong",
        help="Comma-separated tool_status values that count as acceptable, i.e. "
             "decision=ACCEPT (stop the model search). Default: Strong (conservative -- "
             "keep searching until a clean result; the loop's best-so-far handles "
             "the case where nothing is Strong).",
    )

    args = parser.parse_args()
    accept_status = {s.strip() for s in args.accept_status.split(",") if s.strip()}

    fastani = read_fastani(args.fastani)
    clusters = read_clusters(args.clusters)
    metadata = read_labels(args.labels)

    clusters = clusters.merge(metadata[["genome_id"]], on="genome_id", how="inner")

    if clusters.empty:
        sys.exit(
            "ERROR: no overlapping genome IDs between PopPUNK clusters and the labels file "
            "after normalisation (check that the qc_csv genome names match the assembly filenames)."
        )

    genomes = sorted(clusters["genome_id"].unique())

    ani_matrix = make_symmetric_ani_matrix(fastani, genomes)

    # Fail fast if the species group contains genomes too distant to belong to
    # the same species: same-species genomes should be well within FastANI's
    # detection range, so a missing or low-ANI pair signals a mis-assigned genome.
    distant_pairs = check_species_ani(ani_matrix, genomes, args.min_ani)
    if distant_pairs:
        counts = Counter()
        for g1, g2, _ in distant_pairs:
            counts[g1] += 1
            counts[g2] += 1
        print(
            f"ERROR: {len(distant_pairs)} genome pair(s) fall below the within-species "
            f"ANI threshold ({args.min_ani:.2f}); the species group contains genomes too "
            f"distant to belong to the same species. Remove the offending genome(s) and rerun.",
            file=sys.stderr,
        )
        worst = ", ".join(f"{g} ({c} bad pairs)" for g, c in counts.most_common(10))
        print(f"Most-distant genomes: {worst}", file=sys.stderr)
        for g1, g2, ani in distant_pairs[:20]:
            shown = "missing from FastANI (ANI below ~0.80 detection limit)" if ani is None else f"ANI={ani:.3f}"
            print(f"  {g1} vs {g2}: {shown}", file=sys.stderr)
        if len(distant_pairs) > 20:
            print(f"  ... and {len(distant_pairs) - 20} more pair(s)", file=sys.stderr)
        sys.exit(1)

    cluster_metrics = evaluate_clusters(clusters, metadata, ani_matrix, args.min_comparison_cluster_size)
    genome_metrics = evaluate_genomes(clusters, metadata, ani_matrix, args.min_comparison_cluster_size)
    summary_metrics = tool_metrics(clusters, metadata, cluster_metrics, genome_metrics, args.model_name, accept_status)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    cluster_metrics.to_csv(f"{args.out_prefix}.cluster_metrics.tsv", sep="\t", index=False)
    genome_metrics.to_csv(f"{args.out_prefix}.genome_metrics.tsv", sep="\t", index=False)
    summary_metrics.to_csv(f"{args.out_prefix}.tool_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
