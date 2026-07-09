"""
Unit + end-to-end tests for bin/evaluate_poppunk_fastani.py.

Captures the behaviour worked out while hardening the evaluator (weaknesses
B1-B9): the FastANI percent-scale contract, the within-species ANI gate, the
HQ-centric "no structure" logic, the audit-label / binomial-decision split, the
singleton-exclusion rule for separation metrics, and clean handling of
degenerate inputs.

Run with:  python3 -m pytest tests/bin/test_evaluate_poppunk_fastani.py
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "evaluate_poppunk_fastani.py"


@pytest.fixture(scope="session")
def mod():
    """Import the evaluator script as a module (bin/ is not a package)."""
    spec = importlib.util.spec_from_file_location("evaluate_poppunk_fastani", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def write_fastani(path, rows):
    """rows: iterable of (query, reference, ani, frags_mapped, total_frags)."""
    with open(path, "w") as fh:
        for q, r, a, fm, tf in rows:
            fh.write(f"{q}\t{r}\t{a}\t{fm}\t{tf}\n")


def all_vs_all(cluster_of, intra=99.5, inter=96.0, self_ani=100.0, overrides=None):
    """Build all-vs-all FastANI rows (percent scale) from {genome: cluster}.

    overrides: {(a, b): ani} for specific pairs (matched symmetrically).
    """
    overrides = overrides or {}
    genomes = list(cluster_of)
    rows = []
    for a in genomes:
        for b in genomes:
            if (a, b) in overrides:
                ani = overrides[(a, b)]
            elif (b, a) in overrides:
                ani = overrides[(b, a)]
            elif a == b:
                ani = self_ani
            else:
                ani = intra if cluster_of[a] == cluster_of[b] else inter
            rows.append((f"{a}.fna", f"{b}.fna", ani, 950, 1000))
    return rows


def ani_matrix_from(cluster_of, intra=0.995, inter=0.96, self_ani=1.0):
    """0-1 symmetric ANI matrix (as make_symmetric_ani_matrix would produce)."""
    gs = list(cluster_of)
    m = pd.DataFrame(index=gs, columns=gs, dtype=float)
    for a in gs:
        for b in gs:
            m.loc[a, b] = self_ani if a == b else (
                intra if cluster_of[a] == cluster_of[b] else inter
            )
    return m


def write_clusters(path, cluster_of):
    with open(path, "w") as fh:
        fh.write("Taxon,Cluster\n")
        for g, c in cluster_of.items():
            fh.write(f"{g},{c}\n")


def write_labels(path, label_of):
    with open(path, "w") as fh:
        fh.write("genome,label\n")
        for g, lab in label_of.items():
            fh.write(f"{g},{lab}\n")


def run_cli(tmp_path, cluster_of, label_of, ani_rows=None, extra_args=(), tag="run", **ani_kwargs):
    """Run the evaluator CLI end-to-end; return (CompletedProcess, out_prefix)."""
    d = tmp_path / tag
    d.mkdir(exist_ok=True)
    fa, cl, lb, out = d / "ani.tsv", d / "cl.csv", d / "lb.csv", d / "out"
    write_fastani(fa, ani_rows if ani_rows is not None else all_vs_all(cluster_of, **ani_kwargs))
    write_clusters(cl, cluster_of)
    write_labels(lb, label_of)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--fastani", str(fa), "--clusters", str(cl),
         "--labels", str(lb), "--out-prefix", str(out), *extra_args],
        capture_output=True, text=True,
    )
    return proc, out


def read_metrics(out, kind):
    p = Path(f"{out}.{kind}.tsv")
    return pd.read_csv(p, sep="\t") if p.exists() else None


# --------------------------------------------------------------------------- #
# read_fastani -- percent contract + defensive validation (B5)
# --------------------------------------------------------------------------- #

def test_read_fastani_scales_percent_to_fraction(mod, tmp_path):
    fa = tmp_path / "a.tsv"
    write_fastani(fa, [("g1.fna", "g2.fna", 96.0, 950, 1000)])
    df = mod.read_fastani(str(fa))
    assert list(df["ani"]) == [0.96]
    assert "alignment_fraction" not in df.columns  # dead code removed


def test_read_fastani_rejects_fraction_input(mod, tmp_path):
    fa = tmp_path / "a.tsv"
    write_fastani(fa, [("g1.fna", "g2.fna", 0.96, 950, 1000)])
    with pytest.raises(SystemExit):
        mod.read_fastani(str(fa))


def test_read_fastani_rejects_over_100(mod, tmp_path):
    fa = tmp_path / "a.tsv"
    write_fastani(fa, [("g1.fna", "g2.fna", 150.0, 950, 1000)])
    with pytest.raises(SystemExit):
        mod.read_fastani(str(fa))


def test_read_fastani_empty_file(mod, tmp_path):
    fa = tmp_path / "a.tsv"
    fa.write_text("")
    with pytest.raises(SystemExit):
        mod.read_fastani(str(fa))


# --------------------------------------------------------------------------- #
# read_clusters / read_labels parsing + clean errors (B3, B7, B8)
# --------------------------------------------------------------------------- #

def test_read_clusters_standard_taxon_cluster(mod, tmp_path):
    cl = tmp_path / "c.csv"
    write_clusters(cl, {"MGYG1": 1, "MGYG2": 2})  # lineage output is this same format
    df = mod.read_clusters(str(cl))
    assert set(df.columns) == {"genome_id", "cluster_id"}
    assert df.set_index("genome_id")["cluster_id"].to_dict() == {"MGYG1": "1", "MGYG2": "2"}


def test_read_clusters_bad_columns_exits(mod, tmp_path):
    cl = tmp_path / "c.csv"
    cl.write_text("foo,bar\n1,2\n")
    with pytest.raises(SystemExit):
        mod.read_clusters(str(cl))


def test_read_clusters_empty_exits(mod, tmp_path):
    cl = tmp_path / "c.csv"
    cl.write_text("")
    with pytest.raises(SystemExit):
        mod.read_clusters(str(cl))


def test_read_labels_parses_genome_label(mod, tmp_path):
    lb = tmp_path / "l.csv"
    write_labels(lb, {"g1": "HQ", "g2": "MQ"})
    df = mod.read_labels(str(lb))
    assert df.set_index("genome_id")["quality_status"].to_dict() == {"g1": "HQ", "g2": "MQ"}


def test_read_labels_missing_label_column_exits(mod, tmp_path):
    lb = tmp_path / "l.csv"
    lb.write_text("genome,foo\ng1,x\n")
    with pytest.raises(SystemExit):
        mod.read_labels(str(lb))


# --------------------------------------------------------------------------- #
# normalise_genome_id -- reconciles the three ID spaces (B7)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("raw,expected", [
    ("/path/to/genome1.fna.gz", "genome1"),
    ("genome1.fasta", "genome1"),
    ("genome1.fa", "genome1"),
    ("MGYG000001345.fasta", "MGYG000001345"),
    ("already_bare", "already_bare"),
])
def test_normalise_genome_id(mod, raw, expected):
    assert mod.normalise_genome_id(raw) == expected


# --------------------------------------------------------------------------- #
# check_species_ani -- within-species gate (B1)
# --------------------------------------------------------------------------- #

def test_check_species_ani_all_close_passes(mod):
    genomes = ["g1", "g2", "g3"]
    m = pd.DataFrame(0.98, index=genomes, columns=genomes)
    np.fill_diagonal(m.values, 1.0)
    assert mod.check_species_ani(m, genomes, 0.90) == []


def test_check_species_ani_flags_low_pair(mod):
    genomes = ["g1", "g2"]
    m = pd.DataFrame([[1.0, 0.80], [0.80, 1.0]], index=genomes, columns=genomes)
    problems = mod.check_species_ani(m, genomes, 0.90)
    assert len(problems) == 1 and problems[0][2] == pytest.approx(0.80)


def test_check_species_ani_flags_missing_pair(mod):
    genomes = ["g1", "g2"]
    m = pd.DataFrame([[1.0, np.nan], [np.nan, 1.0]], index=genomes, columns=genomes)
    problems = mod.check_species_ani(m, genomes, 0.90)
    assert len(problems) == 1 and problems[0][2] is None  # None == "missing from FastANI"


def test_cli_min_ani_gate_fails_on_distant_genome(tmp_path):
    cluster_of = {"g1": 1, "g2": 1, "g3": 2, "g4": 2}
    label_of = {g: "HQ" for g in cluster_of}
    overrides = {("g4", g): 85.0 for g in ["g1", "g2", "g3"]}  # g4 is too distant
    rows = all_vs_all(cluster_of, overrides=overrides)
    proc, out = run_cli(tmp_path, cluster_of, label_of, ani_rows=rows)
    assert proc.returncode == 1
    assert "below the within-species ANI threshold" in proc.stderr
    assert read_metrics(out, "tool_metrics") is None


# --------------------------------------------------------------------------- #
# confidence_status -- simplified tiers (B9 dead-branch cleanup)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("size,n_hq,ratio,expected", [
    (10, 5, 0.5, "High_confidence"),
    (5, 3, 0.4, "High_confidence"),
    (5, 3, 0.39, "Medium_confidence"),   # hq_ratio just below 0.40
    (4, 3, 0.9, "Medium_confidence"),    # size < 5 -> not High
    (3, 2, 0.9, "Medium_confidence"),
    (2, 2, 1.0, "Low_confidence"),
    (5, 1, 0.2, "Unresolved_low_HQ"),
    (1, 0, 0.0, "Unresolved_low_HQ"),
])
def test_confidence_status(mod, size, n_hq, ratio, expected):
    assert mod.confidence_status(size, n_hq, ratio) == expected


# --------------------------------------------------------------------------- #
# HQ-centric "no structure" -> Weak / TRY_NEXT_MODEL (B2)
# --------------------------------------------------------------------------- #

def test_cli_single_cluster_is_weak(tmp_path):
    cluster_of = {f"g{i}": 1 for i in range(1, 6)}
    label_of = {g: "HQ" for g in cluster_of}
    proc, out = run_cli(tmp_path, cluster_of, label_of)
    tm = read_metrics(out, "tool_metrics").iloc[0]
    assert proc.returncode == 0
    assert tm["tool_status"] == "Weak"
    assert tm["decision"] == "TRY_NEXT_MODEL"
    assert "single cluster" in tm["reason"]


def test_cli_all_singletons_is_weak(tmp_path):
    cluster_of = {f"g{i}": i for i in range(1, 6)}
    label_of = {g: "HQ" for g in cluster_of}
    proc, out = run_cli(tmp_path, cluster_of, label_of)
    tm = read_metrics(out, "tool_metrics").iloc[0]
    assert tm["tool_status"] == "Weak"
    assert "all singletons" in tm["reason"]


def test_cli_mq_only_structure_is_weak(tmp_path):
    # HQ genomes are singletons; only the MQ genomes cluster -> suspect (fragmentation)
    cluster_of = {"g1": 1, "g2": 2, "g3": 3, "g4": 3}
    label_of = {"g1": "HQ", "g2": "HQ", "g3": "MQ", "g4": "MQ"}
    proc, out = run_cli(tmp_path, cluster_of, label_of)
    tm = read_metrics(out, "tool_metrics").iloc[0]
    assert tm["tool_status"] == "Weak"
    assert "only among MQ" in tm["reason"]


# --------------------------------------------------------------------------- #
# Audit label vs binomial decision (B4)
# --------------------------------------------------------------------------- #

def test_cli_records_model_name(tmp_path):
    cluster_of = {f"g{i}": 1 for i in range(1, 4)}
    label_of = {g: "HQ" for g in cluster_of}
    _, out = run_cli(tmp_path, cluster_of, label_of, extra_args=["--model-name", "bgmm"])
    assert read_metrics(out, "tool_metrics").iloc[0]["model"] == "bgmm"


def test_cli_accept_status_controls_decision(tmp_path):
    # Two size-3 clusters: valid but "tiny" -> not Strong; default Strong-only => TRY_NEXT_MODEL.
    cluster_of = {f"g{i}": (1 if i <= 3 else 2) for i in range(1, 7)}
    label_of = {g: "HQ" for g in cluster_of}
    _, out = run_cli(tmp_path, cluster_of, label_of, tag="default")
    tm = read_metrics(out, "tool_metrics").iloc[0]
    status = tm["tool_status"]
    assert tm["decision"] == "TRY_NEXT_MODEL"
    # Widen the acceptance bar to include this status -> decision flips to ACCEPT.
    _, out2 = run_cli(tmp_path, cluster_of, label_of, extra_args=["--accept-status", status], tag="widened")
    assert read_metrics(out2, "tool_metrics").iloc[0]["decision"] == "ACCEPT"


def test_cli_mixed_signal_reason(tmp_path):
    # cluster1 size 7 (not tiny) + cluster2 size 4 (tiny) -> tiny_cluster_rate_HQ ~0.36,
    # in the (0.30, 0.40] band between Moderate and Weak -> tool_status = Mixed.
    cluster_of = {f"g{i}": (1 if i <= 7 else 2) for i in range(1, 12)}
    label_of = {g: "HQ" for g in cluster_of}
    _, out = run_cli(tmp_path, cluster_of, label_of, extra_args=["--model-name", "refine"])
    tm = read_metrics(out, "tool_metrics").iloc[0]
    assert tm["tool_status"] == "Mixed"
    assert tm["reason"].startswith("mixed signal")
    assert "tiny-cluster HQ rate" in tm["reason"]


# --------------------------------------------------------------------------- #
# Degenerate inputs -> clean errors, not tracebacks (B8)
# --------------------------------------------------------------------------- #

def test_cli_zero_overlap_clean_exit(tmp_path):
    cluster_of = {"g1": 1, "g2": 1, "g3": 2, "g4": 2}
    label_of = {"unrelated1": "HQ", "unrelated2": "HQ"}
    proc, _ = run_cli(tmp_path, cluster_of, label_of)
    assert proc.returncode == 1
    assert "no overlapping genome ids" in proc.stderr.lower()
    assert "Traceback" not in proc.stderr


def test_cli_empty_clusters_file_clean_exit(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    fa, cl, lb = d / "ani.tsv", d / "cl.csv", d / "lb.csv"
    write_fastani(fa, all_vs_all({"g1": 1, "g2": 1}))
    cl.write_text("")
    write_labels(lb, {"g1": "HQ", "g2": "HQ"})
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--fastani", str(fa), "--clusters", str(cl),
         "--labels", str(lb), "--out-prefix", str(d / "out")],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1
    assert "empty or unreadable" in proc.stderr
    assert "Traceback" not in proc.stderr


# --------------------------------------------------------------------------- #
# Singletons excluded from separation metrics (B9 methodology fix)
# --------------------------------------------------------------------------- #

def test_min_comparison_size_excludes_singleton_from_nearest_external(tmp_path):
    # Two real clusters (inter 96%) + a singleton g7 sitting at 99% next to cluster 1.
    cluster_of = {"g1": 1, "g2": 1, "g3": 1, "g4": 2, "g5": 2, "g6": 2, "g7": 3}
    label_of = {g: "HQ" for g in cluster_of}
    overrides = {("g7", g): 99.0 for g in ["g1", "g2", "g3"]}
    overrides.update({("g7", g): 95.0 for g in ["g4", "g5", "g6"]})
    rows = all_vs_all(cluster_of, overrides=overrides)

    # Default (min size 2): singleton excluded -> cluster 1's nearest external is the real cluster 2.
    _, out = run_cli(tmp_path, cluster_of, label_of, ani_rows=rows, tag="default")
    cm = read_metrics(out, "cluster_metrics")
    nearest = cm.set_index("cluster_id").loc[1, "nearest_external_cluster"]
    assert int(nearest) == 2

    # min size 1: singleton counts -> cluster 1's nearest external becomes the singleton (3).
    _, out1 = run_cli(tmp_path, cluster_of, label_of, ani_rows=rows,
                      extra_args=["--min-comparison-cluster-size", "1"], tag="minsize1")
    cm1 = read_metrics(out1, "cluster_metrics")
    nearest1 = cm1.set_index("cluster_id").loc[1, "nearest_external_cluster"]
    assert int(nearest1) == 3


def test_is_negative_silhouette_is_nullable_boolean(mod):
    # g5 is a singleton -> its silhouette is NaN, alongside real values from the others.
    cluster_of = {"g1": 1, "g2": 1, "g3": 2, "g4": 2, "g5": 3}
    genomes = list(cluster_of)
    clusters = pd.DataFrame({"genome_id": genomes,
                             "cluster_id": [str(cluster_of[g]) for g in genomes]})
    metadata = pd.DataFrame({"genome_id": genomes, "quality_status": ["HQ"] * len(genomes)})
    matrix = ani_matrix_from(cluster_of)
    gm = mod.evaluate_genomes(clusters, metadata, matrix, 2)
    assert str(gm["is_negative_silhouette"].dtype) == "boolean"
    # the singleton's flag is <NA>, not False
    assert pd.isna(gm.set_index("genome_id").loc["g5", "is_negative_silhouette"])
