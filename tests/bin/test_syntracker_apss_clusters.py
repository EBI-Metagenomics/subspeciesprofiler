"""
Tests for bin/syntracker_apss_clusters.py (SynTracker all-regions APSS -> Leiden -> Taxon,Cluster).

Pins: one reference asserted, pairs on too few regions dropped, low-coverage targets excluded,
edges pruned at min_apss, Leiden communities, targets without edges as singletons, and the
per-region SD / APSS standard-error estimate from a subsampled table.

Needs python-igraph (as in modules/local/syntracker/clusters/environment.yml).

Run with:  python3 -m pytest tests/bin/test_syntracker_apss_clusters.py
"""

import itertools
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("igraph")

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_apss_clusters.py"
COLUMNS = ["Ref_genome", "Sample1", "Sample2", "APSS", "Compared_regions"]


def write_apss(path, scores, regions=300, ref="REF"):
    """scores: {(a, b): apss}; regions: an int or {(a, b): regions}."""
    rows = [(ref, f"{a}.fasta", f"{b}.fasta", s, regions[(a, b)] if isinstance(regions, dict) else regions)
            for (a, b), s in scores.items()]
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, index=False)


def two_groups(within=0.90, between=0.60, bridges=()):
    """Groups a1..a5 and b1..b5; `bridges` lists (a, b) pairs scored like within-group pairs."""
    a = [f"a{i}" for i in range(1, 6)]
    b = [f"b{i}" for i in range(1, 6)]
    scores = {}
    for x, y in itertools.combinations(a + b, 2):
        same = (x[0] == y[0])
        scores[(x, y)] = within if same or (x, y) in bridges or (y, x) in bridges else between
    return scores


def run(tmp_path, scores, regions=300, min_apss="0.75", resolutions="1.0", extra=(), tag="run"):
    d = tmp_path / tag
    d.mkdir(exist_ok=True)
    write_apss(d / "avg_synteny_scores_all_regions.csv", scores, regions)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--apss", str(d / "avg_synteny_scores_all_regions.csv"),
         "--min-apss", min_apss, "--resolutions", resolutions, "--out-prefix", str(d / "sp"), *extra],
        capture_output=True, text=True,
    )
    return proc, d


def clusters(d, model="leiden_apss0.75_r1.0"):
    df = pd.read_csv(d / f"sp_syntracker_{model}_clusters.csv")
    assert list(df.columns) == ["Taxon", "Cluster"]
    return dict(zip(df["Taxon"], df["Cluster"]))


def n_clusters(assignment):
    return len(set(assignment.values()))


# --- clustering ----------------------------------------------------------------------------------

def test_two_groups_are_two_communities(tmp_path):
    proc, d = run(tmp_path, two_groups())
    assert proc.returncode == 0, proc.stderr
    a = clusters(d)
    assert n_clusters(a) == 2
    assert len({a[f"a{i}"] for i in range(1, 6)}) == 1
    assert a["a1"] != a["b1"]


def test_leiden_splits_groups_joined_by_one_surviving_edge(tmp_path):
    proc, d = run(tmp_path, two_groups(bridges=[("a1", "b1")]))
    assert proc.returncode == 0, proc.stderr
    a = clusters(d)
    assert n_clusters(a) == 2
    assert a["a2"] != a["b2"]


def test_sweep_naming_and_full_pruning(tmp_path):
    proc, d = run(tmp_path, two_groups(), min_apss="0.50,0.75,0.95", resolutions="0.5,1.0")
    assert proc.returncode == 0, proc.stderr
    for t in ("0.50", "0.75", "0.95"):
        for r in ("0.5", "1.0"):
            assert (d / f"sp_syntracker_leiden_apss{t}_r{r}_clusters.csv").exists()
    assert n_clusters(clusters(d, "leiden_apss0.95_r1.0")) == 10  # every edge pruned: all singletons


def test_same_seed_same_partition(tmp_path):
    scores = two_groups(within=0.85, between=0.78, bridges=[("a1", "b1"), ("a3", "b4")])
    _, d1 = run(tmp_path, scores, tag="s1")
    _, d2 = run(tmp_path, scores, tag="s2")
    assert clusters(d1) == clusters(d2)


def test_largest_cluster_is_numbered_one(tmp_path):
    scores = {k: v for k, v in two_groups().items() if "b5" not in k}
    proc, d = run(tmp_path, scores)
    assert proc.returncode == 0, proc.stderr
    a = clusters(d)
    assert a["a1"] == 1 and a["b1"] == 2


# --- filters -------------------------------------------------------------------------------------

def test_pairs_on_too_few_regions_are_dropped(tmp_path):
    # every a1 pair is on too few regions, so a1 has no pair left
    scores = two_groups()
    regions = {k: (50 if "a1" in k else 300) for k in scores}
    proc, d = run(tmp_path, scores, regions=regions, extra=["--min-genome-coverage", "0"])
    assert proc.returncode == 0, proc.stderr
    # with no pair left, a1 has no measurable coverage: excluded, whatever the coverage threshold
    assert "a1" not in clusters(d)
    assert "9 below 100 regions dropped" in proc.stderr


def test_low_coverage_target_is_excluded(tmp_path):
    scores = two_groups()
    regions = {k: (120 if "b5" in k else 300) for k in scores}  # b5: 0.4 of the species median
    proc, d = run(tmp_path, scores, regions=regions, extra=["--min-genome-coverage", "0.5"])
    assert proc.returncode == 0, proc.stderr
    assert "b5" not in clusters(d)
    cov = pd.read_csv(d / "sp_syntracker_genome_coverage.tsv", sep="\t").set_index("genome")
    assert list(cov.columns) == ["median_regions", "relative_coverage", "status"]
    assert cov.loc["b5", "status"] == "low_coverage" and cov.loc["a1", "status"] == "ok"
    assert cov.loc["b5", "relative_coverage"] == pytest.approx(0.4)


def test_target_without_pairs_is_low_coverage_and_non_target_ignored(tmp_path):
    d = tmp_path / "tg"
    d.mkdir()
    names = [f"a{i}" for i in range(1, 6)] + [f"b{i}" for i in range(1, 5)] + ["lonely"]  # b5 not a target
    (d / "targets.txt").write_text("".join(f"{n}.fna.gz\n" for n in names))
    proc, d = run(tmp_path, two_groups(), extra=["--targets", str(d / "targets.txt")], tag="tg")
    assert proc.returncode == 0, proc.stderr
    a = clusters(d)
    assert "b5" not in a and "lonely" not in a
    assert "non-target" in proc.stderr


# --- noise ---------------------------------------------------------------------------------------

def test_noise_from_subsampled_table(tmp_path):
    rng = np.random.default_rng(0)
    scores = two_groups()
    sd, n, regions = 0.22, 40, 300
    noisy = {k: v + rng.normal(0, sd * np.sqrt(1 / n - 1 / regions)) for k, v in scores.items()}
    sub = tmp_path / "avg_synteny_scores_40_regions.csv"
    write_apss(sub, noisy, regions=n)
    proc, d = run(tmp_path, scores, regions=regions, extra=["--apss-subsampled", str(sub)])
    assert proc.returncode == 0, proc.stderr
    noise = pd.read_csv(d / "sp_syntracker_noise.tsv", sep="\t").iloc[0]
    assert noise["region_sd"] == pytest.approx(sd, rel=0.3)
    assert noise["apss_se"] == pytest.approx(noise["region_sd"] / np.sqrt(regions), rel=1e-3)  # written to 5 decimals
    assert noise["sd_source"] == sub.name


def test_noise_default_without_subsampled_table(tmp_path):
    proc, d = run(tmp_path, two_groups(), regions=400)
    assert proc.returncode == 0, proc.stderr
    noise = pd.read_csv(d / "sp_syntracker_noise.tsv", sep="\t").iloc[0]
    assert noise["sd_source"] == "default"
    assert noise["apss_se"] == pytest.approx(0.226 / 20)


def test_identical_subsampled_table_falls_back_to_default(tmp_path):
    sub = tmp_path / "avg_synteny_scores_40_regions.csv"
    write_apss(sub, two_groups(), regions=40)  # no noise at all: not a usable estimate
    proc, d = run(tmp_path, two_groups(), extra=["--apss-subsampled", str(sub)])
    assert proc.returncode == 0, proc.stderr
    assert pd.read_csv(d / "sp_syntracker_noise.tsv", sep="\t").iloc[0]["sd_source"] == "default"


# --- input contract ------------------------------------------------------------------------------

def test_more_than_one_reference_is_an_error(tmp_path):
    rows = [("R1", "x.fasta", "y.fasta", 0.9, 300), ("R2", "x.fasta", "z.fasta", 0.9, 300)]
    pd.DataFrame(rows, columns=COLUMNS).to_csv(tmp_path / "avg_synteny_scores_all_regions.csv", index=False)
    proc = subprocess.run([sys.executable, str(SCRIPT), "--apss", str(tmp_path / "avg_synteny_scores_all_regions.csv"),
                           "--out-prefix", str(tmp_path / "sp")], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "reference genomes" in proc.stderr


@pytest.mark.parametrize("content, message", [
    ("Ref_genome,Sample1,Sample2,Compared_regions\nR,x,y,50\n", "lacks columns"),
    ("Ref_genome,Sample1,Sample2,APSS,Compared_regions\n", "no pairs"),
])
def test_bad_apss_fails_clearly(tmp_path, content, message):
    (tmp_path / "avg_synteny_scores_all_regions.csv").write_text(content)
    proc = subprocess.run([sys.executable, str(SCRIPT), "--apss", str(tmp_path / "avg_synteny_scores_all_regions.csv"),
                           "--out-prefix", str(tmp_path / "sp")], capture_output=True, text=True)
    assert proc.returncode != 0
    assert message in proc.stderr
