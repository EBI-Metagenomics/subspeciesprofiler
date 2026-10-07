"""
Tests for bin/syntracker_apss_clusters.py (SynTracker APSS -> Leiden -> Taxon,Cluster tables).

Pins the paper's recipe: one subsampling depth (chosen at the retention cliff), one reference,
Compared_regions >= n asserted on per-n tables, edges pruned at min_apss, Leiden communities,
targets without edges as singletons and non-targets absent.

Needs python-igraph and matplotlib (as in modules/local/syntracker/clusters/environment.yml).

Run with:  python3 -m pytest tests/bin/test_syntracker_apss_clusters.py
"""

import itertools
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("igraph")
pytest.importorskip("matplotlib")

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_apss_clusters.py"
COLUMNS = ["Ref_genome", "Sample1", "Sample2", "APSS", "Compared_regions"]


def write_apss(path, scores, regions=100, ref="REF"):
    """scores: {(a, b): apss}. Written as SynTracker does, with target file names as sample names."""
    rows = [(ref, f"{a}.fasta", f"{b}.fasta", s, regions) for (a, b), s in scores.items()]
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


def call(d, apss_files, *extra):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--apss", *map(str, apss_files), "--out-prefix", str(d / "sp"), *extra],
        capture_output=True, text=True,
    )


def run(tmp_path, scores, depth="40", min_apss="0.75", resolutions="1.0", extra=(), tag="run"):
    """One per-n table at `depth`, clustered at that fixed depth."""
    d = tmp_path / tag
    d.mkdir(exist_ok=True)
    apss = d / f"avg_synteny_scores_{depth}_regions.csv"
    write_apss(apss, scores, regions=100 if depth == "all" else int(depth))
    proc = call(d, [apss], "--depth", depth, "--min-apss", min_apss, "--resolutions", resolutions, *extra)
    return proc, d


def clusters(d, model):
    df = pd.read_csv(d / f"sp_syntracker_{model}_clusters.csv")
    assert list(df.columns) == ["Taxon", "Cluster"]
    return dict(zip(df["Taxon"], df["Cluster"]))


def n_clusters(assignment):
    return len(set(assignment.values()))


# --- clustering ----------------------------------------------------------------------------------

def test_two_groups_are_two_communities(tmp_path):
    proc, d = run(tmp_path, two_groups())
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "leiden_n40_apss0.75_r1.0")
    assert n_clusters(a) == 2
    assert len({a[f"a{i}"] for i in range(1, 6)}) == 1
    assert a["a1"] != a["b1"]


def test_leiden_splits_groups_joined_by_one_surviving_edge(tmp_path):
    # One bridge survives pruning and connects the groups into one component;
    # modularity still separates the two dense groups.
    proc, d = run(tmp_path, two_groups(bridges=[("a1", "b1")]))
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "leiden_n40_apss0.75_r1.0")
    assert n_clusters(a) == 2
    assert a["a2"] != a["b2"]


def test_pruning_floor_sweep_and_naming(tmp_path):
    proc, d = run(tmp_path, two_groups(), min_apss="0.50,0.75,0.95", resolutions="0.5,1.0")
    assert proc.returncode == 0, proc.stderr
    for t in ("0.50", "0.75", "0.95"):
        for r in ("0.5", "1.0"):
            assert (d / f"sp_syntracker_leiden_n40_apss{t}_r{r}_clusters.csv").exists()
            assert (d / f"sp_syntracker_leiden_n40_apss{t}_r{r}_cluster_qc.tsv").exists()
    assert n_clusters(clusters(d, "leiden_n40_apss0.75_r1.0")) == 2
    assert n_clusters(clusters(d, "leiden_n40_apss0.95_r1.0")) == 10  # every edge pruned: all singletons


def test_same_seed_same_partition(tmp_path):
    scores = two_groups(within=0.85, between=0.78, bridges=[("a1", "b1"), ("a3", "b4")])
    _, d1 = run(tmp_path, scores, tag="s1")
    _, d2 = run(tmp_path, scores, tag="s2")
    assert clusters(d1, "leiden_n40_apss0.75_r1.0") == clusters(d2, "leiden_n40_apss0.75_r1.0")


def test_largest_cluster_is_numbered_one(tmp_path):
    scores = {k: v for k, v in two_groups().items() if "b5" not in k}
    proc, d = run(tmp_path, scores)
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "leiden_n40_apss0.75_r1.0")
    assert a["a1"] == 1 and a["b1"] == 2


def test_missing_pair_is_a_missing_edge(tmp_path):
    scores = two_groups()
    del scores[("a1", "a2")]  # too few comparable regions: absent from the APSS table
    proc, d = run(tmp_path, scores)
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "leiden_n40_apss0.75_r1.0")
    assert a["a1"] == a["a2"]  # still together via the other pairs


def test_duplicate_pairs_keep_one(tmp_path):
    d = tmp_path / "dup"
    d.mkdir()
    rows = [("R", "x.fasta", "y.fasta", 0.90, 40), ("R", "y.fasta", "x.fasta", 0.10, 40)]
    pd.DataFrame(rows, columns=COLUMNS).to_csv(d / "avg_synteny_scores_40_regions.csv", index=False)
    proc = call(d, [d / "avg_synteny_scores_40_regions.csv"], "--depth", "40")
    assert proc.returncode == 0, proc.stderr
    assert n_clusters(clusters(d, "leiden_n40_apss0.75_r1.0")) == 1  # first row (0.90) kept


def test_target_without_pairs_is_singleton_and_non_target_absent(tmp_path):
    d = tmp_path / "tg"
    d.mkdir()
    names = [f"a{i}" for i in range(1, 6)] + [f"b{i}" for i in range(1, 6)] + ["lonely"]
    (d / "targets.txt").write_text("\n".join(f"{n}.fasta" for n in names) + "\n")
    proc, d = run(tmp_path, two_groups(), extra=["--targets", str(d / "targets.txt")], tag="tg")
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "leiden_n40_apss0.75_r1.0")
    assert list(a.values()).count(a["lonely"]) == 1

    # a genome that is in the APSS table but not a target is dropped, not clustered
    (d / "targets.txt").write_text("\n".join(names[:9]) + "\n")  # b5 not a target
    proc, d = run(tmp_path, two_groups(), extra=["--targets", str(d / "targets.txt")], tag="tg")
    assert proc.returncode == 0, proc.stderr
    assert "b5" not in clusters(d, "leiden_n40_apss0.75_r1.0")
    assert "non-target" in proc.stderr


def test_cluster_qc(tmp_path):
    scores = {k: v for k, v in two_groups(within=0.90, between=0.60).items() if "b5" not in k}
    d = tmp_path / "qc"
    d.mkdir()
    (d / "targets.txt").write_text("\n".join([f"a{i}" for i in range(1, 6)] + ["b1", "b2", "b3", "b4", "b5"]))
    proc, d = run(tmp_path, scores, extra=["--targets", str(d / "targets.txt")], tag="qc")
    assert proc.returncode == 0, proc.stderr
    qc = pd.read_csv(d / "sp_syntracker_leiden_n40_apss0.75_r1.0_cluster_qc.tsv", sep="\t")
    assert list(qc.columns) == ["cluster", "size", "mean_intra_apss", "max_inter_apss",
                                "low_confidence", "reference_genome", "n"]
    a_row = qc[qc["cluster"] == 1].iloc[0]
    assert a_row["size"] == 5 and a_row["mean_intra_apss"] == pytest.approx(0.90)
    assert a_row["max_inter_apss"] == pytest.approx(0.60)
    assert not a_row["low_confidence"]
    lonely = qc[qc["size"] == 1].iloc[0]  # b5: no pairs at all
    assert lonely["low_confidence"] and pd.isna(lonely["mean_intra_apss"]) and pd.isna(lonely["max_inter_apss"])
    assert set(qc["reference_genome"]) == {"REF"} and set(qc["n"]) == {40}


# --- depth selection -----------------------------------------------------------------------------

def write_depths(d, retained):
    """retained: {n: list of pairs kept at that depth} from two_groups()."""
    scores = two_groups()
    paths = []
    for n, keep in retained.items():
        path = d / f"avg_synteny_scores_{n}_regions.csv"
        write_apss(path, {k: scores[k] for k in keep}, regions=n)
        paths.append(path)
    return paths


def test_auto_depth_stops_before_the_cliff(tmp_path):
    d = tmp_path / "auto"
    d.mkdir()
    pairs = list(two_groups())
    # 45 pairs at 40 and 60, 43 at 80 (> 0.9 x 45), 20 at 100 (cliff), 45 at 200 (after the cliff)
    paths = write_depths(d, {40: pairs, 60: pairs, 80: pairs[:43], 100: pairs[:20], 200: pairs})
    proc = call(d, paths, "--depth", "auto")
    assert proc.returncode == 0, proc.stderr
    assert (d / "sp_syntracker_leiden_n80_apss0.75_r1.0_clusters.csv").exists()
    ret = pd.read_csv(d / "sp_syntracker_depth_retention.tsv", sep="\t")
    assert list(ret.columns) == ["n", "retained_samples", "retained_pairs", "total_targets", "selected"]
    assert ret.loc[ret["selected"], "n"].tolist() == [80]
    assert ret.set_index("n").loc[100, "retained_pairs"] == 20
    assert (d / "sp_syntracker_depth_retention.png").stat().st_size > 0


def test_auto_depth_tolerates_an_empty_deep_table(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    pairs = list(two_groups())
    paths = write_depths(d, {40: pairs, 60: pairs, 200: []})
    proc = call(d, paths, "--depth", "auto")
    assert proc.returncode == 0, proc.stderr
    assert (d / "sp_syntracker_leiden_n60_apss0.75_r1.0_clusters.csv").exists()


def test_auto_depth_ignores_the_all_regions_table(tmp_path):
    d = tmp_path / "all"
    d.mkdir()
    pairs = list(two_groups())
    paths = write_depths(d, {40: pairs})
    write_apss(d / "avg_synteny_scores_all_regions.csv", two_groups(), regions=7)
    proc = call(d, [*paths, d / "avg_synteny_scores_all_regions.csv"], "--depth", "auto")
    assert proc.returncode == 0, proc.stderr
    assert (d / "sp_syntracker_leiden_n40_apss0.75_r1.0_clusters.csv").exists()


def test_fixed_all_depth_has_no_regions_assertion(tmp_path):
    proc, d = run(tmp_path, two_groups(), depth="all")
    assert proc.returncode == 0, proc.stderr
    assert (d / "sp_syntracker_leiden_nall_apss0.75_r1.0_clusters.csv").exists()
    assert not (d / "sp_syntracker_depth_retention.tsv").exists()


# --- input contract ------------------------------------------------------------------------------

def test_compared_regions_below_n_is_an_error(tmp_path):
    d = tmp_path / "short"
    d.mkdir()
    write_apss(d / "avg_synteny_scores_60_regions.csv", two_groups(), regions=59)
    proc = call(d, [d / "avg_synteny_scores_60_regions.csv"], "--depth", "60")
    assert proc.returncode != 0
    assert "Compared_regions < 60" in proc.stderr


def test_more_than_one_reference_is_an_error(tmp_path):
    d = tmp_path / "refs"
    d.mkdir()
    rows = [("R1", "x.fasta", "y.fasta", 0.9, 40), ("R2", "x.fasta", "z.fasta", 0.9, 40)]
    pd.DataFrame(rows, columns=COLUMNS).to_csv(d / "avg_synteny_scores_40_regions.csv", index=False)
    proc = call(d, [d / "avg_synteny_scores_40_regions.csv"], "--depth", "40")
    assert proc.returncode != 0
    assert "reference genomes" in proc.stderr


@pytest.mark.parametrize("content, message", [
    ("Ref_genome,Sample1,Sample2,Compared_regions\nR,x,y,50\n", "lacks columns"),
    ("Ref_genome,Sample1,Sample2,APSS,Compared_regions\n", "no pairs"),
])
def test_bad_apss_fails_clearly(tmp_path, content, message):
    (tmp_path / "avg_synteny_scores_40_regions.csv").write_text(content)
    proc = call(tmp_path, [tmp_path / "avg_synteny_scores_40_regions.csv"], "--depth", "40")
    assert proc.returncode != 0
    assert message in proc.stderr


def test_missing_requested_depth_fails_clearly(tmp_path):
    proc, d = run(tmp_path, two_groups(), depth="40", extra=["--depth", "100"])
    assert proc.returncode != 0
    assert "no avg_synteny_scores_100_regions.csv" in proc.stderr
