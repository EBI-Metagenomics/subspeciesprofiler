"""
Tests for bin/syntracker_apss_clusters.py (SynTracker APSS -> Taxon,Cluster tables).

Pins the behaviour chosen from the real B. longum run: average linkage is the default
because single linkage chains through a few high-APSS "bridging" pairs between groups.

Run with:  python3 -m pytest tests/bin/test_syntracker_apss_clusters.py
"""

import itertools
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_apss_clusters.py"


def write_apss(path, scores):
    """scores: {(a, b): apss}. Written as SynTracker does, with target file names as sample names."""
    rows = [("REF", f"{a}.fasta", f"{b}.fasta", s, 100) for (a, b), s in scores.items()]
    pd.DataFrame(rows, columns=["Ref_genome", "Sample1", "Sample2", "APSS", "Compared_regions"]).to_csv(path, index=False)


def two_groups(within=0.90, between=0.60, bridges=()):
    """Groups a1..a5 and b1..b5; `bridges` lists (a, b) pairs scored like within-group pairs."""
    a = [f"a{i}" for i in range(1, 6)]
    b = [f"b{i}" for i in range(1, 6)]
    scores = {}
    for x, y in itertools.combinations(a + b, 2):
        same = (x[0] == y[0])
        scores[(x, y)] = within if same or (x, y) in bridges or (y, x) in bridges else between
    return scores


def run(tmp_path, scores, thresholds="0.80", extra=(), tag="run"):
    d = tmp_path / tag
    d.mkdir(exist_ok=True)
    write_apss(d / "apss.csv", scores)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--apss", str(d / "apss.csv"), "--thresholds", thresholds,
         "--out-prefix", str(d / "sp"), *extra],
        capture_output=True, text=True,
    )
    return proc, d


def clusters(d, method_tag, threshold):
    df = pd.read_csv(d / f"sp_syntracker_{method_tag}_apss{threshold}_clusters.csv")
    assert list(df.columns) == ["Taxon", "Cluster"]
    return dict(zip(df["Taxon"], df["Cluster"]))


def n_clusters(assignment):
    return len(set(assignment.values()))


def test_bridging_pairs_merge_single_linkage_but_not_average(tmp_path):
    scores = two_groups(bridges=[("a1", "b1"), ("a2", "b2")])
    proc, d = run(tmp_path, scores, extra=["--method", "connected"], tag="cc")
    assert proc.returncode == 0, proc.stderr
    assert n_clusters(clusters(d, "cc", "0.80")) == 1           # chained through the bridges
    proc, d = run(tmp_path, scores, tag="avg")
    assert proc.returncode == 0, proc.stderr
    avg = clusters(d, "avg", "0.80")
    assert n_clusters(avg) == 2                                  # average linkage keeps them apart
    assert len({avg[f"a{i}"] for i in range(1, 6)}) == 1
    assert avg["a1"] != avg["b1"]


def test_threshold_sweep_and_naming(tmp_path):
    proc, d = run(tmp_path, two_groups(), thresholds="0.50,0.80,0.95")
    assert proc.returncode == 0, proc.stderr
    assert n_clusters(clusters(d, "avg", "0.50")) == 1   # everything above 0.50
    assert n_clusters(clusters(d, "avg", "0.80")) == 2   # the two groups
    assert n_clusters(clusters(d, "avg", "0.95")) == 10  # all singletons


def test_largest_cluster_is_numbered_one(tmp_path):
    scores = two_groups()
    scores = {k: v for k, v in scores.items() if not (k[0] == "b5" or k[1] == "b5")}  # b5 absent
    proc, d = run(tmp_path, scores)
    a = clusters(d, "avg", "0.80")
    assert a["a1"] == 1 and a["b1"] == 2


def test_missing_pair_is_ignored_not_dissimilar(tmp_path):
    scores = two_groups()
    del scores[("a1", "a2")]  # too few comparable regions: absent from the APSS table
    proc, d = run(tmp_path, scores)
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "avg", "0.80")
    assert a["a1"] == a["a2"]  # still together via the other pairs


def test_duplicate_pairs_are_averaged(tmp_path):
    d = tmp_path / "dup"
    d.mkdir()
    rows = [("R1", "x.fasta", "y.fasta", 0.90, 50), ("R2", "y.fasta", "x.fasta", 0.70, 50)]  # mean 0.80
    pd.DataFrame(rows, columns=["Ref_genome", "Sample1", "Sample2", "APSS", "Compared_regions"]).to_csv(d / "apss.csv", index=False)
    for t, expected in (("0.79", 1), ("0.81", 2)):
        subprocess.run([sys.executable, str(SCRIPT), "--apss", str(d / "apss.csv"), "--thresholds", t,
                        "--out-prefix", str(d / "sp")], check=True, capture_output=True)
        assert n_clusters(clusters(d, "avg", t)) == expected


def test_labelled_genomes_without_pairs_become_singletons(tmp_path):
    d = tmp_path / "lab"
    d.mkdir()
    pd.DataFrame({"genome": [f"a{i}" for i in range(1, 6)] + [f"b{i}" for i in range(1, 6)] + ["lonely"],
                  "label": "HQ"}).to_csv(d / "labels.csv", index=False)
    proc, d = run(tmp_path, two_groups(), extra=["--labels", str(d / "labels.csv")], tag="lab")
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "avg", "0.80")
    assert "lonely" in a
    assert list(a.values()).count(a["lonely"]) == 1


def test_drep_propagation(tmp_path):
    # SynTracker ran on representatives a1 and b1 only; their dRep clusters carry the label.
    scores = {("a1", "b1"): 0.60}
    cdb = pd.DataFrame({
        "genome": ["a1.fna", "a2.fna", "a3.fna", "b1.fna", "b2.fna", "z1.fna"],
        "secondary_cluster": ["1_1", "1_1", "1_1", "1_2", "1_2", "1_3"],  # 1_3 has no target (dropped by the cap)
    })
    d = tmp_path / "cdb"
    d.mkdir()
    cdb.to_csv(d / "Cdb.csv", index=False)
    proc, d = run(tmp_path, scores, extra=["--drep-cdb", str(d / "Cdb.csv")], tag="cdb")
    assert proc.returncode == 0, proc.stderr
    a = clusters(d, "avg", "0.80")
    assert a["a2"] == a["a3"] == a["a1"]
    assert a["b2"] == a["b1"] != a["a1"]
    assert "z1" not in a  # not represented in SynTracker -> not reported


@pytest.mark.parametrize("content, message", [
    ("Ref_genome,Sample1,Sample2,Compared_regions\nR,x,y,5\n", "lacks columns"),
    ("Ref_genome,Sample1,Sample2,APSS,Compared_regions\n", "no pairs"),
])
def test_bad_apss_fails_clearly(tmp_path, content, message):
    (tmp_path / "apss.csv").write_text(content)
    proc = subprocess.run([sys.executable, str(SCRIPT), "--apss", str(tmp_path / "apss.csv"), "--thresholds", "0.8",
                           "--out-prefix", str(tmp_path / "sp")], capture_output=True, text=True)
    assert proc.returncode != 0
    assert message in proc.stderr
