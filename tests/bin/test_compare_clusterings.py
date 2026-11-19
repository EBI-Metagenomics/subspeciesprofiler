"""
Tests for bin/compare_clusterings.py (agreement between two clusterings).

Pins: a method-2 cluster matches the method-1 cluster that shares most of its genomes (tie: the
larger method-1 cluster, then the first name); fraction_2_in_1 counts genomes in method-2 clusters
of at least 2 genomes, and singletons are counted apart; a split keeps fraction_2_in_1 at 1 and a
merge keeps fraction_1_in_2 at 1; pair agreement and the adjusted Rand index; genomes in one table
only are left out and reported; FASTA extensions do not stop names from matching; the genome
placement labels (agreed / not_agreed / singleton); a full run writes all 6 outputs.

Needs numpy, pandas and matplotlib (the poppunk image); skipped without matplotlib.
Run with:  python3 -m pytest tests/bin/test_compare_clusterings.py
"""

import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "compare_clusterings.py"

_spec = importlib.util.spec_from_file_location("compare_clusterings", SCRIPT)
cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cc)


def series(mapping):
    return pd.Series(mapping).sort_index()


def compare(first, second):
    c1, c2 = series(first), series(second)
    return cc.agreement(c1, c2.loc[c1.index])


def test_identical_partitions_with_renamed_labels():
    first = {"a": "1", "b": "1", "c": "2", "d": "2"}
    second = {"a": "x", "b": "x", "c": "y", "d": "y"}
    _, m, placement = compare(first, second)
    assert m["fraction_2_in_1"] == m["fraction_1_in_2"] == 1.0
    assert m["pair_agreement_2_in_1"] == m["pair_agreement_1_in_2"] == 1.0
    assert m["ARI"] == pytest.approx(1.0)
    assert set(placement["agreement"]) == {"agreed"}


def test_split_keeps_fraction_2_in_1():
    # method 2 splits method-1 cluster 1 into x and y
    first = {"a": "1", "b": "1", "c": "1", "d": "1", "e": "2", "f": "2"}
    second = {"a": "x", "b": "x", "c": "y", "d": "y", "e": "z", "f": "z"}
    _, m, _ = compare(first, second)
    assert m["fraction_2_in_1"] == 1.0
    assert m["fraction_1_in_2"] == pytest.approx(4 / 6)  # a, b (x) agree; c, d (y) do not
    assert m["pair_agreement_2_in_1"] == 1.0 and m["pair_agreement_1_in_2"] < 1.0
    assert m["ARI"] < 1.0


def test_merge_keeps_fraction_1_in_2():
    first = {"a": "1", "b": "1", "c": "2", "d": "2", "e": "3", "f": "3"}
    second = {"a": "x", "b": "x", "c": "x", "d": "x", "e": "y", "f": "y"}
    _, m, placement = compare(first, second)
    assert m["fraction_1_in_2"] == 1.0
    assert m["fraction_2_in_1"] < 1.0
    assert m["pair_agreement_1_in_2"] == 1.0 and m["pair_agreement_2_in_1"] < 1.0
    agreed = (placement["agreement"] == "agreed").sum()
    scored = placement["agreement"].isin(["agreed", "not_agreed"]).sum()
    assert agreed / scored == pytest.approx(m["fraction_2_in_1"])


def test_singletons_are_left_out_and_counted():
    first = {"a": "1", "b": "1", "c": "1", "d": "2"}
    second = {"a": "x", "b": "x", "c": "s1", "d": "s2"}
    _, m, placement = compare(first, second)
    assert m["n_singletons_2"] == 2 and m["n_singletons_1"] == 1
    assert m["fraction_2_in_1"] == 1.0  # only a and b (cluster x) are scored
    labels = dict(zip(placement["genome"], placement["agreement"]))
    assert labels == {"a": "agreed", "b": "agreed", "c": "singleton", "d": "singleton"}


def test_not_agreed_label():
    # x: 2 genomes of cluster 1, 1 of cluster 2 -> matches 1; c (cluster 2) is not_agreed
    first = {"a": "1", "b": "1", "c": "2", "d": "2", "e": "2"}
    second = {"a": "x", "b": "x", "c": "x", "d": "y", "e": "y"}
    _, m, placement = compare(first, second)
    labels = dict(zip(placement["genome"], placement["agreement"]))
    assert labels == {"a": "agreed", "b": "agreed", "c": "not_agreed", "d": "agreed", "e": "agreed"}
    assert m["fraction_2_in_1"] == pytest.approx(4 / 5)


def test_tie_goes_to_the_larger_cluster_then_the_first_name():
    contingency = pd.DataFrame({"x": [1, 1, 0]}, index=["1", "2", "3"])
    sizes_from_other_columns = pd.DataFrame({"y": [0, 5, 0]}, index=["1", "2", "3"])
    assert cc.best_matches(pd.concat([contingency, sizes_from_other_columns], axis=1))["x"] == ("2", 1)
    assert cc.best_matches(contingency)["x"] == ("1", 1)


def test_adjusted_rand_index_hand_example():
    # contingency [[2, 1], [0, 3]]: index 1 + 0 + 0 + 3 = 4; rows 3 + 3 = 6; columns 1 + 6 = 7
    contingency = pd.DataFrame([[2, 1], [0, 3]], index=["1", "2"], columns=["x", "y"])
    expected = 6 * 7 / 15
    assert cc.adjusted_rand_index(contingency) == pytest.approx((4 - expected) / (6.5 - expected))


def test_adjusted_rand_index_of_one_cluster_each_is_one():
    assert cc.adjusted_rand_index(pd.DataFrame([[4]], index=["1"], columns=["x"])) == 1.0


def test_natural_order_of_cluster_names():
    assert sorted(["10", "2", "b", "1"], key=cc.natural_key) == ["1", "2", "10", "b"]


def write(path, rows):
    path.write_text("Taxon,Cluster\n" + "".join(f"{g},{c}\n" for g, c in rows))


def run(d, *extra):
    return subprocess.run([sys.executable, str(SCRIPT), "--clusters1", str(d / "one.csv"),
                           "--clusters2", str(d / "two.csv"), "--label1", "PopPUNK", "--label2", "SynTracker",
                           "--outdir", str(d / "out"), "--prefix", "cmp", *extra],
                          capture_output=True, text=True)


def test_full_run_writes_all_outputs_and_reports_unshared_genomes(tmp_path):
    first = [(f"g{i}.fna.gz", str(1 + i // 10)) for i in range(30)] + [("only_in_1", "1")]
    second = [(f"g{i}", "a" if i < 15 else "b") for i in range(30)] + [("only_in_2", "b")]
    write(tmp_path / "one.csv", first)
    write(tmp_path / "two.csv", second)
    proc = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = tmp_path / "out"
    for suffix in ("_sankey.png", "_sankey.pdf", "_agreement_summary.tsv", "_cluster_matches.tsv",
                   "_contingency.tsv", "_genome_placement.tsv"):
        assert (out / f"cmp{suffix}").exists(), suffix
    summary = pd.read_csv(out / "cmp_agreement_summary.tsv", sep="\t").iloc[0]
    assert (summary["n_shared"], summary["n_only_1"], summary["n_only_2"]) == (30, 1, 1)
    assert summary["label_1"] == "PopPUNK" and summary["n_clusters_1"] == 3 and summary["n_clusters_2"] == 2
    assert not math.isnan(summary["ARI"])
    assert "only in PopPUNK" in proc.stderr and "'only_in_1'" in proc.stderr
    placement = pd.read_csv(out / "cmp_genome_placement.tsv", sep="\t", dtype=str)
    assert list(placement.columns) == ["genome", "cluster_1", "cluster_2", "agreement"]
    assert len(placement) == 30 and placement["genome"].str.endswith(".fna.gz").sum() == 0
    contingency = pd.read_csv(out / "cmp_contingency.tsv", sep="\t", dtype={"cluster_1": str})
    assert contingency["n_genomes"].sum() == 30


def test_duplicated_genome_is_an_error(tmp_path):
    write(tmp_path / "one.csv", [("g1", "1"), ("g1", "2"), ("g2", "1")])
    write(tmp_path / "two.csv", [("g1", "a"), ("g2", "a")])
    proc = run(tmp_path)
    assert proc.returncode != 0 and "more than once" in proc.stderr


def test_no_shared_genomes_is_an_error(tmp_path):
    write(tmp_path / "one.csv", [("g1", "1"), ("g2", "1")])
    write(tmp_path / "two.csv", [("h1", "a"), ("h2", "a")])
    proc = run(tmp_path)
    assert proc.returncode != 0 and "fewer than 2 genomes are in both tables" in proc.stderr
