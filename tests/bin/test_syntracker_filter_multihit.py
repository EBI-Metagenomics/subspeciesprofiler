"""
Tests for bin/syntracker_filter_multihit.py (SynTracker's one-hit-per-genome rule, applied after the run).

Pins: a genome with more than one BLAST hit in a region loses that region from all of its pairs
and only those; APSS is the mean of the kept per-region scores, independent of row order; a pair
left without regions is dropped and counted; two runs that differ only in a third genome (the
Codon P. vulgatus case) give identical values; without multi-hits the table is SynTracker's own.

Run with:  python3 -m pytest tests/bin/test_syntracker_filter_multihit.py
"""

import csv
import itertools
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_filter_multihit.py"
HEADER = '"Ref_genome","Sample1","Sample2","Region","Length1","Length2","Overlap","Blocks","Synteny_score"\n'


def write_run(d, scores, hits):
    """scores: {(g1, g2, region): score}; hits: {region: [genome, ...]} (a genome listed twice = two hits)."""
    d.mkdir(parents=True, exist_ok=True)
    genomes = sorted({g for g1, g2, _ in scores for g in (g1, g2)} | {g for gs in hits.values() for g in gs})
    sample = {g: f"Sample.{i}" for i, g in enumerate(genomes, start=1)}
    (d / "dict.tab").write_text("".join(f"{g}\t{s}\n" for g, s in sample.items()))
    blast = d / "blast_output"
    blast.mkdir(exist_ok=True)
    for region, gs in hits.items():
        (blast / f"{region}.tab").write_text(
            "".join(f"{sample[g]}_contig.{k}\t1\t1000\tplus\n" for k, g in enumerate(gs, start=1)))
    rows = [f'"REF","{g1}","{g2}","{r}",5000,5000,4900,1,{s}\n' for (g1, g2, r), s in scores.items()]
    (d / "per_region.csv").write_text(HEADER + "".join(rows))
    return d


def run(d):
    proc = subprocess.run([sys.executable, str(SCRIPT), "--per-region", "per_region.csv", "--blast-dir", "blast_output",
                           "--dictionary", "dict.tab", "--out-prefix", "job"], capture_output=True, text=True, cwd=d)
    apss, report = {}, {}
    if proc.returncode == 0:
        for row in csv.DictReader(open(d / "avg_synteny_scores_all_regions.csv")):
            apss[(row["Sample1"], row["Sample2"])] = (float(row["APSS"]), int(row["Compared_regions"]), row["APSS"])
        lines = (d / "job_syntracker_multihit.tsv").read_text().splitlines()
        report = dict(zip(lines[0].split("\t"), map(int, lines[1].split("\t"))))
    return proc, apss, report


def three_genomes(tmp_path, b_hits_r2=2):
    """a, b, c on regions R1-R3; b has `b_hits_r2` hits in R2."""
    scores = {}
    for (g1, g2), (r, s) in itertools.product(itertools.combinations("abc", 2), [("R1", 0.9), ("R2", 0.5), ("R3", 0.8)]):
        scores[(g1, g2, r)] = s
    hits = {"R1": list("abc"), "R2": ["a", "c"] + ["b"] * b_hits_r2, "R3": list("abc")}
    return write_run(tmp_path, scores, hits)


def test_multihit_genome_loses_the_region_from_its_pairs_only(tmp_path):
    proc, apss, report = run(three_genomes(tmp_path))
    assert proc.returncode == 0, proc.stderr
    # b is multi-copy in R2: a/b and b/c lose R2, a/c keeps it
    assert apss[("a", "b")][:2] == (pytest.approx((0.9 + 0.8) / 2), 2)
    assert apss[("b", "c")][:2] == (pytest.approx((0.9 + 0.8) / 2), 2)
    assert apss[("a", "c")][:2] == (pytest.approx((0.9 + 0.5 + 0.8) / 3), 3)
    assert report == {"multihit_region_genomes": 1, "rows_read": 9, "rows_removed": 2,
                      "pairs_in": 3, "pairs_changed": 2, "pairs_lost": 0}


def test_without_multihits_the_table_is_the_plain_mean(tmp_path):
    proc, apss, report = run(three_genomes(tmp_path, b_hits_r2=1))
    assert proc.returncode == 0, proc.stderr
    assert all(regions == 3 for _, regions, _ in apss.values())
    assert apss[("a", "b")][0] == pytest.approx((0.9 + 0.5 + 0.8) / 3)
    assert report["rows_removed"] == 0 and report["pairs_changed"] == 0


def test_pair_left_without_regions_is_dropped_and_counted(tmp_path):
    write_run(tmp_path, {("a", "b", "R1"): 0.9, ("a", "c", "R1"): 0.9, ("b", "c", "R1"): 0.9},
              {"R1": ["a", "b", "b", "c"]})
    proc, apss, report = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert set(apss) == {("a", "c")}
    assert report["pairs_lost"] == 2


def test_runs_differing_only_in_a_third_genome_agree(tmp_path):
    # The Codon case: x is multi-copy in region R. SynTracker analysed R (and scored x/y on it) only
    # in the run that also held a single-hit genome z; after the filter both runs agree on x/y.
    with_z = write_run(tmp_path / "with_z",
                       {("x", "y", "R0"): 0.93, ("x", "y", "R"): 0.70, ("x", "z", "R"): 0.6, ("y", "z", "R"): 0.6,
                        ("x", "z", "R0"): 0.9, ("y", "z", "R0"): 0.9},
                       {"R0": list("xyz"), "R": ["x", "x", "y", "z"]})
    without_z = write_run(tmp_path / "without_z", {("x", "y", "R0"): 0.93}, {"R0": list("xy"), "R": ["x", "x", "y"]})
    _, a, _ = run(with_z)
    _, b, _ = run(without_z)
    assert a[("x", "y")] == b[("x", "y")] == (pytest.approx(0.93), 1, a[("x", "y")][2])


def test_apss_does_not_depend_on_row_order(tmp_path):
    scores = {("a", "b", f"R{i}"): v for i, v in enumerate([0.1, 0.7, 0.3, 0.9, 0.2, 0.6, 0.4], start=1)}
    hits = {f"R{i}": ["a", "b"] for i in range(1, 8)}
    _, forward, _ = run(write_run(tmp_path / "fwd", scores, hits))
    _, backward, _ = run(write_run(tmp_path / "rev", dict(reversed(list(scores.items()))), hits))
    assert forward[("a", "b")][2] == backward[("a", "b")][2]  # same text, so the strict merge agrees


def test_blast_sample_missing_from_dictionary_fails(tmp_path):
    d = three_genomes(tmp_path)
    (d / "dict.tab").write_text("a\tSample.1\n")
    proc, _, _ = run(d)
    assert proc.returncode != 0
    assert "not in the names dictionary" in proc.stderr
