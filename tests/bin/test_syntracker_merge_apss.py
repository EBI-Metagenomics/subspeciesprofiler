"""
Tests for bin/syntracker_merge_apss.py (merging SynTracker all-regions APSS from chunked runs).

Pins: tables are concatenated with one row per unordered pair; a pair found in several runs
must agree exactly, otherwise the merge fails; one reference genome; SynTracker's text values
and quoting are kept; target pairs absent from the merge are reported.

Run with:  python3 -m pytest tests/bin/test_syntracker_merge_apss.py
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_merge_apss.py"
HEADER = '"Ref_genome","Sample1","Sample2","APSS","Compared_regions"\n'


def write(path, rows, ref="REF"):
    """rows: [(sample1, sample2, apss_text, regions)]"""
    path.write_text(HEADER + "".join(f'"{ref}","{a}","{b}",{s},{n}\n' for a, b, s, n in rows))
    return path


def run(d, tables, *extra):
    proc = subprocess.run([sys.executable, str(SCRIPT), "--apss", *map(str, tables), "--out-prefix", "sp", *extra],
                          capture_output=True, text=True, cwd=d)
    return proc


def report(d):
    return pd.read_csv(d / "sp_syntracker_merge.tsv", sep="\t").iloc[0]


def test_chunk_tables_merge_into_one_row_per_pair(tmp_path):
    # chunks {a, b} and {c}: runs a+b+c share the pair a/b, written in the other order in the second
    t1 = write(tmp_path / "t1.csv", [("a", "b", "0.861552666869159", 71), ("a", "c", "0.75", 139)])
    t2 = write(tmp_path / "t2.csv", [("b", "a", "0.861552666869159", 71), ("b", "c", "0.7", 200)])
    proc = run(tmp_path, [t1, t2])
    assert proc.returncode == 0, proc.stderr
    text = (tmp_path / "avg_synteny_scores_all_regions.csv").read_text()
    assert text == HEADER + '"REF","a","b",0.861552666869159,71\n"REF","a","c",0.75,139\n"REF","b","c",0.7,200\n'
    r = report(tmp_path)
    assert (r["n_tables"], r["rows_read"], r["duplicate_rows"], r["unique_pairs"]) == (2, 4, 1, 3)


def test_conflicting_duplicate_fails(tmp_path):
    t1 = write(tmp_path / "t1.csv", [("a", "b", "0.80", 100)])
    t2 = write(tmp_path / "t2.csv", [("a", "b", "0.81", 100)])
    proc = run(tmp_path, [t1, t2])
    assert proc.returncode != 0
    assert "not comparable" in proc.stderr


def test_duplicate_with_other_region_count_fails(tmp_path):
    t1 = write(tmp_path / "t1.csv", [("a", "b", "0.80", 100)])
    t2 = write(tmp_path / "t2.csv", [("a", "b", "0.80", 101)])
    assert run(tmp_path, [t1, t2]).returncode != 0


def test_two_references_fail(tmp_path):
    t1 = write(tmp_path / "t1.csv", [("a", "b", "0.80", 100)], ref="R1")
    t2 = write(tmp_path / "t2.csv", [("a", "c", "0.80", 100)], ref="R2")
    proc = run(tmp_path, [t1, t2])
    assert proc.returncode != 0
    assert "reference" in proc.stderr


def test_header_only_tables_fail(tmp_path):
    t1 = write(tmp_path / "t1.csv", [])
    proc = run(tmp_path, [t1])
    assert proc.returncode != 0
    assert "no pairs" in proc.stderr


def test_missing_target_pairs_are_reported(tmp_path):
    t1 = write(tmp_path / "t1.csv", [("a", "b", "0.80", 100), ("a", "c", "0.80", 100)])
    (tmp_path / "targets.txt").write_text("a.fna.gz\nb.fna.gz\nc.fna.gz\n")
    proc = run(tmp_path, [t1], "--targets", "targets.txt")
    assert proc.returncode == 0, proc.stderr
    r = report(tmp_path)
    assert (r["n_targets"], r["expected_pairs"], r["missing_pairs"]) == (3, 3, 1)
    assert "1 of 3 target pairs" in proc.stderr
