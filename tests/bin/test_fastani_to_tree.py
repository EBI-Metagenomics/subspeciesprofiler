"""
Tests for bin/fastani_to_tree.py (neighbour-joining tree from all-vs-all FastANI).

Pins: tips are the file names minus their FASTA extension (as normalise_genome_id, `.contigs`
included); distance = 1 - ANI/100, averaged over A->B and B->A; two genomes give one cherry;
a missing pair, a label collision and a malformed line fail with a clear message.

Run with:  python3 -m pytest tests/bin/test_fastani_to_tree.py  (skipped without Biopython)
"""

import subprocess
import sys
from pathlib import Path

import pytest

Phylo = pytest.importorskip("Bio.Phylo")

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "fastani_to_tree.py"


def write_ani(d, rows):
    """rows: [(query file, reference file, percent ANI)]."""
    (d / "ani.txt").write_text("".join(f"unzipped/{q}\tunzipped/{r}\t{v}\t90\t100\n" for q, r, v in rows))


def all_pairs(names, ani):
    """Every ordered pair, both directions, with the ANI from `ani(a, b)`, plus the self hits."""
    return [(a, b, 100.0 if a == b else ani(a, b)) for a in names for b in names]


def run(d):
    proc = subprocess.run([sys.executable, str(SCRIPT), str(d / "ani.txt"), str(d / "tree.nwk")],
                          capture_output=True, text=True)
    tree = Phylo.read(str(d / "tree.nwk"), "newick") if proc.returncode == 0 else None
    return proc, tree


def test_tips_are_sample_names(tmp_path):
    files = ["a.fna.gz", "b.fasta", "c.contigs", "d.fa"]
    write_ani(tmp_path, all_pairs(files, lambda a, b: 98.0))
    proc, tree = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert sorted(t.name for t in tree.get_terminals()) == ["a", "b", "c", "d"]


def test_distances_follow_ani(tmp_path):
    # a and b are close (99.9), c is far from both (97): a and b are each other's nearest tips
    ani = {frozenset("ab"): 99.9, frozenset("ac"): 97.0, frozenset("bc"): 97.0}
    names = ["a.fna", "b.fna", "c.fna"]
    write_ani(tmp_path, all_pairs(names, lambda x, y: ani[frozenset((x[0], y[0]))]))
    proc, tree = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert tree.distance("a", "b") == pytest.approx(0.001, abs=1e-9)
    assert tree.distance("a", "c") == pytest.approx(0.03, abs=1e-9)


def test_both_directions_are_averaged(tmp_path):
    write_ani(tmp_path, [("a.fna", "b.fna", 98.0), ("b.fna", "a.fna", 99.0)])
    proc, tree = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert tree.distance("a", "b") == pytest.approx(0.015, abs=1e-9)


def test_missing_pair_is_an_error(tmp_path):
    rows = [r for r in all_pairs(["a.fna", "b.fna", "c.fna"], lambda a, b: 98.0) if {r[0], r[1]} != {"a.fna", "c.fna"}]
    write_ani(tmp_path, rows)
    proc, _ = run(tmp_path)
    assert proc.returncode == 1
    assert "1 genome pairs are missing" in proc.stderr and "a <-> c" in proc.stderr


def test_label_collision_is_an_error(tmp_path):
    write_ani(tmp_path, all_pairs(["a.fna", "a.fasta"], lambda a, b: 98.0))
    proc, _ = run(tmp_path)
    assert proc.returncode == 1
    assert "collision" in proc.stderr


def test_malformed_line_is_an_error(tmp_path):
    (tmp_path / "ani.txt").write_text("a.fna\tb.fna\t98.0\n")
    proc, _ = run(tmp_path)
    assert proc.returncode == 1
    assert "expected five columns" in proc.stderr
