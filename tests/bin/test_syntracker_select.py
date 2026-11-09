"""
Tests for bin/syntracker_select.py (SynTracker reference and target choice).

Pins: the reference is the most central HQ genome (mean FastANI to the other HQ genomes) among
those with N50 above the floor; with none above it, the highest-N50 HQ genome; an MQ genome is never
the reference; targets are every gated HQ and MQ genome, the highest-N50 ones first when capped
(0 = no cap), the reference always kept.

Run with:  python3 -m pytest tests/bin/test_syntracker_select.py
"""

import itertools
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_select.py"


def write_inputs(d, n50, ani, mq=()):
    """n50: {genome: N50}; ani: {(a, b): percent ANI}, written in both directions; mq: the MQ genomes."""
    pd.DataFrame({"file": [f"{g}.fna.gz" for g in n50], "format": "FASTA", "N50": list(n50.values()),
                  "label": ["MQ" if g in mq else "HQ" for g in n50]}) \
        .to_csv(d / "stats.tsv", sep="\t", index=False)
    rows = [(f"/x/{a}.fna.gz", f"/x/{b}.fna.gz", v, 100, 100) for (a, b), v in ani.items()]
    rows += [(f"/x/{b}.fna.gz", f"/x/{a}.fna.gz", v, 100, 100) for (a, b), v in ani.items()]
    rows += [(f"/x/{g}.fna.gz", f"/x/{g}.fna.gz", 100, 100, 100) for g in n50]
    pd.DataFrame(rows).to_csv(d / "ani.txt", sep="\t", header=False, index=False)


def run(d, *extra):
    proc = subprocess.run([sys.executable, str(SCRIPT), "--ani", str(d / "ani.txt"), "--stats", str(d / "stats.tsv"),
                           "--out-prefix", str(d / "sp"), *extra], capture_output=True, text=True)
    out = {}
    if proc.returncode == 0:
        out["reference"] = (d / "sp_syntracker_reference.txt").read_text().strip()
        out["targets"] = (d / "sp_syntracker_targets.txt").read_text().split()
        out["table"] = pd.read_csv(d / "sp_syntracker_selection.tsv", sep="\t")
    return proc, out


def star(center, others, near=99.0, far=97.0):
    """`center` is close to every genome; the others are far from each other."""
    ani = {}
    for a, b in itertools.combinations([center] + others, 2):
        ani[(a, b)] = near if center in (a, b) else far
    return ani


def test_reference_is_most_central_above_the_n50_floor(tmp_path):
    n50 = {"hub": 200_000, "x": 3_000_000, "y": 150_000, "z": 120_000}
    write_inputs(tmp_path, n50, star("hub", ["x", "y", "z"]))
    proc, out = run(tmp_path, "--min-n50", "100000")
    assert proc.returncode == 0, proc.stderr
    assert out["reference"] == "hub.fna.gz"        # most central, although x is more contiguous
    assert out["targets"][0] == "hub.fna.gz"
    assert sorted(out["targets"]) == sorted(f"{g}.fna.gz" for g in n50)


def test_fragmented_central_genome_is_not_the_reference(tmp_path):
    n50 = {"hub": 20_000, "x": 3_000_000, "y": 150_000, "z": 120_000}
    write_inputs(tmp_path, n50, star("hub", ["x", "y", "z"]))
    proc, out = run(tmp_path, "--min-n50", "100000")
    assert proc.returncode == 0, proc.stderr
    assert out["reference"] != "hub.fna.gz"
    table = out["table"].set_index("genome")
    assert table.loc["hub", "role"] == "target"


def test_no_genome_above_floor_falls_back_to_highest_n50(tmp_path):
    n50 = {"hub": 20_000, "x": 60_000, "y": 30_000}
    write_inputs(tmp_path, n50, star("hub", ["x", "y"]))
    proc, out = run(tmp_path, "--min-n50", "100000")
    assert proc.returncode == 0, proc.stderr
    assert out["reference"] == "x.fna.gz"
    assert "highest N50" in proc.stderr


def test_cap_keeps_reference_and_highest_n50(tmp_path):
    n50 = {"hub": 150_000, "a": 900_000, "b": 800_000, "c": 700_000, "d": 10_000}
    write_inputs(tmp_path, n50, star("hub", ["a", "b", "c", "d"]))
    proc, out = run(tmp_path, "--min-n50", "100000", "--max-targets", "3")
    assert proc.returncode == 0, proc.stderr
    assert out["targets"] == ["hub.fna.gz", "a.fna.gz", "b.fna.gz"]
    roles = out["table"].set_index("genome")["role"]
    assert roles["hub"] == "reference" and roles["c"] == "capped" and roles["d"] == "capped"
    assert list(out["table"]["role"]) == ["reference", "target", "target", "capped", "capped"]
    assert "leaving 2 out" in proc.stderr


def test_max_targets_zero_keeps_every_genome(tmp_path):
    n50 = {"hub": 150_000, "a": 900_000, "b": 800_000, "c": 700_000, "d": 10_000}
    write_inputs(tmp_path, n50, star("hub", ["a", "b", "c", "d"]))
    proc, out = run(tmp_path, "--min-n50", "100000", "--max-targets", "0")
    assert proc.returncode == 0, proc.stderr
    # reference first, then by N50; nobody capped and no warning
    assert out["targets"] == ["hub.fna.gz", "a.fna.gz", "b.fna.gz", "c.fna.gz", "d.fna.gz"]
    assert "capped" not in set(out["table"]["role"])
    assert "WARNING" not in proc.stderr


def test_max_targets_one_is_an_error(tmp_path):
    write_inputs(tmp_path, {"a": 500_000, "b": 400_000}, {("a", "b"): 99.0})
    proc, _ = run(tmp_path, "--max-targets", "1")
    assert proc.returncode != 0
    assert "--max-targets" in proc.stderr


def test_fewer_than_two_genomes_is_an_error(tmp_path):
    write_inputs(tmp_path, {"only": 500_000}, {})
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "at least 2" in proc.stderr


def test_stats_without_n50_column_fails_clearly(tmp_path):
    pd.DataFrame({"file": ["a.fna", "b.fna"], "num_seqs": [1, 2]}).to_csv(tmp_path / "stats.tsv", sep="\t", index=False)
    (tmp_path / "ani.txt").write_text("a.fna\tb.fna\t99\t1\t1\n")
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "N50" in proc.stderr


def test_mq_genome_is_a_target_but_never_the_reference(tmp_path):
    # mq is the most central and most contiguous genome, but only HQ genomes can be the reference.
    n50 = {"mq": 3_000_000, "hub": 200_000, "y": 150_000, "z": 120_000}
    write_inputs(tmp_path, n50, star("mq", ["hub", "y", "z"]), mq={"mq"})
    proc, out = run(tmp_path, "--min-n50", "100000")
    assert proc.returncode == 0, proc.stderr
    assert out["reference"] != "mq.fna.gz"
    assert "mq.fna.gz" in out["targets"]
    table = out["table"].set_index("genome")
    assert table.loc["mq", "label"] == "MQ" and table.loc["mq", "role"] == "target"
    assert pd.isna(table.loc["mq", "centrality"])  # centrality is over the HQ genomes only


def test_mq_genome_is_not_the_fallback_reference(tmp_path):
    n50 = {"mq": 3_000_000, "hub": 20_000, "y": 30_000}
    write_inputs(tmp_path, n50, star("hub", ["mq", "y"]), mq={"mq"})
    proc, out = run(tmp_path, "--min-n50", "100000")
    assert proc.returncode == 0, proc.stderr
    assert out["reference"] == "y.fna.gz"


def test_no_hq_genome_is_an_error(tmp_path):
    write_inputs(tmp_path, {"a": 500_000, "b": 400_000}, {("a", "b"): 99.0}, mq={"a", "b"})
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "HQ" in proc.stderr
