"""
Tests for bin/syntracker_gate.py (contiguity gate for SynTracker targets).

Pins: a genome passes when its N50 (default) or N90 (--metric N90) is at least --min-length bp,
the boundary included; HQ and MQ genomes are gated alike; the passed table keeps the seqkit
columns and adds `label`; missing columns and unlabelled genomes fail clearly.

Run with:  python3 -m pytest tests/bin/test_syntracker_gate.py
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "syntracker_gate.py"

# genome: (label, N50, N90)
GENOMES = {
    "solid": ("HQ", 250_000, 40_000),
    "edge": ("HQ", 5_000, 5_000),       # exactly at the threshold: passes both
    "n50only": ("HQ", 30_000, 2_000),   # passes N50, fails N90
    "frag": ("HQ", 3_000, 900),         # fails both
    "mqok": ("MQ", 80_000, 9_000),
    "mqfrag": ("MQ", 4_000, 1_000),
}


def write_inputs(d, genomes=GENOMES, columns=("N50", "N90")):
    stats = pd.DataFrame({"file": [f"genomes/{g}.fna.gz" for g in genomes], "format": "FASTA",
                          "num_seqs": 10, "N50": [v[1] for v in genomes.values()],
                          "N50_num": 3, "N90": [v[2] for v in genomes.values()]})
    stats[["file", "format", "num_seqs", *columns]].to_csv(d / "stats.tsv", sep="\t", index=False)
    pd.DataFrame({"genome": list(genomes), "label": [v[0] for v in genomes.values()]}) \
        .to_csv(d / "labels.csv", index=False)


def run(d, *extra):
    proc = subprocess.run([sys.executable, str(SCRIPT), "--stats", str(d / "stats.tsv"),
                           "--labels", str(d / "labels.csv"), "--out-prefix", str(d / "sp"), *extra],
                          capture_output=True, text=True)
    out = {}
    if proc.returncode == 0:
        out["report"] = pd.read_csv(d / "sp_syntracker_gate.tsv", sep="\t").set_index("genome")
        out["passed"] = pd.read_csv(d / "sp_syntracker_gate_stats.tsv", sep="\t")
    return proc, out


def passed(out):
    return sorted(out["report"].index[out["report"]["status"] == "pass"])


def test_n50_gate_is_the_default(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert passed(out) == ["edge", "mqok", "n50only", "solid"]
    assert set(out["report"]["metric"]) == {"N50"} and set(out["report"]["min_length"]) == {5000}


def test_n90_gate_is_stricter(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path, "--metric", "N90")
    assert proc.returncode == 0, proc.stderr
    assert passed(out) == ["edge", "mqok", "solid"]
    assert out["report"].loc["n50only", "status"] == "fail"


def test_min_length_is_a_parameter(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path, "--min-length", "50000")
    assert proc.returncode == 0, proc.stderr
    assert passed(out) == ["mqok", "solid"]


def test_passed_table_keeps_seqkit_columns_and_adds_label(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    table = out["passed"]
    assert list(table.columns) == ["file", "format", "num_seqs", "N50", "N90", "label"]
    assert sorted(table["file"]) == ["edge.fna.gz", "mqok.fna.gz", "n50only.fna.gz", "solid.fna.gz"]
    assert dict(zip(table["file"], table["label"]))["mqok.fna.gz"] == "MQ"
    assert "4 of 6 genomes pass" in proc.stderr and "MQ 1/2" in proc.stderr


def test_missing_n90_column_fails_clearly(tmp_path):
    write_inputs(tmp_path, columns=("N50",))
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "N90" in proc.stderr and "-N 90" in proc.stderr


def test_unlabelled_genome_is_an_error(tmp_path):
    write_inputs(tmp_path)
    pd.DataFrame({"genome": ["solid"], "label": ["HQ"]}).to_csv(tmp_path / "labels.csv", index=False)
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "no label" in proc.stderr
