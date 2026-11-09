"""
Tests for bin/spp_eligibility_from_qc.py (per-genome QC and species eligibility).

Pins: HQ/MQ/DISCARDED from completeness/contamination; with seqkit stats, an HQ/MQ genome whose
N50 (default) or N90 (--contiguity_metric N90) is below --min_contiguity is DISCARDED (reason
`contiguity`), the threshold itself passing; counts and the species label use only the genomes that
pass both; without stats there is no contiguity filter; a genome missing from the stats and a
stats table without N90 fail clearly.

Run with:  python3 -m pytest tests/bin/test_spp_eligibility_from_qc.py
"""

import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "spp_eligibility_from_qc.py"

# genome: (completeness, contamination, N50, N90)
GENOMES = {
    "hq_ok": (99.0, 0.5, 250_000, 40_000),
    "hq_edge": (99.0, 0.5, 5_000, 5_000),       # exactly at the threshold: passes both
    "hq_n50only": (99.0, 0.5, 30_000, 2_000),   # passes N50, fails N90
    "hq_frag": (99.0, 0.5, 3_000, 900),         # fails both
    "mq_ok": (85.0, 3.0, 80_000, 9_000),
    "mq_frag": (85.0, 3.0, 4_000, 1_000),
    "low": (50.0, 0.5, 900_000, 100_000),       # fails quality, whatever its contiguity
}


def write_inputs(d, genomes=GENOMES, stats_columns=("N50", "N90")):
    with open(d / "qc.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genome", "completeness", "contamination"])
        for g, (comp, cont, _, _) in genomes.items():
            w.writerow([f"{g}.fna.gz", comp, cont])
    with open(d / "stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        cols = {"N50": 2, "N90": 3}
        w.writerow(["file", "format", "num_seqs", *stats_columns])
        for g, v in genomes.items():
            w.writerow([f"genomes/{g}.fna.gz", "FASTA", 10, *(v[cols[c]] for c in stats_columns)])


def run(d, *extra, stats=True):
    cmd = [sys.executable, str(SCRIPT), "--species_name", "sp", "--qc_csv", str(d / "qc.csv"),
           "--output", str(d / "report.tsv"), "--rfile_output", str(d / "rfile.tsv"),
           "--labels_output", str(d / "labels.csv"), "--genome_qc_output", str(d / "genome_qc.tsv"), *extra]
    if stats:
        cmd += ["--stats", str(d / "stats.tsv")]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    out = {}
    if proc.returncode == 0:
        with open(d / "report.tsv") as fh:
            out["report"] = next(csv.DictReader(fh, delimiter="\t"))
        with open(d / "genome_qc.tsv") as fh:
            out["genome_qc"] = {r["genome"]: r for r in csv.DictReader(fh, delimiter="\t")}
        out["rfile"] = [line.split("\t")[0] for line in (d / "rfile.tsv").read_text().splitlines()]
        with open(d / "labels.csv") as fh:
            out["labels"] = {r["genome"]: r["label"] for r in csv.DictReader(fh)}
    return proc, out


def test_without_stats_only_quality_filters(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path, stats=False)
    assert proc.returncode == 0, proc.stderr
    assert sorted(out["rfile"]) == sorted(g for g in GENOMES if g != "low")
    assert out["report"]["n_hq"] == "4" and out["report"]["n_mq"] == "2"
    assert out["report"]["n_discarded"] == "1" and out["report"]["n_discarded_contiguity"] == "0"
    assert out["genome_qc"]["hq_frag"]["reason"] == "pass" and out["genome_qc"]["hq_frag"]["N50"] == ""


def test_n50_gate_is_the_default(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert sorted(out["rfile"]) == ["hq_edge", "hq_n50only", "hq_ok", "mq_ok"]
    assert out["labels"] == {"hq_ok": "HQ", "hq_edge": "HQ", "hq_n50only": "HQ", "mq_ok": "MQ"}
    rep = out["report"]
    assert (rep["n_hq"], rep["n_mq"], rep["n_total_passing_qc"]) == ("3", "1", "4")
    assert (rep["n_discarded"], rep["n_discarded_contiguity"]) == ("3", "2")
    qc = out["genome_qc"]
    assert (qc["hq_frag"]["quality"], qc["hq_frag"]["label"], qc["hq_frag"]["reason"]) == ("HQ", "DISCARDED", "contiguity")
    assert (qc["low"]["label"], qc["low"]["reason"]) == ("DISCARDED", "quality")
    assert qc["hq_edge"]["reason"] == "pass"
    assert "2 contiguity: N50 < 5000 bp" in proc.stderr


def test_n90_gate_is_stricter(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path, "--contiguity_metric", "N90")
    assert proc.returncode == 0, proc.stderr
    assert sorted(out["rfile"]) == ["hq_edge", "hq_ok", "mq_ok"]
    assert out["genome_qc"]["hq_n50only"]["reason"] == "contiguity"


def test_min_contiguity_is_a_parameter(tmp_path):
    write_inputs(tmp_path)
    proc, out = run(tmp_path, "--min_contiguity", "0")
    assert proc.returncode == 0, proc.stderr
    assert len(out["rfile"]) == 6


def test_label_uses_post_filter_counts(tmp_path):
    # 110 HQ genomes would be STRONG; 70 of them are fragmented, leaving 40 HQ -> not eligible
    genomes = {f"g{i}": (99.0, 0.5, 100_000 if i < 40 else 1_000, 10_000) for i in range(110)}
    write_inputs(tmp_path, genomes)
    proc, out = run(tmp_path, stats=False)
    assert proc.returncode == 0, proc.stderr
    assert out["report"]["spp_label"].startswith("STRONG")
    proc, out = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert out["report"]["n_hq"] == "40"
    assert out["report"]["spp_label"].startswith("DISCARDED")


def test_genome_missing_from_stats_is_an_error(tmp_path):
    write_inputs(tmp_path)
    lines = (tmp_path / "stats.tsv").read_text().splitlines()
    (tmp_path / "stats.tsv").write_text("\n".join(line for line in lines if "mq_ok" not in line) + "\n")
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "no row in the seqkit stats" in proc.stderr and "mq_ok" in proc.stderr


def test_stats_without_n90_fails_clearly(tmp_path):
    write_inputs(tmp_path, stats_columns=("N50",))
    proc, _ = run(tmp_path)
    assert proc.returncode != 0
    assert "N90" in proc.stderr and "-N 90" in proc.stderr
