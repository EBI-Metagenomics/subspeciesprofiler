"""
Tests for bin/ani_qc.py (the FastANI-based genome filter).

Pins: a missing pair and a pair below the minimum ANI are bad pairs; the rule discards the fewest
genomes with its greedy order; ties go against MQ, then the lower N50, then the first name; clean
input is unchanged; every output table agrees on the kept genomes.

Run with:  python3 -m pytest tests/bin/test_ani_qc.py
"""

import csv
import itertools
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "ani_qc.py"


def run_qc(tmp_path, genomes, missing=(), ani=None, labels=None, n50=None):
    """Write the inputs, run the script and return the discarded genomes in order."""
    labels = labels or {}
    n50 = n50 or {}
    ani = ani or {}
    absent = {frozenset(pair) for pair in missing}
    rows = [f"genomes/{g}.fna.gz\tgenomes/{g}.fna.gz\t100\t100\t100\n" for g in genomes]
    for a, b in itertools.permutations(genomes, 2):
        pair = frozenset((a, b))
        if pair not in absent:
            rows.append(f"genomes/{a}.fna.gz\tgenomes/{b}.fna.gz\t{ani.get(pair, 98.0)}\t90\t100\n")
    (tmp_path / "ani.txt").write_text("".join(rows))
    (tmp_path / "rfile.tsv").write_text("".join(f"{g}\t./genomes/{g}.fna.gz\n" for g in genomes))
    (tmp_path / "labels.csv").write_text(
        "genome,label\n" + "".join(f"{g},{labels.get(g, 'HQ')}\n" for g in genomes)
    )
    header = "genome\tfile\tN50\tquality\tlabel\treason\n"
    qc_rows = "".join(
        f"{g}\t{g}.fna.gz\t{n50.get(g, 100000)}\t{labels.get(g, 'HQ')}\t{labels.get(g, 'HQ')}\tpass\n" for g in genomes
    )
    (tmp_path / "genome_qc.tsv").write_text(header + qc_rows + "old\told.fna.gz\t10\tHQ\tDISCARDED\tcontiguity\n")
    subprocess.run(
        [
            sys.executable, str(SCRIPT), "--species-name", "sp", "--ani", "ani.txt", "--rfile", "rfile.tsv",
            "--labels", "labels.csv", "--genome-qc", "genome_qc.tsv", "--out-prefix", "sp",
        ],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    return [row["genome"] for row in read_tsv(tmp_path / "sp_ani_qc.tsv")]


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_clean_input_is_unchanged(tmp_path):
    assert run_qc(tmp_path, ["a", "b", "c"]) == []
    assert (tmp_path / "sp.qc.ani.txt").read_text() == (tmp_path / "ani.txt").read_text()
    assert (tmp_path / "sp_ani_qc_rfile.tsv").read_text() == (tmp_path / "rfile.tsv").read_text()
    assert (tmp_path / "sp_ani_qc_labels.csv").read_text() == (tmp_path / "labels.csv").read_text()


def test_missing_pair_discards_one_genome_and_updates_all_outputs(tmp_path):
    discarded = run_qc(tmp_path, ["a", "b", "c"], missing=[("b", "c")], n50={"b": 500, "c": 900})
    assert discarded == ["b"]
    assert "b" not in (tmp_path / "sp_ani_qc_rfile.tsv").read_text()
    assert "b," not in (tmp_path / "sp_ani_qc_labels.csv").read_text()
    assert "genomes/b" not in (tmp_path / "sp.qc.ani.txt").read_text()
    qc = {row["genome"]: row for row in read_tsv(tmp_path / "sp_ani_qc_genome_qc.tsv")}
    assert (qc["b"]["label"], qc["b"]["reason"]) == ("DISCARDED", "ani")
    assert (qc["c"]["label"], qc["c"]["reason"]) == ("HQ", "pass")
    assert read_tsv(tmp_path / "sp_ani_qc.tsv")[0]["bad_pairs"] == "c:missing"
    report = read_tsv(tmp_path / "sp_ani_qc_eligibility.tsv")[0]
    assert report["n_total_passing_qc"] == "2"
    assert report["n_discarded"] == "2"
    assert report["n_discarded_contiguity"] == "1"
    assert report["n_discarded_ani"] == "1"


def test_low_ani_is_a_bad_pair_and_the_limit_passes(tmp_path):
    low = {frozenset(("a", "b")): 89.9}
    assert len(run_qc(tmp_path, ["a", "b", "c"], ani=low)) == 1
    limit = {frozenset(("a", "b")): 90.0}
    assert run_qc(tmp_path, ["a", "b", "c"], ani=limit) == []


def test_mean_of_both_directions_is_used(tmp_path):
    run_qc(tmp_path, ["a", "b"])
    text = (tmp_path / "ani.txt").read_text().replace(
        "genomes/a.fna.gz\tgenomes/b.fna.gz\t98.0", "genomes/a.fna.gz\tgenomes/b.fna.gz\t81.0"
    )
    (tmp_path / "ani.txt").write_text(text)
    subprocess.run(
        [
            sys.executable, str(SCRIPT), "--species-name", "sp", "--ani", "ani.txt", "--rfile", "rfile.tsv",
            "--labels", "labels.csv", "--genome-qc", "genome_qc.tsv", "--out-prefix", "sp",
        ],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    assert len(read_tsv(tmp_path / "sp_ani_qc.tsv")) == 1  # (81 + 98) / 2 = 89.5


def test_tie_discards_mq_before_hq(tmp_path):
    discarded = run_qc(tmp_path, ["a", "b"], missing=[("a", "b")], labels={"b": "MQ"}, n50={"a": 1, "b": 999})
    assert discarded == ["b"]


def test_tie_discards_lower_n50(tmp_path):
    assert run_qc(tmp_path, ["a", "b"], missing=[("a", "b")], n50={"a": 999, "b": 1}) == ["b"]


def test_tie_discards_first_name(tmp_path):
    assert run_qc(tmp_path, ["a", "b"], missing=[("a", "b")]) == ["a"]


def test_b_fragilis_pattern_needs_two_genomes(tmp_path):
    genomes = ["g106561", "g244638", "g825227", "g852776", "g896483", "g902452", "other"]
    missing = [
        ("g106561", "g244638"), ("g106561", "g896483"), ("g244638", "g825227"), ("g244638", "g852776"),
        ("g244638", "g902452"), ("g852776", "g896483"), ("g896483", "g902452"),
    ]
    assert run_qc(tmp_path, genomes, missing=missing) == ["g244638", "g896483"]
