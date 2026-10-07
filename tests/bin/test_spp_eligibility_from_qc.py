"""
Tests for the dRep --genomeInfo output of bin/spp_eligibility_from_qc.py.

Run with:  python3 -m pytest tests/bin/test_spp_eligibility_from_qc.py
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "spp_eligibility_from_qc.py"


def test_genomeinfo_lists_hq_genomes_with_file_names(tmp_path):
    # CheckM2-style headers, matched case-insensitively; one HQ, one MQ, one discarded genome
    (tmp_path / "qc.csv").write_text(
        "Name,Completeness_CheckM2,Contamination_CheckM2\n"
        "g1.fna.gz,99.5,0.5\n"
        "g2.fa,85.0,3.0\n"
        "g3.fasta,50.0,10.0\n"
    )
    subprocess.run(
        [sys.executable, str(SCRIPT), "--species_name", "sp", "--qc_csv", str(tmp_path / "qc.csv"),
         "--output", str(tmp_path / "report.tsv"), "--genomeinfo_output", str(tmp_path / "genomeinfo.csv")],
        check=True, capture_output=True,
    )
    info = pd.read_csv(tmp_path / "genomeinfo.csv")
    assert list(info.columns) == ["genome", "completeness", "contamination"]
    assert info.to_dict("records") == [{"genome": "g1.fna.gz", "completeness": 99.5, "contamination": 0.5}]


def test_genomeinfo_is_optional(tmp_path):
    (tmp_path / "qc.csv").write_text("genome,completeness,contamination\ng1.fna,99.5,0.5\n")
    subprocess.run(
        [sys.executable, str(SCRIPT), "--species_name", "sp", "--qc_csv", str(tmp_path / "qc.csv"),
         "--output", str(tmp_path / "report.tsv")],
        check=True, capture_output=True,
    )
    assert not list(tmp_path.glob("*genomeinfo*"))
