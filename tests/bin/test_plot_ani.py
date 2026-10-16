"""
Tests for bin/plot_ani.py (per-species ANI histogram + clustered ANI heatmap).

Pins: one comparison per unordered pair (the two FastANI directions averaged, self-hits
dropped); pairs FastANI did not report are counted as missing and stay NaN in the matrix;
UPGMA keeps each ANI group contiguous in the heatmap order; PNG, PDF and summary are written,
also when there is nothing to plot.

Needs matplotlib and scipy (both in the poppunk image the module runs on).

Run with:  python3 -m pytest tests/bin/test_plot_ani.py
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("scipy")

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "plot_ani.py"
sys.path.insert(0, str(REPO / "bin"))

import plot_ani  # noqa: E402
from evaluate_poppunk_fastani import read_fastani  # noqa: E402


def write_fastani(path, genomes, ani_of, missing=()):
    """ani_of(a, b) -> percent ANI; pairs in `missing` (frozensets) are left out, both directions."""
    rows = []
    for a in genomes:
        for b in genomes:
            if frozenset((a, b)) in missing:
                continue
            rows.append((f"genomes/{a}.fna", f"genomes/{b}.fna", 100.0 if a == b else ani_of(a, b), 900, 1000))
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)
    return path


def two_groups(a, b):
    return 99.5 if a[0] == b[0] else 96.0


GENOMES = ["A1", "B1", "A2", "B2", "A3", "B3"]


def test_pairs_are_unordered_and_directions_averaged(tmp_path):
    path = tmp_path / "ani.txt"
    pd.DataFrame([("x.fna", "y.fna", 99.0, 1, 1), ("y.fna", "x.fna", 98.0, 1, 1), ("x.fna", "x.fna", 100.0, 1, 1)]) \
        .to_csv(path, sep="\t", header=False, index=False)
    pairs = plot_ani.pairwise(read_fastani(str(path)))
    assert len(pairs) == 1
    assert pairs.iloc[0][["g1", "g2"]].tolist() == ["x", "y"]
    assert pairs.iloc[0]["ani"] == pytest.approx(98.5)


def test_missing_pair_stays_nan_and_groups_cluster_together(tmp_path):
    path = write_fastani(tmp_path / "ani.txt", GENOMES, two_groups, missing={frozenset(("A1", "B3"))})
    fastani = read_fastani(str(path))
    genomes = sorted(set(fastani["query"]))
    matrix = plot_ani.ani_matrix(plot_ani.pairwise(fastani), genomes)
    assert np.isnan(matrix.at["A1", "B3"]) and np.isnan(matrix.at["B3", "A1"])
    assert matrix.at["A1", "A1"] == 100.0
    _, order = plot_ani.cluster_order(matrix)
    groups = [genomes[i][0] for i in order]
    # each group is one contiguous block in the heatmap order
    assert groups in (["A"] * 3 + ["B"] * 3, ["B"] * 3 + ["A"] * 3)


def test_figure_and_summary_are_written(tmp_path):
    write_fastani(tmp_path / "ani.txt", GENOMES, two_groups, missing={frozenset(("A1", "B3"))})
    (tmp_path / "labels.csv").write_text("genome,label\n" + "".join(f"{g},{'MQ' if g == 'B2' else 'HQ'}\n" for g in GENOMES))
    proc = subprocess.run([sys.executable, str(SCRIPT), "--ani", "ani.txt", "--labels", "labels.csv",
                           "--title", "sp", "--out-prefix", "sp"], capture_output=True, text=True, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    for ext in ("png", "pdf"):
        assert (tmp_path / f"sp_ani_overview.{ext}").stat().st_size > 1000
    summary = pd.read_csv(tmp_path / "sp_ani_overview.tsv", sep="\t").iloc[0]
    assert (summary["n_genomes"], summary["n_pairs"], summary["n_missing_pairs"]) == (6, 14, 1)
    assert summary["ani_min"] == pytest.approx(96.0)
    assert summary["colour_min"] == pytest.approx(96.0)


def test_single_genome_writes_a_placeholder_figure(tmp_path):
    write_fastani(tmp_path / "ani.txt", ["A1"], two_groups)
    proc = subprocess.run([sys.executable, str(SCRIPT), "--ani", "ani.txt", "--out-prefix", "sp"],
                          capture_output=True, text=True, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "sp_ani_overview.png").exists()
    summary = pd.read_csv(tmp_path / "sp_ani_overview.tsv", sep="\t").iloc[0]
    assert (summary["n_genomes"], summary["n_pairs"]) == (1, 0)
