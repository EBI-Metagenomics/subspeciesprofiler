"""
Tests for the iTOL colour file of bin/ordination_enrichment.py.

Pins: one TAB-separated `range` line per genome, in the ordination plot's cluster colours
(singletons get no line, so iTOL leaves them uncoloured) with a `<cluster> (n=<size>)` legend label; node IDs use the same character
replacement as bin/fastani_to_tree.py ([^A-Za-z0-9_.-] -> '_'); with --tree, genomes that are not
leaves are left out and reported, leaves without a cluster are reported, and no match at all is
an error; Newick leaf parsing skips internal labels, branch lengths and comments.
With --metadata: one plot and one iTOL DATASET_COLORSTRIP file per column, with every genome,
missing values as grey 'NA', the same node IDs and tree check as the cluster file, and a
duplicated metadata genome as an error.

Needs numpy, pandas, scipy and matplotlib (the poppunk image); skipped without scipy.
Run with:  python3 -m pytest tests/bin/test_ordination_enrichment.py
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "ordination_enrichment.py"

_spec = importlib.util.spec_from_file_location("ordination_enrichment", SCRIPT)
oe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(oe)

# genome -> cluster; 'GCF_7+8' needs the Newick-safe rename; 'solo' is a singleton
CLUSTERS = {"g1": "A", "g2": "A", "g3": "A", "g4": "B", "g5": "B", "g6": "B", "GCF_7+8": "B", "solo": "C"}


def write_inputs(d, clusters=CLUSTERS):
    genomes = list(clusters)
    # two gene blocks, so the ordination has at least two axes
    genes = {f"gene{i}": [int((j + i) % 3 == 0 or (clusters[g] == "A") == (i % 2 == 0))
                          for j, g in enumerate(genomes)] for i in range(12)}
    with open(d / "gpa.Rtab", "w") as fh:
        fh.write("Gene\t" + "\t".join(genomes) + "\n")
        for gene, row in genes.items():
            fh.write(gene + "\t" + "\t".join(map(str, row)) + "\n")
    with open(d / "clusters.csv", "w") as fh:
        fh.write("Taxon,Cluster\n" + "".join(f"{g},{c}\n" for g, c in clusters.items()))
    (d / "anno.tsv").write_text("gene\tannotation\n" + "".join(f"gene{i}\thypothetical\n" for i in range(12)))


def run(d, *extra):
    return subprocess.run([sys.executable, str(SCRIPT), "--rtab", str(d / "gpa.Rtab"),
                           "--clusters", str(d / "clusters.csv"), "--annotations", str(d / "anno.tsv"),
                           "--outdir", str(d), "--prefix", "sp", "--anosim-perm", "9", *extra],
                          capture_output=True, text=True)


def itol_rows(d):
    lines = (d / "sp_itol_colors.txt").read_text().splitlines()
    assert lines[0] == "TREE_COLORS" and "SEPARATOR TAB" in lines and "DATA" in lines
    data = lines[lines.index("DATA") + 1:]
    return {f[0]: f[1:] for f in (line.split("\t") for line in data if not line.startswith("#"))}


def test_newick_leaf_names_skip_internal_labels_and_lengths(tmp_path):
    (tmp_path / "t.nwk").write_text("((g1:0.1,'g 2':0.2)Inner1:0.0[&c],g3:0.3)Inner2:0.0;\n")
    assert oe.newick_leaf_names(tmp_path / "t.nwk") == ["g1", "g 2", "g3"]


def test_one_range_per_genome_with_plot_colours(tmp_path):
    write_inputs(tmp_path)
    proc = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    rows = itol_rows(tmp_path)
    assert set(rows) == {"g1", "g2", "g3", "g4", "g5", "g6", "GCF_7_8"}  # the singleton 'solo' is left out
    assert all(len(r) == 3 and r[0] == "range" for r in rows.values())
    assert rows["g1"][2] == "A (n=3)" and rows["g4"][2] == "B (n=4)"
    # same colour within a cluster, different between clusters, and the plot's tab10 colours
    assert rows["g1"][1] == rows["g3"][1] != rows["g4"][1] == rows["GCF_7_8"][1]
    assert rows["g1"][1] == oe.to_hex(oe.plt.get_cmap("tab10")(0))
    assert "names not checked" in proc.stdout


def test_tree_names_match_after_newick_safe_rename(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "t.nwk").write_text("(((g1:1,g2:1):1,g3:1):1,((g4:1,g5:1):1,(g6:1,GCF_7_8:1):1):1,solo:1);\n")
    proc = run(tmp_path, "--tree", str(tmp_path / "t.nwk"))
    assert proc.returncode == 0, proc.stderr
    assert len(itol_rows(tmp_path)) == 7
    assert "not leaves of the tree" not in proc.stderr and "have no cluster" not in proc.stderr
    assert "7 genome[s] coloured, 0 not in the tree" in proc.stdout


def test_mismatches_are_reported_and_left_out(tmp_path):
    write_inputs(tmp_path)
    # g6 is missing from the tree; 'extra' is a leaf without a cluster
    (tmp_path / "t.nwk").write_text("((g1:1,g2:1,g3:1):1,(g4:1,g5:1,GCF_7_8:1):1,solo:1,extra:1);\n")
    proc = run(tmp_path, "--tree", str(tmp_path / "t.nwk"))
    assert proc.returncode == 0, proc.stderr
    assert "g6" not in itol_rows(tmp_path) and len(itol_rows(tmp_path)) == 6
    assert "1 genome(s) are not leaves of the tree" in proc.stderr and "'g6'" in proc.stderr
    assert "1 tree leaves have no cluster" in proc.stderr and "'extra'" in proc.stderr


def test_no_genome_in_the_tree_is_an_error(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "t.nwk").write_text("(x1:1,x2:1,x3:1);\n")
    proc = run(tmp_path, "--tree", str(tmp_path / "t.nwk"))
    assert proc.returncode == 1
    assert "no genome of the clustering is a leaf of the tree" in proc.stderr


METADATA = (
    "genome\tdisease_status\tgenome_type\n"
    "g1\tCRC\tMAG\ng2\tCRC\tisolate\ng3\tcontrol\tMAG\ng4\tadenoma\tMAG\n"
    "g5\tunknown\tisolate\ng6\t\tMAG\nGCF_7+8\tCRC\tMAG\n"  # g6: empty cell; 'solo' not listed
    "not_analysed\tcontrol\tMAG\n"
)


def strip_file(d, column):
    """Return the header lines and the data rows {node: (colour, category)} of a ring file."""
    lines = (d / f"sp_itol_{column}.txt").read_text().splitlines()
    data = lines[lines.index("DATA") + 1:]
    header = {f[0]: f[1:] for f in (line.split("\t") for line in lines[:lines.index("DATA")])}
    rows = {f[0]: tuple(f[1:]) for f in (line.split("\t") for line in data if not line.startswith("#"))}
    return lines, header, rows


def test_metadata_gives_one_plot_and_one_ring_per_column(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "meta.tsv").write_text(METADATA)
    proc = run(tmp_path, "--metadata", str(tmp_path / "meta.tsv"))
    assert proc.returncode == 0, proc.stderr
    for column in ("disease_status", "genome_type"):
        assert (tmp_path / f"sp_ordination_{column}.png").exists()
        assert (tmp_path / f"sp_ordination_{column}.pdf").exists()
        lines, header, rows = strip_file(tmp_path, column)
        assert lines[0] == "DATASET_COLORSTRIP" and "SEPARATOR TAB" in lines
        assert header["DATASET_LABEL"] == [column] and header["LEGEND_TITLE"] == [column]
        assert len(rows) == 8  # every genome, singleton included
        # the legend lists each category once, in the colour of its rows
        legend = dict(zip(header["LEGEND_LABELS"], header["LEGEND_COLORS"]))
        assert len(header["LEGEND_SHAPES"]) == len(legend)
        assert all(legend[category] == colour for colour, category in rows.values())
    _, _, rows = strip_file(tmp_path, "disease_status")
    assert rows["g1"] == rows["g2"] != rows["g3"]
    assert rows["GCF_7_8"][1] == "CRC"


def test_missing_metadata_is_na_in_grey(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "meta.tsv").write_text(METADATA)
    proc = run(tmp_path, "--metadata", str(tmp_path / "meta.tsv"))
    assert proc.returncode == 0, proc.stderr
    _, _, rows = strip_file(tmp_path, "disease_status")
    assert rows["g6"] == (oe.GREY, "NA")  # empty cell
    assert rows["solo"] == (oe.GREY, "NA")  # not in the table
    assert rows["g5"] == (oe.GREY, "unknown")
    assert "1 genome(s) are not in the metadata table" in proc.stderr and "'solo'" in proc.stderr
    assert "1 metadata row(s) are not analysed genomes" in proc.stderr
    coords = (tmp_path / "sp_ordination_coords.csv").read_text().splitlines()
    assert coords[0].endswith(",disease_status,genome_type")


def test_metadata_rings_use_tree_node_ids_and_skip_non_leaves(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "meta.tsv").write_text(METADATA)
    # g6 is not a leaf of the tree
    (tmp_path / "t.nwk").write_text("((g1:1,g2:1,g3:1):1,(g4:1,g5:1,GCF_7_8:1):1,solo:1);\n")
    proc = run(tmp_path, "--metadata", str(tmp_path / "meta.tsv"), "--tree", str(tmp_path / "t.nwk"))
    assert proc.returncode == 0, proc.stderr
    _, _, rows = strip_file(tmp_path, "genome_type")
    assert set(rows) == {"g1", "g2", "g3", "g4", "g5", "GCF_7_8", "solo"}
    assert proc.stderr.count("not leaves of the tree") == 1  # reported once, not once per file


def test_duplicated_metadata_genome_is_an_error(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "meta.tsv").write_text("genome,disease_status\ng1,CRC\ng1,control\n")
    proc = run(tmp_path, "--metadata", str(tmp_path / "meta.tsv"))
    assert proc.returncode == 1
    assert "lists genomes more than once" in proc.stderr


def test_no_metadata_files_without_metadata(tmp_path):
    write_inputs(tmp_path)
    proc = run(tmp_path)
    assert proc.returncode == 0, proc.stderr
    names = {p.name for p in tmp_path.iterdir()}
    assert not any(n.startswith("sp_ordination_") and not n.startswith("sp_ordination.")
                   and n != "sp_ordination_coords.csv" for n in names)
    assert names & {"sp_itol_colors.txt"} and not any(n.startswith("sp_itol_") and n != "sp_itol_colors.txt"
                                                      for n in names)


def test_category_colours_grey_for_missing_and_unknown():
    colours = oe.category_colours(["CRC", "NA", "control", "Unknown", "CRC"])
    assert list(colours) == ["CRC", "NA", "Unknown", "control"]
    assert colours["NA"] == colours["Unknown"] == oe.GREY
    assert colours["CRC"] != colours["control"] and oe.GREY not in (colours["CRC"], colours["control"])


def test_node_id_collision_is_an_error(tmp_path):
    write_inputs(tmp_path, {**CLUSTERS, "GCF_7_8": "A"})  # 'GCF_7+8' and 'GCF_7_8' -> same ID
    proc = run(tmp_path)
    assert proc.returncode == 1
    assert "collide on the same iTOL node ID" in proc.stderr
