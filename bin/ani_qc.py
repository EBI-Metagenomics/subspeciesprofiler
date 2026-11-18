#!/usr/bin/env python3
"""
ANI QC for one species
======================
Makes the final genome list of a species from its all-vs-all FastANI. Every branch of the
pipeline uses this list.

A bad pair is a genome pair with ANI below --min-ani, or with no FastANI value in either
direction. FastANI omits a pair when the ANI is too low or the two genomes share too little
sequence. The ANI of a pair is the mean of its two directions, as in the evaluator.

The script discards the fewest genomes it can with a greedy rule. It discards the genome with
the most bad pairs and repeats until no bad pair remains. On a tie, it discards MQ before HQ,
then the lower N50, then the first name in sort order.

Outputs (`<p>` = --out-prefix):
  <p>_ani_qc_rfile.tsv        r-file of the kept genomes
  <p>_ani_qc_labels.csv       genome,label of the kept genomes
  <p>_ani_qc_genome_qc.tsv    the genome_qc table; discarded genomes get label DISCARDED, reason ani
  <p>.qc.ani.txt              FastANI rows of the kept genomes
  <p>_ani_qc.tsv              the discarded genomes, in discard order, with their bad pairs
  <p>_ani_qc_eligibility.tsv  species label and genome counts after this step

Usage:
    ani_qc.py --species-name B_fragilis --ani B_fragilis.ani.txt --rfile B_fragilis_rfile.tsv \
        --labels B_fragilis_labels.csv --genome-qc B_fragilis_genome_qc.tsv \
        --min-ani 0.90 --out-prefix B_fragilis
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple

from evaluate_poppunk_fastani import (
    check_species_ani,
    make_symmetric_ani_matrix,
    normalise_genome_id,
    read_fastani,
)
from spp_eligibility_from_qc import compute_label

logger = logging.getLogger(__name__)


class BadPair(NamedTuple):
    """A genome pair that fails the ANI check. ``ani`` is a fraction, or None if missing."""

    genome_a: str
    genome_b: str
    ani: float | None


class Discard(NamedTuple):
    """A discarded genome and the bad pairs that it still had when it was discarded."""

    genome: str
    label: str
    n50: float
    partners: list[tuple[str, float | None]]


def select_discards(bad_pairs: list[BadPair], labels: dict[str, str], n50: dict[str, float]) -> list[Discard]:
    """Choose the genomes to discard so that no bad pair remains.

    :param bad_pairs: The bad pairs of the species.
    :param labels: The HQ or MQ label of each genome.
    :param n50: The N50 of each genome.
    :return: The discarded genomes, in discard order.
    """
    partners: dict[str, dict[str, float | None]] = defaultdict(dict)
    for genome_a, genome_b, ani in bad_pairs:
        partners[genome_a][genome_b] = ani
        partners[genome_b][genome_a] = ani

    discards = []
    while partners:
        genome = min(
            partners,
            key=lambda name: (-len(partners[name]), labels[name] != "MQ", n50.get(name, 0.0), name),
        )
        removed = partners.pop(genome)
        discards.append(Discard(genome, labels[genome], n50.get(genome, 0.0), sorted(removed.items())))
        for partner in removed:
            del partners[partner][genome]
            if not partners[partner]:
                del partners[partner]
    return discards


def read_rfile(path: Path) -> list[list[str]]:
    """Read the r-file.

    :param path: The r-file (sample, path; no header).
    :return: The rows of the r-file.
    """
    with path.open(newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t") if row]


def read_table(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a table with a header.

    :param path: The table.
    :param delimiter: The column delimiter.
    :return: The column names and the rows.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), list(reader)


def write_table(path: Path, columns: list[str], rows: list[dict[str, object]], delimiter: str) -> None:
    """Write a table with a header.

    :param path: The output table.
    :param columns: The column names.
    :param rows: The rows.
    :param delimiter: The column delimiter.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_filtered_ani(source: Path, output: Path, discarded: set[str]) -> None:
    """Copy the FastANI rows that do not include a discarded genome.

    :param source: The FastANI table.
    :param output: The filtered FastANI table.
    :param discarded: The discarded genomes.
    """
    names: dict[str, str] = {}

    def genome_id(path: str) -> str:
        if path not in names:
            names[path] = normalise_genome_id(path)
        return names[path]

    with source.open() as lines, output.open("w") as out:
        for line in lines:
            query, reference = line.split("\t", 2)[:2]
            if genome_id(query) not in discarded and genome_id(reference) not in discarded:
                out.write(line)


def format_partners(partners: list[tuple[str, float | None]]) -> str:
    """Format the bad pairs of a discarded genome.

    :param partners: The other genome of each pair, and the ANI fraction or None.
    :return: Items such as ``genome2:missing`` or ``genome3:89.50``, separated by semicolons.
    """
    return ";".join(f"{name}:{'missing' if ani is None else f'{ani * 100:.2f}'}" for name, ani in partners)


def main() -> None:
    """Run the ANI QC of one species."""
    parser = argparse.ArgumentParser(description="Discard the genomes behind missing or low-ANI FastANI pairs.")
    parser.add_argument("--species-name", required=True, help="Species name for the eligibility table.")
    parser.add_argument("--ani", required=True, type=Path, help="All-vs-all FastANI output.")
    parser.add_argument("--rfile", required=True, type=Path, help="r-file of SPECIESQC.")
    parser.add_argument("--labels", required=True, type=Path, help="Labels CSV of SPECIESQC (genome,label).")
    parser.add_argument("--genome-qc", required=True, type=Path, help="genome_qc table of SPECIESQC.")
    parser.add_argument("--min-ani", type=float, default=0.90, help="Minimum pair ANI, 0-1 scale (default 0.90).")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    rfile = read_rfile(args.rfile)
    genomes = sorted({normalise_genome_id(row[0]) for row in rfile})
    _, label_rows = read_table(args.labels, ",")
    labels = {normalise_genome_id(row["genome"]): row["label"] for row in label_rows}
    qc_columns, qc_rows = read_table(args.genome_qc, "\t")
    qc_ids = {normalise_genome_id(row["genome"]): row for row in qc_rows}

    missing = [genome for genome in genomes if genome not in labels or genome not in qc_ids]
    if missing:
        sys.exit(f"ERROR: {len(missing)} r-file genomes are not in the labels or genome_qc table: {missing[:5]}")

    n50 = {genome: float(qc_ids[genome].get("N50") or 0) for genome in genomes}
    fastani = read_fastani(str(args.ani))
    matrix = make_symmetric_ani_matrix(fastani, genomes)
    bad_pairs = [BadPair(*pair) for pair in check_species_ani(matrix, genomes, args.min_ani)]
    discards = select_discards(bad_pairs, labels, n50)
    discarded = {item.genome for item in discards}

    prefix = args.out_prefix
    kept_rfile = [row for row in rfile if normalise_genome_id(row[0]) not in discarded]
    with Path(f"{prefix}_ani_qc_rfile.tsv").open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(kept_rfile)

    kept_labels = [row for row in label_rows if normalise_genome_id(row["genome"]) not in discarded]
    write_table(Path(f"{prefix}_ani_qc_labels.csv"), ["genome", "label"], kept_labels, ",")

    for genome in discarded:
        qc_ids[genome].update(label="DISCARDED", reason="ani")
    write_table(Path(f"{prefix}_ani_qc_genome_qc.tsv"), qc_columns, qc_rows, "\t")

    write_filtered_ani(args.ani, Path(f"{prefix}.qc.ani.txt"), discarded)

    write_table(
        Path(f"{prefix}_ani_qc.tsv"),
        ["order", "genome", "label", "N50", "n_bad_pairs", "bad_pairs"],
        [
            {
                "order": order,
                "genome": item.genome,
                "label": item.label,
                "N50": f"{item.n50:g}",
                "n_bad_pairs": len(item.partners),
                "bad_pairs": format_partners(item.partners),
            }
            for order, item in enumerate(discards, start=1)
        ],
        "\t",
    )

    kept_label_values = [row["label"] for row in kept_labels]
    n_hq, n_mq = kept_label_values.count("HQ"), kept_label_values.count("MQ")
    spp_label, n_eff, hq_ratio = compute_label(n_hq, n_mq)
    reasons = [row.get("reason") for row in qc_rows if row.get("label") == "DISCARDED"]
    write_table(
        Path(f"{prefix}_ani_qc_eligibility.tsv"),
        [
            "species_name", "spp_label", "n_hq", "n_mq", "n_eff", "hq_ratio", "n_total_passing_qc",
            "n_discarded", "n_discarded_contiguity", "n_discarded_ani",
        ],
        [
            {
                "species_name": args.species_name,
                "spp_label": spp_label,
                "n_hq": n_hq,
                "n_mq": n_mq,
                "n_eff": n_eff,
                "hq_ratio": hq_ratio,
                "n_total_passing_qc": n_hq + n_mq,
                "n_discarded": len(reasons),
                "n_discarded_contiguity": reasons.count("contiguity"),
                "n_discarded_ani": reasons.count("ani"),
            }
        ],
        "\t",
    )

    logger.info(
        "%s: %d bad pairs, %d genomes discarded, %d kept -> %s",
        args.species_name, len(bad_pairs), len(discards), n_hq + n_mq, spp_label,
    )


if __name__ == "__main__":
    main()
