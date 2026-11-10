#!/usr/bin/env python3
"""
Neighbour-joining tree from all-vs-all FastANI
==============================================
Turns one species' all-vs-all FastANI output (five columns: query, reference, ANI in
percent, matched and total fragments) into a Newick tree. The distance of a pair is
1 - ANI/100, averaged over the A->B and B->A estimates; the tree is built by
neighbour joining (Biopython), or by hand for two genomes. Every pair must be present:
FastANI omits pairs that are too divergent, which fails the run with examples.

Tip names are the file names minus their FASTA extension, as in
normalise_genome_id() (bin/evaluate_poppunk_fastani.py), with any character outside
[A-Za-z0-9_.-] replaced by `_` to keep the Newick valid.

Usage:
    fastani_to_tree.py B_longum.ani.txt B_longum_ani_nj_tree.nwk
"""

import argparse
import re
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean

from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def genome_label(path):
    """Convert a FASTA path into a Newick-safe genome name."""
    name = Path(path).name

    extensions = (
        ".fasta.gz", ".fna.gz", ".fa.gz", ".fas.gz",
        ".fasta", ".fna", ".fa", ".fas", ".contigs",
    )

    for extension in extensions:
        if name.lower().endswith(extension):
            name = name[:-len(extension)]
            break

    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def parse_fastani(filename):
    pairwise_ani = defaultdict(list)
    label_to_path = {}
    genomes = set()

    with open(filename, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) != 5:
                raise ValueError(
                    f"Line {line_number}: expected five columns, "
                    f"found {len(fields)}"
                )

            query_path, reference_path, ani_text, matched, total = fields

            try:
                ani = float(ani_text)
                matched = int(matched)
                total = int(total)
            except ValueError as error:
                raise ValueError(
                    f"Line {line_number}: invalid numeric value"
                ) from error

            if not 0 <= ani <= 100:
                raise ValueError(
                    f"Line {line_number}: ANI must be between 0 and 100"
                )

            if total <= 0 or not 0 <= matched <= total:
                raise ValueError(
                    f"Line {line_number}: invalid fragment counts"
                )

            query = genome_label(query_path)
            reference = genome_label(reference_path)

            for label, path in (
                (query, query_path),
                (reference, reference_path),
            ):
                previous_path = label_to_path.get(label)

                if previous_path is not None and previous_path != path:
                    raise ValueError(
                        f"Genome label collision: {previous_path!r} and "
                        f"{path!r} both become {label!r}"
                    )

                label_to_path[label] = path

            genomes.update((query, reference))

            if query != reference:
                pair = tuple(sorted((query, reference)))
                pairwise_ani[pair].append(ani)

    return sorted(genomes), pairwise_ani


def create_distance_matrix(genomes, pairwise_ani):
    missing_pairs = [
        (genome_a, genome_b)
        for genome_a, genome_b in combinations(genomes, 2)
        if (genome_a, genome_b) not in pairwise_ani
    ]

    if missing_pairs:
        examples = ", ".join(
            f"{genome_a} <-> {genome_b}"
            for genome_a, genome_b in missing_pairs[:10]
        )

        raise ValueError(
            f"{len(missing_pairs)} genome pairs are missing. "
            f"Examples: {examples}. fastANI may omit pairs that are "
            f"too divergent."
        )

    # Biopython requires the lower triangular matrix, including the diagonal.
    lower_triangle = []

    for row_index, genome_a in enumerate(genomes):
        row = []

        for column_index in range(row_index + 1):
            genome_b = genomes[column_index]

            if genome_a == genome_b:
                row.append(0.0)
            else:
                pair = tuple(sorted((genome_a, genome_b)))

                # Average A→B and B→A estimates when both are present.
                ani = mean(pairwise_ani[pair])
                distance = 1.0 - ani / 100.0
                row.append(distance)

        lower_triangle.append(row)

    return DistanceMatrix(
        names=genomes,
        matrix=lower_triangle,
    )


def build_tree(distance_matrix):
    if len(distance_matrix.names) == 2:
        genome_a, genome_b = distance_matrix.names
        distance = distance_matrix[genome_a, genome_b]

        root = Clade()
        root.clades = [
            Clade(name=genome_a, branch_length=distance / 2),
            Clade(name=genome_b, branch_length=distance / 2),
        ]

        return Tree(root=root, rooted=False)

    constructor = DistanceTreeConstructor()
    return constructor.nj(distance_matrix)


def main():
    parser = argparse.ArgumentParser(
        description="Convert all-versus-all fastANI output into a Newick tree."
    )

    parser.add_argument(
        "fastani_output",
        help="Five-column fastANI output file",
    )

    parser.add_argument(
        "newick_output",
        help="Output Newick filename",
    )

    args = parser.parse_args()

    genomes, pairwise_ani = parse_fastani(args.fastani_output)

    if len(genomes) < 2:
        raise ValueError("At least two genomes are required")

    distance_matrix = create_distance_matrix(genomes, pairwise_ani)
    tree = build_tree(distance_matrix)

    negative_lengths = [
        clade.branch_length
        for clade in tree.find_clades()
        if clade.branch_length is not None and clade.branch_length < 0
    ]

    if negative_lengths:
        print(
            "Warning: the tree contains negative branch lengths; "
            "the ANI distances may not be perfectly tree-like.",
            file=sys.stderr,
        )

    Phylo.write(tree, args.newick_output, "newick")

    print(
        f"Wrote neighbour-joining tree for {len(genomes)} genomes "
        f"to {args.newick_output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
