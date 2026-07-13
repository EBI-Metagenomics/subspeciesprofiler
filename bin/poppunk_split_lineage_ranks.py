#!/usr/bin/env python3
"""Split a PopPUNK lineage `_lineages.csv` into one clustering per rank.

PopPUNK lineage mode fits several rank resolutions in one run and writes them as columns
(`Rank_1_Lineage`, `Rank_2_Lineage`, ...) of `<prefix>_lineages.csv`, plus an `overall_Lineage`
column that equals the finest rank (already the content of `<prefix>_clusters.csv`). Each rank is
a genuinely different partition, so to score them all we emit one `Taxon,Cluster` clustering per
rank into `lineage_rank<N>/lineage_rank<N>_clusters.csv`, ready for the ANI evaluator. The
`overall_Lineage` column is skipped (redundant with the finest rank).
"""
import argparse
import csv
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lineages", required=True, help="PopPUNK <prefix>_lineages.csv")
    parser.add_argument("--outdir", default=".", help="Directory to write per-rank clustering dirs into.")
    args = parser.parse_args()

    with open(args.lineages, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        id_col = fieldnames[0] if fieldnames else None
        rank_cols = [c for c in fieldnames if re.fullmatch(r"Rank_\d+_Lineage", c)]
        if id_col is None or not rank_cols:
            sys.exit(
                f"ERROR: expected an id column plus Rank_<N>_Lineage columns in '{args.lineages}'; "
                f"found {fieldnames}"
            )
        rows = list(reader)

    for col in rank_cols:
        rank = re.search(r"Rank_(\d+)_Lineage", col).group(1)
        name = f"lineage_rank{rank}"
        rank_dir = os.path.join(args.outdir, name)
        os.makedirs(rank_dir, exist_ok=True)
        with open(os.path.join(rank_dir, f"{name}_clusters.csv"), "w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["Taxon", "Cluster"])
            for row in rows:
                writer.writerow([row[id_col], row[col]])
        print(f"wrote {name}_clusters.csv ({len(rows)} genomes)", file=sys.stderr)


if __name__ == "__main__":
    main()
