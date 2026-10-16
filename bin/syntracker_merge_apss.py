#!/usr/bin/env python3
"""
Merge SynTracker all-regions APSS tables from chunked runs
==========================================================
The pipeline splits a species' SynTracker targets into chunks and runs SynTracker once per
chunk pair (Ci + Cj), all against the same reference. A pair's all-regions APSS depends only on
the two genomes and the reference (BLAST hits are filtered per genome, synteny is compared pair
by pair), so the tables of the chunk runs can be concatenated. Verified on E. lenta: 20 genomes
run alone gave the same APSS and Compared_regions as in the 199-genome run for 190/190 pairs.

Only the all-regions tables (`avg_synteny_scores_all_regions.csv`) can be merged. The
subsampled `avg_synteny_scores_<n>_regions.csv` tables cannot: SynTracker draws each pair's n
regions with one fixed seed over the whole sorted table, so a pair's draw depends on every
pair before it, i.e. on which other targets were in the run.

Steps:
  1. Read every table (columns Ref_genome, Sample1, Sample2, APSS, Compared_regions).
  2. Assert a single reference genome across all tables.
  3. Key each pair without order (a pair may appear as A,B in one run and B,A in another).
  4. A pair found in more than one table must have exactly the same APSS and
     Compared_regions; otherwise the merge fails (the runs are not comparable).
  5. Write one row per pair, sorted, keeping SynTracker's text values and quoting.

With --targets, the number of target pairs missing from the merge is reported (SynTracker
leaves out pairs that share no region, so this is a warning, not an error).

Outputs (`<p>` = --out-prefix):
  avg_synteny_scores_all_regions.csv   the merged table, in SynTracker's format
  <p>_syntracker_merge.tsv             n_tables, rows_read, duplicate_rows, unique_pairs,
                                       n_targets, expected_pairs, missing_pairs

Usage:
    syntracker_merge_apss.py --apss chunk_01/avg_synteny_scores_all_regions.csv ... \
        --targets targets.txt --out-prefix B_fragilis
"""

import argparse
import csv
import sys
from pathlib import Path

from evaluate_poppunk_fastani import normalise_genome_id

COLUMNS = ["Ref_genome", "Sample1", "Sample2", "APSS", "Compared_regions"]
OUTPUT = "avg_synteny_scores_all_regions.csv"


def read_table(path: str) -> list:
    """Rows of one all-regions table, as dicts of the raw text values."""
    try:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle)
            missing = set(COLUMNS) - set(reader.fieldnames or [])
            if missing:
                sys.exit(f"ERROR: APSS table '{path}' lacks columns {sorted(missing)}; found {reader.fieldnames}")
            return [{c: row[c] for c in COLUMNS} for row in reader]
    except FileNotFoundError as err:
        sys.exit(f"ERROR: cannot read APSS table '{path}': {err}")


def merge(paths: list) -> tuple:
    """(merged rows sorted by pair, number of rows read, number of duplicate rows)."""
    merged, refs, rows_read, duplicates = {}, set(), 0, 0
    for path in paths:
        for row in read_table(path):
            rows_read += 1
            refs.add(row["Ref_genome"])
            if row["Sample1"] == row["Sample2"]:
                continue
            key = tuple(sorted((row["Sample1"], row["Sample2"])))
            seen = merged.get(key)
            if seen is None:
                merged[key] = (row, path)
                continue
            duplicates += 1
            first, first_path = seen
            if float(first["APSS"]) != float(row["APSS"]) or int(first["Compared_regions"]) != int(row["Compared_regions"]):
                sys.exit(f"ERROR: pair {key[0]} / {key[1]} differs between '{first_path}' "
                         f"(APSS {first['APSS']}, {first['Compared_regions']} regions) and '{path}' "
                         f"(APSS {row['APSS']}, {row['Compared_regions']} regions); the chunk runs are not comparable.")
    if len(refs) > 1:
        sys.exit(f"ERROR: the APSS tables use {len(refs)} reference genomes ({sorted(refs)}); "
                 "chunk runs must share one reference.")
    if not merged:
        sys.exit("ERROR: the APSS tables hold no pairs.")
    return [merged[key][0] for key in sorted(merged)], rows_read, duplicates


def write_table(rows: list, path: str) -> None:
    """SynTracker's format: quoted names, unquoted numbers."""
    with open(path, "w") as handle:
        handle.write(",".join(f'"{c}"' for c in COLUMNS) + "\n")
        for row in rows:
            handle.write(f'"{row["Ref_genome"]}","{row["Sample1"]}","{row["Sample2"]}",{row["APSS"]},{row["Compared_regions"]}\n')


def main():
    parser = argparse.ArgumentParser(description="Merge SynTracker all-regions APSS tables from chunked runs.")
    parser.add_argument("--apss", nargs="+", required=True, help="avg_synteny_scores_all_regions.csv of each run.")
    parser.add_argument("--targets", default=None, help="All SynTracker targets (file or sample names), one per line.")
    parser.add_argument("--out-prefix", required=True, help="Prefix of the merge report.")
    args = parser.parse_args()

    rows, rows_read, duplicates = merge(args.apss)
    write_table(rows, OUTPUT)

    report = {"n_tables": len(args.apss), "rows_read": rows_read, "duplicate_rows": duplicates,
              "unique_pairs": len(rows), "n_targets": "", "expected_pairs": "", "missing_pairs": ""}
    if args.targets:
        targets = {normalise_genome_id(t.strip()) for t in Path(args.targets).read_text().splitlines() if t.strip()}
        found = {tuple(sorted((normalise_genome_id(r["Sample1"]), normalise_genome_id(r["Sample2"])))) for r in rows}
        expected = len(targets) * (len(targets) - 1) // 2
        missing = expected - len({p for p in found if p[0] in targets and p[1] in targets})
        report.update(n_targets=len(targets), expected_pairs=expected, missing_pairs=missing)
        if missing:
            print(f"WARNING: {missing} of {expected} target pairs are absent from the merged APSS "
                  "(SynTracker leaves out pairs that share no region).", file=sys.stderr)
    with open(f"{args.out_prefix}_syntracker_merge.tsv", "w") as handle:
        handle.write("\t".join(report) + "\n" + "\t".join(str(v) for v in report.values()) + "\n")
    print(f"merged {len(args.apss)} tables: {rows_read} rows, {duplicates} duplicates, {len(rows)} pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
