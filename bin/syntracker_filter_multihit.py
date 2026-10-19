#!/usr/bin/env python3
"""
Apply SynTracker's one-hit-per-genome rule and recompute the all-regions APSS
=============================================================================
SynTracker (1.4.0) BLASTs every reference region against the targets and marks a genome with
more than one hit in a region as invalid there (the region is multi-copy in that genome); a
region is analysed only if at least two genomes are valid. But its extraction loop
(`blast_per_region_process` in syntracker_first_stage/blast.py) then extracts *every* genome,
valid or not, taking an invalid genome's first hit. So whether the pairs of a multi-copy genome
get that region depends on whether some *other* single-hit genome is in the run: the result of
a pair depends on the other targets. Found on P. vulgatus (Codon run, 2026-10): region
DBGDQY010000030.1_120000_121000, where GCA_944320055 has two hits, was analysed in the chunk
jobs that also held a single-hit genome (620 regions for GCA_944320055 / GCF_021891395) and
skipped in the others (619).

This script applies the rule SynTracker meant: every per-region row in which either genome had
more than one BLAST hit for that region is dropped, and each pair's APSS is recomputed as the
mean of its remaining per-region scores (SynTracker's all-regions APSS is exactly that mean).
The result no longer depends on which other genomes were in the run, so chunk jobs agree.

Inputs (all from one SynTracker run):
  --per-region   summary_output/synteny_scores_per_region.csv (streamed; ~600 MB per 150 targets)
  --blast-dir    <ref>/blast_output/ (one <region>.tab per region: sseqid sstart send sstrand)
  --dictionary   combined_targets/sample_names_dictionary.tab (genome <TAB> Sample.N)

Outputs (`<p>` = --out-prefix):
  avg_synteny_scores_all_regions.csv   SynTracker's format: Ref_genome, Sample1, Sample2, APSS,
                                       Compared_regions (APSS summed with math.fsum, so the value
                                       does not depend on row order, written with 15 significant
                                       digits like R)
  <p>_syntracker_multihit.tsv          multihit_region_genomes, rows_read, rows_removed,
                                       pairs_in, pairs_changed, pairs_lost

Usage:
    syntracker_filter_multihit.py --per-region synteny_scores_per_region.csv \
        --blast-dir syntracker/GCA_048344735/blast_output \
        --dictionary syntracker/combined_targets/sample_names_dictionary.tab --out-prefix job
"""

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

OUTPUT = "avg_synteny_scores_all_regions.csv"
SAMPLE = re.compile(r"^(Sample\.\d+)_contig")  # as SynTracker parses its renamed contigs


def read_dictionary(path) -> dict:
    """Sample.N -> genome name, from SynTracker's sample_names_dictionary.tab (genome <TAB> Sample.N)."""
    names = {}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            genome, sample = line.rstrip("\n").split("\t")[:2]
            names[sample.strip()] = genome.strip()
    if not names:
        sys.exit(f"ERROR: sample names dictionary '{path}' is empty")
    return names


def multihit(blast_dir, names: dict) -> set:
    """{(region, genome)} for every genome with more than one BLAST hit in a region."""
    files = sorted(Path(blast_dir).glob("*.tab"))
    if not files:
        sys.exit(f"ERROR: no BLAST outputs (*.tab) in '{blast_dir}'")
    found = set()
    for f in files:
        counts = Counter()
        for line in f.read_text().splitlines():
            m = SAMPLE.match(line)
            if m:
                counts[m.group(1)] += 1
        region = f.name[: -len(".tab")]
        for sample, n in counts.items():
            if n > 1:
                if sample not in names:
                    sys.exit(f"ERROR: BLAST sample '{sample}' ({f.name}) is not in the names dictionary")
                found.add((region, names[sample]))
    return found


def filter_rows(per_region, drop: set):
    """Stream the per-region table: kept scores and the original region count per pair."""
    kept, original, read, removed = defaultdict(list), Counter(), 0, 0
    with open(per_region, newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"Ref_genome", "Sample1", "Sample2", "Region", "Synteny_score"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: per-region table '{per_region}' lacks columns {sorted(missing)}")
        for row in reader:
            read += 1
            pair = (row["Ref_genome"], row["Sample1"], row["Sample2"])
            original[pair] += 1
            region = row["Region"]
            if (region, row["Sample1"]) in drop or (region, row["Sample2"]) in drop:
                removed += 1
                continue
            kept[pair].append(float(row["Synteny_score"]))
    return kept, original, read, removed


def main():
    parser = argparse.ArgumentParser(description="Drop multi-hit genome/region rows from SynTracker's per-region "
                                                 "scores and recompute the all-regions APSS.")
    parser.add_argument("--per-region", required=True, help="synteny_scores_per_region.csv of one SynTracker run.")
    parser.add_argument("--blast-dir", required=True, help="The run's <ref>/blast_output directory.")
    parser.add_argument("--dictionary", required=True, help="The run's combined_targets/sample_names_dictionary.tab.")
    parser.add_argument("--out-prefix", required=True, help="Prefix of the report.")
    args = parser.parse_args()

    drop = multihit(args.blast_dir, read_dictionary(args.dictionary))
    kept, original, read, removed = filter_rows(args.per_region, drop)
    if read == 0:
        sys.exit(f"ERROR: per-region table '{args.per_region}' has no rows")

    with open(OUTPUT, "w") as handle:
        handle.write('"Ref_genome","Sample1","Sample2","APSS","Compared_regions"\n')
        for (ref, s1, s2) in sorted(kept):
            scores = kept[(ref, s1, s2)]
            handle.write(f'"{ref}","{s1}","{s2}",{math.fsum(scores) / len(scores):.15g},{len(scores)}\n')

    changed = sum(1 for pair, n in original.items() if len(kept.get(pair, [])) != n)
    lost = sum(1 for pair in original if pair not in kept)
    report = {"multihit_region_genomes": len(drop), "rows_read": read, "rows_removed": removed,
              "pairs_in": len(original), "pairs_changed": changed, "pairs_lost": lost}
    with open(f"{args.out_prefix}_syntracker_multihit.tsv", "w") as handle:
        handle.write("\t".join(report) + "\n" + "\t".join(str(v) for v in report.values()) + "\n")
    print(f"{len(drop)} multi-hit region/genome combinations: removed {removed} of {read} per-region rows; "
          f"{changed} of {len(original)} pairs changed, {lost} lost", file=sys.stderr)


if __name__ == "__main__":
    main()
