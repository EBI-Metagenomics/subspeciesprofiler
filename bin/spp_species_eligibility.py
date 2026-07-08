#!/usr/bin/env python3
"""
SPP Species Group Eligibility Assessment
=========================================
Joins GTDB metadata (bacteria + archaea) to RefSeq assembly summary on
accession, extracts CheckM2 QC values, groups by GTDB species epithet,
and classifies each target species group against SPP inclusion criteria.

Usage:
    python spp_species_eligibility.py \
        --species   target_species.tsv \
        --gtdb_bac  bac120_metadata.tsv \
        --gtdb_arc  ar53_metadata.tsv \
        --refseq    assembly_summary_refseq.txt \
        --output    spp_eligibility_report.tsv

Species list TSV format (tab-separated):
    query_name                   synonyms                      notes
    Escherichia coli                                           Proteobacteria
    Phocaeicola vulgatus         Bacteroides vulgatus          Bacteroidota; reclassified
    Ruminococcus gnavus          Mediterraneibacter gnavus     Firmicutes

    - query_name is required; synonyms and notes are optional.
    - synonyms: semicolon-separated list of alternate names to match in GTDB taxonomy.

Optional:
    --genbank      assembly_summary_genbank.txt
    --min_completeness   80   (default)
    --hq_completeness    90   (default)
    --hq_contamination    1   (default; max contamination for HQ)
    --max_contamination   5   (default; max contamination for MQ)
"""

import argparse
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: pip install pandas")


# ---------------------------------------------------------------------------
# Species list loader — replaces hardcoded TARGET_SPECIES
# ---------------------------------------------------------------------------

def load_species_list(path: str):
    """
    Load target species from a TSV file.
    Required column:  query_name   (binomial, e.g. 'Escherichia coli')
    Optional columns: synonyms     (semicolon-separated, e.g. 'Bacteroides vulgatus')
                      notes        (free text)

    Returns a list of (query_name, [synonyms], notes) tuples,
    matching the structure previously used by TARGET_SPECIES.
    """
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "query_name" not in df.columns:
        sys.exit(
            f"ERROR: species list file '{path}' must have a 'query_name' column. "
            f"Found: {list(df.columns)}"
        )
    result = []
    for _, row in df.iterrows():
        query    = row["query_name"].strip()
        synonyms = [s.strip() for s in row.get("synonyms", "").split(";") if s.strip()]
        notes    = row.get("notes", "").strip()
        if query:
            result.append((query, synonyms, notes))
    if not result:
        sys.exit(f"ERROR: no valid species found in '{path}'")
    return result


# ---------------------------------------------------------------------------
# SPP thresholds
# ---------------------------------------------------------------------------
MIN_COMP   = 80.0
HQ_COMP    = 90.0
HQ_CONT    = 1.0   # max contamination allowed for HQ; MAX_CONT (5.0) is the MQ ceiling
MAX_CONT   = 5.0

STRONG_N_HQ   = 75
STRONG_N_EFF  = 100
ACCEPT_N_HQ   = 50
ACCEPT_N_EFF  = 80

BOUNDARY_HQ_MARGIN  = 5
BOUNDARY_EFF_MARGIN = 8
NEAR_THRESH_N_HQ    = 40
NEAR_THRESH_N_EFF   = 64


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_gtdb(path: str, source_label: str) -> pd.DataFrame:
    """
    Load GTDB metadata TSV (bac120 or ar53).
    Returns accession, completeness, contamination, gtdb_taxonomy, gtdb_species.
    Accession format in GTDB: GB_GCA_... or RS_GCF_...
    We strip the two-letter prefix to get a plain GCA_/GCF_ accession.
    """
    needed = {
        "accession",
        "checkm2_completeness",
        "checkm2_contamination",
        "gtdb_taxonomy",
    }
    df = pd.read_csv(
        path, sep="\t", low_memory=False,
        usecols=lambda c: c in needed,
    )
    df = df.rename(columns={
        "checkm2_completeness": "completeness",
        "checkm2_contamination": "contamination",
    })

    # Strip GTDB prefix and version:
    # "RS_GCF_000001405.40" -> "GCF_000001405"
    # "GB_GCA_000007325.1"  -> "GCA_000007325"
    df["accession_clean"] = (
        df["accession"]
        .str.replace(r"^[A-Z]{2}_", "", regex=True)   # remove RS_ / GB_
        .str.replace(r"\.\d+$",     "", regex=True)   # remove version suffix
    )

    # Extract species epithet from gtdb_taxonomy: last s__ field
    # e.g. "d__Bacteria;...;s__Escherichia coli" -> "Escherichia coli"
    df["gtdb_species"] = df["gtdb_taxonomy"].str.extract(r"s__([^;]+)$")[0].str.strip()

    df["domain"] = source_label  # "bacteria" or "archaea"
    return df[["accession_clean", "completeness", "contamination",
               "gtdb_species", "gtdb_taxonomy", "domain"]]


def load_assembly_summary(path: str, source: str) -> pd.DataFrame:
    """
    Load NCBI assembly_summary (RefSeq or GenBank).
    Keeps only latest assemblies.

    RefSeq:  filters out excluded_from_refseq entries.
    GenBank: filters out assemblies where paired_asm_comp == 'identical',
             which flags GenBank entries that are exact duplicates of a
             RefSeq record — keeping them would cause double-counting when
             both tables are loaded.

    Note: in the GenBank summary, 'excluded_from_refseq' is concatenated
    with the next column name (missing tab) and is therefore not parseable
    as a clean column — the paired_asm_comp filter handles deduplication
    instead.
    """
    wanted = {
        "#assembly_accession", "version_status",
        "excluded_from_refseq", "ftp_path", "paired_asm_comp",
        "biosample", "bioproject",
    }
    df = pd.read_csv(
        path, sep="\t", skiprows=1, low_memory=False,
        usecols=lambda c: c in wanted,
    )
    df = df.rename(columns={"#assembly_accession": "assembly_accession"})

    # Keep only latest assemblies
    df = df[df["version_status"] == "latest"].copy()

    if source == "refseq":
        # Filter suppressed/excluded assemblies (RefSeq only)
        if "excluded_from_refseq" in df.columns:
            df = df[
                df["excluded_from_refseq"].isna() |
                (df["excluded_from_refseq"].astype(str) == "na")
            ].copy()
    else:
        # GenBank: drop assemblies that are identical to a RefSeq record
        # to avoid double-counting when both tables are provided
        if "paired_asm_comp" in df.columns:
            df = df[
                df["paired_asm_comp"].astype(str).str.lower() != "identical"
            ].copy()

    df["accession_clean"] = df["assembly_accession"].str.replace(r"\.\d+$", "", regex=True)
    df["ncbi_source"] = source

    # Ensure biosample and bioproject columns exist even if absent in the file
    for col in ("biosample", "bioproject"):
        if col not in df.columns:
            df[col] = "na"

    return df[["accession_clean", "ncbi_source", "ftp_path", "biosample", "bioproject"]]


def deduplicate_species_group(grp, biosample_threshold=1000, bioproject_cap=0.20):
    """
    Remove technical replicates from a species group.
    Both deduplication levels are applied only when n_total > biosample_threshold
    to avoid inadvertently reducing small species groups below inclusion thresholds.

    Level 1 — biosample deduplication:
      Multiple assemblies sharing the same biosample accession are the same
      biological sample sequenced more than once. Keep one per biosample,
      preferring HQ over MQ, then highest completeness within the same tier.

    Level 2 — bioproject cap:
      A single large surveillance project can dominate the species group.
      Cap each bioproject to at most bioproject_cap * n_total genomes,
      keeping the highest-completeness genomes within each project.
      Genomes with missing/na bioproject accession are not capped.

    Returns the deduplicated DataFrame and a dict of counts for reporting.
    """
    n_before = len(grp)
    stats = {"n_before": n_before, "n_after_biosample": n_before,
             "n_after_bioproject": n_before, "bioprojects_capped": []}

    # Both levels only apply when the species group is large enough
    if n_before <= biosample_threshold:
        return grp, stats

    # --- Level 1: biosample deduplication ---
    has_biosample = (
        grp["biosample"].notna() &
        (grp["biosample"].astype(str).str.strip() != "na") &
        (grp["biosample"].astype(str).str.strip() != "")
    )
    with_bs    = grp[has_biosample].copy()
    without_bs = grp[~has_biosample].copy()

    if not with_bs.empty:
        # Sort so HQ comes before MQ, then highest completeness first
        # qc_label: HQ < MQ alphabetically, which is the order we want
        with_bs = with_bs.sort_values(
            ["biosample", "qc_label", "completeness"],
            ascending=[True, True, False],
        )
        with_bs = with_bs.drop_duplicates(subset=["biosample"], keep="first")

    grp = pd.concat([with_bs, without_bs], ignore_index=True)
    stats["n_after_biosample"] = len(grp)

    # --- Level 2: bioproject cap ---
    n_current = len(grp)
    max_per_project = max(1, int(n_current * bioproject_cap))

    has_bp = (
        grp["bioproject"].notna() &
        (grp["bioproject"].astype(str).str.strip() != "na") &
        (grp["bioproject"].astype(str).str.strip() != "")
    )
    with_bp    = grp[has_bp].copy()
    without_bp = grp[~has_bp].copy()

    if not with_bp.empty:
        # Sort by completeness descending to keep best genomes per project
        with_bp = with_bp.sort_values("completeness", ascending=False)
        project_counts = with_bp.groupby("bioproject").size()
        capped = project_counts[project_counts > max_per_project].index.tolist()
        stats["bioprojects_capped"] = capped

        with_bp = with_bp.groupby("bioproject", group_keys=False).apply(
            lambda x: x.head(max_per_project)
        )

    grp = pd.concat([with_bp, without_bp], ignore_index=True)
    stats["n_after_bioproject"] = len(grp)

    return grp, stats


def make_short_name(query: str) -> str:
    """
    Build a 6-character short name from a species binomial.
    Rule: first letter of genus + first 5 letters of species epithet,
    all lowercase. Spaces/hyphens in the epithet are ignored.
    Examples:
        'Leptospira interrogans' -> 'linter'
        'Escherichia coli'       -> 'ecoli'
        'Fusobacterium nucleatum'-> 'fnucle'
    """
    parts = query.strip().split()
    genus   = parts[0][0].lower() if parts else "x"
    epithet = "".join(parts[1:]).lower()[:5] if len(parts) > 1 else "xxxxx"
    return genus + epithet


# ---------------------------------------------------------------------------
# QC and inclusion logic
# ---------------------------------------------------------------------------

def classify_genome(comp, cont, min_comp, hq_comp, hq_cont, max_cont):
    try:
        comp = float(comp)
        cont = float(cont)
    except (ValueError, TypeError):
        return "DISCARDED"
    if comp >= hq_comp and cont <= hq_cont:      # HQ: >=90 completeness, <=1% contamination
        return "HQ"
    elif comp >= min_comp and cont <= max_cont:  # MQ: >=80 completeness, <=5% contamination
        return "MQ"
    return "DISCARDED"


def compute_label(n_hq, n_mq):
    n_eff    = n_hq + 0.5 * n_mq
    hq_ratio = n_hq / (n_hq + n_mq) if (n_hq + n_mq) > 0 else 0.0

    if n_hq >= STRONG_N_HQ and n_eff >= STRONG_N_EFF:
        base = "STRONG"
    elif n_hq >= ACCEPT_N_HQ and n_eff >= ACCEPT_N_EFF:
        base = "ACCEPTABLE"
    elif n_hq >= ACCEPT_N_HQ:
        base = "NOT_ELIGIBLE"
    else:
        base = "DISCARDED"

    suffix = ""
    if base in ("STRONG", "ACCEPTABLE"):
        thr_hq  = STRONG_N_HQ  if base == "STRONG" else ACCEPT_N_HQ
        thr_eff = STRONG_N_EFF if base == "STRONG" else ACCEPT_N_EFF
        if (n_hq - thr_hq) < BOUNDARY_HQ_MARGIN or (n_eff - thr_eff) < BOUNDARY_EFF_MARGIN:
            suffix = "_BOUNDARY"
    elif base == "DISCARDED":
        if n_hq >= NEAR_THRESH_N_HQ and n_eff >= NEAR_THRESH_N_EFF:
            suffix = "_NEAR_THRESHOLD"

    return base + suffix, round(n_eff, 1), round(hq_ratio, 3)


# ---------------------------------------------------------------------------
# Species matching
# ---------------------------------------------------------------------------

def matches_target(gtdb_species: str, query: str, synonyms: list) -> bool:
    """
    Match query name against GTDB species epithet.

    Two-step logic:
      1. Exact match: 'Desulfovibrio desulfuricans' matches only itself,
         not 'Desulfovibrio desulfuricans_B' or '_C' etc.
      2. Space-prefix match: 'Escherichia coli' matches
         'Escherichia coli K-12' or 'Escherichia coli str. O157'
         because the suffix follows a space, not an underscore.

    GTDB alphabetic suffixes (_B, _C ...) are separated by underscore,
    not space, so they are correctly excluded from unsuffixed queries.
    To target a specific suffixed lineage, include the full suffix in
    query_name (e.g. 'Desulfovibrio desulfuricans_B').
    """
    if not isinstance(gtdb_species, str):
        return False
    s = gtdb_species.strip().lower()
    candidates = [query.strip().lower()] + [syn.strip().lower() for syn in synonyms]
    for cand in candidates:
        if s == cand:
            return True
        if s.startswith(cand + " "):
            return True
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="SPP species eligibility assessment")
    parser.add_argument("--species",   required=True,  help=(
        "TSV file with target species. Required column: query_name. "
        "Optional columns: synonyms (semicolon-separated), notes."
    ))
    parser.add_argument("--gtdb_bac", required=True,  help="bac120_metadata.tsv")
    parser.add_argument("--gtdb_arc", required=True,  help="ar53_metadata.tsv")
    parser.add_argument("--refseq",   required=True,  help="assembly_summary_refseq.txt")
    parser.add_argument("--genbank",  default=None,   help="assembly_summary_genbank.txt (optional)")
    parser.add_argument("--output",   default="spp_eligibility_report.tsv")
    parser.add_argument("--min_completeness",  type=float, default=MIN_COMP)
    parser.add_argument("--hq_completeness",   type=float, default=HQ_COMP)
    parser.add_argument("--hq_contamination",  type=float, default=HQ_CONT)
    parser.add_argument("--max_contamination", type=float, default=MAX_CONT)
    args = parser.parse_args()

    min_comp = args.min_completeness
    hq_comp  = args.hq_completeness
    hq_cont  = args.hq_contamination
    max_cont = args.max_contamination

    print(f"Loading species list from: {args.species}", file=sys.stderr)
    target_species = load_species_list(args.species)
    print(f"  {len(target_species)} species loaded.", file=sys.stderr)

    # ------------------------------------------------------------------
    # Load and concatenate GTDB tables
    # ------------------------------------------------------------------
    print("Loading GTDB bacteria...", file=sys.stderr)
    gtdb_bac = load_gtdb(args.gtdb_bac, "bacteria")

    print("Loading GTDB archaea...", file=sys.stderr)
    gtdb_arc = load_gtdb(args.gtdb_arc, "archaea")

    gtdb = pd.concat([gtdb_bac, gtdb_arc], ignore_index=True)
    print(f"  Total GTDB genomes loaded: {len(gtdb):,}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Load NCBI assembly summaries
    # ------------------------------------------------------------------
    print("Loading RefSeq assembly summary...", file=sys.stderr)
    refseq = load_assembly_summary(args.refseq, "refseq")

    ncbi = refseq.copy()
    if args.genbank:
        print("Loading GenBank assembly summary...", file=sys.stderr)
        genbank = load_assembly_summary(args.genbank, "genbank")
        # Prefer RefSeq label when accession appears in both
        ncbi = pd.concat([refseq, genbank], ignore_index=True)
        ncbi = ncbi.sort_values("ncbi_source").drop_duplicates(
            subset="accession_clean", keep="first"
        )

    print(f"  Total NCBI assemblies loaded: {len(ncbi):,}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Join GTDB <-> NCBI on accession
    # Note: GTDB accession_clean is GCA_/GCF_ without version.
    # NCBI accession_clean is also stripped of version above.
    # For RefSeq (GCF_) genomes GTDB stores the GCF_ directly.
    # For GenBank-only (GCA_) genomes GTDB stores GCA_.
    # ------------------------------------------------------------------
    merged = gtdb.merge(ncbi, on="accession_clean", how="left")

    # Drop genomes not present in the NCBI summary (suppressed, withdrawn,
    # or otherwise unavailable). These must be excluded before QC counting
    # to avoid inflating N_HQ / N_MQ with genomes that cannot be downloaded.
    n_before = len(merged)
    unavailable = merged["ftp_path"].isna() | (merged["ftp_path"].astype(str) == "na")
    n_unavailable = unavailable.sum()
    merged = merged[~unavailable].copy()
    print(
        f"  Dropped {n_unavailable:,} GTDB genomes not found in NCBI summary "
        f"(suppressed/withdrawn). {len(merged):,} remain.",
        file=sys.stderr,
    )

    merged["ncbi_source"] = merged["ncbi_source"].fillna("not_in_ncbi_summary")

    # Build ready-to-download FTP URL:
    # NCBI ftp_path is the assembly directory, e.g.:
    #   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/325/GCF_000007325.1_ASM732v1
    # The genomic FASTA filename is the directory basename + _genomic.fna.gz
    def build_fasta_url(ftp_dir: str) -> str:
        ftp_dir = str(ftp_dir).rstrip("/")
        if not ftp_dir or ftp_dir == "na":
            return "na"
        basename = ftp_dir.split("/")[-1]
        return f"{ftp_dir}/{basename}_genomic.fna.gz"

    merged["ftp_url"] = merged["ftp_path"].apply(build_fasta_url)

    # Apply QC classification
    merged["qc_label"] = merged.apply(
        lambda r: classify_genome(
            r["completeness"], r["contamination"], min_comp, hq_comp, hq_cont, max_cont
        ),
        axis=1,
    )

    # ------------------------------------------------------------------
    # Match each target species against GTDB species epithet
    # ------------------------------------------------------------------
    results    = []
    ftp_rows   = []   # genome-level FTP output

    for query, synonyms, notes in target_species:
        mask = merged["gtdb_species"].apply(
            lambda s: matches_target(s, query, synonyms)
        )
        grp = merged[mask]

        if grp.empty:
            results.append({
                "query_name":           query,
                "gtdb_species_matched": "NO_MATCH",
                "n_gtdb_matched":       0,
                "n_refseq":             0,
                "n_genbank_only":       0,
                "n_hq":                 0,
                "n_mq":                 0,
                "n_discarded":          0,
                "n_total_passing_qc":   0,
                "n_eff":                0.0,
                "hq_ratio":             0.0,
                "spp_label":            "NO_MATCH",
                "domain":               "unknown",
                "synonyms_checked":     "; ".join(synonyms),
                "notes":                notes,
            })
            print(f"  {query}: NO MATCH in GTDB", file=sys.stderr)
            continue

        # Report which GTDB species names were matched (may be >1 after reclassification)
        matched_names = grp["gtdb_species"].dropna().unique()
        matched_str   = "; ".join(sorted(matched_names))
        domain        = grp["domain"].mode()[0]

        # Apply deduplication before counting and FTP collection
        grp, dedup_stats = deduplicate_species_group(
            grp, biosample_threshold=1000, bioproject_cap=0.20
        )

        n_hq   = (grp["qc_label"] == "HQ").sum()
        n_mq   = (grp["qc_label"] == "MQ").sum()
        n_disc = (grp["qc_label"] == "DISCARDED").sum()

        n_refseq  = (grp["ncbi_source"] == "refseq").sum()
        n_gb_only = (grp["ncbi_source"] == "genbank").sum()

        label, n_eff, hq_ratio = compute_label(n_hq, n_mq)

        # Build deduplication summary string for the report
        dedup_note = ""
        if dedup_stats["n_before"] != dedup_stats["n_after_bioproject"]:
            n_removed = dedup_stats["n_before"] - dedup_stats["n_after_bioproject"]
            n_bs_removed = dedup_stats["n_before"] - dedup_stats["n_after_biosample"]
            n_bp_removed = dedup_stats["n_after_biosample"] - dedup_stats["n_after_bioproject"]
            dedup_note = (
                "biosample_dedup={} bioproject_cap={}({} projects capped)".format(
                    n_bs_removed, n_bp_removed,
                    len(dedup_stats["bioprojects_capped"])
                )
            )

        results.append({
            "query_name":           query,
            "gtdb_species_matched": matched_str,
            "n_gtdb_matched":       dedup_stats["n_before"],
            "n_after_dedup":        dedup_stats["n_after_bioproject"],
            "n_refseq":             n_refseq,
            "n_genbank_only":       n_gb_only,
            "n_hq":                 n_hq,
            "n_mq":                 n_mq,
            "n_discarded":          n_disc,
            "n_total_passing_qc":   n_hq + n_mq,
            "n_eff":                n_eff,
            "hq_ratio":             hq_ratio,
            "spp_label":            label,
            "domain":               domain,
            "dedup_summary":        dedup_note,
            "synonyms_checked":     "; ".join(synonyms),
            "notes":                notes,
        })

        print(
            f"  {query}: {dedup_stats['n_before']} matched -> "
            f"{dedup_stats['n_after_bioproject']} after dedup "
            f"({n_hq} HQ, {n_mq} MQ, {n_disc} discarded) -> {label}"
            + (f" [{dedup_note}]" if dedup_note else ""),
            file=sys.stderr,
        )

        # Collect genome-level FTP rows for genomes in the final deduplicated set
        short = make_short_name(query)
        passing = grp[grp["qc_label"].isin(["HQ", "MQ"])].copy()
        for _, genome in passing.iterrows():
            ftp_rows.append({
                "short_name":    short,
                "gtdb_species":  genome["gtdb_species"],
                "accession":     genome["accession_clean"],
                "qc_label":      genome["qc_label"],
                "completeness":  genome["completeness"],
                "contamination": genome["contamination"],
                "ncbi_source":   genome["ncbi_source"],
                "ftp_url":       genome["ftp_url"],
            })

    # ------------------------------------------------------------------
    # Write TSV output
    # ------------------------------------------------------------------
    col_order = [
        "query_name", "spp_label", "domain",
        "n_hq", "n_mq", "n_eff", "hq_ratio",
        "n_total_passing_qc", "n_discarded",
        "n_gtdb_matched", "n_after_dedup",
        "n_refseq", "n_genbank_only",
        "dedup_summary", "gtdb_species_matched", "synonyms_checked", "notes",
    ]
    out = pd.DataFrame(results)[col_order]
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nReport written to: {args.output}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Write genome-level FTP output
    # Columns: short_name | gtdb_species | accession | qc_label |
    #          completeness | contamination | ncbi_source | ftp_path
    # Only genomes passing QC (HQ or MQ) are included.
    # Genomes with ftp_path == 'na' are included but flagged — they are
    # in GTDB but absent from the NCBI summary (suppressed/withdrawn).
    # ------------------------------------------------------------------
    ftp_output = args.output.replace(".tsv", "_genome_ftp.tsv")
    if ftp_rows:
        ftp_df = pd.DataFrame(ftp_rows)[[
            "short_name", "gtdb_species", "accession",
            "qc_label", "completeness", "contamination",
            "ncbi_source", "ftp_url",
        ]]
        ftp_df = ftp_df.sort_values(["short_name", "qc_label", "accession"])
        ftp_df.to_csv(ftp_output, sep="\t", index=False)
        print(f"Genome FTP list written to: {ftp_output}", file=sys.stderr)
    else:
        print("No passing-QC genomes found; FTP file not written.", file=sys.stderr)

    # ------------------------------------------------------------------
    # Print summary to stdout
    # ------------------------------------------------------------------
    width = 36
    print("\n" + "=" * 90)
    print(
        f"{'QUERY NAME':<{width}} {'LABEL':<28} "
        f"{'N_HQ':>5} {'N_MQ':>5} {'N_EFF':>7} {'HQ_RATIO':>9} {'DOMAIN':<12}"
    )
    print("=" * 90)
    for _, row in out.iterrows():
        print(
            f"{row['query_name']:<{width}} {row['spp_label']:<28} "
            f"{row['n_hq']:>5} {row['n_mq']:>5} {row['n_eff']:>7.1f} "
            f"{row['hq_ratio']:>9.3f} {row['domain']:<12}"
        )
    print("=" * 90)

    label_counts = out["spp_label"].value_counts()
    print("\nLabel distribution:")
    for lbl, cnt in label_counts.items():
        print(f"  {lbl:<35} n = {cnt}")


if __name__ == "__main__":
    main()
