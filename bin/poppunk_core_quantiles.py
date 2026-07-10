#!/usr/bin/env python3
"""
Core-distance quantiles from a PopPUNK database
===============================================
Derives the core-distance threshold sweep values for `poppunk --fit-model
threshold` from the *observed* distribution of a (QC'd) PopPUNK database, rather
than guessing arbitrary thresholds.

Input is the database's `<name>.dists.npy` — a `(n_pairs, 2)` float32 array where
column 0 is the core distance and column 1 the accessory distance (verified
against a real poppunk 2.7.8 database).

Writes a CSV `quantile,core_threshold`, one row per *usable* quantile. Quantiles
whose threshold is <= 0 (common at the low end: many genome pairs have a core
distance of exactly 0) are dropped, and duplicate thresholds are collapsed, so
the sweep only runs on the informative (>0, distinct) range.

Usage:
    poppunk_core_quantiles.py \
        --dists   pp_db_qc.dists.npy \
        --output  core_quantiles.csv \
        [--quantiles 0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.15,0.20]
"""

import argparse
import sys

import numpy as np

DEFAULT_QUANTILES = "0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.15,0.20"
# Round thresholds to this many decimals to stabilise deduplication and collapse
# near-zero values (core distances here span ~1e-5 to 1e-2, so 9 dp is ample).
ROUND_DP = 9


def main():
    parser = argparse.ArgumentParser(
        description="Compute core-distance quantiles from a PopPUNK <db>.dists.npy."
    )
    parser.add_argument("--dists", required=True,
                        help="PopPUNK <db>.dists.npy (n_pairs x 2; column 0 = core distance).")
    parser.add_argument("--output", required=True,
                        help="Output CSV: quantile,core_threshold.")
    parser.add_argument("--quantiles", default=DEFAULT_QUANTILES,
                        help=f"Comma-separated quantiles in [0,1] (default: {DEFAULT_QUANTILES}).")
    args = parser.parse_args()

    try:
        dists = np.load(args.dists)
    except FileNotFoundError:
        sys.exit(f"ERROR: distances file not found: '{args.dists}'.")
    except Exception as exc:  # noqa: BLE001 - surface any load error cleanly
        sys.exit(f"ERROR: could not load distances '{args.dists}': {exc}")

    if dists.ndim != 2 or dists.shape[0] == 0 or dists.shape[1] < 1:
        sys.exit(f"ERROR: expected a (n_pairs, 2) distances array in '{args.dists}', got shape {dists.shape}.")

    core = dists[:, 0].astype(float)
    core = core[~np.isnan(core)]
    if core.size == 0:
        sys.exit(f"ERROR: no (non-NaN) core distances in '{args.dists}'.")

    try:
        quantiles = [float(q) for q in args.quantiles.split(",") if q.strip()]
    except ValueError:
        sys.exit(f"ERROR: --quantiles must be comma-separated numbers, got '{args.quantiles}'.")
    for q in quantiles:
        if not 0.0 <= q <= 1.0:
            sys.exit(f"ERROR: quantile {q} is outside [0, 1].")

    # Compute thresholds; drop <= 0 (degenerate) and duplicates, keep sorted.
    rows = []
    seen = set()
    for q in sorted(quantiles):
        thr = round(float(np.quantile(core, q)), ROUND_DP)
        if thr <= 0.0 or thr in seen:
            continue
        seen.add(thr)
        rows.append((q, thr))

    if not rows:
        sys.exit(
            "ERROR: no usable (>0, distinct) core-distance thresholds from the requested quantiles; "
            "the core-distance distribution may be degenerate (many exact zeros). Try higher quantiles."
        )

    with open(args.output, "w") as out:
        out.write("quantile,core_threshold\n")
        for q, thr in rows:
            out.write(f"{q},{thr:.9g}\n")

    kept = ", ".join(f"{q}->{thr:.3g}" for q, thr in rows)
    print(f"Core-distance quantiles ({len(rows)} usable of {len(quantiles)} requested): {kept}",
          file=sys.stderr)
    print(f"Written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
