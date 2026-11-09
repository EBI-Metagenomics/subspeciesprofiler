# ebi-metagenomics/subspeciesprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of ebi-metagenomics/subspeciesprofiler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- SynTracker contiguity gate (`SYNTRACKER_GATE`): HQ and MQ genomes go to SynTracker only if their N50, or N90 with `--syntracker_strict`, is at least `--syntracker_min_nx` (default 5000 bp, the SynTracker window). MQ genomes are now SynTracker targets; the reference is still HQ. Report: `<species>/syntracker/selection/<species>_syntracker_gate.tsv`.

### `Fixed`

### `Dependencies`

### `Deprecated`
