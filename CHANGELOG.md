# ebi-metagenomics/subspeciesprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of ebi-metagenomics/subspeciesprofiler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Contiguity in the species QC (`SPECIESQC`, with `SEQKIT_STATS` on every genome): an HQ or MQ genome is DISCARDED if its N50, or N90 with `--qc_strict_contiguity`, is below `--qc_min_nx` (default 5000 bp, the SynTracker window). Every branch uses the genomes that pass, and the eligibility label is computed from them. New outputs: `<species>/speciesqc/<species>_genome_qc.tsv`, `<species>_seqkit_stats.tsv`, and the `n_discarded_contiguity` report column.
- SynTracker now runs on the HQ and MQ genomes (before: HQ only); the reference is still HQ.
- Pangenome branch (`PANGENOME`, skip with `--skip_pangenome`): ggCaller 1.5.2 gene calls (`GGCALLER`, with its database from `--ggcaller_db` or downloaded by `GGCALLER_DOWNLOADDB`) clustered by Panaroo (nf-core `panaroo/run`, `--clean-mode sensitive --refind-mode off`) into `<species>/pangenome/panaroo/gene_presence_absence.Rtab`, over the genomes that passed the species QC. The gene-call proteins are published in `<species>/pangenome/ggcaller/`.

### `Fixed`

### `Dependencies`

### `Deprecated`
