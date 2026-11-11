# ebi-metagenomics/subspeciesprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [date]

The input samplesheet is unchanged. The outputs change: new pangenome and ANI-tree outputs, a new species QC table, FastANI in a new folder, and a contiguity filter that can remove genomes from every branch.

### `Added`

- Pangenome branch (`PANGENOME`, skip with `--skip_pangenome`), over the genomes that passed the species QC:
  - ggCaller 1.5.2 gene calls (`GGCALLER`). Its database comes from `--ggcaller_db`, or `GGCALLER_DOWNLOADDB` downloads it once per run (`--save_ggcaller_db` publishes it to `databases/ggcaller_db/`).
  - Panaroo (nf-core `panaroo/run`, `--clean-mode sensitive --refind-mode off --remove-invalid-genes`) makes `<species>/pangenome/panaroo/gene_presence_absence.Rtab` and the other Panaroo results.
  - `<species>/pangenome/ggcaller/` holds one gzipped GFF (`gff/<sample>.gff.gz`), protein (`faa/<sample>.faa.gz`) and CDS (`ffn/<sample>.ffn.gz`) file per genome. ggCaller writes one file for all genomes; the pipeline splits it.
  - Genome names in all these files are the PopPUNK sample names.
- ANI tree (`ANI_TREE`, `bin/fastani_to_tree.py`): a neighbour-joining Newick tree of each species' FastANI distances (1 − ANI/100), `<species>/<species>_ani_nj_tree.nwk`.
- Species QC outputs: `<species>/speciesqc/<species>_genome_qc.tsv` (every genome with its completeness, contamination, N50, N90, quality, label and reason) and `<species>_seqkit_stats.tsv`. New column `n_discarded_contiguity` in the eligibility reports.
- `bin/ordination_enrichment.py`, a post-pipeline script (not run by the pipeline) for the gene presence/absence table and a clustering: a PCoA plot by cluster (singletons hollow), pairwise ANOSIM, Fisher gene enrichment, and an iTOL `TREE_COLORS` file for the ANI tree (`--tree` checks the node names against the tree).
- New parameters: `--qc_min_nx`, `--qc_strict_contiguity`, `--skip_pangenome`, `--ggcaller_db`, `--save_ggcaller_db`, `--panaroo_clean_mode`.

### `Changed`

- Contiguity filter in the species QC (`SPECIESQC`, with `SEQKIT_STATS` on every genome): an HQ or MQ genome is DISCARDED (reason `contiguity`) if its N50, or N90 with `--qc_strict_contiguity`, is below `--qc_min_nx` (default 5000 bp, the SynTracker window). It is the only genome filter: PopPUNK, FastANI, SynTracker and the pangenome all use the genomes that pass. The species eligibility label is computed after this filter, so a species can now fall below `--qc_filter`. Set `--qc_min_nx 0` to get the v1.0.0 genome set.
- SynTracker runs on the HQ and MQ genomes (v1.0.0: HQ only). The reference is still an HQ genome.
- Curation guide: read the SynTracker rows by `ani_separated_fraction` too, not only by `tool_status`. Kateryna Pantiukh found on her data that the row with the highest `ani_separated_fraction` reproduced her curated subspecies, although it was rated `Weak` or `Moderate`.
- FastANI output moved to `<species>/fastani/` (was `<species>/poppunk/fastani/`), because it does not depend on PopPUNK.

### `Fixed`

### `Dependencies`

| Tool      | New version |
| --------- | ----------- |
| ggCaller  | 1.5.2       |
| Panaroo   | 1.6.0       |
| Biopython | 1.84        |

- The nf-core `panaroo/run` module is patched (`modules/nf-core/panaroo/run/panaroo-run.diff`) to decompress gzipped GFF inputs, because Panaroo reads only uncompressed `.gff`.

### `Deprecated`

## v1.0.0

Initial release of ebi-metagenomics/subspeciesprofiler, created with the [nf-core](https://nf-co.re/) template.
