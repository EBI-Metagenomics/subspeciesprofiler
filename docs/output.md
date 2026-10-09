# ebi-metagenomics/subspeciesprofiler: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Species eligibility](#species-eligibility) - Per-species QC classification and a combined cross-species report
- [Subspecies clustering](#subspecies-clustering) - PopPUNK database, model sweep, FastANI evaluation and per-model report
- [SynTracker](#syntracker) - Synteny-based clustering, an independent signal scored like the PopPUNK fits
- [MultiQC](#multiqc) - Aggregate subspecies dashboard for the whole run
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Results are organised **species-first**: each species in the samplesheet gets its own top-level directory, so multi-species runs stay easy to navigate.

```
<outdir>/
├── species_eligibility_report.tsv        # combined eligibility, one row per species
├── <species>/
│   ├── speciesqc/                         # this species' eligibility report + QC classification
│   ├── poppunk/
│   │   ├── <species>_model_report.tsv     # every PopPUNK model tried + its verdict
│   │   ├── createdb/  qcdb/  quantiles/
│   │   ├── fastani/
│   │   ├── fitmodel/<model>/              # incl. PopPUNK diagnostic plots (*.png)
│   │   ├── evaluate/<model>/
│   │   └── microreact/                    # Strong + Moderate models, one file set each
│   └── syntracker/
│       ├── <species>_syntracker_model_report.tsv  # every SynTracker clustering + its verdict
│       ├── selection/                     # reference and target choice
│       ├── apss/                          # SynTracker APSS tables + its log
│       ├── clusters/                      # one Taxon,Cluster table per min_apss x resolution, coverage, noise
│       └── evaluate/<model>/              # APSS (+ ANI concordance) evaluation of each clustering
├── multiqc/
└── pipeline_info/
```

### Species eligibility

<details markdown="1">
<summary>Output files</summary>

- `species_eligibility_report.tsv` (top level): one row per species — `species_name`, `spp_label`, `n_hq`, `n_mq`, `n_eff`, `hq_ratio` — appended across every species in the samplesheet.
- `<species>/speciesqc/`: the per-species eligibility report and the derived PopPUNK inputs (r-file + HQ/MQ labels).

</details>

Each species' genome QC table is classified into HQ / MQ genomes and assigned a subspecies-clustering eligibility label. Only species whose label is kept by `--qc_filter` continue to clustering.

### Subspecies clustering

<details markdown="1">
<summary>Output files</summary>

- `<species>/poppunk/<species>_model_report.tsv`: the "profiler history" — every PopPUNK model fitted for this species with:
  - PopPUNK's own network score (`poppunk_network_score`, the plain score from the fit's network summary). This is context only: it peaks for the most fragmented fits, so it does not feed the verdict;
  - its evaluation verdict (`tool_status`);
  - the cluster counts (`n_clusters`, `n_nonsingleton_clusters`, `largest_cluster_fraction`);
  - the structure scores and a plain-text `eval_summary`.

  The binary `decision` is kept only in the per-fit `evaluate/<model>/*.tool_metrics.tsv`. See the [curation guide](curation.md) for how to read this table and choose a model.

- `<species>/poppunk/createdb/`, `qcdb/`, `quantiles/`: the PopPUNK database, its distance-QC output, and the core-distance quantiles that seed the threshold sweep.
- `<species>/poppunk/fastani/`: the all-vs-all FastANI distances used to score each model externally.
- `<species>/poppunk/fitmodel/<model>/`: per-model fit output, including PopPUNK diagnostic plots (`*.png`) and the run log (`*_fit.log`).
- `<species>/poppunk/evaluate/<model>/`: per-model evaluation metrics (`*.tool_metrics.tsv`, `*.cluster_metrics.tsv`, `*.genome_metrics.tsv`).
- `<species>/poppunk/microreact/`: [Microreact](https://microreact.org) visualisations for every model the evaluator rated **Strong or Moderate**, so the credible candidates can be inspected side by side and the best one picked by eye. All models share this one directory and each file carries the model name (family plus its swept parameters), e.g.:
  - `<species>_<model>.microreact` — upload this single file to [microreact.org/upload](https://microreact.org/upload) to open the visualisation.
  - `<species>_<model>_microreact_clusters.csv` — cluster assignments with Microreact headers.
  - `<species>_<model>_core_NJ.nwk` — neighbour-joining tree from the core distances.
  - `<species>_<model>_perplexity<P>_accessory_mandrake.dot` — 2D embedding of the accessory distances.

  For example `E_coli_refine_from_dbscan_D5_mcp0.01.microreact` is the refined DBSCAN fit with `D=5`, `min-cluster-prop=0.01`.

</details>

[PopPUNK](https://poppunk.readthedocs.io/) builds a genome database and fits candidate clustering models, which are scored against all-vs-all [FastANI](https://github.com/ParBLiSS/FastANI) distances to pick a subspecies partition.

### SynTracker

<details markdown="1">
<summary>Output files</summary>

- `<species>/syntracker/<species>_syntracker_model_report.tsv`: one row per SynTracker clustering (`syntracker_leiden_apss<t>_r<res>`): its verdict in APSS space (`tool_status`, structure scores, silhouettes, singleton rates), the APSS standard error used as the gap margin (`apss_se`), the same partition scored in ANI space (`median_silhouette_ANI_HQ`, `negative_silhouette_ANI_HQ_fraction`, `ani_separated_fraction`) and the resulting `signal` (`synteny_and_ani`, `synteny_only` or `none`). It has its own columns, separate from the PopPUNK model report.
- `<species>/syntracker/selection/<species>_syntracker_selection.tsv`: every HQ genome with its N50, its centrality (mean FastANI to the other HQ genomes) and its role: `reference`, `target`, or `capped` (left out by `--syntracker_max_targets`). The `*_reference.txt` and `*_targets.txt` files list the chosen genomes.
- `<species>/syntracker/apss/`:
  - `avg_synteny_scores_all_regions.csv`: the average pairwise synteny score (APSS) of every pair of targets over all the regions they share (`Ref_genome, Sample1, Sample2, APSS, Compared_regions`). This is the table that is clustered and evaluated.
  - `avg_synteny_scores_<N>_regions.csv` (N = 40, 60, 80, 100, 200): APSS from N regions subsampled per pair; pairs with fewer than N regions are absent. SynTracker fixes the subsampling seed, so these are reproducible, but they are noisy (sampling noise about 0.23/sqrt(N)); only the 40-region table is used, to estimate that noise.
  - `SynTracker_log.txt`: SynTracker's own log, with each region's outcome.
- `<species>/syntracker/clusters/`:
  - `<species>_syntracker_leiden_apss<t>_r<res>_clusters.csv`: one `Taxon,Cluster` table per APSS floor `t` (`--syntracker_min_apss`) and Leiden resolution `res` (`--syntracker_resolutions`); cluster 1 is the largest. A target with no pair at or above `t` is a singleton.
  - `<species>_syntracker_genome_coverage.tsv`: each target's median number of regions per pair, relative to the species median, and its status; `low_coverage` targets (below `--syntracker_min_genome_coverage`) are left out of the clusterings.
  - `<species>_syntracker_noise.tsv`: the per-region synteny-score SD, the median regions per pair and the resulting APSS standard error.
- `<species>/syntracker/evaluate/<model>/`: the evaluation of each clustering (`*.tool_metrics.tsv`, `*.cluster_metrics.tsv` with `*_APSS` columns and Cliff's delta, `*.genome_metrics.tsv` with `silhouette_APSS`).

</details>

[SynTracker](https://github.com/leylabmpi/SynTracker) compares genomes by **synteny**, the arrangement of genes, in regions of a single reference genome, which makes it independent of PopPUNK's k-mer distances and able to separate genomes of near-identical ANI that differ by rearrangements. APSS is close to 1 minus the share of 5-kb windows that a pair has broken into more than one synteny block.

The reference is the most central HQ genome (highest mean FastANI to the other HQ genomes) among those with N50 of at least `--syntracker_ref_min_n50`: a target loses regions the further it is from the reference, and a central reference spreads that loss evenly across lineages. Every HQ genome is a target; SynTracker's cost grows with the square of their number, so above `--syntracker_max_targets` only the genomes with the highest N50 are kept. The all-regions APSS is turned into clusters as in the SynTracker paper's network analysis: an APSS-weighted graph of the targets, pruned at the APSS floor, partitioned into Leiden (modularity) communities. Each clustering is scored in APSS space, with a gap margin of twice the APSS standard error, and in ANI space for concordance. See the [curation guide](curation.md) for reading these rows.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

For this pipeline the MultiQC report is a run-level **subspecies dashboard**: it collects a cross-species eligibility table and the selected model per species (as custom-content tables), alongside the software versions used, for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
