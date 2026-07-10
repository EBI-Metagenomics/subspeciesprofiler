<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-subspeciesprofiler_logo_dark.png">
    <img alt="nf-core/subspeciesprofiler" src="docs/images/nf-core-subspeciesprofiler_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/subspeciesprofiler)
[![GitHub Actions CI Status](https://github.com/nf-core/subspeciesprofiler/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/subspeciesprofiler/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/subspeciesprofiler/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/subspeciesprofiler/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/subspeciesprofiler/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.6.0.dev0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.6.0.dev0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/subspeciesprofiler)

## Introduction

**ebi-metagenomics/subspeciesprofiler** is a bioinformatics pipeline that generates subspecies clusters from a group of genomes of the same species.

Given, per species, a directory of assembled genomes and a per-genome completeness/contamination table, for each eligible species the pipeline:

1. **Classifies** each genome as HQ / MQ / DISCARDED from its completeness/contamination and decides whether the species group is eligible for subspecies clustering (`speciesqc`).
2. Builds and QCs a **[PopPUNK](https://poppunk.bacpop.org)** database from the HQ/MQ genomes and derives a **core-distance threshold sweep** from the observed distance distribution.
3. Fits PopPUNK clustering models and **scores each against all-vs-all [FastANI](https://github.com/ParBLiSS/FastANI)**, then selects the model(s) that best recover the ANI-defined subspecies structure.

Outputs are a per-species set of PopPUNK databases, fitted models (with diagnostic plots), ANI-based scores, and a [MultiQC](http://multiqc.info/) report.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet describing the species to cluster — one row per species:

`samplesheet.csv`:

```csv
species,genomes_dir,qc_csv
bacteroides_xylanisolvens,/path/to/bxylanisolvens/genomes,/path/to/bxylanisolvens/qc.csv
```

| Column        | Description                                                                                                 |
| ------------- | ----------------------------------------------------------------------------------------------------------- |
| `species`     | Unique species name (no spaces).                                                                            |
| `genomes_dir` | Directory of assembled genomes for that species. FASTA (`.fasta`/`.fa`/`.fna`), optionally gzipped (`.gz`). |
| `qc_csv`      | Per-genome completeness/contamination table for that species (see below).                                   |

The `qc_csv` has a header row and one row per genome. The `genome` column is the **assembly filename with its extension** (it is stripped to the PopPUNK sample name), and the completeness/contamination columns are matched case-insensitively by name (so CheckM/CheckM2-style headers work directly):

```csv
genome,completeness,contamination
MGYG000001345.fasta,99.5,0.5
MGYG000005313.fasta,98.2,1.1
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/subspeciesprofiler \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!NOTE]
> On x86_64 Linux clusters, use `-profile singularity` (or `docker`); no special settings are needed. The PopPUNK, FastANI and pandas steps all ship pinned containers.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/subspeciesprofiler/usage) and the [parameter documentation](https://nf-co.re/subspeciesprofiler/parameters).

## How it works

### Genome classification and species eligibility

Each genome is classified from its completeness (`comp`) and contamination (`cont`):

| Class         | Rule                                     |
| ------------- | ---------------------------------------- |
| **HQ**        | `comp ≥ 90` and `cont ≤ 1%`              |
| **MQ**        | `comp ≥ 80` and `cont ≤ 5%`              |
| **DISCARDED** | otherwise (dropped; not sent to PopPUNK) |

Only HQ + MQ genomes are used to build the PopPUNK database. Each species is then given an eligibility label from its HQ / effective-genome counts (`STRONG`, `ACCEPTABLE`, `NOT_ELIGIBLE`, `DISCARDED`, with `_BOUNDARY` / `_NEAR_THRESHOLD` variants). Only species whose label is in `--qc_filter` (default `STRONG,STRONG_BOUNDARY,ACCEPTABLE,ACCEPTABLE_BOUNDARY`) proceed to clustering.

### Model scoring and selection

Candidate PopPUNK models are scored against all-vs-all FastANI:

- **Within-species ANI gate** — every genome pair in a species group must be ≥ `--min-ani` (0–1 scale, default `0.90`) similar. A pair below the threshold, or one FastANI cannot align, **fails the species**: it flags a genome too distant to belong (e.g. a mis-assignment during data preparation).
- **`tool_status`** ∈ `Strong` / `Moderate` / `Weak` / `Mixed` — the clustering quality, weighted towards the trustworthy HQ genomes (cluster cohesion + separation, penalising over-splitting).
- **`decision`** ∈ `ACCEPT` / `TRY_NEXT_MODEL` — `ACCEPT` when `tool_status` is in `--accept-status` (default `Strong`).

**All accepted models are reported** (e.g. several `Strong` fits from one sweep), together with the single best-ranked model as a fallback.

## Pipeline output

Results are organised per species under `<outdir>/`:

| Path                                  | Contents                                                                                                                                    |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `speciesqc/<species>/`                | Eligibility report + the PopPUNK r-file and per-genome HQ/MQ labels.                                                                        |
| `poppunk/<species>/createdb/`         | The PopPUNK sketch database + core/accessory distances (with diagnostic plots).                                                             |
| `poppunk/<species>/qcdb/`             | The QC'd (pruned) database (with plots).                                                                                                    |
| `poppunk/<species>/quantiles/`        | The data-derived core-distance thresholds (`*_core_quantiles.csv`).                                                                         |
| `poppunk/<species>/fastani/`          | The all-vs-all ANI (`*.ani.txt`).                                                                                                           |
| `poppunk/<species>/fitmodel/<model>/` | Each fitted PopPUNK model — cluster assignments (`*_clusters.csv`) and fit plots — one dir per model.                                       |
| `poppunk/<species>/evaluate/<model>/` | The per-model verdict: `*.tool_metrics.tsv` (`model`, `tool_status`, `decision`, `reason`, scores) plus per-cluster and per-genome metrics. |

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/subspeciesprofiler/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/subspeciesprofiler/output).

## Credits

nf-core/subspeciesprofiler was originally written by Alejandra Escobar.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#subspeciesprofiler` channel](https://nfcore.slack.com/channels/subspeciesprofiler) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/subspeciesprofiler for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
