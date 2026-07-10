# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`nf-core/subspeciesprofiler` is a Nextflow (DSL2) bioinformatics pipeline built from the nf-core template (v3.6.0.dev0). Its intended purpose is to **generate subspecies clusters from a group of genomes of the same species**. The repository is currently a fresh template scaffold: the only real analysis step wired up so far is FastQC → MultiQC. Most files still contain `TODO nf-core:` markers where pipeline-specific logic is meant to go.

Requires Nextflow `>=25.04.0`.

## Common commands

Run the pipeline (test profile, needs a container engine):

```bash
nextflow run . -profile test,docker --outdir results
```

Run on real data:

```bash
nextflow run . -profile docker --input samplesheet.csv --outdir results
```

Available engine profiles: `docker`, `singularity`, `conda`, `mamba`, `podman`, `shifter`, `charliecloud`, `apptainer`, `wave`, `gpu`. Data profiles: `test` (minimal), `test_full`. Combine as `-profile test,docker`. **Parameters must be passed via CLI or `-params-file`, never via `-c`** (custom configs are for infrastructure config only).

Testing (nf-test, config in `nf-test.config`):

```bash
nf-test test tests/default.nf.test              # the pipeline-level snapshot test
nf-test test tests/default.nf.test --profile docker
nf-test test --update-snapshot tests/default.nf.test   # regenerate snapshots after intended output changes
```

`nf-test.config` ignores tests under `modules/nf-core/**` and `subworkflows/nf-core/**`, and auto-triggers a full run when core config files change.

Linting / formatting:

```bash
nf-core pipelines lint     # nf-core structural + template linting
pre-commit run --all-files # prettier + whitespace/EOF fixers
```

## Architecture

Nextflow entrypoints form a layered call chain — read them in this order:

1. **`main.nf`** — top-level `workflow {}`. Calls `PIPELINE_INITIALISATION` (validates params, parses samplesheet), then `NFCORE_SUBSPECIESPROFILER`, then `PIPELINE_COMPLETION` (email/Slack notifications). Also sets `params.fasta` from igenomes via `getGenomeAttribute('fasta')`.
2. **`workflows/subspeciesprofiler.nf`** — the `SUBSPECIESPROFILER` workflow: the actual science. Takes the samplesheet channel, runs modules, and aggregates software versions + MultiQC report. **New analysis steps go here.**
3. **`subworkflows/local/utils_nfcore_subspeciesprofiler_pipeline/main.nf`** — pipeline-specific glue: `PIPELINE_INITIALISATION`, `PIPELINE_COMPLETION`, samplesheet parsing (`samplesheetToList` against `assets/schema_input.json`), `validateInputSamplesheet`, and `methodsDescriptionText`.
4. **`subworkflows/nf-core/*`** and **`modules/nf-core/*`** — vendored, unmodified nf-core components. Do **not** hand-edit these; manage them with `nf-core modules`/`nf-core subworkflows` commands (tracked in `modules.json`).

### Samplesheet channel shape

`PIPELINE_INITIALISATION` emits `ch_samplesheet` as tuples of `[ meta, [fastqs] ]` where `meta` carries `id` and `single_end`. `validateInputSamplesheet` enforces that all runs of a sample share the same single-end/paired-end datatype. The input schema lives in `assets/schema_input.json`.

### Adding a module and configuring it

- Install with `nf-core modules install <tool>`, then `include { ... }` it in `workflows/subspeciesprofiler.nf` and wire its channels.
- Feed each module's `.out.versions` into the `ch_versions` mix so it reaches the software-versions YAML.
- To collect a module's output into the MultiQC report, mix it into `ch_multiqc_files`.
- Per-process CLI args (`ext.args`), output filename prefixes (`ext.prefix`), and `publishDir` paths are set in **`conf/modules.config`** via `withName:` blocks — not inside the module files. Default publish path derives from the process name.

### Config layering

`nextflow.config` is the root; it pulls in `conf/base.config` (resource labels/retries), `conf/modules.config` (per-module options), and the selected profile from `conf/test.config` / `conf/test_full.config`. Parameters are declared and JSON-schema-validated against `nextflow_schema.json` (edit it with `nf-core pipelines schema build`). igenomes reference data is defined in `conf/igenomes.config`.

## Conventions

- This is a template scaffold — search for `TODO nf-core:` to find the spots intended for pipeline-specific customization (README intro, test data paths in `conf/test.config`, citations, etc.).
- Test data currently points at borrowed `viralrecon` samplesheets in `conf/test.config`; replace with real subspeciesprofiler test data before relying on the test profile.
- After changing anything that affects pipeline outputs, regenerate the nf-test snapshot (`--update-snapshot`) — the pipeline test asserts against a stored snapshot of all output file names and stable file contents.
- Keep `modules/nf-core/` and `subworkflows/nf-core/` in sync via nf-core tooling rather than manual edits; prettier/whitespace hooks explicitly exclude them.
