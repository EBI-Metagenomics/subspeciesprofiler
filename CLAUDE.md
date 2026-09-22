# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

`ebi-metagenomics/subspeciesprofiler` — a Nextflow DSL2 pipeline (built from the nf-core template 3.6.0.dev0, but **not** an nf-core pipeline: `.nf-core.yml` sets `is_nfcore: false`) that finds **subspecies clusters within a species**. Input is one row per _species_ — a directory of assembled genomes plus a per-genome completeness/contamination table. For each species it builds a PopPUNK database, fits **many** clustering models, scores every fit against all-vs-all FastANI, and reports which model(s) best recover ANI-defined subspecies structure.

The core idea: this is a **model profiler, not a single-model clustering run**. There is no "the" model — every family is swept, every fit is scored, and a species that no model resolves is itself a result.

Two deliberate exclusions, both from PI review: PopPUNK's **lineage** model is absent (it sub-subclusters _within_ a strain, a different question), and **refinement is applied to dbscan only** — the non-standard `--multi-boundary` / `--unconstrained` variants are not used. Do not re-add either without checking first.

## Commands

```bash
# Run on the bundled two-species fixture
nextflow run . -profile test,docker --outdir results

# nf-test: everything, one file, one test, or by tag.
# A container profile is REQUIRED even for -stub runs: modules capture versions with
# eval() in the output block, which executes in stub mode too ('poppunk: command not found').
nf-test test --profile test,docker
nf-test test modules/local/poppunk/fitmodel/tests/main.nf.test --profile docker
nf-test test --tag poppunk/fitmodel --profile docker
nf-test test tests/default.nf.test --profile test,docker --update-snapshot

# Python unit tests for the evaluator (36 tests; not run by CI)
python3 -m pytest tests/bin/test_evaluate_poppunk_fastani.py

# Lint (both run in CI on every PR)
prek run --all-files          # or: pre-commit run --all-files
nf-core pipelines lint    # smaller rule set now that is_nfcore: false

# After adding/removing a param in nextflow.config
nf-core pipelines schema build
```

### Local-environment gotchas (hard-won; see the comments in `modules/local/poppunk/*/tests/nextflow.config`)

- **`-profile conda`/`mamba` hangs**: `mamba env create` blocks on a stdin prompt when spawned by Nextflow. The PopPUNK module test configs set `conda.useMicromamba = true` to work around it. Any new conda-based module test needs the same.
- **Apple Silicon + PopPUNK**: bioconda's `osx-arm64` graph-tool crashes every `poppunk --fit-model` right after the fit (`No static implementation was found…`). Export `CONDA_SUBDIR=osx-64` **in the shell** before `nf-test`/`nextflow` — `conda.createOptions` is _not_ forwarded by Nextflow, so it cannot be fixed from a config file. Do not set this on Linux/x86 CI.
- **Apple Silicon + Docker**: the `poppunk` biocontainer is `linux/amd64` only. Under emulation on an arm64 Mac, `poppunk --create-db` **segfaults** (exit 139) inside pp-sketchlib partway through distance calculation. So the real-data PopPUNK module tests and the full `tests/default.nf.test` cannot pass locally on Apple Silicon under any engine — stub tests are fine. Verify those on a Linux/x86 runner.
- **PopPUNK 2.7.8 flags to avoid**: `--qc-keep` triggers an `UnboundLocalError` in `remove_qc_fail` whenever any sample fails QC; `--auto-max-dists` does not exist. Both are documented at the `POPPUNK_QCDB` selector in `conf/modules.config`.

## Architecture

### Data contract

`assets/schema_input.json` defines the samplesheet: `species` (→ `meta.id`), `genomes_dir`, `qc_csv`. The `qc_csv` is `genome,completeness,contamination` where `genome` is the **assembly filename with extension**; stripping that extension yields the PopPUNK sample name. `bin/evaluate_poppunk_fastani.py:normalise_genome_id` is the single place that reconciles FastANI paths with PopPUNK sample names — change it there, nowhere else.

### Flow

`main.nf` (boilerplate) → `workflows/subspeciesprofiler.nf` → `SPECIESQC` gate → `POPPUNK_METHODS` subworkflow → MultiQC.

1. **`SPECIESQC`** (`modules/local/speciesqc`, runs `bin/spp_eligibility_from_qc.py`) classifies genomes HQ/MQ/DISCARDED and labels the _species group_ (`STRONG`, `ACCEPTABLE`, …). It emits three things the rest of the pipeline depends on: the eligibility report, an **r-file** (`name<TAB>./genomes/<file>`), and a **labels CSV** (`genome,label`).
2. `workflows/subspeciesprofiler.nf` branches on `params.qc_filter`; species whose label isn't in the list are dropped with a `log.warn` and never reach PopPUNK.
3. **`POPPUNK_METHODS`** (`subworkflows/local/poppunk_methods/main.nf`) does the real work — see below. Its `take:` contract expects `meta.spp_label` and `meta.n_genomes`, both folded on from the `SPECIESQC` report.

**r-file staging contract**: the r-file's second column holds _relative_ `./genomes/<file>` paths, so any process consuming it must stage genomes with `path(genomes, stageAs: 'genomes/*')`. `POPPUNK_CREATEDB` and `FASTANI_ALLVSALL` both do. This is why `FASTANI_ALLVSALL` is a local module rather than the nf-core one (that one needs absolute paths).

### The model sweep (`POPPUNK_METHODS`)

Stage 0 runs once per species: `POPPUNK_CREATEDB` → `POPPUNK_QCDB` → `POPPUNK_QUANTILES` (derives the threshold sweep from the _observed_ core-distance distribution via `bin/poppunk_core_quantiles.py`), plus `FASTANI_ALLVSALL` in parallel (model-independent).

Then **three families sweep**: `threshold` (one fit per quantile), `bgmm` (sweeps `--K`), `dbscan` (sweeps the `D` × `min-cluster-prop` grid). Grids come from the `params.poppunk_*` settings.

**BGMM is size-gated.** It needs clearly separated distance components, which only holds for small collections, so the whole family is skipped for species with `meta.n_genomes >= params.poppunk_bgmm_max_genomes` (default 200) and dbscan carries the fit. `meta.n_genomes` is the post-QC count (`n_hq + n_mq`, i.e. `n_total_passing_qc`) folded onto meta in `workflows/subspeciesprofiler.nf` — **not** `n_eff`, which is a weighted score rather than a genome count. A caller that supplies no `n_genomes` runs BGMM rather than silently losing the family.

**Every dbscan fit is then refined** with PopPUNK's standard `--fit-model refine`, seeded from that fit's directory. This is ungated — it runs whether or not anything already scored `ACCEPT`. Refined fits are evaluated and ranked alongside the swept ones; nothing else is refined.

**Microreact visuals** are generated for every fit the evaluator rates **Strong or Moderate** (deliberately _not_ `decision == ACCEPT`, which is Strong-only — the point is to compare all credible models by eye). `POPPUNK_VISUALISE` wraps `poppunk_visualise --microreact`, on the same poppunk container (`rapidnj` and `mandrake` are bundled in the bioconda package).

Two traps in that module: PopPUNK resolves a model as `<model-dir>/<basename(model-dir)>_fit.pkl`, so **the fit and db directories must never be renamed via `stageAs`**; and it writes its output into a `<prefix>/` directory whose files are all prefixed with that basename, which the module then lifts into the task dir so everything publishes flat into one per-species `poppunk/microreact/` folder with the model name carried in each filename.

**One module, aliased.** `POPPUNK_REFINE_DBSCAN` is `POPPUNK_FITMODEL` aliased, and `POPPUNK_EVALUATE_DBSCAN_REFINE` is `POPPUNK_EVALUATE` aliased. The model spec is a _per-task value input_ (`fit_args`, e.g. `"bgmm --K 4"`, `"refine"`), **not** `ext.args` — that's what lets one process cover the whole sweep. `meta.model` labels each fit and drives both the join-back and the publish path. Follow this pattern rather than adding a new process.

### The evaluator is the pipeline's judgement

`bin/evaluate_poppunk_fastani.py` (wrapped by `modules/local/poppunk/evaluate`) is where all scientific decisions live. It emits `tool_metrics.tsv`, `cluster_metrics.tsv`, `genome_metrics.tsv`. Two fields drive control flow:

- **`tool_status`** — audit label: `Strong` / `Moderate` / `Mixed` / `Weak`, from HQ and total structure scores plus singleton / tiny-cluster / negative-silhouette rates.
- **`decision`** — binary: `ACCEPT` if `tool_status` ∈ `--accept-status` (default `Strong`), else `TRY_NEXT_MODEL`. `ch_accepted` reports every ACCEPT row; `ch_best` is the single top-ranked fit and is the fallback when nothing accepts.

Scores are **HQ-centric**: `tool_structure_score_HQ` is `NaN` when no HQ genome sits in a scorable cluster, which is forced to `Weak` so it ranks last. The 36 pytest cases in `tests/bin/` pin this behaviour (FastANI percent scale, within-species ANI gate, singleton exclusion, degenerate inputs) — run them after touching the evaluator.

### Failure policy

Fit and evaluate processes use `errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'ignore' }` in `conf/modules.config`. A degenerate grid point (bad `--K`/`--D`) fails _deterministically_, so it's ignored — the fit simply drops out of the ranking instead of killing the species. Only transient/resource failures retry (escalating via `task.attempt` in `conf/base.config`). Preserve this distinction.

### Config layering

`nextflow.config` holds params and engine profiles → `conf/base.config` (resource labels) → institutional configs → `conf/modules.config` (publish paths, `ext.args`, error strategies). Adding a param means editing `nextflow.config` **and** `nextflow_schema.json` (`nf-core pipelines schema build`). Sweep grids are comma-separated strings, tokenised in the subworkflow.

Output is organised **species-first**: `<outdir>/<species>/poppunk/<stage>/`, with `fitmodel/` and `evaluate/` further split by `meta.model`. `publishDir` publishes whole directories so PopPUNK's diagnostic plots come along.

Versions use **both** `ch_versions` and `Channel.topic("versions")`; local modules here emit topic tuples via `eval(...)` in the output block. Both merge into the MultiQC versions YAML.

## Conventions

- New tools go in `modules/local/<tool>/` with `main.nf` + `environment.yml` + `meta.yml` + `tests/`, a real-data test and a `-stub` test, and a stub block whose output filenames match the real ones (the subworkflow's stub test depends on this).
- **Stubs must produce varied, realistic output where downstream code branches on it.** `poppunk/evaluate`'s stub returns `Strong` for `refine_from_*`, `Moderate` for `bgmm_*` and `Weak` otherwise, precisely so stub runs exercise the ACCEPT and Microreact branches. A stub returning one constant verdict leaves those paths untested.
- Never hand-edit `modules/nf-core/` or `subworkflows/nf-core/` — vendored and SHA-pinned in `modules.json`; use `nf-core modules install/update`.
- PRs target `dev`, not `master`.
- New tool → add to `CITATIONS.md` and `docs/output.md`.

## Planned work (do NOT "clean these up")

- **`checkm2/predict`** is installed in `modules.json` but unwired _on purpose_. It will be wired behind a flag that generates completeness/contamination when the samplesheet's `qc_csv` is absent. Separate task, not yet started.
- **`drep/dereplicate`** is likewise installed and unwired _on purpose_. It will pick the representative genome per cluster to use as the SynTracker reference.
- **SynTracker** (non-PopPUNK synteny signal) is the **next task to pick up**. Planned as an independent branch, with `drep/dereplicate` selecting the representative genome per cluster as its reference. It is not on bioconda and not in nf-core/modules, so it needs a custom `microbiome-informatics/syntracker:1.4.0` container (BLAST+ plus python, r-base=4.0.5, bioconductor-decipher, r-tidyverse, and the GitHub source tree). Full detail in the plan file `~/.claude/plans/we-are-going-to-abstract-hartmanis.md`.

## Known loose ends

- Test data is a local fixture (`tests/fixtures/`, `conf/test.config`). Now that this isn't an nf-core pipeline it no longer _needs_ to move to nf-core/test-datasets, but `conf/test.config` still carries a TODO saying it should, and `tests/nextflow.config` points `pipelines_testdata_base_path` at a `subspeciesprofiler` branch the fixtures don't use.
- `ro-crate-metadata.json` is generated — `nf-core pipelines lint` refreshes its README snapshot automatically. Never hand-edit it.
- Logo assets are still named `nf-core-subspeciesprofiler_logo_*.png` (in `assets/` and `docs/images/`) and referenced by that name from `README.md`. Renaming needs new artwork, so it was left alone.
- Template `TODO nf-core:` markers remain in `README.md`, `conf/base.config`, `conf/test_full.config`, `docs/usage.md`, `assets/methods_description_template.yml` and the `toolCitationText`/`toolBibliographyText` stubs in the local utils subworkflow. They need real content, not deletion.
- `bin/spp_species_eligibility.py` (GTDB/RefSeq-wide species screening) is not called by any module; it's an upstream/offline companion to `spp_eligibility_from_qc.py`.
