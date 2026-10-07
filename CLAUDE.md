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

# Python unit tests for the evaluator, the QC script and the APSS clustering (not run by CI;
# the SynTracker tests need python-igraph + matplotlib and are skipped without them)
python3 -m pytest tests/bin

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

`tool_metrics.tsv` column order is `model`, `poppunk_network_score`, `tool_status`, `n_clusters`, `n_nonsingleton_clusters`, `largest_cluster_fraction`, the scores, `decision`, `eval_summary`. **Never feed `poppunk_network_score` into `tool_status`.** It is `transitivity × (1 − density)`, so it peaks for over-fragmented fits (hundreds of singletons plus tiny near-identical cliques) that the evaluator must reject. A `Weak` fit that fails only on the singleton or tiny-cluster rules while its HQ structure score is ≥ 0.80 gets an `over-fragmentation:` prefix in `eval_summary`. That is a fit-level statement: "no structure" is a species-level conclusion (every fit fails). Don't reword it as strain-level; strains aren't defined here. The curation logic is documented in `docs/curation.md`. `poppunk_network_score` is PopPUNK's plain network `Score` (not the betweenness variants), parsed from the `<prefix>_fit.log` that `poppunk/fitmodel` saves into every fit dir. `*_model_report.tsv` (built in `workflows/subspeciesprofiler.nf`) is the concatenation of every fit's `tool_metrics.tsv` **minus `decision`**. `decision` stays in `tool_metrics.tsv` because control flow reads it.

Scores are **HQ-centric**: `tool_structure_score_HQ` is `NaN` when no HQ genome sits in a scorable cluster, which is forced to `Weak` so it ranks last. The evaluator's 44 pytest cases in `tests/bin/` pin this behaviour (FastANI percent scale, within-species ANI gate, singleton exclusion, degenerate inputs) — run them after touching the evaluator.

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

### The SynTracker branch (`SYNTRACKER_METHODS`)

`subworkflows/local/syntracker_methods/main.nf`, an independent non-PopPUNK signal, skipped with `--skip_syntracker`. SynTracker separates genomes by gene **arrangement**, so it can split genomes of near-identical ANI; the design avoids collapsing or labelling genomes by ANI before SynTracker sees them.

1. `DREP_DEREPLICATE` (**local** `modules/local/drep/dereplicate`, a copy of the nf-core module with a `genome_info` input; `ext.args '-sa 0.95'`) on the **HQ** genomes, with `--genomeInfo` from `SPECIESQC.out.genomeinfo` (`genome,completeness,contamination`, filename with extension). dRep only picks the **reference**: normally one representative per species; if several, the highest N50 from `data_tables/genomeInformation.csv` (`genome,completeness,contamination,length,N50,centrality`), with a `log.warn`.
2. Targets = **every HQ genome** (no dereplication, no propagation), capped at `--syntracker_max_targets` by N50 (reference always kept, `log.warn`). SynTracker's cost grows with the **square** of the targets (~570 core-hours for 588 B. longum genomes).
3. `SYNTRACKER_RUN` (`modules/local/syntracker/run`, container-only `quay.io/microbiome-informatics/syntracker:1.4.0_patch1`, **linux/amd64 only**). Always `-mode new`; the task fails if the all-regions APSS table has no pairs (SynTracker reports R failures only in its log).
4. `SYNTRACKER_CLUSTERS` (`bin/syntracker_apss_clusters.py`), the SynTracker paper's recipe: one depth (`--syntracker_regions auto` = highest N before targets/pairs retained drop below `--syntracker_min_retention` of the lowest N, with a retention TSV/PNG), one reference (asserted), `Compared_regions >= N` asserted on per-N tables (not a filter), an igraph graph of the targets pruned at each `--syntracker_min_apss`, Leiden (modularity, weights = raw APSS, seeded) at each `--syntracker_resolutions`. A target without edges is a singleton; non-targets (MQ, capped) are absent. Plus a per-cluster QC table.
5. `SYNTRACKER_EVALUATE` = `POPPUNK_EVALUATE` aliased; rows `syntracker_leiden_n<N>_apss<t>_r<res>` join the model report. FastANI comes from `POPPUNK_METHODS.out.ani`.

Hard-won SynTracker facts (do not relearn them):

- **The R stack must stay pinned** (R 4.0.5, DECIPHER 2.18.1, RSQLite). An unpinned current DECIPHER fails `Seqs2DB` on every region (`N function calls resulted in an error`) and SynTracker then crashes in `left_join()`. The image is built from the EBI containers repo, `syntracker/1.4.0_patch1`.
- **APSS is not on the ANI scale.** On B. longum (earlier method: average linkage on all-regions APSS of dRep 0.99 representatives), subspecies separated at APSS 0.72-0.80 and 0.90 fragmented the species; that reproduced dbscan's 3 groups (ARI 0.98). Single linkage chained through a few high-APSS bridging pairs. The Leiden method has not been validated on real data yet.
- **Never use SynTracker's `-mode continue` in the pipeline.** It loads every finished region into the R parent before forking workers, so continue runs OOM where fresh runs don't; it also renames a region `_done` _before_ saving its result, so a kill can silently lose regions.
- Sample names: SynTracker names a target by its file basename minus extension; the module stages every target as `<sample>.fasta`, so names reconcile with `normalise_genome_id`.
- `normalise_genome_id` is imported by `bin/syntracker_apss_clusters.py` from `bin/evaluate_poppunk_fastani.py` (still the single place). The clusters env pins **pandas < 3**: the evaluator breaks under pandas 3's copy-on-write (`underlying array is read-only` in `test_check_species_ani_all_close_passes`).
- `procps-ng` is linux-only on conda-forge, so the clusters `environment.yml` cannot be solved on macOS; tests there need a local env with pandas, python-igraph and matplotlib.

Follow-up (documented, not implemented): after pruning, every edge weight lies in `[min_apss, 1]`, so modularity is driven mostly by which edges exist. Rescaling to `(APSS − min_apss)/(1 − min_apss)` restores the dynamic range; planned as a flag, raw APSS stays the default.

## Planned work (do NOT "clean these up")

- **`checkm2/predict`** is installed in `modules.json` but unwired _on purpose_. It will be wired behind a flag that generates completeness/contamination when the samplesheet's `qc_csv` is absent. Separate task, not yet started.

## Known loose ends

- Test data is a local fixture (`tests/fixtures/`, `conf/test.config`). Now that this isn't an nf-core pipeline it no longer _needs_ to move to nf-core/test-datasets, but `conf/test.config` still carries a TODO saying it should, and `tests/nextflow.config` points `pipelines_testdata_base_path` at a `subspeciesprofiler` branch the fixtures don't use.
- `ro-crate-metadata.json` is generated — `nf-core pipelines lint` refreshes its README snapshot automatically. Never hand-edit it.
- Logo assets are still named `nf-core-subspeciesprofiler_logo_*.png` (in `assets/` and `docs/images/`) and referenced by that name from `README.md`. Renaming needs new artwork, so it was left alone.
- Template `TODO nf-core:` markers remain in `README.md`, `conf/base.config`, `conf/test_full.config`, `docs/usage.md`, `assets/methods_description_template.yml` and the `toolCitationText`/`toolBibliographyText` stubs in the local utils subworkflow. They need real content, not deletion.
- `bin/spp_species_eligibility.py` (GTDB/RefSeq-wide species screening) is not called by any module; it's an upstream/offline companion to `spp_eligibility_from_qc.py`.
