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
# the SynTracker clustering tests need python-igraph and are skipped without it)
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
4. **`ANI_PLOT`**, fed by `POPPUNK_METHODS`'s FastANI (`modules/local/ani/plot`, `bin/plot_ani.py`): one figure per species, `<species>/<species>_ani_overview.{png,pdf,tsv}`, with a pairwise ANI histogram plus KDE next to a UPGMA-clustered ANI heatmap with dendrograms and an HQ/MQ strip. It runs on the **poppunk image on purpose**: that is the only image in the pipeline with matplotlib + scipy (no seaborn, so the clustered heatmap is built from scipy + matplotlib). Its tests (`tests/bin/test_plot_ani.py`) skip without scipy, so run them in that image.

**r-file staging contract**: the r-file's second column holds _relative_ `./genomes/<file>` paths, so any process consuming it must stage genomes with `path(genomes, stageAs: 'genomes/*')`. `POPPUNK_CREATEDB` and `FASTANI_ALLVSALL` both do. This is why `FASTANI_ALLVSALL` is a local module rather than the nf-core one (that one needs absolute paths).

### The model sweep (`POPPUNK_METHODS`)

Stage 0 runs once per species: `POPPUNK_CREATEDB` → `POPPUNK_QCDB` → `POPPUNK_QUANTILES` (derives the threshold sweep from the _observed_ core-distance distribution via `bin/poppunk_core_quantiles.py`), plus `FASTANI_ALLVSALL` in parallel (model-independent).

Then **three families sweep**: `threshold` (one fit per quantile), `bgmm` (sweeps `--K`), `dbscan` (sweeps the `D` × `min-cluster-prop` grid). Grids come from the `params.poppunk_*` settings. `min-cluster-prop` 0.02 was dropped from the default (it failed for every species on the 2026-10 Codon run).

**BGMM is size-gated.** It needs clearly separated distance components, which only holds for small collections, so the whole family is skipped for species with `meta.n_genomes >= params.poppunk_bgmm_max_genomes` (default 200) and dbscan carries the fit. `meta.n_genomes` is the post-QC count (`n_hq + n_mq`, i.e. `n_total_passing_qc`) folded onto meta in `workflows/subspeciesprofiler.nf` — **not** `n_eff`, which is a weighted score rather than a genome count. A caller that supplies no `n_genomes` runs BGMM rather than silently losing the family.

**Every dbscan fit is then refined** with PopPUNK's standard `--fit-model refine`, seeded from that fit's directory. This is ungated — it runs whether or not anything already scored `ACCEPT`. Refined fits are evaluated and ranked alongside the swept ones; nothing else is refined.

**Microreact visuals** are generated for every fit the evaluator rates **Strong or Moderate** (deliberately _not_ `decision == ACCEPT`, which is Strong-only — the point is to compare all credible models by eye). Each is triggered **per fit, as soon as that fit is evaluated** (keyed join on species + model), not from the ranking: `ch_ranked` uses `groupTuple()` without a size (failed fits are dropped, so the count is unknown) and therefore waits for every species' evaluations, which on the 2026-10 Codon run held back all Microreact behind the slowest species. Only `accepted_models` / `best_model` (MultiQC) still come from the ranking. `POPPUNK_VISUALISE` wraps `poppunk_visualise --microreact`, on the same poppunk container (`rapidnj` and `mandrake` are bundled in the bioconda package).

Two traps in that module: PopPUNK resolves a model as `<model-dir>/<basename(model-dir)>_fit.pkl`, so **the fit and db directories must never be renamed via `stageAs`**; and it writes its output into a `<prefix>/` directory whose files are all prefixed with that basename, which the module then lifts into the task dir so everything publishes flat into one per-species `poppunk/microreact/` folder with the model name carried in each filename.

**One module, aliased.** `POPPUNK_REFINE_DBSCAN` is `POPPUNK_FITMODEL` aliased, and `POPPUNK_EVALUATE_DBSCAN_REFINE` is `POPPUNK_EVALUATE` aliased. The model spec is a _per-task value input_ (`fit_args`, e.g. `"bgmm --K 4"`, `"refine"`), **not** `ext.args` — that's what lets one process cover the whole sweep. `meta.model` labels each fit and drives both the join-back and the publish path. Follow this pattern rather than adding a new process.

### The evaluator is the pipeline's judgement

`bin/evaluate_poppunk_fastani.py` (wrapped by `modules/local/poppunk/evaluate`) is where all scientific decisions live. It emits `tool_metrics.tsv`, `cluster_metrics.tsv`, `genome_metrics.tsv`. Two fields drive control flow:

- **`tool_status`** — audit label: `Strong` / `Moderate` / `Mixed` / `Weak`, from HQ and total structure scores plus singleton / tiny-cluster / negative-silhouette rates.
- **`decision`** — binary: `ACCEPT` if `tool_status` ∈ `--accept-status` (default `Strong`), else `TRY_NEXT_MODEL`. `ch_accepted` reports every ACCEPT row; `ch_best` is the single top-ranked fit and is the fallback when nothing accepts.

`tool_metrics.tsv` column order is `model`, `poppunk_network_score`, `tool_status`, `n_clusters`, `n_nonsingleton_clusters`, `largest_cluster_fraction`, the scores, `decision`, `eval_summary`. **Never feed `poppunk_network_score` into `tool_status`.** It is `transitivity × (1 − density)`, so it peaks for over-fragmented fits (hundreds of singletons plus tiny near-identical cliques) that the evaluator must reject. A `Weak` fit that fails only on the singleton or tiny-cluster rules while its HQ structure score is ≥ 0.80 gets an `over-fragmentation:` prefix in `eval_summary`. That is a fit-level statement: "no structure" is a species-level conclusion (every fit fails). Don't reword it as strain-level; strains aren't defined here. The curation logic is documented in `docs/curation.md`. `poppunk_network_score` is PopPUNK's plain network `Score` (not the betweenness variants), parsed from the `<prefix>_fit.log` that `poppunk/fitmodel` saves into every fit dir. `<species>/poppunk/<species>_model_report.tsv` (built in `workflows/subspeciesprofiler.nf`) is the concatenation of every PopPUNK fit's `tool_metrics.tsv` **minus `decision`**; SynTracker rows have different columns and go to their own `<species>/syntracker/<species>_syntracker_model_report.tsv`, built the same way. `decision` stays in `tool_metrics.tsv` because control flow reads it.

Scores are **HQ-centric**: `tool_structure_score_HQ` is `NaN` when no HQ genome sits in a scorable cluster, which is forced to `Weak` so it ranks last. The evaluator's pytest cases in `tests/bin/` pin this behaviour (FastANI percent scale, within-species ANI gate, singleton exclusion, degenerate inputs, and the APSS mode below) — run them after touching the evaluator. ANI mode must stay byte-for-byte what PopPUNK fits have always got.

**APSS mode** (`--distance apss --apss <all-regions table> [--apss-noise <noise.tsv>]`, set automatically by `poppunk/evaluate` when its 5th input `apss` is non-empty; PopPUNK callers pass `[]`) scores SynTracker clusterings with the same code: the gap floor becomes 2 × the APSS standard error (instead of ANI's 0.002), cohesion has weight 0 (`APSS_WEIGHTS`; its ANI anchors are meaningless and the APSS 5th percentile does not discriminate), columns are suffixed `_APSS`, Cliff's delta is reported per cluster, `poppunk_network_score` is dropped, and the same partition is scored in ANI space → `ani_separated_fraction` (HQ share in clusters whose ANI Cliff's delta, within-cluster vs nearest-cluster, is ≥ 0.474; not the tail gap, which a few distant genomes inside an under-split cluster pull below the floor — that gave a false `synteny_only` on B. longum) and `signal` (`synteny_and_ani` / `synteny_only` / `none`). Because Leiden builds the clusters on the APSS graph, a clean APSS score alone is not evidence of subspecies; `signal` is.

### Failure policy

Fit and evaluate processes use `errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'ignore' }` in `conf/modules.config`. A degenerate grid point (bad `--K`/`--D`) fails _deterministically_, so it's ignored — the fit simply drops out of the ranking instead of killing the species. Only transient/resource failures retry (escalating via `task.attempt` in `conf/base.config`). Preserve this distinction.

`SYNTRACKER_RUN` is the exception, because every chunk job is required: resource/kill exits retry twice (memory and time escalate with `task.attempt`), any other exit retries once, then `ignore`. The errorStrategy checks `task.attempt` itself, because returning `retry` past `maxRetries` stops the whole pipeline. A species with any job still missing gets no SynTracker rows and a `log.warn` naming the failed chunks, never a partial merge: missing pairs would make genomes look unconnected (false singletons). `-resume` re-runs only the failed jobs. PopPUNK results are kept either way.

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

`subworkflows/local/syntracker_methods/main.nf`, an independent non-PopPUNK signal, skipped with `--skip_syntracker`. SynTracker separates genomes by gene **arrangement**, so it can split genomes of near-identical ANI; nothing collapses or labels genomes by ANI before SynTracker sees them.

1. `SEQKIT_STATS` (nf-core, `--all`) gives the HQ genomes' N50; `SYNTRACKER_SELECT` (`bin/syntracker_select.py`) picks the **reference** = the most central HQ genome (mean FastANI to the other HQ genomes, from `POPPUNK_METHODS.out.ani`) with N50 ≥ `--syntracker_ref_min_n50`, and the **targets** = every HQ genome, reference first then by N50. `--syntracker_max_targets` defaults to 0 = no cap; a positive value keeps the highest N50 (reference always kept). Waiting for FastANI costs ≤ 1 h against a SynTracker run of ≥ 16 h.
2. `SYNTRACKER_RUN` (`modules/local/syntracker/run`, container-only `quay.io/microbiome-informatics/syntracker:1.4.0_patch1`, **linux/amd64 only**). Always `-mode new`; the task fails if the all-regions APSS table has no pairs (SynTracker reports R failures only in its log). **Chunked:** with more than 2 × `--syntracker_chunk_size` (75) targets, they are dealt round-robin (in the select order) into k = ⌈n/75⌉ chunks and one `SYNTRACKER_RUN` runs per chunk pair Cᵢ ∪ Cⱼ (k(k−1)/2 jobs, ≤ 150 targets each, same reference); otherwise one job (`meta.chunk = 'all'`). `meta.chunk` (`c<i>_c<j>`) sets the publish folder `apss/<chunk>/`; `meta.n_chunk_jobs` feeds a `groupKey` so a species is released as soon as its own jobs finish. A job that still fails after its retries leaves its species incomplete (see "Failure policy"). The subworkflow joins the expected jobs to the finished ones, warns naming the missing chunks, and gives that species no SynTracker rows rather than a partial merge. Pairs within a chunk are computed in k − 1 jobs; compute ≈ 2(k−1)/k × one big run. B. fragilis (435 HQ): 6 chunks, 15 jobs.
3. `SYNTRACKER_FILTER` (`bin/syntracker_filter_multihit.py`), one task per job, applies SynTracker's **one-hit-per-genome rule** (see the bug below): it drops every per-region row in which either genome had > 1 BLAST hit for that region (from the job's `<ref>/blast_output/*.tab` + `combined_targets/sample_names_dictionary.tab`, emitted by `SYNTRACKER_RUN` as `blast`) and recomputes the all-regions APSS as the `math.fsum` mean of the kept per-region scores (order-independent, written `%.15g` like R; 95% of untouched pairs match SynTracker's text exactly, the rest within 1e-15). Report: `apss/<chunk>/*_syntracker_multihit.tsv`. Failure other than resource exits **terminates** the run.
4. `SYNTRACKER_MERGE` (`bin/syntracker_merge_apss.py`) concatenates the jobs' **filtered all-regions** tables (always, even for one job): one row per unordered pair, sorted, SynTracker's text values kept; it **fails** if a pair from two jobs differs in APSS or `Compared_regions`, or if the tables mix references, and that failure **stops the run on purpose** (a misbehaviour must not be skipped). `*_syntracker_merge.tsv` counts rows, duplicates and target pairs absent (SynTracker drops pairs sharing no region). The 40-region table for the noise estimate comes from the first job in chunk order only.
5. `SYNTRACKER_CLUSTERS` (`bin/syntracker_apss_clusters.py`) on the **all-regions** table: pairs with `Compared_regions < --syntracker_min_regions` dropped; targets whose median regions per pair is below `--syntracker_min_genome_coverage` × the species median excluded (`*_genome_coverage.tsv`); Leiden (modularity, weights = raw APSS, seeded) on the graph pruned at each floor, at each `--syntracker_resolutions`. Floors: `--syntracker_min_apss auto` (default) = `--syntracker_apss_steps` (8) values **evenly spaced** between the `--syntracker_apss_range` (1st–99th) percentiles of the species' cleaned APSS (`*_floors.tsv`); or a fixed list. Not quantiles: the separations sought are gaps, which hold almost no pairs (B. longum's 0.70–0.78 subspecies gap holds 0.2%), and species differ in where APSS sits (B. fragilis median 0.918, above any fixed grid suited to B. longum). It also writes `*_noise.tsv`: the per-region SD estimated from the 40-region table vs all-regions (var(APSS_40 − APSS_all) = SD²(1/40 − 1/R)) and the APSS SE = SD/√(median regions).
6. `SYNTRACKER_EVALUATE` = `POPPUNK_EVALUATE` aliased, given `[all-regions table, noise.tsv]` → APSS mode (above); rows `syntracker_leiden_apss<t>_r<res>` go to the SynTracker model report.

Hard-won SynTracker facts (do not relearn them):

- **The R stack must stay pinned** (R 4.0.5, DECIPHER 2.18.1, RSQLite). An unpinned current DECIPHER fails `Seqs2DB` on every region (`N function calls resulted in an error`) and SynTracker then crashes in `left_join()`. The image is built from the EBI containers repo, `syntracker/1.4.0_patch1`.
- **What APSS is.** Each region's score is 1 + log10(overlap/length / blocks), so it is set by the number of synteny blocks (1 ≈ 0.99, 2 ≈ 0.69, 3 ≈ 0.51); a pair's APSS ≈ 1 − its share of 5-kb windows broken into > 1 block (ρ = −0.98 on E. lenta). Near-identical genomes still break ~3% of windows. All-regions APSS is exactly the mean of the per-region scores.
- **SynTracker bug: its one-hit rule is not applied at extraction (found 2026-10-19).** `blast_per_region_process` (`/opt/syntracker/syntracker_first_stage/blast.py`) marks a genome with > 1 BLAST hit in a region `valid = False` and analyses the region only if ≥ 2 genomes are valid, but its extraction loop then extracts **every** genome, valid or not (an invalid one via its first hit). So whether a multi-copy genome's pairs get that region depends on whether another single-hit genome is in the run. Found on the Codon run (`spp_dev/syntracker/spp_final_results`): P. vulgatus `GCA_944320055 / GCF_021891395` had 620 regions in 7 chunk jobs and 619 in `c4_c8`, `c5_c8`, `c7_c8` (ΔAPSS 0.00023); the one region, `DBGDQY010000030.1_120000_121000`, holds two hits of GCA_944320055 and was analysed only where the reference (a single-hit target, chunk 1) or another single-hit genome was present (work dirs `82151b…` = c1_c8, `b42fa5…` = c4_c8, copied into that results folder). It was the only disagreement in ~278k duplicate pair computations over 4 chunked species. It affects unchunked runs too. `SYNTRACKER_FILTER` fixes it after the run (0.04–0.05% of per-region rows removed, ≤ 0.002 APSS change; both jobs then give 619 / 0.927080 and merge strictly). Patching the image instead would re-run every job. Worth reporting upstream: `if not hits_by_sample_dict[sample]["valid"]: continue` in the extraction loop.
- **Chunking is valid once the one-hit rule is applied (verified 2026-10-16, refined 2026-10-19).** A pair's per-region scores otherwise depend only on the two genomes and the reference (one DECIPHER `FindSynteny` per region, compared pair by pair). Check: 20 E. lenta genomes run alone (`spp_dev/codon_results/syntracker_results_sub`) vs the 199-genome run (`spp_dev/syntracker/e_lenta/syntracker_results`), same reference: 190/190 pairs identical in APSS and `Compared_regions`, all 54,622 per-region rows identical; merging the two tables reproduces the 199-run table exactly. The 40-region subsampled table agreed for only 19/190 pairs (max Δ 0.125), so **never concatenate the subsampled `avg_synteny_scores_<n>_regions.csv` tables** (SynTracker's `set.seed(1)` + `sample_n` over the sorted table makes a pair's draw depend on every pair before it).
- **Use all-regions, not the subsampled tables.** Subsampling adds noise of SD ≈ 0.226/√n (same in E. lenta and B. longum); at n = 40 that equalled the real spread between E. lenta pairs (r = 0.71 with all-regions). SynTracker fixes the subsampling seed (1), so the tables are reproducible, but a pair's draw depends on every pair sorted before it.
- **Missing regions are not random.** A target keeps fewer regions when it is fragmented (E. lenta: ρ = 0.87 with N50; the most fragmented fifth also get APSS +0.012 above their ANI trend) or far from the reference (B. longum minority subspecies: ~0.66-0.71 of the species median, hence the 0.4 coverage default — 0.5 started removing them).
- **First pipeline run (2026-10-10, `spp_dev/syntracker/spp_results_synt`), re-clustered offline with `auto` floors.** B. longum: 13 Strong / `synteny_and_ani`, a 3-group plateau at 0.710–0.807 identical to dbscan on the targets (ARI 1.0). E. lenta: 24/24 Weak / `none` (PopPUNK: all Weak). B. fragilis: clonal groups at 0.945–0.962 (Mixed, ~40 groups, 10–13% singletons), else Weak — matching PopPUNK's Moderate fits (157 clusters, largest 7%): strain-level clonal complexes, no subspecies. Two-thirds of B. fragilis's singletons come from the 150-target cap (21/150 targets lose all clonal partners to it); that is why the cap is now off and large species run chunked. A weight-shuffle permutation null is **not** a usable test: E. lenta's ANI-correlated gradient beats it by many SDs.
- **No Microreact for SynTracker, on purpose.** `poppunk_visualise --external-clustering --include-files` could draw SynTracker clusters, but only on PopPUNK's k-mer core-distance NJ tree; colouring synteny clusters on a non-synteny tree was judged misleading. Do not add it back without a synteny-based tree.
- **dRep at 95% is not a species check here.** On E. lenta it split 10 fragmented genomes (N50 ~7-11 kb) off at 96.8-99.3% ANI because of its 0.6 aligned-fraction threshold; that is why dRep was dropped for FastANI centrality.
- **Never use SynTracker's `-mode continue` in the pipeline.** It loads every finished region into the R parent before forking workers, so continue runs OOM where fresh runs don't; it also renames a region `_done` _before_ saving its result, so a kill can silently lose regions.
- Sample names: SynTracker names a target by its file basename minus extension; the module stages every target as `<sample>.fasta`, so names reconcile with `normalise_genome_id`.
- `normalise_genome_id` (and `read_fastani`) are imported by `bin/syntracker_*.py` from `bin/evaluate_poppunk_fastani.py` (still the single place). The select/filter/merge/clusters env pins **pandas < 3**: the evaluator breaks under pandas 3's copy-on-write (`underlying array is read-only` in `test_check_species_ani_all_close_passes`). pandas' Spearman needs scipy, which is not in any image.
- `procps-ng` is linux-only on conda-forge, so the select/filter/merge/clusters `environment.yml` cannot be solved on macOS; tests there need a local env with pandas and python-igraph. `matplotlib-base` is in that image but unused.

Follow-ups (documented, not implemented):

- Edge weights: after pruning, every edge weight lies in `[min_apss, 1]`, so modularity is driven mostly by which edges exist. Rescaling to `(APSS − min_apss)/(1 − min_apss)` restores the dynamic range; planned as a flag, raw APSS stays the default.

## Planned work (do NOT "clean these up")

- **Finish the Codon run** (`spp_dev/syntracker/spp_final_results`) with `-resume`: every `SYNTRACKER_RUN` should be cached (output declarations are not in the task hash; checked locally), only FILTER, MERGE and downstream re-run, and the P. vulgatus merge must pass strictly. Then check whether B. fragilis's clonal-level fits move from Mixed to Moderate, like PopPUNK's.

- **`checkm2/predict`** is installed in `modules.json` but unwired _on purpose_. It will be wired behind a flag that generates completeness/contamination when the samplesheet's `qc_csv` is absent. Separate task, not yet started.

## Known loose ends

- Test data is a local fixture (`tests/fixtures/`, `conf/test.config`). Now that this isn't an nf-core pipeline it does not need to move to nf-core/test-datasets, but `tests/nextflow.config` still points `pipelines_testdata_base_path` at a `subspeciesprofiler` branch the fixtures don't use.
- `ro-crate-metadata.json` is generated — `nf-core pipelines lint` refreshes its README snapshot automatically. Never hand-edit it.
- The nf-core logo images and the template `TODO nf-core:` markers were removed on purpose (2026-10-16); the completion email has no inline logo. Don't restore them. `conf/test_full.config` still holds the template's placeholder input (viralrecon) and `genome`, not real full-size test data.
- `bin/spp_species_eligibility.py` (GTDB/RefSeq-wide species screening) is not called by any module; it's an upstream/offline companion to `spp_eligibility_from_qc.py`.
