# ebi-metagenomics/subspeciesprofiler: Curating the model report

## Introduction

The pipeline fits many PopPUNK models per species and scores every fit against all-vs-all FastANI. The result is `<species>/poppunk/<species>_model_report.tsv`, one row per fit. This page explains how to read that table and pick a model (or conclude there is none), using worked examples from a real run.

No single column decides. Read each row as a **pattern** across several columns.

## Key columns

| Column                                                     | What it tells you                                                                                                                                                |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `poppunk_network_score`                                    | PopPUNK's own network score, reported for context only. It is **not** a verdict (see below).                                                                     |
| `tool_status`                                              | The evaluator's verdict: `Strong` / `Moderate` / `Mixed` / `Weak`. It asks whether the clusters are subspecies-scale groups separated by ANI.                    |
| `n_clusters`, `n_nonsingleton_clusters`                    | How many clusters the fit produced, and how many have more than one genome. Hundreds of clusters of which few have more than one genome means fragmentation.     |
| `largest_cluster_fraction`                                 | Share of genomes in the biggest cluster. Near 1 means one dominant group; near 0 means no large groups.                                                          |
| `tool_structure_score_*`, `*_silhouette_*`                 | Cluster quality (cohesion and ANI separation). **Computed over non-singleton genomes only**: they describe the clustered genomes, however few of them there are. |
| `singleton_rate_HQ`, `tiny_cluster_rate_HQ`                | Fraction of HQ genomes that are singletons, or in clusters of fewer than 5. **Always read these together with the scores above.**                                |
| `defective_HQ_fraction`, `negative_silhouette_HQ_fraction` | Fraction of HQ genomes in clusters that overlap on ANI with their nearest neighbour, or that sit closer to another cluster than to their own.                    |
| `eval_summary`                                             | Plain-text reason for the verdict, naming the thresholds that failed.                                                                                            |

## Why the network score is not a verdict

PopPUNK's network score is `transitivity × (1 − density)`:

- **Transitivity** is 1 when every connected component is a clique (every member linked to every other).
- **Density** is close to 0 when most genomes are unlinked.

A very fine clustering therefore scores close to 1: hundreds of singletons plus a few pairs and triples of near-identical genomes. A large, genuine subspecies cluster is never a perfect clique, because its members are not all within the threshold of each other. That lowers its transitivity, and so its score.

The highest network score thus often marks the **most fragmented** fit, not the best one. For the same reason, PopPUNK's `refine` step (which maximises this score) can gain score by splitting a few genomes off an otherwise clean solution. Use the network score only as a tie-break between fits that are equally good on everything else.

## Patterns

| Pattern                            | Signature                                                                                                                                                                                                             | Meaning                                                                 | Action                                                                                      |
| ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| **Over-fragmentation** (one fit)   | `Weak`, `eval_summary` starts with `over-fragmentation`; high network score; `singleton_rate_HQ` > 0.3; `tiny_cluster_rate_HQ` high; `largest_cluster_fraction` near 0; structure score and silhouette look excellent | The fit groups only near-identical genomes and leaves the rest isolated | Reject this fit. It says nothing about whether the species has structure.                   |
| **Over-split continuum** (one fit) | `Weak`; high `defective_HQ_fraction` and `negative_silhouette_HQ_fraction`                                                                                                                                            | The cutoff slices through a group that is continuous in ANI             | Reject this fit                                                                             |
| **Subspecies structure**           | `Strong` (or `Moderate`); singleton and tiny-cluster rates near 0; median silhouette around 0.4–0.6; the network score can be low                                                                                     | Groups separated by a clear ANI gap                                     | Candidate: compare the candidates (below) and inspect them in Microreact                    |
| **No structure** (the species)     | **Every** fit is `Weak`: over-fragmented at fine cutoffs, while at coarser cutoffs the large clusters are Defective and overlap on ANI                                                                                | The genomes form an ANI continuum; no cutoff yields separated groups    | Report "no subspecies structure" for this species at this sampling. This is a valid result. |

**Over-fragmentation is not the same as "no structure".** Over-fragmentation is a verdict on **one fit**. "No structure" is a conclusion about the **species**, and is warranted only when **no** fit in the sweep produces separated groups.

## Choosing among several good fits

When several fits are `Strong` or `Moderate`:

1. Prefer the highest `accepted_HQ_fraction` with the lowest singleton and tiny-cluster rates.
2. Prefer **stable** solutions: many grid points, or several model families, that give the same clustering (identical metrics across rows) are strong evidence the structure is real.
3. Only then use a higher `poppunk_network_score` as a tie-break.
4. Open the candidates' `<species>_<model>.microreact` files and check that the groups make sense (tree, metadata, known subspecies).

The pipeline's own `best_model` follows the same order: `tool_status` first, then the HQ structure score.

## SynTracker rows

Rows named `syntracker_leiden_n<N>_apss<t>_r<res>` come from the SynTracker branch: HQ genomes clustered by **synteny** (average pairwise synteny score, APSS, at subsampling depth `N`) into Leiden communities of the APSS graph pruned at `t`, with resolution `res`, then scored by the same FastANI evaluator. Read them like any other row, with three differences:

- **`poppunk_network_score` is empty.** There is no PopPUNK network behind them.
- **APSS is not on the ANI scale.** Pairs within a species typically sit around APSS 0.80-0.90; on _B. longum_, subspecies separated at **0.72-0.80** (with the earlier average-linkage method). A high floor leaves targets without edges, which become singletons, so the same over-fragmentation flag applies; higher resolutions split communities further.
- **Coverage can be partial.** Only HQ genomes are targets, capped at `--syntracker_max_targets` by N50; MQ genomes and genomes dropped by the cap are absent from the tables and not scored. The run log warns when the cap applies. The per-cluster QC table flags clusters of two genomes or fewer as `low_confidence`.

SynTracker is an **independent** signal: it does not use PopPUNK's k-mer distances. When a SynTracker plateau and a PopPUNK solution agree on the same groups, that is strong evidence the structure is real.

## Worked examples

These come from a test run on _Bifidobacterium longum_ (568 genomes after QC) and _Eggerthella lenta_ (199 genomes).

### _B. longum_: structure found

| Fit                             | Network   | Status | Clusters (non-singleton) | Largest | Reading                                                                                                                                                                               |
| ------------------------------- | --------- | ------ | ------------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `threshold_q0.001`              | 0.969     | Weak   | 473 (57)                 | 0.01    | **Over-fragmentation.** 71% of HQ genomes are singletons; the rest sit in pairs and triples at 99.98%+ ANI. The structure score (0.93) and silhouette (0.98) describe only those few. |
| `threshold_q0.005`              | 0.053     | Weak   | 154 (22)                 | 0.65    | **Over-split continuum.** 65% of HQ genomes are in Defective clusters, and 83% have a negative silhouette.                                                                            |
| `dbscan_*` (all 12 grid points) | 0.11      | Strong | 4 (3)                    | 0.94    | **Subspecies structure.** Clusters of 535, 18 and 14 genomes, 98.2–99.4% ANI within and ~96.6% between (a 1.4–2.2% gap with no overlap). Identical across the whole grid.             |
| `threshold_q0.1`–`q0.2`         | 0.34–0.41 | Strong | —                        | —       | Converge on the dbscan solution.                                                                                                                                                      |
| `refine_from_dbscan_*`          | 0.40      | Strong | —                        | —       | Same solution, but splits off a few genomes to raise the network score (HQ accepted 0.966 vs 1.0).                                                                                    |

SynTracker on the same species, with the **earlier method** (dRep 99% representatives, average linkage or connected components on all-regions APSS; 588 genomes, scored against the same FastANI). The current Leiden method has not been run on this species yet:

| Fit                               | Status   | Clusters | Reading                                                                                                                    |
| --------------------------------- | -------- | -------- | -------------------------------------------------------------------------------------------------------------------------- |
| `syntracker_avg_apss0.72`–`0.80`  | Strong   | 4–9      | **The same three groups as dbscan** (544, 22 and 21 genomes; agreement with dbscan ARI 0.98), stable over five thresholds. |
| `syntracker_avg_apss0.82`, `0.85` | Weak     | 23–144   | **Over-split continuum**: the largest group starts to break up.                                                            |
| `syntracker_avg_apss0.90`         | Weak     | 426      | **Over-fragmentation.**                                                                                                    |
| `syntracker_cc_apss0.74`–`0.80`   | Moderate | 3        | Connected components merge two of the groups through a few bridging pairs; only 0.82–0.85 recovers all three.              |

**Conclusion:** choose the dbscan solution, confirmed independently by the SynTracker plateau. It has the cleanest metrics and is stable across the grid, even though its network score is among the lowest in the table.

### _E. lenta_: no subspecies structure

| Fit                                     | Network    | Status | Clusters (non-singleton) | Largest   | Reading                                                                                                                                                 |
| --------------------------------------- | ---------- | ------ | ------------------------ | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bgmm_K10`, `threshold_q0.002`          | 0.85–0.997 | Weak   | 159–176 (12–25)          | 0.04      | **Over-fragmentation**, as above.                                                                                                                       |
| `threshold_q0.02`, `dbscan_D3_mcp0.005` | 0.03–0.72  | Weak   | 83–89 (6–25)             | 0.18–0.40 | Coarser cuts, but the large clusters are **Defective**: intra-cluster p5 ANI ~98.6% against nearest-cluster p95 ANI ~98.9–99.1%, overlapping by 65–97%. |

**Conclusion:** no fit separates the genomes. Fine cutoffs over-fragment, and coarse cutoffs merge overlapping groups, because the species is an ANI continuum (~98.2–99.1%). Report _E. lenta_ as having no subspecies structure at this sampling.
