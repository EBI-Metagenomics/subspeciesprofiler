//
// SYNTRACKER_METHODS: synteny-based clustering for one eligible species, independent of PopPUNK.
//
// 1. SEQKIT_STATS gives each HQ genome's N50; SYNTRACKER_SELECT picks the reference (the most
//    central HQ genome by mean FastANI, among those with N50 >= --syntracker_ref_min_n50) and the
//    targets (every HQ genome, reference first and then by N50; --syntracker_max_targets > 0 caps
//    them, keeping the highest N50).
// 2. SynTracker scores every target pair (APSS) over the reference's regions. Its cost grows with
//    the square of the targets, so above 2 x --syntracker_chunk_size they are split round-robin
//    (in N50 order) into k chunks and SynTracker runs once per chunk pair Ci + Cj, k(k-1)/2 jobs
//    against the same reference. A pair's all-regions APSS depends only on the two genomes and
//    the reference (checked on E. lenta: 190/190 pairs identical run alone vs in the 199-genome
//    run), so SYNTRACKER_MERGE concatenates the all-regions tables, failing if a pair computed in
//    two jobs (every pair within a chunk) differs. The subsampled tables cannot be merged (a pair's
//    subsample depends on the other targets): the noise estimate uses one job's 40-region table.
// 3. SYNTRACKER_CLUSTERS clusters the all-regions APSS: pairs on fewer than
//    --syntracker_min_regions regions are dropped, low-coverage targets excluded, and Leiden runs
//    on the APSS graph pruned at each floor (`--syntracker_min_apss auto`: evenly spaced over the
//    species' own APSS range), at each --syntracker_resolutions.
// 4. Every table is scored in APSS space by the evaluator (gap margin = 2x the APSS standard error,
//    estimated from the 40-region table), together with its ANI concordance; the rows go to the
//    species' SynTracker model report (model `syntracker_leiden_apss<t>_r<res>`).
//

include { SEQKIT_STATS                            } from '../../../modules/nf-core/seqkit/stats/main'
include { SYNTRACKER_SELECT                       } from '../../../modules/local/syntracker/select/main'
include { SYNTRACKER_RUN                          } from '../../../modules/local/syntracker/run/main'
include { SYNTRACKER_MERGE                        } from '../../../modules/local/syntracker/merge/main'
include { SYNTRACKER_CLUSTERS                     } from '../../../modules/local/syntracker/clusters/main'
include { POPPUNK_EVALUATE as SYNTRACKER_EVALUATE } from '../../../modules/local/poppunk/evaluate/main'

workflow SYNTRACKER_METHODS {

    take:
    ch_input // channel: [ val(meta), [ genome files ], path(rfile), path(labels) ]  (as POPPUNK_METHODS)
    ch_ani   // channel: [ val(meta), path(ani) ]  all-vs-all FastANI, from POPPUNK_METHODS

    main:

    // HQ genomes only: the r-file maps sample -> ./genomes/<file>, the labels map sample -> HQ/MQ,
    // so no FASTA-suffix logic is needed here.
    ch_hq = ch_input
        .map { meta, genomes, rfile, labels ->
            def label_of = labels.splitCsv(header: true).collectEntries { row -> [ row.genome, row.label ] }
            def hq_files = rfile.splitCsv(sep: '\t')
                .findAll { row -> label_of[row[0]] == 'HQ' }
                .collect { row -> file(row[1]).name } as Set
            [ meta, [ genomes ].flatten().findAll { genome -> genome.name in hq_files } ]
        }
        .filter { meta, hq ->
            if ( hq.size() < 2 ) {
                log.warn("Species '${meta.id}': fewer than 2 HQ genomes; skipping SynTracker.")
            }
            hq.size() >= 2
        }

    SEQKIT_STATS( ch_hq )
    SYNTRACKER_SELECT(
        ch_ani
            .map { meta, ani -> [ meta.id, ani ] }
            .join( SEQKIT_STATS.out.stats.map { meta, stats -> [ meta.id, meta, stats ] } )
            .map { id, ani, meta, stats -> [ meta, ani, stats ] }
    )

    // Back from the selected file names to the staged genome files.
    ch_selected = ch_hq
        .map { meta, hq -> [ meta.id, meta, hq ] }
        .join( SYNTRACKER_SELECT.out.reference.map { meta, ref -> [ meta.id, ref.text.trim() ] } )
        .join( SYNTRACKER_SELECT.out.targets.map { meta, targets -> [ meta.id, targets.readLines()*.trim().findAll() ] } )
        .map { id, meta, hq, reference, targets ->
            def file_of = [ hq ].flatten().collectEntries { [ it.name, it ] }
            [ meta, file_of[reference], targets.collect { file_of[it] } ]
        }

    // One SynTracker job per chunk pair (k > 2), else a single job on every target (chunk 'all').
    // meta.chunk names the job (publish path); meta.n_chunk_jobs lets the gather below release a
    // species as soon as its own jobs are done.
    def chunk_size = params.syntracker_chunk_size as int
    ch_chunk_runs = ch_selected
        .flatMap { meta, reference, targets ->
            def k = ( targets.size() + chunk_size - 1 ).intdiv( chunk_size )
            if ( k <= 2 ) {
                return [ [ meta + [ chunk: 'all', n_chunk_jobs: 1 ], reference, targets ] ]
            }
            def chunks = ( 0..<k ).collect { i -> targets.withIndex().findAll { t, j -> j % k == i }.collect { t, j -> t } }
            def pairs  = ( 0..<k ).collectMany { i -> ( (i + 1)..<k ).collect { j -> [ i, j ] } }
            log.info("Species '${meta.id}': ${targets.size()} SynTracker targets in ${k} chunks, ${pairs.size()} jobs.")
            pairs.collect { i, j ->
                [ meta + [ chunk: "c${i + 1}_c${j + 1}", n_chunk_jobs: pairs.size() ], reference, chunks[i] + chunks[j] ]
            }
        }

    SYNTRACKER_RUN( ch_chunk_runs )

    // Gather each species' jobs (in chunk order, for -resume) and check that every job came back.
    // A missing job (failed after its retries and ignored, see conf/modules.config) would leave all
    // of its chunk pair's genome pairs out of the merge, and those genomes would then look
    // unconnected (false singletons). So an incomplete species gets no SynTracker rows, with a
    // warning naming the failed jobs; `-resume` re-runs only those. Its PopPUNK results are kept.
    ch_jobs_expected = ch_chunk_runs
        .map { meta, reference, targets -> [ groupKey( meta.id, meta.n_chunk_jobs ), meta.chunk ] }
        .groupTuple()
        .map { id, chunks -> [ id.toString(), chunks ] }
    ch_jobs_done = SYNTRACKER_RUN.out.apss
        .map { meta, files ->
            def all = [ files ].flatten()
            [ groupKey( meta.id, meta.n_chunk_jobs ), meta.chunk,
              all.find { it.name == 'avg_synteny_scores_all_regions.csv' },
              all.find { it.name == 'avg_synteny_scores_40_regions.csv' } ?: [] ]
        }
        .groupTuple( remainder: true )
        .map { id, chunks, all_regions, sub40 -> [ id.toString(), [ chunks, all_regions, sub40 ] ] }
    ch_runs = ch_jobs_expected
        .join( ch_jobs_done, remainder: true )
        .filter { id, expected, done ->
            def missing = expected - ( done ? done[0] : [] )
            if ( missing ) {
                log.warn("Species '${id}': ${missing.size()} of ${expected.size()} SynTracker jobs failed after retries " +
                    "(${missing.sort().join(', ')}); no SynTracker results for this species. " +
                    "Its PopPUNK results are kept. Fix the cause and re-run with -resume: only the failed jobs run again.")
            }
            !missing
        }
        .map { id, expected, done ->
            def ( chunks, all_regions, sub40 ) = done
            def order = chunks.withIndex().sort { it[0] }*.get(1)
            [ id, order.collect { all_regions[it] }, sub40[order[0]] ]
        }

    SYNTRACKER_MERGE(
        ch_selected
            .map { meta, reference, targets -> [ meta.id, meta ] }
            .join( ch_runs )
            .join( SYNTRACKER_SELECT.out.targets.map { meta, targets -> [ meta.id, targets ] } )
            .map { id, meta, all_regions, sub40, targets -> [ meta, all_regions, targets ] }
    )

    // Cluster the merged all-regions table; the 40-region table (of the first job only) feeds the
    // noise estimate.
    ch_apss = SYNTRACKER_MERGE.out.apss
        .map { meta, all_regions -> [ meta.id, meta, all_regions ] }
        .join( ch_runs.map { id, all_regions, sub40 -> [ id, sub40 ] } )
    SYNTRACKER_CLUSTERS(
        ch_apss
            .join( SYNTRACKER_SELECT.out.targets.map { meta, targets -> [ meta.id, targets ] } )
            .map { id, meta, all_regions, sub40, targets -> [ meta, all_regions, sub40, targets ] }
    )

    // One evaluation per table: meta.model is the table name minus '<id>_' and '_clusters.csv'.
    ch_eval_in = SYNTRACKER_CLUSTERS.out.clusters
        .flatMap { meta, tables ->
            [ tables ].flatten().collect { table ->
                def model = table.name - "${meta.id}_" - '_clusters.csv'
                [ meta.id, meta + [ model: model ], table ]
            }
        }
        .combine( ch_ani.map { meta, ani -> [ meta.id, ani ] }, by: 0 )
        .combine( ch_input.map { meta, genomes, rfile, labels -> [ meta.id, labels ] }, by: 0 )
        .combine( ch_apss.map { id, meta, all_regions, sub40 -> [ id, all_regions ] }, by: 0 )
        .combine( SYNTRACKER_CLUSTERS.out.noise.map { meta, noise -> [ meta.id, noise ] }, by: 0 )
        .map { id, meta, table, ani, labels, all_regions, noise -> [ meta, table, ani, labels, [ all_regions, noise ] ] }
    SYNTRACKER_EVALUATE( ch_eval_in )

    emit:
    tool_metrics = SYNTRACKER_EVALUATE.out.tool_metrics                          // channel: [ val(meta), path(tool_metrics.tsv) ]
    apss         = SYNTRACKER_MERGE.out.apss                                     // channel: [ val(meta), path(avg_synteny_scores_all_regions.csv) ] merged over the jobs
    merge_report = SYNTRACKER_MERGE.out.report                                   // channel: [ val(meta), path(*_syntracker_merge.tsv) ]
    runs         = SYNTRACKER_RUN.out.apss                                       // channel: [ val(meta + chunk), [ avg_synteny_scores_*.csv ] ] one per job
    clusters     = SYNTRACKER_CLUSTERS.out.clusters                              // channel: [ val(meta), [ *_clusters.csv ] ]
    coverage     = SYNTRACKER_CLUSTERS.out.coverage                              // channel: [ val(meta), path(*_genome_coverage.tsv) ]
    noise        = SYNTRACKER_CLUSTERS.out.noise                                 // channel: [ val(meta), path(*_noise.tsv) ]
    floors       = SYNTRACKER_CLUSTERS.out.floors                                // channel: [ val(meta), path(*_floors.tsv) ]
    selection    = SYNTRACKER_SELECT.out.selection                               // channel: [ val(meta), path(*_selection.tsv) ]
    reference    = ch_selected.map { meta, reference, targets -> [ meta, reference ] } // channel: [ val(meta), path(reference fasta) ]
}
