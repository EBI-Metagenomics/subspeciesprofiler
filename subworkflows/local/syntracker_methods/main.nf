//
// SYNTRACKER_METHODS: synteny-based clustering for one eligible species, independent of PopPUNK.
//
// 1. SEQKIT_STATS gives each HQ genome's N50; SYNTRACKER_SELECT picks the reference (the most
//    central HQ genome by mean FastANI, among those with N50 >= --syntracker_ref_min_n50) and the
//    targets (every HQ genome, the highest-N50 ones first when capped at --syntracker_max_targets:
//    SynTracker's cost grows with the square of the target count).
// 2. SynTracker scores every target pair (APSS) over the reference's regions.
// 3. SYNTRACKER_CLUSTERS clusters the all-regions APSS: pairs on fewer than
//    --syntracker_min_regions regions are dropped, low-coverage targets excluded, and Leiden runs
//    on the APSS graph pruned at each --syntracker_min_apss, at each --syntracker_resolutions.
// 4. Every table is scored in APSS space by the evaluator (gap margin = 2x the APSS standard error,
//    estimated from the 40-region table), together with its ANI concordance; the rows go to the
//    species' SynTracker model report (model `syntracker_leiden_apss<t>_r<res>`).
//

include { SEQKIT_STATS                            } from '../../../modules/nf-core/seqkit/stats/main'
include { SYNTRACKER_SELECT                       } from '../../../modules/local/syntracker/select/main'
include { SYNTRACKER_RUN                          } from '../../../modules/local/syntracker/run/main'
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

    SYNTRACKER_RUN( ch_selected )

    // Cluster the all-regions table; the 40-region table only feeds the noise estimate.
    ch_apss = SYNTRACKER_RUN.out.apss
        .map { meta, files ->
            def all = [ files ].flatten()
            [ meta.id, meta,
              all.find { it.name == 'avg_synteny_scores_all_regions.csv' },
              all.find { it.name == 'avg_synteny_scores_40_regions.csv' } ?: [] ]
        }
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
    apss         = SYNTRACKER_RUN.out.apss                                       // channel: [ val(meta), [ avg_synteny_scores_*.csv ] ]
    clusters     = SYNTRACKER_CLUSTERS.out.clusters                              // channel: [ val(meta), [ *_clusters.csv ] ]
    coverage     = SYNTRACKER_CLUSTERS.out.coverage                              // channel: [ val(meta), path(*_genome_coverage.tsv) ]
    noise        = SYNTRACKER_CLUSTERS.out.noise                                 // channel: [ val(meta), path(*_noise.tsv) ]
    selection    = SYNTRACKER_SELECT.out.selection                               // channel: [ val(meta), path(*_selection.tsv) ]
    reference    = ch_selected.map { meta, reference, targets -> [ meta, reference ] } // channel: [ val(meta), path(reference fasta) ]
}
