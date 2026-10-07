//
// SYNTRACKER_METHODS: synteny-based clustering for one eligible species, independent of PopPUNK.
//
// 1. dRep dereplicates the species' HQ genomes at a strict secondary ANI (ext.args, -sa 0.99), so
//    near-identical genomes collapse to one representative.
// 2. The representatives are SynTracker's targets, capped at --syntracker_max_targets (largest
//    dRep clusters first): SynTracker's cost grows with the square of the target count. The
//    winner of the largest dRep cluster is the single reference genome.
// 3. SynTracker scores every target pair (APSS); SYNTRACKER_CLUSTERS turns the APSS into one
//    Taxon,Cluster table per threshold and propagates each representative's cluster to its
//    dRep cluster members.
// 4. Every table is scored against FastANI by the same evaluator as the PopPUNK fits, so the
//    SynTracker rows land in the species' model report (model `syntracker_<avg|cc>_apss<t>`).
//

include { DREP_DEREPLICATE                        } from '../../../modules/nf-core/drep/dereplicate/main'
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
    ch_hq = ch_input.map { meta, genomes, rfile, labels ->
        def label_of = labels.splitCsv(header: true).collectEntries { row -> [ row.genome, row.label ] }
        def hq_files = rfile.splitCsv(sep: '\t')
            .findAll { row -> label_of[row[0]] == 'HQ' }
            .collect { row -> file(row[1]).name } as Set
        [ meta, [ genomes ].flatten().findAll { genome -> genome.name in hq_files } ]
    }
    DREP_DEREPLICATE( ch_hq, [ [], [] ] )

    // Reference and targets from dRep's tables: Wdb = winner per secondary cluster,
    // Cdb = cluster membership (also passed on for propagation).
    def max_targets = params.syntracker_max_targets as int
    ch_selected = DREP_DEREPLICATE.out.summary_tables
        .join( DREP_DEREPLICATE.out.fastas )
        .map { meta, tables, fastas ->
            def table    = { name -> [ tables ].flatten().find { it.name == name } }
            def rows     = { name -> table(name)?.size() ? table(name).splitCsv(header: true) : [] }
            def fasta_of = [ fastas ].flatten().collectEntries { [ it.name, it ] }
            def size     = rows('Cdb.csv').countBy { it.secondary_cluster }
            def winners  = rows('Wdb.csv').sort { a, b ->
                (size[b.cluster] <=> size[a.cluster]) ?: ((b.score as double) <=> (a.score as double)) ?: (a.genome <=> b.genome)
            }*.genome
            if ( !winners ) {
                // dRep tables without rows (e.g. the module's -stub): take its FASTAs in name order.
                winners = fasta_of.keySet().sort()
            }
            def kept = winners.take(max_targets)
            if ( winners.size() > max_targets ) {
                def wdb_cluster = rows('Wdb.csv').collectEntries { [ it.genome, it.cluster ] }
                def dropped = winners.drop(max_targets).sum { size[wdb_cluster[it]] ?: 1 }
                log.warn("Species '${meta.id}': ${winners.size()} dRep representatives > --syntracker_max_targets (${max_targets}); " +
                    "SynTracker runs on the ${max_targets} largest clusters, leaving ${dropped} genome(s) unclustered.")
            }
            def cdb = size ? table('Cdb.csv') : []  // only propagate from a real cluster table
            [ meta, fasta_of[kept[0]], kept.collect { fasta_of[it] }, cdb ]
        }
        .filter { meta, reference, targets, cdb ->
            def ok = targets.size() >= 2
            if ( !ok ) {
                log.warn("Species '${meta.id}': fewer than 2 dRep representatives; skipping SynTracker.")
            }
            ok
        }

    SYNTRACKER_RUN( ch_selected.map { meta, reference, targets, cdb -> [ meta, reference, targets ] } )

    // The APSS table to cluster: all regions (deterministic) or one subsampled level.
    def apss_name = "avg_synteny_scores_${params.syntracker_regions}_regions.csv"
    ch_clusters_in = SYNTRACKER_RUN.out.apss
        .map { meta, files -> [ meta.id, meta, [ files ].flatten().find { it.name == apss_name } ] }
        .filter { id, meta, apss ->
            if ( apss == null ) {
                log.warn("Species '${meta.id}': SynTracker produced no ${apss_name}; skipping its clustering.")
            }
            apss != null
        }
        .join( ch_input.map { meta, genomes, rfile, labels -> [ meta.id, labels ] } )
        .join( ch_selected.map { meta, reference, targets, cdb -> [ meta.id, cdb ] } )
        .map { id, meta, apss, labels, cdb -> [ meta, apss, labels, cdb ] }
    SYNTRACKER_CLUSTERS( ch_clusters_in )

    // One evaluation per threshold: meta.model is the table name minus '<id>_' and '_clusters.csv'.
    ch_eval_in = SYNTRACKER_CLUSTERS.out.clusters
        .flatMap { meta, tables ->
            [ tables ].flatten().collect { table ->
                def model = table.name - "${meta.id}_" - '_clusters.csv'
                [ meta.id, meta + [ model: model ], table ]
            }
        }
        .combine( ch_ani.map { meta, ani -> [ meta.id, ani ] }, by: 0 )
        .combine( ch_input.map { meta, genomes, rfile, labels -> [ meta.id, labels ] }, by: 0 )
        .map { id, meta, table, ani, labels -> [ meta, table, ani, labels ] }
    SYNTRACKER_EVALUATE( ch_eval_in )

    emit:
    tool_metrics = SYNTRACKER_EVALUATE.out.tool_metrics                                          // channel: [ val(meta), path(tool_metrics.tsv) ]
    apss         = SYNTRACKER_RUN.out.apss                                                       // channel: [ val(meta), [ avg_synteny_scores_*.csv ] ]
    clusters     = SYNTRACKER_CLUSTERS.out.clusters                                              // channel: [ val(meta), [ *_clusters.csv ] ]
    reference    = ch_selected.map { meta, reference, targets, cdb -> [ meta, reference ] }      // channel: [ val(meta), path(reference fasta) ]
}
