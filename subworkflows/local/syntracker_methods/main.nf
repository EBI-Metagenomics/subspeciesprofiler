//
// SYNTRACKER_METHODS: synteny-based clustering for one eligible species, independent of PopPUNK.
//
// 1. dRep (with the HQ genomes' completeness/contamination, -sa 0.95) picks the SynTracker
//    reference: normally one representative per species; if the species splits at 95% ANI, the
//    representative with the highest N50 (least fragmented) is used, with a warning.
// 2. Every HQ genome is a SynTracker target (no dereplication), capped at
//    --syntracker_max_targets by N50 (least fragmented first): SynTracker's cost grows with the
//    square of the target count, and contig breaks look like synteny breaks.
// 3. SynTracker scores every target pair (APSS). SYNTRACKER_CLUSTERS picks one subsampling depth
//    (--syntracker_regions, `auto` = the highest before retention falls off a cliff) and runs
//    Leiden on the APSS graph pruned at each --syntracker_min_apss, at each
//    --syntracker_resolutions: one Taxon,Cluster table per combination.
// 4. Every table is scored against FastANI by the same evaluator as the PopPUNK fits, so the
//    SynTracker rows land in the species' model report
//    (model `syntracker_leiden_n<n>_apss<t>_r<res>`).
//

include { DREP_DEREPLICATE                        } from '../../../modules/local/drep/dereplicate/main'
include { SYNTRACKER_RUN                          } from '../../../modules/local/syntracker/run/main'
include { SYNTRACKER_CLUSTERS                     } from '../../../modules/local/syntracker/clusters/main'
include { POPPUNK_EVALUATE as SYNTRACKER_EVALUATE } from '../../../modules/local/poppunk/evaluate/main'

workflow SYNTRACKER_METHODS {

    take:
    ch_input      // channel: [ val(meta), [ genome files ], path(rfile), path(labels) ]  (as POPPUNK_METHODS)
    ch_ani        // channel: [ val(meta), path(ani) ]  all-vs-all FastANI, from POPPUNK_METHODS
    ch_genomeinfo // channel: [ val(meta), path(genomeinfo) ]  dRep --genomeInfo for the HQ genomes, from SPECIESQC

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
    DREP_DEREPLICATE(
        ch_hq
            .map { meta, hq -> [ meta.id, meta, hq ] }
            .join( ch_genomeinfo.map { meta, info -> [ meta.id, info ] } )
            .map { id, meta, hq, info -> [ meta, hq, info ] }
    )

    // Reference = the dRep winner with the highest N50; targets = every HQ genome, least
    // fragmented first when capped. N50 comes from dRep's genomeInformation.csv
    // (genome,completeness,contamination,length,N50,centrality), keyed by file name.
    def max_targets = params.syntracker_max_targets as int
    ch_selected = DREP_DEREPLICATE.out.summary_tables
        .map { meta, tables -> [ meta.id, meta, tables ] }
        .join( ch_hq.map { meta, hq -> [ meta.id, hq ] } )
        .map { id, meta, tables, hq ->
            def table   = { name -> [ tables ].flatten().find { it.name == name } }
            def rows    = { name -> table.call(name)?.size() ? table.call(name).splitCsv(header: true) : [] }
            def n50     = rows.call('genomeInformation.csv').collectEntries { [ it.genome, (it.N50 ?: 0) as double ] }
            def score   = rows.call('Wdb.csv').collectEntries { [ it.genome, (it.score ?: 0) as double ] }
            def by_frag = { a, b -> (n50[b] ?: 0d) <=> (n50[a] ?: 0d) ?: (score[b] ?: 0d) <=> (score[a] ?: 0d) ?: a <=> b }
            def hq_of   = [ hq ].flatten().collectEntries { [ it.name, it ] }

            def winners = score.keySet().findAll { it in hq_of }.sort(false, by_frag)
            if ( !winners ) {
                // dRep tables without rows: fall back to the least fragmented HQ genome.
                winners = hq_of.keySet().sort(false, by_frag).take(1)
            }
            if ( winners.size() > 1 ) {
                log.warn("Species '${meta.id}': ${winners.size()} dRep representatives at 95% ANI; " +
                    "SynTracker reference = ${winners[0]} (highest N50).")
            }
            def reference = winners[0]

            def targets = hq_of.keySet().sort(false, by_frag)
            if ( targets.size() > max_targets ) {
                // keep the reference among the targets
                def kept = ( [ reference ] + targets.findAll { it != reference } ).take(max_targets)
                log.warn("Species '${meta.id}': ${targets.size()} HQ genomes > --syntracker_max_targets (${max_targets}); " +
                    "SynTracker runs on the ${max_targets} with the highest N50, leaving ${targets.size() - max_targets} genome(s) unclustered.")
                targets = kept
            }
            [ meta, hq_of[reference], targets.collect { hq_of[it] } ]
        }
        .filter { meta, reference, targets ->
            def ok = targets.size() >= 2
            if ( !ok ) {
                log.warn("Species '${meta.id}': fewer than 2 HQ genomes; skipping SynTracker.")
            }
            ok
        }

    SYNTRACKER_RUN( ch_selected )

    // The APSS table(s) to cluster: every per-n table for `auto`, else the one requested depth.
    def regions = params.syntracker_regions.toString()
    def wanted  = { name ->
        regions == 'auto' ? ( name ==~ /avg_synteny_scores_\d+_regions\.csv/ ) : name == "avg_synteny_scores_${regions}_regions.csv"
    }
    // Target names for the clustering graph: SynTracker names a target by its file basename minus
    // extension (the module stages every target as <sample>.fasta), as normalise_genome_id does.
    ch_targets = ch_selected
        .collectFile { item ->
            def (meta, reference, targets) = item
            def names = targets.collect { it.name.replaceFirst(/\.gz$/, '').replaceFirst(/\.[^.]+$/, '') }
            [ "${meta.id}_targets.txt", names.join('\n') + '\n' ]
        }
        .map { f -> [ f.name - '_targets.txt', f ] }
    ch_clusters_in = SYNTRACKER_RUN.out.apss
        .map { meta, files -> [ meta.id, meta, [ files ].flatten().findAll { wanted.call(it.name) } ] }
        .filter { id, meta, apss ->
            if ( !apss ) {
                log.warn("Species '${meta.id}': SynTracker produced no APSS table for --syntracker_regions ${regions}; skipping its clustering.")
            }
            apss as boolean
        }
        .join( ch_targets )
        .map { id, meta, apss, targets -> [ meta, apss, targets ] }
    SYNTRACKER_CLUSTERS( ch_clusters_in )

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
        .map { id, meta, table, ani, labels -> [ meta, table, ani, labels ] }
    SYNTRACKER_EVALUATE( ch_eval_in )

    emit:
    tool_metrics = SYNTRACKER_EVALUATE.out.tool_metrics                                          // channel: [ val(meta), path(tool_metrics.tsv) ]
    apss         = SYNTRACKER_RUN.out.apss                                                       // channel: [ val(meta), [ avg_synteny_scores_*.csv ] ]
    clusters     = SYNTRACKER_CLUSTERS.out.clusters                                              // channel: [ val(meta), [ *_clusters.csv ] ]
    qc           = SYNTRACKER_CLUSTERS.out.qc                                                    // channel: [ val(meta), [ *_cluster_qc.tsv ] ]
    retention    = SYNTRACKER_CLUSTERS.out.retention                                             // channel: [ val(meta), [ *_depth_retention.{tsv,png} ] ]
    reference    = ch_selected.map { meta, reference, targets -> [ meta, reference ] }           // channel: [ val(meta), path(reference fasta) ]
}
