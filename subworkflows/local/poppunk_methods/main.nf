//
// POPPUNK_METHODS: staged PopPUNK model-selection for one eligible species.
//
// Stage 0 (one-time): build + QC the database, derive core-distance quantiles,
// and compute all-vs-all FastANI (model-independent). Stage 1A: fan a
// core-distance threshold sweep out over the quantiles, fit each with PopPUNK,
// score every fit against FastANI, then rank and select the best model.
//
// Cheap models are run in parallel and scored; escalation to the expensive
// refine stages (gated on no acceptable model) is a later phase.
//

include { POPPUNK_CREATEDB  } from '../../../modules/local/poppunk/createdb/main'
include { POPPUNK_QCDB      } from '../../../modules/local/poppunk/qcdb/main'
include { POPPUNK_QUANTILES } from '../../../modules/local/poppunk/quantiles/main'
include { POPPUNK_FITMODEL  } from '../../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_EVALUATE  } from '../../../modules/local/poppunk/evaluate/main'
include { FASTANI_ALLVSALL  } from '../../../modules/local/fastani_allvsall/main'

workflow POPPUNK_METHODS {

    take:
    ch_input // channel: [ val(meta), [ genome files ], path(rfile), path(labels) ]

    main:

    // Stage 0: database + QC + core-distance quantiles
    ch_db_in = ch_input.map { meta, genomes, rfile, labels -> [ meta, genomes, rfile ] }
    POPPUNK_CREATEDB( ch_db_in )
    POPPUNK_QCDB( POPPUNK_CREATEDB.out.db )
    POPPUNK_QUANTILES( POPPUNK_QCDB.out.qc_db )

    // All-vs-all ANI: model-independent, computed once per species
    FASTANI_ALLVSALL( ch_db_in )

    // Stage 1A: core-distance threshold sweep. Fan out one fit per quantile,
    // pairing each threshold with its species' QC'd database (join on meta.id).
    ch_qcdb = POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, meta, db ] }
    ch_thresholds = POPPUNK_QUANTILES.out.quantiles
        .splitCsv( header: true )
        .map { meta, row -> [ meta.id, row.quantile, row.core_threshold ] }

    ch_fit_in = ch_qcdb
        .combine( ch_thresholds, by: 0 )
        .map { id, meta, db, quantile, thr ->
            // meta.id stays the species id (for the join back); meta.model labels the fit.
            [ meta + [ model: "threshold_q${quantile}" ], db, "threshold --threshold ${thr}" ]
        }
    POPPUNK_FITMODEL( ch_fit_in )

    // Score every fit against the species' FastANI + labels (join back on meta.id)
    ch_ani    = FASTANI_ALLVSALL.out.ani.map { meta, ani    -> [ meta.id, ani ] }
    ch_labels = ch_input.map                 { meta, g, r, labels -> [ meta.id, labels ] }
    ch_eval_in = POPPUNK_FITMODEL.out.model
        .map { meta, model -> [ meta.id, meta, model ] }
        .combine( ch_ani, by: 0 )
        .combine( ch_labels, by: 0 )
        .map { id, meta, model, ani, labels -> [ meta, model, ani, labels ] }
    POPPUNK_EVALUATE( ch_eval_in )

    // Select: gather each species' fits and rank them.
    ch_ranked = POPPUNK_EVALUATE.out.tool_metrics
        .splitCsv( header: true, sep: '\t' )
        .map { meta, row -> [ meta.id, row ] }
        .groupTuple()

    // Report ALL accepted models (decision == ACCEPT) -- if several fits are Strong in one
    // sweep, all of them are reported, not just one.
    ch_accepted = ch_ranked.map { id, rows -> [ id, rows.findAll { it.decision == 'ACCEPT' } ] }

    // The single best-ranked fit (tool_status, then HQ structure score) -- the fallback when
    // nothing is accepted, and the seed for the later refine escalation.
    ch_best = ch_ranked.map { id, rows ->
        def order = [ 'Strong': 3, 'Moderate': 2, 'Mixed': 1, 'Weak': 0 ]
        def score = { r -> (r.tool_structure_score_HQ ?: '').isNumber() ? r.tool_structure_score_HQ.toDouble() : -1d }
        def best = rows.sort { a, b ->
            ((order[b.tool_status] ?: -1) <=> (order[a.tool_status] ?: -1)) ?: (score(b) <=> score(a))
        }.first()
        [ id, best ]
    }

    emit:
    tool_metrics    = POPPUNK_EVALUATE.out.tool_metrics // channel: [ val(meta), path(tool_metrics.tsv) ]  (one per fit)
    accepted_models = ch_accepted                       // channel: [ val(species_id), [ all ACCEPT rows ] ]
    best_model      = ch_best                           // channel: [ val(species_id), map(best row) ]
}
