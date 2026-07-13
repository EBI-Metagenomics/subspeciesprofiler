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
include { POPPUNK_FITMODEL      } from '../../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_LINEAGE_RANKS } from '../../../modules/local/poppunk/lineage_ranks/main'
include { POPPUNK_EVALUATE      } from '../../../modules/local/poppunk/evaluate/main'
include { FASTANI_ALLVSALL      } from '../../../modules/local/fastani_allvsall/main'

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

    // Stage 1: fan out one fit per model family, each paired with its species' QC'd database.
    // meta.id stays the species id (for the join back); meta.model labels the fit. This is a
    // profiler: all four families always run, sweeping the grids set by params. A fit that dies on a
    // degenerate grid point is dropped from the ranking (errorStrategy in conf/modules.config), and a
    // species that no default model resolves is itself informative.

    // Threshold sweep: quantile-derived core-distance cutoffs.
    ch_thr_fits = POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, meta, db ] }
        .combine(
            POPPUNK_QUANTILES.out.quantiles.splitCsv( header: true ).map { meta, row -> [ meta.id, row.quantile, row.core_threshold ] },
            by: 0
        )
        .map { id, meta, db, quantile, thr -> [ meta + [ model: "threshold_q${quantile}" ], db, "threshold --threshold ${thr}" ] }

    // Lineage: a single fit per species holds every rank; split + scored per rank below.
    ch_lin_fits = POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta + [ model: 'lineage' ], db, "lineage --ranks ${params.poppunk_lineage_ranks}" ] }

    // BGMM: sweep the number of mixture components K.
    def bgmm_k = params.poppunk_bgmm_k.toString().tokenize(',')*.trim()
    ch_bgmm_fits = POPPUNK_QCDB.out.qc_db.flatMap { meta, db -> bgmm_k.collect { k -> [ meta + [ model: "bgmm_K${k}" ], db, "bgmm --K ${k}" ] } }

    // DBSCAN: sweep the D x min-cluster-prop grid.
    def dbscan_d   = params.poppunk_dbscan_d.toString().tokenize(',')*.trim()
    def dbscan_mcp = params.poppunk_dbscan_min_cluster_prop.toString().tokenize(',')*.trim()
    ch_dbscan_fits = POPPUNK_QCDB.out.qc_db.flatMap { meta, db ->
        [ dbscan_d, dbscan_mcp ].combinations().collect { d, mcp ->
            [ meta + [ model: "dbscan_D${d}_mcp${mcp}" ], db, "dbscan --D ${d} --min-cluster-prop ${mcp}" ]
        }
    }

    ch_fit_in = ch_thr_fits.mix( ch_lin_fits, ch_bgmm_fits, ch_dbscan_fits )
    POPPUNK_FITMODEL( ch_fit_in )

    // A lineage fit is split into one clustering per rank (lineage_rank1, ...); every other family
    // evaluates its fit directory directly.
    POPPUNK_FITMODEL.out.model
        .branch { meta, model ->
            lineage: meta.model == 'lineage'
            other:   true
        }
        .set { ch_fitted }

    POPPUNK_LINEAGE_RANKS( ch_fitted.lineage )
    ch_lineage_evals = POPPUNK_LINEAGE_RANKS.out.ranks
        .transpose()
        .map { meta, rank_dir -> [ meta + [ model: rank_dir.name ], rank_dir ] }

    ch_to_eval = ch_fitted.other.mix( ch_lineage_evals )

    // Score every fit against the species' FastANI + labels (join back on meta.id)
    ch_ani    = FASTANI_ALLVSALL.out.ani.map { meta, ani    -> [ meta.id, ani ] }
    ch_labels = ch_input.map                 { meta, g, r, labels -> [ meta.id, labels ] }
    ch_eval_in = ch_to_eval
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
        def best = rows.sort { a, b ->
            def sa = (a.tool_structure_score_HQ ?: '').isNumber() ? a.tool_structure_score_HQ.toDouble() : -1d
            def sb = (b.tool_structure_score_HQ ?: '').isNumber() ? b.tool_structure_score_HQ.toDouble() : -1d
            ((order[b.tool_status] ?: -1) <=> (order[a.tool_status] ?: -1)) ?: (sb <=> sa)
        }.first()
        [ id, best ]
    }

    emit:
    tool_metrics    = POPPUNK_EVALUATE.out.tool_metrics // channel: [ val(meta), path(tool_metrics.tsv) ]  (one per fit)
    accepted_models = ch_accepted                       // channel: [ val(species_id), [ all ACCEPT rows ] ]
    best_model      = ch_best                           // channel: [ val(species_id), map(best row) ]
}
