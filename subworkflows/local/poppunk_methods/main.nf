//
// POPPUNK_METHODS: PopPUNK model profiling for one eligible species.
//
// Stage 0 (one-time): build + QC the database, derive core-distance quantiles,
// and compute all-vs-all FastANI (model-independent).
//
// Stage 1: fan out a sweep across three model families (threshold / bgmm /
// dbscan), then refine every dbscan fit with PopPUNK's standard refinement.
// Every fit -- swept or refined -- is scored against FastANI, then ranked.
//
// Every fit the evaluator rates Strong or Moderate additionally gets Microreact
// visualisations, for manual inspection and selection among the worthy models.
//
// Lineage is deliberately absent: it clusters *within* a strain (sub-sub-
// clustering), which is a different question from subspecies structure.
// Refinement is applied to dbscan only, and is never gated: the non-standard
// refine variants (--multi-boundary, --unconstrained) are not used, and a
// species that no model resolves is itself an informative result.
//

include { POPPUNK_CREATEDB  } from '../../../modules/local/poppunk/createdb/main'
include { POPPUNK_QCDB      } from '../../../modules/local/poppunk/qcdb/main'
include { POPPUNK_QUANTILES } from '../../../modules/local/poppunk/quantiles/main'
include { POPPUNK_FITMODEL  } from '../../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_FITMODEL as POPPUNK_REFINE_DBSCAN } from '../../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_EVALUATE  } from '../../../modules/local/poppunk/evaluate/main'
include { POPPUNK_EVALUATE as POPPUNK_EVALUATE_DBSCAN_REFINE } from '../../../modules/local/poppunk/evaluate/main'
include { POPPUNK_VISUALISE } from '../../../modules/local/poppunk/visualise/main'
include { FASTANI_ALLVSALL  } from '../../../modules/local/fastani_allvsall/main'

workflow POPPUNK_METHODS {

    take:
    ch_input // channel: [ val(meta), [ genome files ], path(rfile), path(labels) ]
             // meta.n_genomes (post-QC genome count) gates the BGMM family -- see below.

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
    // profiler: all three families always run, sweeping the grids set by params. A fit that dies on
    // a degenerate grid point is dropped from the ranking (errorStrategy in conf/modules.config).

    // Threshold sweep: quantile-derived core-distance cutoffs.
    ch_thr_fits = POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, meta, db ] }
        .combine(
            POPPUNK_QUANTILES.out.quantiles.splitCsv( header: true ).map { meta, row -> [ meta.id, row.quantile, row.core_threshold ] },
            by: 0
        )
        .map { id, meta, db, quantile, thr -> [ meta + [ model: "threshold_q${quantile}" ], db, "threshold --threshold ${thr}", [] ] }

    // BGMM: sweep the number of mixture components K.
    // BGMM needs the core/accessory distance components to be clearly separated, which in practice
    // only holds for small collections -- PopPUNK recommends it for small sample sets and dbscan for
    // larger ones. So the whole family is skipped for species at or above --poppunk_bgmm_max_genomes,
    // where dbscan (plus its refinement) carries the fit. meta.n_genomes is the post-QC count
    // (n_hq + n_mq, set in workflows/subspeciesprofiler.nf); if a caller supplies no count, run BGMM
    // rather than silently dropping the family.
    def bgmm_k = params.poppunk_bgmm_k.toString().tokenize(',')*.trim()
    ch_bgmm_fits = POPPUNK_QCDB.out.qc_db
        .filter { meta, db ->
            def keep = meta.n_genomes == null || (meta.n_genomes as int) < (params.poppunk_bgmm_max_genomes as int)
            if ( !keep ) {
                log.info("Species '${meta.id}': ${meta.n_genomes} genomes >= --poppunk_bgmm_max_genomes (${params.poppunk_bgmm_max_genomes}); skipping the BGMM sweep.")
            }
            keep
        }
        .flatMap { meta, db -> bgmm_k.collect { k -> [ meta + [ model: "bgmm_K${k}" ], db, "bgmm --K ${k}", [] ] } }

    // DBSCAN: sweep the D x min-cluster-prop grid.
    def dbscan_d   = params.poppunk_dbscan_d.toString().tokenize(',')*.trim()
    def dbscan_mcp = params.poppunk_dbscan_min_cluster_prop.toString().tokenize(',')*.trim()
    ch_dbscan_fits = POPPUNK_QCDB.out.qc_db.flatMap { meta, db ->
        [ dbscan_d, dbscan_mcp ].combinations().collect { d, mcp ->
            [ meta + [ model: "dbscan_D${d}_mcp${mcp}" ], db, "dbscan --D ${d} --min-cluster-prop ${mcp}", [] ]
        }
    }

    ch_fit_in = ch_thr_fits.mix( ch_bgmm_fits, ch_dbscan_fits )
    POPPUNK_FITMODEL( ch_fit_in )

    // Score every swept fit against the species' FastANI + labels (join back on meta.id)
    ch_eval_in = POPPUNK_FITMODEL.out.model
        .map { meta, model -> [ meta.id, meta, model ] }
        .combine( FASTANI_ALLVSALL.out.ani.map { meta, ani -> [ meta.id, ani ] }, by: 0 )
        .combine( ch_input.map { meta, g, r, labels -> [ meta.id, labels ] }, by: 0 )
        .map { id, meta, model, ani, labels -> [ meta, model, ani, labels ] }
    POPPUNK_EVALUATE( ch_eval_in )

    // ---- Standard refinement of every dbscan fit (ungated) ----
    // PopPUNK's dbscan always designates some points as noise, so the boundary is worth refining;
    // the other families are left as fitted. `refine` seeds from the dbscan fit directory, which the
    // fitmodel module stages as `seed_model`.
    ch_dbscan_refine_in = POPPUNK_FITMODEL.out.model
        .filter { meta, dir -> meta.model.startsWith('dbscan_') }
        .map { meta, dir -> [ meta.id, meta, dir ] }
        .combine( POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, db ] }, by: 0 )
        .map { id, meta, dir, db -> [ meta + [ model: "refine_from_${meta.model}" ], db, 'refine', dir ] }
    POPPUNK_REFINE_DBSCAN( ch_dbscan_refine_in )

    ch_refine_eval_in = POPPUNK_REFINE_DBSCAN.out.model
        .map { meta, model -> [ meta.id, meta, model ] }
        .combine( FASTANI_ALLVSALL.out.ani.map { meta, ani -> [ meta.id, ani ] }, by: 0 )
        .combine( ch_input.map { meta, g, r, labels -> [ meta.id, labels ] }, by: 0 )
        .map { id, meta, model, ani, labels -> [ meta, model, ani, labels ] }
    POPPUNK_EVALUATE_DBSCAN_REFINE( ch_refine_eval_in )

    // ---- Final selection over the swept + dbscan-refined fits ----
    ch_tool_metrics = POPPUNK_EVALUATE.out.tool_metrics
        .mix( POPPUNK_EVALUATE_DBSCAN_REFINE.out.tool_metrics )
    ch_ranked = ch_tool_metrics
        .splitCsv( header: true, sep: '\t' )
        .map { meta, row -> [ meta.id, row ] }
        .groupTuple()

    // Report ALL accepted models (decision == ACCEPT) -- if several fits are Strong, all are kept.
    ch_accepted = ch_ranked.map { id, rows -> [ id, rows.findAll { it.decision == 'ACCEPT' } ] }

    // The single best-ranked fit (tool_status, then HQ structure score) -- the fallback when nothing
    // is accepted.
    ch_best = ch_ranked.map { id, rows ->
        def order = [ 'Strong': 3, 'Moderate': 2, 'Mixed': 1, 'Weak': 0 ]
        def best = rows.sort { a, b ->
            def sa = (a.tool_structure_score_HQ ?: '').isNumber() ? a.tool_structure_score_HQ.toDouble() : -1d
            def sb = (b.tool_structure_score_HQ ?: '').isNumber() ? b.tool_structure_score_HQ.toDouble() : -1d
            ((order[b.tool_status] ?: -1) <=> (order[a.tool_status] ?: -1)) ?: (sb <=> sa)
        }.first()
        [ id, best ]
    }

    // ---- Microreact visuals for every model worth inspecting ----
    // The evaluator's `decision` is deliberately not used here: it is ACCEPT only for Strong,
    // whereas the point of these visuals is to eyeball and choose among ALL the credible models.
    // So the gate is `tool_status` in (Strong, Moderate).
    ch_viz_wanted = ch_ranked
        .flatMap { id, rows ->
            rows.findAll { it.tool_status in [ 'Strong', 'Moderate' ] }.collect { [ id, it.model ] }
        }

    // Pair each wanted model name with its own fit directory (swept or dbscan-refined) and the
    // species' QC'd db. Same combine + filter idiom used elsewhere to rejoin a model name to its dir.
    ch_fit_dirs = POPPUNK_FITMODEL.out.model
        .mix( POPPUNK_REFINE_DBSCAN.out.model )
        .map { meta, dir -> [ meta.id, meta.model, dir ] }

    ch_viz_in = ch_viz_wanted
        .combine( ch_fit_dirs, by: 0 )
        .filter { id, wanted, candidate, dir -> wanted == candidate }
        .combine( POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, meta, db ] }, by: 0 )
        .map { id, wanted, candidate, dir, meta, db -> [ meta + [ model: wanted ], db, dir ] }
    POPPUNK_VISUALISE( ch_viz_in )


    emit:
    tool_metrics    = ch_tool_metrics                   // channel: [ val(meta), path(tool_metrics.tsv) ]  (swept + dbscan-refined)
    accepted_models = ch_accepted                       // channel: [ val(species_id), [ all ACCEPT rows ] ]
    best_model      = ch_best                           // channel: [ val(species_id), map(best row) ]
    microreact      = POPPUNK_VISUALISE.out.microreact  // channel: [ val(meta), path(.microreact) ]  (Strong/Moderate fits)
}
