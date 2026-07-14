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
include { POPPUNK_FITMODEL as POPPUNK_REFINE } from '../../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_FITMODEL as POPPUNK_MULTIBOUNDARY } from '../../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_LINEAGE_RANKS } from '../../../modules/local/poppunk/lineage_ranks/main'
include { POPPUNK_EVALUATE      } from '../../../modules/local/poppunk/evaluate/main'
include { POPPUNK_EVALUATE as POPPUNK_EVALUATE_REFINE } from '../../../modules/local/poppunk/evaluate/main'
include { POPPUNK_EVALUATE as POPPUNK_EVALUATE_MULTIBOUNDARY } from '../../../modules/local/poppunk/evaluate/main'
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
        .map { id, meta, db, quantile, thr -> [ meta + [ model: "threshold_q${quantile}" ], db, "threshold --threshold ${thr}", [] ] }

    // Lineage: a single fit per species holds every rank; split + scored per rank below.
    ch_lin_fits = POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta + [ model: 'lineage' ], db, "lineage --ranks ${params.poppunk_lineage_ranks}", [] ] }

    // BGMM: sweep the number of mixture components K.
    def bgmm_k = params.poppunk_bgmm_k.toString().tokenize(',')*.trim()
    ch_bgmm_fits = POPPUNK_QCDB.out.qc_db.flatMap { meta, db -> bgmm_k.collect { k -> [ meta + [ model: "bgmm_K${k}" ], db, "bgmm --K ${k}", [] ] } }

    // DBSCAN: sweep the D x min-cluster-prop grid.
    def dbscan_d   = params.poppunk_dbscan_d.toString().tokenize(',')*.trim()
    def dbscan_mcp = params.poppunk_dbscan_min_cluster_prop.toString().tokenize(',')*.trim()
    ch_dbscan_fits = POPPUNK_QCDB.out.qc_db.flatMap { meta, db ->
        [ dbscan_d, dbscan_mcp ].combinations().collect { d, mcp ->
            [ meta + [ model: "dbscan_D${d}_mcp${mcp}" ], db, "dbscan --D ${d} --min-cluster-prop ${mcp}", [] ]
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

    // Cheap-sweep ranking -- drives the refine gate below.
    ch_ranked_cheap = POPPUNK_EVALUATE.out.tool_metrics
        .splitCsv( header: true, sep: '\t' )
        .map { meta, row -> [ meta.id, row ] }
        .groupTuple()

    // ---- Stage 3: gated refine escalation ----
    // Ladder gate: a species escalates only if NO cheap fit was accepted (Strong). Then refine the
    // top-N candidates by HQ structure score, POOLED across the non-lineage families (lineage ranks
    // aren't boundary models, so they can't seed refine). An empty seed list = no escalation.
    ch_refine_seeds = ch_ranked_cheap.map { id, rows ->
        def accepted = rows.any { it.decision == 'ACCEPT' }
        def seeds = []
        if ( !accepted ) {
            seeds = rows
                .findAll { r -> !r.model.startsWith('lineage') && (r.tool_structure_score_HQ ?: '').isNumber() && r.tool_structure_score_HQ.toDouble() > 0 }
                .sort { a, b -> b.tool_structure_score_HQ.toDouble() <=> a.tool_structure_score_HQ.toDouble() }
                .take( params.poppunk_refine_top_n as int )
                .collect { it.model }
        }
        [ id, seeds ]
    }

    // Pair each seed model-name with its fit directory and the species' QC'd db, then refine it.
    ch_refine_in = ch_refine_seeds
        .flatMap { id, names -> names.collect { name -> [ id, name ] } }
        .combine( POPPUNK_FITMODEL.out.model.map { meta, dir -> [ meta.id, meta.model, dir ] }, by: 0 )
        .filter { it[2] == it[1] }
        .combine( POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, meta, db ] }, by: 0 )
        .map { id, seed_name, cand_model, cand_dir, meta, db ->
            [ meta + [ model: "refine_from_${seed_name}" ], db, 'refine', cand_dir ]
        }
    POPPUNK_REFINE( ch_refine_in )

    // Score the refine fits (same evaluator; join back on species id).
    ch_refine_eval_in = POPPUNK_REFINE.out.model
        .map { meta, model -> [ meta.id, meta, model ] }
        .combine( FASTANI_ALLVSALL.out.ani.map { meta, ani -> [ meta.id, ani ] }, by: 0 )
        .combine( ch_input.map { meta, g, r, labels -> [ meta.id, labels ] }, by: 0 )
        .map { id, meta, model, ani, labels -> [ meta, model, ani, labels ] }
    POPPUNK_EVALUATE_REFINE( ch_refine_eval_in )

    // ---- Stage 4: gated multi-boundary refine ----
    // Gate on the post-refine ranking: escalate only if still no ACCEPT. `--multi-boundary` sweeps
    // several boundary positions around each seed (each position is its own clustering), seeded from
    // the best 1-2 boundary-model candidates so far (cheap or refine; lineage excluded).
    ch_ranked_refined = POPPUNK_EVALUATE.out.tool_metrics
        .mix( POPPUNK_EVALUATE_REFINE.out.tool_metrics )
        .splitCsv( header: true, sep: '\t' )
        .map { meta, row -> [ meta.id, row ] }
        .groupTuple()

    ch_mb_seeds = ch_ranked_refined.map { id, rows ->
        def accepted = rows.any { it.decision == 'ACCEPT' }
        def seeds = []
        if ( !accepted ) {
            seeds = rows
                .findAll { r -> !r.model.startsWith('lineage') && (r.tool_structure_score_HQ ?: '').isNumber() && r.tool_structure_score_HQ.toDouble() > 0 }
                .sort { a, b -> b.tool_structure_score_HQ.toDouble() <=> a.tool_structure_score_HQ.toDouble() }
                .take( params.poppunk_multiboundary_top_n as int )
                .collect { it.model }
        }
        [ id, seeds ]
    }

    // Candidate fit dirs that can seed a refine: cheap fits + Stage-3 refine fits.
    ch_seed_dirs = POPPUNK_FITMODEL.out.model
        .mix( POPPUNK_REFINE.out.model )
        .map { meta, dir -> [ meta.id, meta.model, dir ] }

    ch_mb_in = ch_mb_seeds
        .flatMap { id, names -> names.collect { name -> [ id, name ] } }
        .combine( ch_seed_dirs, by: 0 )
        .filter { it[2] == it[1] }
        .combine( POPPUNK_QCDB.out.qc_db.map { meta, db -> [ meta.id, meta, db ] }, by: 0 )
        .map { id, seed_name, cand_model, cand_dir, meta, db ->
            [ meta + [ model: "multiboundary_from_${seed_name}" ], db, "refine --multi-boundary ${params.poppunk_multiboundary_n}", cand_dir ]
        }
    POPPUNK_MULTIBOUNDARY( ch_mb_in )

    // Each boundary position writes its own <prefix>_boundary<K>_clusters.csv; fan them all to eval.
    ch_mb_eval_in = POPPUNK_MULTIBOUNDARY.out.model
        .flatMap { meta, dir ->
            files( "${dir}/*_boundary*_clusters.csv" ).collect { f ->
                def k = (f.name =~ /_boundary(\d+)_clusters\.csv/)[0][1]
                [ meta + [ model: "${meta.model}_b${k}" ], f ]
            }
        }
        .map { meta, clusters -> [ meta.id, meta, clusters ] }
        .combine( FASTANI_ALLVSALL.out.ani.map { meta, ani -> [ meta.id, ani ] }, by: 0 )
        .combine( ch_input.map { meta, g, r, labels -> [ meta.id, labels ] }, by: 0 )
        .map { id, meta, clusters, ani, labels -> [ meta, clusters, ani, labels ] }
    POPPUNK_EVALUATE_MULTIBOUNDARY( ch_mb_eval_in )

    // ---- Final selection over cheap + refine + multi-boundary fits ----
    ch_tool_metrics = POPPUNK_EVALUATE.out.tool_metrics
        .mix( POPPUNK_EVALUATE_REFINE.out.tool_metrics )
        .mix( POPPUNK_EVALUATE_MULTIBOUNDARY.out.tool_metrics )
    ch_ranked = ch_tool_metrics
        .splitCsv( header: true, sep: '\t' )
        .map { meta, row -> [ meta.id, row ] }
        .groupTuple()

    // Report ALL accepted models (decision == ACCEPT) -- if several fits are Strong, all are kept.
    ch_accepted = ch_ranked.map { id, rows -> [ id, rows.findAll { it.decision == 'ACCEPT' } ] }

    // The single best-ranked fit (tool_status, then HQ structure score) -- the fallback when nothing
    // is accepted, and the seed for the later refine rungs.
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
    tool_metrics    = ch_tool_metrics                   // channel: [ val(meta), path(tool_metrics.tsv) ]  (cheap + refine)
    accepted_models = ch_accepted                       // channel: [ val(species_id), [ all ACCEPT rows ] ]
    best_model      = ch_best                           // channel: [ val(species_id), map(best row) ]
}
