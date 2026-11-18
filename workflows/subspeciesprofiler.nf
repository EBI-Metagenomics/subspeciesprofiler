/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SEQKIT_STATS           } from '../modules/nf-core/seqkit/stats/main'
include { SPECIESQC              } from '../modules/local/speciesqc/main'
include { FASTANI_ALLVSALL       } from '../modules/local/fastani_allvsall/main'
include { ANI_QC                 } from '../modules/local/ani/qc/main'
include { POPPUNK_METHODS        } from '../subworkflows/local/poppunk_methods/main'
include { SYNTRACKER_METHODS     } from '../subworkflows/local/syntracker_methods/main'
include { PANGENOME              } from '../subworkflows/local/pangenome/main'
include { ANI_PLOT               } from '../modules/local/ani/plot/main'
include { ANI_TREE               } from '../modules/local/ani/tree/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_subspeciesprofiler_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SUBSPECIESPROFILER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: N50/N90 of every genome in each species' genomes_dir, for the contiguity gate
    //
    SEQKIT_STATS (
        ch_samplesheet.map { meta, genomes_dir, qc_csv -> [ meta, files("${genomes_dir}/*", checkIfExists: true) ] }
    )

    //
    // MODULE: per-genome QC (HQ/MQ by completeness/contamination, then N50 or N90 >= --qc_min_nx)
    // and the species' eligibility label. ANI_QC below makes the final genome list from its r-file.
    //
    SPECIESQC (
        ch_samplesheet
            .map { meta, genomes_dir, qc_csv -> [ meta.id, meta, qc_csv ] }
            .join( SEQKIT_STATS.out.stats.map { meta, stats -> [ meta.id, stats ] } )
            .map { id, meta, qc_csv, stats -> [ meta, qc_csv, stats ] }
    )

    // Combined eligibility report: append every species' single-row report into one table,
    // published at the top of the outdir (alongside the per-species directories).
    def ch_eligibility_report = SPECIESQC.out.report
        .map { meta, report -> report }
        .collectFile(
            name: 'species_eligibility_report.tsv',
            storeDir: params.outdir,
            keepHeader: true,
            skip: 1,
            sort: true,
        )

    //
    // Keep only species whose eligibility label is in --qc_filter, folding the label and the
    // post-QC genome count onto meta so they travel downstream. Join is keyed on meta.id
    // because adding fields changes the meta map.
    //
    def qc_keep = params.qc_filter.tokenize(',')

    // n_total_passing_qc is n_hq + n_mq: exactly the genomes written to the r-file and handed to
    // POPPUNK_CREATEDB, so it is the right count for the BGMM size gate in POPPUNK_METHODS.
    // (n_eff is a weighted HQ/MQ score, not a genome count -- don't use it here.)
    def ch_spp_info = SPECIESQC.out.report
        .splitCsv(header: true, sep: '\t')
        .map { meta, row -> [ meta.id, [ spp_label: row.spp_label, n_genomes: row.n_total_passing_qc as Integer ] ] }
    def ch_rfile  = SPECIESQC.out.rfile.map  { meta, rf -> [ meta.id, rf ] }
    def ch_labels = SPECIESQC.out.labels.map { meta, lb -> [ meta.id, lb ] }

    ch_samplesheet
        .map { meta, genomes_dir, qc_csv -> [ meta.id, meta, genomes_dir ] }
        .join(ch_spp_info)
        .join(ch_rfile)
        .join(ch_labels)
        .map { id, meta, genomes_dir, spp_info, rfile, labels ->
            [ meta + spp_info, genomes_dir, rfile, labels ]
        }
        .branch { meta, genomes_dir, rfile, labels ->
            pass: meta.spp_label in qc_keep
            fail: true
        }
        .set { ch_eligibility }

    // Species filtered out by --qc_filter (available for reporting / logging)
    ch_eligibility.fail
        .subscribe { meta, genomes_dir, rfile, labels -> log.warn("Species '${meta.id}' dropped: spp_label='${meta.spp_label}' not in --qc_filter (${params.qc_filter})") }

    //
    // MODULE: all-vs-all FastANI of the genomes that passed SPECIESQC, once per eligible species.
    //
    FASTANI_ALLVSALL(
        ch_eligibility.pass.map { meta, genomes_dir, rfile, labels ->
            [ meta, files("${genomes_dir}/*", checkIfExists: true), rfile ]
        }
    )

    //
    // MODULE: discard the fewest genomes so that every remaining pair has a FastANI value at or
    // above the minimum ANI. Its outputs are the one set of genomes that every branch below uses.
    //
    ANI_QC(
        FASTANI_ALLVSALL.out.ani
            .map { meta, ani -> [ meta.id, meta, ani ] }
            .join( ch_eligibility.pass.map { meta, genomes_dir, rfile, labels -> [ meta.id, rfile, labels ] } )
            .join( SPECIESQC.out.genome_qc.map { meta, genome_qc -> [ meta.id, genome_qc ] } )
            .map { id, meta, ani, rfile, labels, genome_qc -> [ meta, ani, rfile, labels, genome_qc ] }
    )

    // ANI_QC can lower the genome counts, so apply --qc_filter again with its label and count.
    ch_eligibility.pass
        .map { meta, genomes_dir, rfile, labels -> [ meta.id, meta, genomes_dir ] }
        .join(
            ANI_QC.out.report
                .splitCsv(header: true, sep: '\t')
                .map { meta, row -> [ meta.id, [ spp_label: row.spp_label, n_genomes: row.n_total_passing_qc as Integer ] ] }
        )
        .join( ANI_QC.out.rfile.map { meta, rfile -> [ meta.id, rfile ] } )
        .join( ANI_QC.out.labels.map { meta, labels -> [ meta.id, labels ] } )
        .join( ANI_QC.out.ani.map { meta, ani -> [ meta.id, ani ] } )
        .join( ANI_QC.out.genome_qc.map { meta, genome_qc -> [ meta.id, genome_qc ] } )
        .map { id, meta, genomes_dir, spp_info, rfile, labels, ani, genome_qc ->
            [ meta + spp_info, genomes_dir, rfile, labels, ani, genome_qc ]
        }
        .branch { meta, genomes_dir, rfile, labels, ani, genome_qc ->
            pass: meta.spp_label in qc_keep
            fail: true
        }
        .set { ch_ani_qc }

    ch_ani_qc.fail
        .subscribe { meta, genomes_dir, rfile, labels, ani, genome_qc -> log.warn("Species '${meta.id}' dropped after ANI QC: spp_label='${meta.spp_label}' not in --qc_filter (${params.qc_filter})") }

    def ch_ani       = ch_ani_qc.pass.map { meta, genomes_dir, rfile, labels, ani, genome_qc -> [ meta, ani ] }
    def ch_genome_qc = ch_ani_qc.pass.map { meta, genomes_dir, rfile, labels, ani, genome_qc -> [ meta, genome_qc ] }

    //
    // SUBWORKFLOW: staged PopPUNK model selection for each eligible species.
    // Expand genomes_dir into the assembly files (staged as genomes/); carry the
    // r-file + HQ/MQ labels through. The other branches use the same input channel.
    //
    ch_poppunk_in = ch_ani_qc.pass.map { meta, genomes_dir, rfile, labels, ani, genome_qc ->
        [ meta, files("${genomes_dir}/*", checkIfExists: true), rfile, labels ]
    }
    if ( !params.skip_poppunk ) {
        POPPUNK_METHODS( ch_poppunk_in, ch_ani )
    }

    //
    // MODULE: per-species ANI overview (pairwise ANI histogram + UPGMA-clustered ANI heatmap with
    // an HQ/MQ strip), from the FastANI table of the final genome list.
    //
    ANI_PLOT(
        ch_ani
            .map { meta, ani -> [ meta.id, meta, ani ] }
            .join( ch_poppunk_in.map { meta, genomes, rfile, labels -> [ meta.id, labels ] } )
            .map { id, meta, ani, labels -> [ meta, ani, labels ] }
    )

    //
    // MODULE: per-species neighbour-joining tree (Newick) of the FastANI distances (1 - ANI/100).
    //
    ANI_TREE( ch_ani )

    //
    // SUBWORKFLOW: SynTracker (synteny), an independent non-PopPUNK signal. Its clusterings are
    // scored by the same evaluator, so they join the PopPUNK fits in the model report.
    //
    // Per-species model reports (Phase-4 "profiler history"): every fitted model's verdict, one
    // table per branch, since PopPUNK fits (scored in ANI space) and SynTracker clusterings (scored
    // in APSS space, with ANI concordance) have different columns. collectFile reads each per-fit
    // tool_metrics by content, so identical filenames across fits don't collide. The evaluator's
    // `decision` column is dropped: it only drives control flow (ch_accepted, the MultiQC
    // best-model table) and would duplicate tool_status in the report.
    def model_report = { ch, branch, suffix ->
        ch.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: params.outdir) { meta, tsv ->
            def rows = tsv.readLines().collect { it.split('\t', -1) as List }
            def drop = rows[0].indexOf('decision')
            if ( drop >= 0 ) {
                rows.each { it.remove(drop as int) }
            }
            [ "${meta.id}/${branch}/${meta.id}_${suffix}.tsv", rows*.join('\t').join('\n') + '\n' ]
        }
    }
    if ( !params.skip_poppunk ) {
        model_report.call( POPPUNK_METHODS.out.tool_metrics, 'poppunk', 'model_report' )
    }
    //
    // SUBWORKFLOW: gene presence/absence (ggCaller gene calls -> Panaroo), over the same genomes.
    //
    if ( !params.skip_pangenome ) {
        PANGENOME( ch_poppunk_in )
    }

    if ( !params.skip_syntracker ) {
        SYNTRACKER_METHODS( ch_poppunk_in, ch_ani, ch_genome_qc )
        model_report.call( SYNTRACKER_METHODS.out.tool_metrics, 'syntracker', 'syntracker_model_report' )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'subspeciesprofiler_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    //
    // Subspecies dashboard: cross-species eligibility + selected-model tables as MultiQC
    // custom content (this replaces feeding the raw eligibility TSVs to MultiQC).
    //
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_eligibility_report
            .map { report ->
                "# id: 'subspecies_eligibility'\n" +
                "# section_name: 'Subspecies eligibility'\n" +
                "# description: 'Per-species QC eligibility for subspecies clustering (HQ/MQ counts and label).'\n" +
                "# plot_type: 'table'\n" +
                report.text
            }
            .collectFile(name: 'subspecies_eligibility_mqc.tsv')
    )
    if ( !params.skip_poppunk ) {
        ch_multiqc_files = ch_multiqc_files.mix(
            POPPUNK_METHODS.out.best_model
                .map { species_id, row ->
                    [ species_id, row.model, row.tool_status, row.decision, row.tool_structure_score_HQ ].join('\t') + '\n'
                }
                .collectFile(
                    name: 'subspecies_best_model_mqc.tsv',
                    sort: true,
                    seed: "# id: 'subspecies_best_model'\n# section_name: 'Selected model per species'\n# description: 'Best PopPUNK model chosen per species with its evaluation verdict.'\n# plot_type: 'table'\nspecies\tmodel\ttool_status\tdecision\ttool_structure_score_HQ\n",
                )
        )
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
