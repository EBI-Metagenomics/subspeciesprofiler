/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SPECIESQC              } from '../modules/local/speciesqc/main'
include { POPPUNK_METHODS        } from '../subworkflows/local/poppunk_methods/main'
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
    // MODULE: Assess subspecies-clustering eligibility from each species' genome QC table
    //
    SPECIESQC (
        ch_samplesheet.map { meta, genomes_dir, qc_csv -> [ meta, qc_csv ] }
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
    // Keep only species whose eligibility label is in --qc_filter, folding the
    // label onto meta so it travels downstream. Join is keyed on meta.id because
    // adding spp_label changes the meta map.
    //
    def qc_keep = params.qc_filter.tokenize(',')

    def ch_spp_label = SPECIESQC.out.report
        .splitCsv(header: true, sep: '\t')
        .map { meta, row -> [ meta.id, row.spp_label ] }
    def ch_rfile  = SPECIESQC.out.rfile.map  { meta, rf -> [ meta.id, rf ] }
    def ch_labels = SPECIESQC.out.labels.map { meta, lb -> [ meta.id, lb ] }

    ch_samplesheet
        .map { meta, genomes_dir, qc_csv -> [ meta.id, meta, genomes_dir ] }
        .join(ch_spp_label)
        .join(ch_rfile)
        .join(ch_labels)
        .map { id, meta, genomes_dir, label, rfile, labels ->
            [ meta + [ spp_label: label ], genomes_dir, rfile, labels ]
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
    // SUBWORKFLOW: staged PopPUNK model selection for each eligible species.
    // Expand genomes_dir into the assembly files (staged as genomes/); carry the
    // r-file + HQ/MQ labels through.
    //
    ch_poppunk_in = ch_eligibility.pass.map { meta, genomes_dir, rfile, labels ->
        [ meta, files("${genomes_dir}/*", checkIfExists: true), rfile, labels ]
    }
    POPPUNK_METHODS( ch_poppunk_in )

    // Per-species model report (Phase-4 "profiler history"): append every fitted model's
    // verdict into one table under that species' poppunk directory. collectFile reads each
    // per-fit tool_metrics by content, so the identical filenames across fits don't collide.
    POPPUNK_METHODS.out.tool_metrics
        .collectFile(keepHeader: true, skip: 1, sort: true, storeDir: params.outdir) { meta, tsv ->
            [ "${meta.id}/poppunk/${meta.id}_model_report.tsv", tsv ]
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
