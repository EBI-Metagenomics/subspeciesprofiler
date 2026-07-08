process SPECIESQC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(qc_tsv)

    output:
    tuple val(meta), path("*_spp_eligibility_report.tsv"), emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args can override QC thresholds, e.g.:
    //   --min_completeness 80 --hq_completeness 90 --max_contamination 5
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    spp_eligibility_from_qc.py \\
        --species_name "${meta.id}" \\
        --qc_tsv ${qc_tsv} \\
        --output ${prefix}_spp_eligibility_report.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'species_name\\tspp_label\\tn_hq\\tn_mq\\tn_eff\\thq_ratio\\tn_total_passing_qc\\tn_discarded\\n' > ${prefix}_spp_eligibility_report.tsv
    printf '${meta.id}\\tSTRONG\\t80\\t40\\t100.0\\t0.667\\t120\\t5\\n' >> ${prefix}_spp_eligibility_report.tsv
    """
}
