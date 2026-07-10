process POPPUNK_QUANTILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(qc_db)

    output:
    tuple val(meta), path("*_core_quantiles.csv"), emit: quantiles
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args can override the quantile list, e.g. --quantiles 0.05,0.1,0.2
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    poppunk_core_quantiles.py \\
        --dists ${qc_db}/*.dists.npy \\
        --output ${prefix}_core_quantiles.csv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'quantile,core_threshold\\n0.05,0.000608\\n0.1,0.001198\\n0.2,0.009238\\n' > ${prefix}_core_quantiles.csv
    """
}
