process POPPUNK_QCDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("${qc_prefix}")       , emit: qc_db
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version"), topic: versions, emit: versions_poppunk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    qc_prefix = "${prefix}_qc"
    """
    poppunk \\
        --qc-db \\
        --ref-db ${db} \\
        --output ${qc_prefix} \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    qc_prefix = "${prefix}_qc"
    """
    echo $args

    mkdir -p ${qc_prefix}
    touch ${qc_prefix}/${qc_prefix}.h5
    touch ${qc_prefix}/${qc_prefix}.dists.npy
    touch ${qc_prefix}/${qc_prefix}.dists.pkl
    touch ${qc_prefix}/${qc_prefix}_qcreport.txt
    """
}
