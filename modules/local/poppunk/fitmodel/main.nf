process POPPUNK_FITMODEL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("${fit_prefix}")               , emit: model
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version"), topic: versions, emit: versions_poppunk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'bgmm'
    prefix     = task.ext.prefix ?: "${meta.id}"
    fit_prefix = "${prefix}_fitmodel"
    """
    poppunk \\
        --fit-model $args \\
        --ref-db ${db} \\
        --output ${fit_prefix} \\
        --threads $task.cpus
    """

    stub:
    def args = task.ext.args ?: 'bgmm'
    prefix     = task.ext.prefix ?: "${meta.id}"
    fit_prefix = "${prefix}_fitmodel"
    """
    echo $args

    mkdir -p ${fit_prefix}
    touch ${fit_prefix}/${fit_prefix}_fit.npz
    touch ${fit_prefix}/${fit_prefix}_fit.pkl
    touch ${fit_prefix}/${fit_prefix}_DPGMM_fit.png
    touch ${fit_prefix}/${fit_prefix}_DPGMM_fit_contours.pdf
    touch ${fit_prefix}/${fit_prefix}_graph.gt
    touch ${fit_prefix}/${fit_prefix}_clusters.csv
    touch ${fit_prefix}/${fit_prefix}_unword_clusters.csv
    touch ${fit_prefix}/${fit_prefix}.refs
    touch ${fit_prefix}/${fit_prefix}.refs_graph.gt
    touch ${fit_prefix}/${fit_prefix}.refs.dists.npy
    touch ${fit_prefix}/${fit_prefix}.refs.dists.pkl
    touch ${fit_prefix}/${fit_prefix}.refs.h5
    """
}
