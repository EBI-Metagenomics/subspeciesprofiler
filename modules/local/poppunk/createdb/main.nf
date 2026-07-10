process POPPUNK_CREATEDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    // `genomes` are the HQ/MQ assemblies, staged into `genomes/` so the r-file's
    // relative `./genomes/<file>` paths (written by the speciesqc module) resolve.
    tuple val(meta), path(genomes, stageAs: 'genomes/*'), path(rfile)

    output:
    tuple val(meta), path("${prefix}"), emit: db
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version"), topic: versions, emit: versions_poppunk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    poppunk \\
        --create-db \\
        --r-files ${rfile} \\
        --output ${prefix} \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.h5
    touch ${prefix}/${prefix}.dists.npy
    touch ${prefix}/${prefix}.dists.pkl
    """
}
