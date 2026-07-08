process POPPUNK_CREATEDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("${prefix}")          , emit: db
    tuple val(meta), path("${prefix}_rlist.txt"), emit: rlist
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version"), topic: versions, emit: versions_poppunk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    for fasta in ${assemblies}; do
        sample=\$(basename "\$fasta")
        sample="\${sample%.*}"
        printf "%s\\t%s\\n" "\$sample" "\$fasta" >> ${prefix}_rlist.txt
    done

    poppunk \\
        --create-db \\
        --r-files ${prefix}_rlist.txt \\
        --output ${prefix} \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_rlist.txt
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.h5
    touch ${prefix}/${prefix}.dists.npy
    touch ${prefix}/${prefix}.dists.pkl
    """
}
