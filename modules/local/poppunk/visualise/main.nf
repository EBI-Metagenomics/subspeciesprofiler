process POPPUNK_VISUALISE {
    tag "${meta.id} - ${meta.model}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    // db  = the species' QC'd PopPUNK database (`<id>_qc`)
    // fit = the fit directory of ONE model (`<id>_fitmodel`).
    // NEITHER may be renamed on staging: poppunk_visualise resolves a model as
    // `<model-dir>/<basename(model-dir)>_fit.pkl`, so a `stageAs` that changes the
    // directory name makes it fail with "Unable to locate previous model fit".
    tuple val(meta), path(db), path(fit)

    output:
    tuple val(meta), path("${prefix}_microreact_clusters.csv"), emit: clusters
    tuple val(meta), path("${prefix}.microreact")             , emit: microreact
    tuple val(meta), path("${prefix}_core_NJ.nwk")            , emit: tree     , optional: true
    tuple val(meta), path("${prefix}*_mandrake.dot")          , emit: embedding, optional: true
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version"), topic: versions, emit: versions_poppunk

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args can add e.g. `--perplexity 15` (mandrake needs perplexity < n genomes),
    // `--info-csv`, or `--api-key <key>` to upload the instance to microreact.org automatically.
    // Pass an API key as a Nextflow secret or env var -- never as a plain pipeline param.
    def args = task.ext.args ?: ''
    // The model name (family + its swept params, e.g. `dbscan_D5_mcp0.01`) goes in the prefix, so
    // every file is self-describing once published flat into <outdir>/<species>/poppunk/microreact.
    prefix   = task.ext.prefix ?: "${meta.id}_${meta.model}"
    """
    poppunk_visualise \\
        --ref-db ${db} \\
        --model-dir ${fit} \\
        --output ${prefix} \\
        --microreact \\
        --rapidnj rapidnj \\
        --threads $task.cpus \\
        $args

    # poppunk_visualise writes <prefix>/<prefix>_*; lift those up so they publish as flat,
    # uniquely-named files into the shared per-species microreact directory.
    mv ${prefix}/* ./
    rm -rf ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_${meta.model}"
    """
    echo "$args"

    printf 'id,Cluster__autocolour\\ngenome_a,1\\ngenome_b,1\\ngenome_c,2\\n' > ${prefix}_microreact_clusters.csv
    printf '(genome_a:0.01,genome_b:0.01,genome_c:0.02);\\n' > ${prefix}_core_NJ.nwk
    printf '{}\\n' > ${prefix}.microreact
    printf 'graph {\\n}\\n' > ${prefix}_perplexity20.0_accessory_mandrake.dot
    """
}
