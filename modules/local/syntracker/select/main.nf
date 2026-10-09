process SYNTRACKER_SELECT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same Seqera Containers image as syntracker/clusters.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:b51112881a2c3d71' :
        'community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:9b8cabfc8e8a2ced' }"

    input:
    // ani   = all-vs-all FastANI of the species (percent ANI)
    // stats = seqkit stats --all --tabular of the HQ genomes (one row per FASTA, with N50)
    tuple val(meta), path(ani), path(stats)

    output:
    tuple val(meta), path("*_syntracker_reference.txt"), emit: reference
    tuple val(meta), path("*_syntracker_targets.txt")  , emit: targets
    tuple val(meta), path("*_syntracker_selection.tsv"), emit: selection
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args: `--min-n50 100000 --max-targets 150`
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    syntracker_select.py \\
        --ani ${ani} \\
        --stats ${stats} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Reference = the first genome of the stats table, every genome a target.
    """
    tail -n +2 ${stats} | cut -f1 | sed 's#.*/##' > files.txt
    if [ ! -s files.txt ]; then printf 'genome1.fna.gz\\ngenome2.fna.gz\\n' >| files.txt; fi
    head -n 1 files.txt > ${prefix}_syntracker_reference.txt
    cp files.txt ${prefix}_syntracker_targets.txt
    printf 'genome\\tfile\\tN50\\tcentrality\\trole\\n' > ${prefix}_syntracker_selection.tsv
    """
}
