process POPPUNK_LINEAGE_RANKS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    // The PopPUNK lineage fit directory (contains <prefix>_lineages.csv with one column per rank).
    tuple val(meta), path(fit)

    output:
    // One directory per rank, each holding a Taxon,Cluster clustering for the ANI evaluator.
    tuple val(meta), path("lineage_rank*", type: 'dir'), emit: ranks
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    poppunk_split_lineage_ranks.py \\
        --lineages \$(ls ${fit}/*_lineages.csv | head -n 1) \\
        --outdir .
    """

    stub:
    """
    for r in 1 2; do
        mkdir -p lineage_rank\${r}
        printf 'Taxon,Cluster\\ngenome_a,1\\ngenome_b,1\\ngenome_c,2\\n' > lineage_rank\${r}/lineage_rank\${r}_clusters.csv
    done
    """
}
