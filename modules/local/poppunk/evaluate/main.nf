process POPPUNK_EVALUATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(model), path(fastani), path(labels)

    output:
    tuple val(meta), path("*.tool_metrics.tsv")   , emit: tool_metrics
    tuple val(meta), path("*.cluster_metrics.tsv"), emit: cluster_metrics
    tuple val(meta), path("*.genome_metrics.tsv") , emit: genome_metrics
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args can set evaluation thresholds, e.g.:
    //   --min-ani 0.90 --accept-status Strong --min-comparison-cluster-size 2
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def model_name = meta.model ?: meta.id
    """
    # PopPUNK writes both <p>_clusters.csv (Taxon,Cluster) and <p>_unword_clusters.csv;
    # evaluate the former.
    clusters=\$(ls ${model}/*_clusters.csv | grep -v unword | head -n 1)

    evaluate_poppunk_fastani.py \\
        --fastani ${fastani} \\
        --clusters "\$clusters" \\
        --labels ${labels} \\
        --model-name '${model_name}' \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def model_name = meta.model ?: meta.id
    """
    printf 'model\\ttool_status\\tdecision\\treason\\ttool_structure_score_HQ\\n' > ${prefix}.tool_metrics.tsv
    printf '${model_name}\\tWeak\\tTRY_NEXT_MODEL\\tstub\\t0.5\\n' >> ${prefix}.tool_metrics.tsv
    touch ${prefix}.cluster_metrics.tsv
    touch ${prefix}.genome_metrics.tsv
    """
}
