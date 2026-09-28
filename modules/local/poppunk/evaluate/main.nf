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
    # `model` is either a fit directory (PopPUNK writes <p>_clusters.csv + <p>_unword_clusters.csv;
    # evaluate the former) or a single Taxon,Cluster clusters CSV (e.g. one multi-boundary position).
    # A fit directory also carries the fitmodel log (<p>_fit.log), the source of the PopPUNK
    # network score; without one the score is reported as NaN.
    fit_log_arg=""
    if [ -d "${model}" ]; then
        clusters=\$(ls ${model}/*_clusters.csv | grep -v unword | head -n 1)
        fit_log=\$(find ${model}/ -maxdepth 1 -name '*_fit.log' | head -n 1)
        [ -n "\$fit_log" ] && fit_log_arg="--fit-log \$fit_log"
    else
        clusters="${model}"
    fi

    evaluate_poppunk_fastani.py \\
        --fastani ${fastani} \\
        --clusters "\$clusters" \\
        --labels ${labels} \\
        --model-name '${model_name}' \\
        --out-prefix ${prefix} \\
        \$fit_log_arg \\
        ${args}
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def model_name = meta.model ?: meta.id
    // The verdict is varied by model family so that stub runs actually exercise the downstream
    // branches: Strong -> ACCEPT + Microreact, Moderate -> Microreact only, Weak -> neither.
    // A stub that always returned one status would leave those paths untested.
    def status     = model_name.startsWith('refine_from_') ? 'Strong'
                   : model_name.startsWith('bgmm')         ? 'Moderate'
                   : 'Weak'
    def decision   = status == 'Strong' ? 'ACCEPT' : 'TRY_NEXT_MODEL'
    def score      = status == 'Strong' ? '0.85' : ( status == 'Moderate' ? '0.70' : '0.50' )
    """
    printf 'model\\tpoppunk_network_score\\ttool_status\\tn_clusters\\tn_nonsingleton_clusters\\tlargest_cluster_fraction\\ttool_structure_score_HQ\\tdecision\\teval_summary\\n' > ${prefix}.tool_metrics.tsv
    printf '${model_name}\\t0.9\\t${status}\\t3\\t3\\t0.8\\t${score}\\t${decision}\\tstub\\n' >> ${prefix}.tool_metrics.tsv
    touch ${prefix}.cluster_metrics.tsv
    touch ${prefix}.genome_metrics.tsv
    """
}
