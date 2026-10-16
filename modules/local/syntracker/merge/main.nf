process SYNTRACKER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same Seqera Containers image as syntracker/select and syntracker/clusters.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:b51112881a2c3d71' :
        'community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:9b8cabfc8e8a2ced' }"

    input:
    // apss    = the avg_synteny_scores_all_regions.csv of every SynTracker run of the species (one
    //           per chunk pair, all against the same reference); never the subsampled tables
    // targets = every SynTracker target of the species, one per line (for the missing-pair count)
    tuple val(meta), path(apss, stageAs: 'chunk_??/*'), path(targets)

    output:
    tuple val(meta), path("avg_synteny_scores_all_regions.csv"), emit: apss
    tuple val(meta), path("*_syntracker_merge.tsv")            , emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    syntracker_merge_apss.py \\
        --apss ${apss} \\
        --targets ${targets} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // One row per unordered pair, first occurrence kept: the stub runs write identical values.
    """
    awk -F, 'FNR == 1 { if (NR == 1) print; next }
             { key = (\$2 < \$3) ? \$2 FS \$3 : \$3 FS \$2; if (!(key in seen)) { seen[key] = 1; print } }' ${apss} \\
        > avg_synteny_scores_all_regions.csv
    printf 'n_tables\\trows_read\\tduplicate_rows\\tunique_pairs\\tn_targets\\texpected_pairs\\tmissing_pairs\\n' > ${prefix}_syntracker_merge.tsv
    """
}
