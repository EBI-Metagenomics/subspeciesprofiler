process SYNTRACKER_FILTER {
    tag "${meta.id} - ${meta.chunk}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Same Seqera Containers image as syntracker/select, merge and clusters.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:b51112881a2c3d71' :
        'community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:9b8cabfc8e8a2ced' }"

    input:
    // One SynTracker run (one chunk job): its per-region scores, the per-region BLAST hits and
    // the Sample.N -> genome dictionary. Staged under fixed names (the BLAST dir is named after
    // the reference in the run).
    tuple val(meta), path(per_region), path(blast_output, stageAs: 'blast_output'), path(dictionary)

    output:
    tuple val(meta), path("avg_synteny_scores_all_regions.csv"), emit: apss
    tuple val(meta), path("*_syntracker_multihit.tsv")         , emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.chunk}"
    """
    syntracker_filter_multihit.py \\
        --per-region ${per_region} \\
        --blast-dir blast_output \\
        --dictionary ${dictionary} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.chunk}"
    // One APSS row per pair of the per-region table (the stub SynTracker run has no multi-hits).
    """
    echo '"Ref_genome","Sample1","Sample2","APSS","Compared_regions"' > avg_synteny_scores_all_regions.csv
    awk -F, 'NR > 1 && !seen[\$1 FS \$2 FS \$3]++ { printf "%s,%s,%s,0.95,1\\n", \$1, \$2, \$3 }' ${per_region} >> avg_synteny_scores_all_regions.csv
    printf 'multihit_region_genomes\\trows_read\\trows_removed\\tpairs_in\\tpairs_changed\\tpairs_lost\\n' > ${prefix}_syntracker_multihit.tsv
    """
}
