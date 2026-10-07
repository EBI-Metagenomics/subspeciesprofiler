process SYNTRACKER_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    // apss     = one SynTracker avg_synteny_scores_*.csv (all-regions by default)
    // labels   = per-genome HQ/MQ labels (genome,label) from speciesqc
    // drep_cdb = dRep Cdb.csv when SynTracker ran on dRep representatives (propagates each
    //            representative's cluster to its dRep cluster); pass [] to skip
    tuple val(meta), path(apss), path(labels), path(drep_cdb)

    output:
    tuple val(meta), path("*_syntracker_*_clusters.csv"), emit: clusters
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args sets the sweep, e.g. `--thresholds 0.70,0.72,0.74 --method average`.
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cdb    = drep_cdb ? "--drep-cdb ${drep_cdb}" : ''
    """
    syntracker_apss_clusters.py \\
        --apss ${apss} \\
        --labels ${labels} \\
        ${cdb} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Same file names as the real run: one table per threshold, tagged by method.
    def thresholds = (args =~ /--thresholds\s+(\S+)/) ? (args =~ /--thresholds\s+(\S+)/)[0][1].tokenize(',') : [ '0.80' ]
    def tag        = args.contains('--method connected') ? 'cc' : 'avg'
    def files      = thresholds.collect { t -> "${prefix}_syntracker_${tag}_apss${t}_clusters.csv" }.join(' ')
    """
    for f in ${files}; do
        printf 'Taxon,Cluster\\n' > "\$f"
    done
    """
}
