process SYNTRACKER_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Seqera Containers image built from environment.yml (pandas, python-igraph, matplotlib-base,
    // procps-ng); the pandas-only biocontainer has no igraph.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:b51112881a2c3d71' :
        'community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:9b8cabfc8e8a2ced' }"

    input:
    // apss    = SynTracker avg_synteny_scores_<n|all>_regions.csv table(s) from one run: every
    //           per-n table for `--depth auto`, or the one table of a fixed depth
    // targets = SynTracker target sample names, one per line (each target is a graph node)
    tuple val(meta), path(apss), path(targets)

    output:
    tuple val(meta), path("*_syntracker_*_clusters.csv")              , emit: clusters
    tuple val(meta), path("*_syntracker_*_cluster_qc.tsv")            , emit: qc
    tuple val(meta), path("*_syntracker_depth_retention.{tsv,png}")   , emit: retention, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('python-igraph'), eval("python -c 'import igraph; print(igraph.__version__)'"), topic: versions, emit: versions_igraph

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args sets the depth and the sweep, e.g.
    //   `--depth auto --min-retention 0.9 --min-apss 0.70,0.75 --resolutions 0.5,1.0 --seed 42`
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    syntracker_apss_clusters.py \\
        --apss ${apss} \\
        --targets ${targets} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Same file names as the real run: one table (and QC table) per (min_apss, resolution), tagged
    // with the depth (n40 when it is chosen automatically). Targets are spread over two clusters.
    def tokens      = args.tokenize(' ')
    def opt         = { name, fallback -> def i = tokens.indexOf("--${name}".toString()); i >= 0 && i + 1 < tokens.size() ? tokens[i + 1] : fallback }
    def depth       = opt.call('depth', 'auto')
    def n           = depth == 'auto' ? '40' : depth
    def models      = opt.call('min-apss', '0.75').tokenize(',').collectMany { t ->
        opt.call('resolutions', '1.0').tokenize(',').collect { r -> "${prefix}_syntracker_leiden_n${n}_apss${t}_r${r}" }
    }.join(' ')
    """
    for m in ${models}; do
        printf 'Taxon,Cluster\\n' > "\$m"_clusters.csv
        sort ${targets} | awk 'NF { printf "%s,%d\\n", \$0, (NR % 2) + 1 }' >> "\$m"_clusters.csv
        printf 'cluster\\tsize\\tmean_intra_apss\\tmax_inter_apss\\tlow_confidence\\treference_genome\\tn\\n' > "\$m"_cluster_qc.tsv
    done
    if [ "${depth}" = auto ]; then
        printf 'n\\tretained_samples\\tretained_pairs\\ttotal_targets\\tselected\\n40\\t0\\t0\\t0\\tTrue\\n' > ${prefix}_syntracker_depth_retention.tsv
        touch ${prefix}_syntracker_depth_retention.png
    fi
    """
}
