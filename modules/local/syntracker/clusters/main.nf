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
    // apss            = SynTracker avg_synteny_scores_all_regions.csv (what is clustered)
    // apss_subsampled = one subsampled table (avg_synteny_scores_40_regions.csv), used only to
    //                   estimate the per-region score SD; [] to use the default
    // targets         = SynTracker target sample names, one per line (each target is a graph node)
    tuple val(meta), path(apss), path(apss_subsampled), path(targets)

    output:
    tuple val(meta), path("*_syntracker_*_clusters.csv")         , emit: clusters
    tuple val(meta), path("*_syntracker_genome_coverage.tsv")    , emit: coverage
    tuple val(meta), path("*_syntracker_noise.tsv")              , emit: noise
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('python-igraph'), eval("python -c 'import igraph; print(igraph.__version__)'"), topic: versions, emit: versions_igraph

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args sets the filters and the sweep, e.g.
    //   `--min-regions 100 --min-genome-coverage 0.4 --min-apss 0.70,0.75 --resolutions 0.25,0.5,1.0 --seed 42`
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def subsampled = apss_subsampled ? "--apss-subsampled ${apss_subsampled}" : ''
    """
    syntracker_apss_clusters.py \\
        --apss ${apss} \\
        ${subsampled} \\
        --targets ${targets} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Same file names as the real run: one table per (min_apss, resolution). Targets are spread
    // over two clusters.
    def tokens = args.tokenize(' ')
    def opt    = { name, fallback -> def i = tokens.indexOf("--${name}".toString()); i >= 0 && i + 1 < tokens.size() ? tokens[i + 1] : fallback }
    def models = opt.call('min-apss', '0.75').tokenize(',').collectMany { t ->
        opt.call('resolutions', '1.0').tokenize(',').collect { r -> "${prefix}_syntracker_leiden_apss${t}_r${r}" }
    }.join(' ')
    """
    for m in ${models}; do
        printf 'Taxon,Cluster\\n' > "\$m"_clusters.csv
        # sample names, as normalise_genome_id derives them from the target file names
        sed -e 's/\\.gz\$//' -e 's/\\.\\(fasta\\|fna\\|fa\\|fas\\)\$//' ${targets} | sort | awk 'NF { printf "%s,%d\\n", \$0, (NR % 2) + 1 }' >> "\$m"_clusters.csv
    done
    printf 'genome\\tmedian_regions\\trelative_coverage\\tstatus\\n' > ${prefix}_syntracker_genome_coverage.tsv
    sort ${targets} | awk 'NF { printf "%s\\t300\\t1.0\\tok\\n", \$0 }' >> ${prefix}_syntracker_genome_coverage.tsv
    printf 'region_sd\\tmedian_regions\\tapss_se\\tsd_source\\tsubsampled_n\\tn_pairs\\tn_targets\\tn_low_coverage\\n0.226\\t300\\t0.013\\tdefault\\t\\t1\\t2\\t0\\n' > ${prefix}_syntracker_noise.tsv
    """
}
