process ANI_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // PopPUNK's image, reused for its matplotlib + scipy + pandas (no other image has all three).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    // ani    = all-vs-all FastANI of the species (percent ANI)
    // labels = SPECIESQC `genome,label` (HQ/MQ) table, for the strip beside the heatmap
    tuple val(meta), path(ani), path(labels)

    output:
    tuple val(meta), path("*_ani_overview.png"), emit: png
    tuple val(meta), path("*_ani_overview.pdf"), emit: pdf
    tuple val(meta), path("*_ani_overview.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('matplotlib'), eval("python -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib
    tuple val("${task.process}"), val('scipy'), eval("python -c 'import scipy; print(scipy.__version__)'"), topic: versions, emit: versions_scipy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def labels_arg = labels ? "--labels ${labels}" : ''
    """
    # matplotlib and fontconfig want a writable cache dir; the container home may not be one
    export MPLCONFIGDIR="\$PWD/.mplconfig" XDG_CACHE_HOME="\$PWD/.cache"

    plot_ani.py \\
        --ani ${ani} \\
        ${labels_arg} \\
        --title ${meta.id} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ani_overview.png ${prefix}_ani_overview.pdf
    printf 'n_genomes\\tn_pairs\\tn_missing_pairs\\tani_min\\tani_median\\tani_max\\tcolour_min\\n' > ${prefix}_ani_overview.tsv
    """
}
