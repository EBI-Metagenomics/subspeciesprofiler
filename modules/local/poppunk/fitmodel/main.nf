process POPPUNK_FITMODEL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    // fit_args = the model type + its params for THIS task, e.g. "threshold --threshold 0.0012",
    // "bgmm --K 4", "dbscan --D 5 --min-cluster-prop 0.01",
    // or "refine" (standard refinement of a prior dbscan fit). Supplied per-task so one module
    // covers the whole model sweep; task.ext.args adds any common flags.
    // model_dir = optional starting-model directory for `refine` (a prior fit's output). Pass `[]`
    // for the from-scratch families (threshold/bgmm/dbscan). It must NOT be renamed on staging:
    // PopPUNK resolves the seed as `<model-dir>/<basename(model-dir)>_fit.pkl`, so a renaming
    // `stageAs` makes every refine fail on start-up. Staging it under `seed/` keeps its name while
    // avoiding a clash with this task's own `<prefix>_fitmodel` output directory.
    tuple val(meta), path(db), val(fit_args), path(model_dir, stageAs: 'seed/*')

    output:
    tuple val(meta), path("${fit_prefix}")               , emit: model
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version"), topic: versions, emit: versions_poppunk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seed = model_dir ? "--model-dir ${model_dir}" : ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    fit_prefix = "${prefix}_fitmodel"
    // PopPUNK's stderr is kept as <fit_prefix>_fit.log inside the fit directory: its "Network
    // summary" (network score) is read from there by poppunk/evaluate. On failure the log is
    // echoed and PopPUNK's own exit status is kept, which the errorStrategy relies on.
    """
    poppunk \\
        --fit-model ${fit_args} \\
        --ref-db ${db} \\
        ${seed} \\
        --output ${fit_prefix} \\
        --threads $task.cpus \\
        $args \\
        2> poppunk_fit.log || { rc=\$?; cat poppunk_fit.log >&2; exit \$rc; }

    cat poppunk_fit.log >&2
    mv poppunk_fit.log ${fit_prefix}/${fit_prefix}_fit.log
    """

    stub:
    def args = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    fit_prefix = "${prefix}_fitmodel"
    def seed_check = model_dir ? "test -f ${model_dir}/\$(basename ${model_dir})_fit.pkl" : ''
    """
    echo "${fit_args} ${args}"
    # Fail like PopPUNK would if the refine seed was renamed on staging
    ${seed_check}

    mkdir -p ${fit_prefix}
    printf 'Network summary:\\n\\tComponents\\t\\t\\t\\t2\\n\\tScore\\t\\t\\t\\t\\t0.9000\\n\\tScore (w/ betweenness)\\t\\t\\t0.8000\\n' > ${fit_prefix}/${fit_prefix}_fit.log
    touch ${fit_prefix}/${fit_prefix}_fit.npz
    touch ${fit_prefix}/${fit_prefix}_fit.pkl
    touch ${fit_prefix}/${fit_prefix}_DPGMM_fit.png
    touch ${fit_prefix}/${fit_prefix}_DPGMM_fit_contours.pdf
    touch ${fit_prefix}/${fit_prefix}_graph.gt
    touch ${fit_prefix}/${fit_prefix}_clusters.csv
    touch ${fit_prefix}/${fit_prefix}_unword_clusters.csv
    touch ${fit_prefix}/${fit_prefix}.refs
    touch ${fit_prefix}/${fit_prefix}.refs_graph.gt
    touch ${fit_prefix}/${fit_prefix}.refs.dists.npy
    touch ${fit_prefix}/${fit_prefix}.refs.dists.pkl
    touch ${fit_prefix}/${fit_prefix}.refs.h5
    """
}
