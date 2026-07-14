process POPPUNK_FITMODEL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"

    input:
    // fit_args = the model type + its params for THIS task, e.g. "threshold --threshold 0.0012",
    // "bgmm --K 4", "lineage --ranks 1,2,3", "dbscan --D 5 --min-cluster-prop 0.01", "refine",
    // "refine --multi-boundary 20", "refine --unconstrained". Supplied per-task so one module covers
    // the whole model sweep; task.ext.args adds any common flags.
    // model_dir = optional starting-model directory for `refine` (a prior fit's output). Pass `[]`
    // for the from-scratch families (threshold/lineage/bgmm/dbscan). Staged as `seed_model` so it
    // never collides with this task's own `<prefix>_fitmodel` output directory.
    tuple val(meta), path(db), val(fit_args), path(model_dir, stageAs: 'seed_model')

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
    """
    poppunk \\
        --fit-model ${fit_args} \\
        --ref-db ${db} \\
        ${seed} \\
        --output ${fit_prefix} \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    fit_prefix = "${prefix}_fitmodel"
    """
    echo "${fit_args} ${args}"

    mkdir -p ${fit_prefix}
    touch ${fit_prefix}/${fit_prefix}_fit.npz
    touch ${fit_prefix}/${fit_prefix}_fit.pkl
    touch ${fit_prefix}/${fit_prefix}_DPGMM_fit.png
    touch ${fit_prefix}/${fit_prefix}_DPGMM_fit_contours.pdf
    touch ${fit_prefix}/${fit_prefix}_graph.gt
    touch ${fit_prefix}/${fit_prefix}_clusters.csv
    touch ${fit_prefix}/${fit_prefix}_unword_clusters.csv
    # lineage fits also write a per-rank table (harmless placeholder for other model families)
    printf 'id,Rank_1_Lineage,Rank_2_Lineage,Rank_3_Lineage,overall_Lineage\\ngenome_a,1,1,1,1-1-1\\ngenome_b,2,1,1,2-1-1\\ngenome_c,3,2,1,3-2-1\\n' > ${fit_prefix}/${fit_prefix}_lineages.csv
    touch ${fit_prefix}/${fit_prefix}.refs
    touch ${fit_prefix}/${fit_prefix}.refs_graph.gt
    touch ${fit_prefix}/${fit_prefix}.refs.dists.npy
    touch ${fit_prefix}/${fit_prefix}.refs.dists.pkl
    touch ${fit_prefix}/${fit_prefix}.refs.h5
    """
}
