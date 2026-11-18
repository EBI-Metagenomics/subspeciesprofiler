process ANI_QC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same Seqera Containers image as syntracker/select: bin/ani_qc.py imports the evaluator (pandas).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:b51112881a2c3d71' :
        'community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:9b8cabfc8e8a2ced' }"

    input:
    // ani       = all-vs-all FastANI of the species (percent ANI)
    // rfile, labels, genome_qc = the SPECIESQC outputs of the species
    tuple val(meta), path(ani), path(rfile), path(labels), path(genome_qc)

    output:
    tuple val(meta), path("*_ani_qc_rfile.tsv")      , emit: rfile
    tuple val(meta), path("*_ani_qc_labels.csv")     , emit: labels
    tuple val(meta), path("*_ani_qc_genome_qc.tsv")  , emit: genome_qc
    tuple val(meta), path("*.qc.ani.txt")            , emit: ani
    tuple val(meta), path("*_ani_qc.tsv")            , emit: discarded
    tuple val(meta), path("*_ani_qc_eligibility.tsv"), emit: report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args: `--min-ani 0.90`
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ani_qc.py \\
        --species-name "${meta.id}" \\
        --ani ${ani} \\
        --rfile ${rfile} \\
        --labels ${labels} \\
        --genome-qc ${genome_qc} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    // Keep the species values from SPECIESQC, so that stub runs pass the second --qc_filter gate.
    def spp_label = meta.spp_label ?: 'STRONG'
    def n_genomes = meta.n_genomes ?: 120
    """
    cp ${rfile} ${prefix}_ani_qc_rfile.tsv
    cp ${labels} ${prefix}_ani_qc_labels.csv
    cp ${genome_qc} ${prefix}_ani_qc_genome_qc.tsv
    cp ${ani} ${prefix}.qc.ani.txt
    printf 'order\\tgenome\\tlabel\\tN50\\tn_bad_pairs\\tbad_pairs\\n' > ${prefix}_ani_qc.tsv
    printf 'species_name\\tspp_label\\tn_hq\\tn_mq\\tn_eff\\thq_ratio\\tn_total_passing_qc\\tn_discarded\\tn_discarded_contiguity\\tn_discarded_ani\\n' > ${prefix}_ani_qc_eligibility.tsv
    printf '${meta.id}\\t${spp_label}\\t80\\t40\\t100.0\\t0.667\\t${n_genomes}\\t5\\t2\\t0\\n' >> ${prefix}_ani_qc_eligibility.tsv
    """
}
