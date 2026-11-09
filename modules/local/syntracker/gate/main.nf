process SYNTRACKER_GATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same Seqera Containers image as syntracker/select.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:b51112881a2c3d71' :
        'community.wave.seqera.io/library/matplotlib-base_pandas_procps-ng_python-igraph_python:9b8cabfc8e8a2ced' }"

    input:
    // stats  = seqkit stats --all -N 90 --tabular of the HQ/MQ genomes (one row per FASTA, with N50 and N90)
    // rfile  = PopPUNK r-file (sample<TAB>./genomes/<file>); read by the stub only
    // labels = per-genome labels CSV (genome,label)
    tuple val(meta), path(stats), path(rfile), path(labels)

    output:
    tuple val(meta), path("*_syntracker_gate.tsv")      , emit: report
    tuple val(meta), path("*_syntracker_gate_stats.tsv"), emit: stats
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args: `--metric N50 --min-length 5000`
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    syntracker_gate.py \\
        --stats ${stats} \\
        --labels ${labels} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Every r-file genome passes, with its label and fixed N50/N90. The seqkit stub writes an
    // empty table, so the genomes come from the r-file instead.
    """
    awk -F, 'NR > 1 { print \$1 "\\t" \$2 }' ${labels} > labels.tsv
    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tN50\\tN90\\tlabel\\n' > ${prefix}_syntracker_gate_stats.tsv
    awk -F'\\t' 'NR == FNR { label[\$1] = \$2; next } (\$1 in label) { n = split(\$2, p, "/"); print p[n] "\\tFASTA\\tDNA\\t1\\t5000000\\t1000000\\t100000\\t" label[\$1] }' \\
        labels.tsv ${rfile} >> ${prefix}_syntracker_gate_stats.tsv
    printf 'genome\\tfile\\tlabel\\tN50\\tN90\\tmetric\\tmin_length\\tstatus\\n' > ${prefix}_syntracker_gate.tsv
    """
}
