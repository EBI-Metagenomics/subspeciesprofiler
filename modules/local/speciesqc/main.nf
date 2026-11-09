process SPECIESQC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    // qc_csv = per-genome completeness/contamination (genome,completeness,contamination)
    // stats  = seqkit stats --all -N 90 --tabular of the species' genome files (N50 and N90)
    tuple val(meta), path(qc_csv), path(stats)

    output:
    tuple val(meta), path("*_spp_eligibility_report.tsv"), emit: report
    tuple val(meta), path("*_rfile.tsv")                 , emit: rfile
    tuple val(meta), path("*_labels.csv")                , emit: labels
    tuple val(meta), path("*_genome_qc.tsv")             , emit: genome_qc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // task.ext.args sets the contiguity gate and can override the QC thresholds, e.g.:
    //   --contiguity_metric N50 --min_contiguity 5000
    //   --min_completeness 80 --hq_completeness 90 --hq_contamination 1 --max_contamination 5
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    spp_eligibility_from_qc.py \\
        --species_name "${meta.id}" \\
        --qc_csv ${qc_csv} \\
        --stats ${stats} \\
        --output ${prefix}_spp_eligibility_report.tsv \\
        --rfile_output ${prefix}_rfile.tsv \\
        --labels_output ${prefix}_labels.csv \\
        --genome_qc_output ${prefix}_genome_qc.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'species_name\\tspp_label\\tn_hq\\tn_mq\\tn_eff\\thq_ratio\\tn_total_passing_qc\\tn_discarded\\tn_discarded_contiguity\\n' > ${prefix}_spp_eligibility_report.tsv
    printf '${meta.id}\\tSTRONG\\t80\\t40\\t100.0\\t0.667\\t120\\t5\\t2\\n' >> ${prefix}_spp_eligibility_report.tsv
    printf 'genome1\\t./genomes/genome1.fna.gz\\n' > ${prefix}_rfile.tsv
    printf 'genome2\\t./genomes/genome2.fna.gz\\n' >> ${prefix}_rfile.tsv
    # both genomes HQ, so stub runs reach the SynTracker branch (it needs at least 2 HQ genomes)
    printf 'genome,label\\ngenome1,HQ\\ngenome2,HQ\\n' > ${prefix}_labels.csv
    printf 'genome\\tfile\\tcompleteness\\tcontamination\\tN50\\tN90\\tquality\\tlabel\\treason\\n' > ${prefix}_genome_qc.tsv
    printf 'genome1\\tgenome1.fna.gz\\t99.5\\t0.5\\t1000000\\t100000\\tHQ\\tHQ\\tpass\\n' >> ${prefix}_genome_qc.tsv
    printf 'genome2\\tgenome2.fna.gz\\t98.2\\t0.8\\t1000000\\t100000\\tHQ\\tHQ\\tpass\\n' >> ${prefix}_genome_qc.tsv
    """
}
