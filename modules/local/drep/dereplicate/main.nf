// Local copy of nf-core drep/dereplicate with a `genome_info` input: the nf-core module cannot pass
// dRep's --genomeInfo (completeness/contamination), so dRep would otherwise run CheckM or need
// --ignoreGenomeQuality, which drops completeness/contamination from the winner score.
process DREP_DEREPLICATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/drep:3.6.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/drep:3.6.2--pyhdfd78af_0'}"

    input:
    // genome_info = dRep --genomeInfo CSV (genome,completeness,contamination), `genome` being each
    //               FASTA's file name; it must list every input genome
    tuple val(meta), path(fastas, stageAs: 'input_fastas/*'), path(genome_info)

    output:
    tuple val(meta), path("dereplicated_genomes/*"), emit: fastas
    tuple val(meta), path("data_tables/*.csv")     , emit: summary_tables
    tuple val(meta), path("figures/*pdf")          , emit: figures, optional: true
    tuple val(meta), path("logger.log")            , emit: log
    tuple val("${task.process}"), val("drep"), eval("dRep | sed '2!d;s/.*v//g;s/ .*//g'"), emit: versions_drep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    find -L input_fastas/ -type f > fastas_paths.txt

    dRep \\
        dereplicate \\
        drep_work/ \\
        -p ${task.cpus} \\
        -g fastas_paths.txt \\
        --genomeInfo ${genome_info} \\
        ${args}

    mkdir dereplicated_genomes/ figures/ data_tables
    cp drep_work/dereplicated_genomes/* dereplicated_genomes/
    cp drep_work/figures/* figures/ 2>/dev/null || true
    cp drep_work/data_tables/* data_tables/
    cp drep_work/log/logger.log logger.log
    """

    stub:
    // One winner (the first FASTA by name), and genomeInformation.csv rows for every genome with
    // N50 rising by name, so stub runs exercise the reference choice and the N50-ranked target cap.
    """
    mkdir -p dereplicated_genomes/ figures/ data_tables/
    ls input_fastas | sort > names.txt
    first=\$(head -n 1 names.txt)
    cp -L "input_fastas/\$first" dereplicated_genomes/
    echo 'genome,completeness,contamination,length,N50,centrality' > data_tables/genomeInformation.csv
    awk '{ printf "%s,99.0,0.5,2000000,%d,0.99\\n", \$0, 10000 * NR }' names.txt >> data_tables/genomeInformation.csv
    echo 'genome,secondary_cluster,primary_cluster,threshold,cluster_method,comparison_algorithm' > data_tables/Cdb.csv
    awk '{ printf "%s,1_1,1,0.05,average,fastANI\\n", \$0 }' names.txt >> data_tables/Cdb.csv
    echo 'genome,cluster,score' > data_tables/Wdb.csv
    echo "\$first,1_1,100.0" >> data_tables/Wdb.csv
    touch figures/Winning_genomes.pdf
    touch logger.log
    """
}
