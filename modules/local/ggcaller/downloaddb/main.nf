process GGCALLER_DOWNLOADDB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ggcaller%3A1.5.2--hd2be7a0_0' :
        'biocontainers/ggcaller:1.5.2--hd2be7a0_0' }"

    output:
    path "ggcaller_db"  , emit: db
    tuple val("${task.process}"), val('ggcaller'), eval("ggcaller --version | sed 's/^ggcaller //'"), topic: versions, emit: versions_ggcaller

    when:
    task.ext.when == null || task.ext.when

    script:
    // ggCaller's download-only branch (--balrog-db without inputs) fetches ggCallerdb.tar.bz2
    // (~120 MB, ftp.ebi.ac.uk) into the directory, which must already exist (otherwise the
    // tarball is saved as a file of that name), and unpacks it. ggCaller unpacks the BALROG
    // models only on first use, which would write into a read-only input; unpack them here so
    // GGCALLER never writes to the database.
    """
    mkdir ggcaller_db
    ggcaller --balrog-db ggcaller_db
    tar -xzf ggcaller_db/ggCallerdb/balrog_models.tar.gz -C ggcaller_db/ggCallerdb
    test -d ggcaller_db/ggCallerdb/balrog_models
    rm -rf ggCaller_output
    """

    stub:
    """
    mkdir -p ggcaller_db/ggCallerdb/balrog_models
    touch ggcaller_db/ggCallerdb.tar.bz2
    """
}
