process FASTANI_ALLVSALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.34--hb66fcc3_7' :
        'quay.io/biocontainers/fastani:1.34--hb66fcc3_7' }"

    input:
    // Same genome staging as POPPUNK_CREATEDB: the HQ/MQ genomes go into `genomes/`
    // and the r-file's 2nd column gives their `./genomes/<file>` paths. Unlike the
    // nf-core FASTANI module (which needs absolute paths and stages no genomes), this
    // keeps the all-vs-all portable.
    tuple val(meta), path(genomes, stageAs: 'genomes/*'), path(rfile)

    output:
    tuple val(meta), path("*.ani.txt"), emit: ani
    tuple val("${task.process}"), val('fastani'), eval('fastANI --version 2>&1 | head -1 | sed "s/version //"'), topic: versions, emit: versions_fastani

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Build the FastANI genome list from the r-file paths, decompressing gzipped
    # assemblies (fastANI reads uncompressed FASTA only).
    mkdir -p unzipped
    : > genome_paths.txt
    cut -f2 ${rfile} | while IFS= read -r p; do
        [ -z "\$p" ] && continue
        case "\$p" in
            *.gz) out="unzipped/\$(basename "\${p%.gz}")"; gunzip -c "\$p" > "\$out"; echo "\$out" >> genome_paths.txt ;;
            *)    echo "\$p" >> genome_paths.txt ;;
        esac
    done

    fastANI \\
        --ql genome_paths.txt \\
        --rl genome_paths.txt \\
        --threads $task.cpus \\
        -o ${prefix}.ani.txt \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ani.txt
    """
}
