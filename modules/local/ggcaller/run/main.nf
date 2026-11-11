process GGCALLER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ggcaller%3A1.5.2--hd2be7a0_0' :
        'biocontainers/ggcaller:1.5.2--hd2be7a0_0' }"

    input:
    // rfile   = SPECIESQC r-file (sample<TAB>./genomes/<file>): the genomes to call genes in
    // genomes = those genome files, staged under genomes/ to match the r-file paths
    // db      = ggCaller database directory (from GGCALLER_DOWNLOADDB or --ggcaller_db)
    tuple val(meta), path(genomes, stageAs: 'genomes/*'), path(rfile)
    path db

    output:
    tuple val(meta), path("gff/*.gff.gz")       , emit: gff
    tuple val(meta), path("faa/*.faa.gz")       , emit: proteins
    tuple val(meta), path("ffn/*.ffn.gz")       , emit: genes
    tuple val("${task.process}"), val('ggcaller'), eval("ggcaller --version | sed 's/^ggcaller //'"), topic: versions, emit: versions_ggcaller

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    // Bifrost refuses more threads than the node has cores and ggCaller then exits 0 with no genes,
    // so cap the threads at nproc and check below that every genome got a GFF.
    // ggCaller names each genome by its file name minus the last extension (x.fa.gz -> x.fa), in
    // the GFF file names, the GFF IDs and the FASTA headers; rename them to the r-file sample
    // names, the names the rest of the pipeline uses (normalise_genome_id).
    // ggCaller writes the gene calls of all genomes into one gene_calls.{faa,ffn} and has no option
    // to split them, so they are split here into one file per sample; every output is gzipped.
    """
    for p in ggCallerdb.tar.bz2 ggCallerdb/balrog_models; do
        if [ ! -e "${db}/\$p" ]; then
            echo "ERROR: ggCaller database '${db}' lacks \$p (use the GGCALLER_DOWNLOADDB output, see --ggcaller_db)" >&2
            exit 1
        fi
    done

    threads=\$(( ${task.cpus} < \$(nproc) ? ${task.cpus} : \$(nproc) ))
    cut -f2 ${rfile} > refs.txt

    ggcaller \\
        --refs refs.txt \\
        --balrog-db ${db} \\
        --out ggcaller \\
        --threads \$threads \\
        ${args}

    mkdir gff faa ffn
    : > stems.tsv
    while IFS=\$'\\t' read -r sample path; do
        stem=\$(basename "\$path")
        stem=\${stem%.*}
        if [ ! -s "ggcaller/GFF/\${stem}.gff" ]; then
            echo "ERROR: ggCaller wrote no GFF for '\$sample' (\$path); see the Bifrost messages above" >&2
            exit 1
        fi
        sed "s/ID=\${stem}_/ID=\${sample}_/" "ggcaller/GFF/\${stem}.gff" > "gff/\${sample}.gff"
        printf '%s\\t%s\\n' "\$stem" "\$sample" >> stems.tsv
    done < ${rfile}
    # FASTA headers are >\${stem}_<ORF id> (the id may be negative): one pass per file, renaming each
    # header and writing each record to <ext>/<sample>.<ext>. The records of a genome are contiguous,
    # so only one output file is open at a time.
    for ext in faa ffn; do
        awk -F'\\t' -v ext=\$ext 'NR == FNR { name[\$1] = \$2; next }
            /^>/ { id = substr(\$0, 2); n = match(id, /_-?[0-9]+(\$| )/); stem = substr(id, 1, n - 1)
                   if (!(stem in name)) { print "ERROR: gene call " id " matches no r-file genome" > "/dev/stderr"; exit 1 }
                   if (name[stem] != cur) { if (out != "") close(out); cur = name[stem]; out = ext "/" cur "." ext }
                   \$0 = ">" cur substr(id, n) }
            { print >> out }' stems.tsv ggcaller/gene_calls.\$ext
    done
    while IFS=\$'\\t' read -r sample path; do
        for ext in faa ffn; do
            if [ ! -s "\$ext/\$sample.\$ext" ]; then
                echo "ERROR: ggCaller called no genes for '\$sample' (no \$ext/\$sample.\$ext)" >&2
                exit 1
            fi
        done
    done < ${rfile}
    ls gff/*.gff faa/*.faa ffn/*.ffn | xargs -P \$threads -n 20 gzip -n
    """

    stub:
    // One minimal Panaroo-readable GFF3 (one CDS + the FASTA) and one gene per r-file sample.
    """
    mkdir gff faa ffn
    cut -f1 ${rfile} | while read -r sample; do
        printf '##gff-version 3\\n##sequence-region contig1 1 9\\ncontig1\\tggCaller:1.5.2\\tCDS\\t1\\t9\\t0\\t+\\t0\\tID=%s_1;annotation=hypothetical protein;inference=prediction\\n##FASTA\\n>contig1\\nATGAAATAA\\n' "\$sample" | gzip -n > gff/\$sample.gff.gz
        printf '>%s_1\\nMK\\n' "\$sample" | gzip -n > faa/\$sample.faa.gz
        printf '>%s_1\\nATGAAATAA\\n' "\$sample" | gzip -n > ffn/\$sample.ffn.gz
    done
    """
}
