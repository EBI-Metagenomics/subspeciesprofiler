process ANI_TREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // biopython 1.84 (Python 3.12); the 1.65 image is Python 2.7 and cannot run the script.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'biocontainers/biopython:1.84' }"

    input:
    // ani = all-vs-all FastANI of the species (percent ANI, five columns)
    tuple val(meta), path(ani)

    output:
    tuple val(meta), path("*_ani_nj_tree.nwk"), emit: tree
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Neighbour joining is O(n^3) in pure Python: ~20 s for 400 genomes, ~30 min for 1,800.
    """
    fastani_to_tree.py ${ani} ${prefix}_ani_nj_tree.nwk
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '(genome1:0.005,genome2:0.005);' > ${prefix}_ani_nj_tree.nwk
    """
}
