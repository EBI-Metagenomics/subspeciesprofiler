process SYNTRACKER_RUN {
    tag "$meta.id"
    label 'process_high'

    // Container-only: SynTracker is not on bioconda, and its R stack must stay pinned to the
    // versions it was written for (R 4.0.5, DECIPHER 2.18.1); a current DECIPHER fails on every
    // region. The image is built from the EBI-Metagenomics containers repo (syntracker/1.4.0_patch1).
    container "quay.io/microbiome-informatics/syntracker:1.4.0_patch1"

    input:
    // reference = the single reference genome (the most central HQ genome, from syntracker/select)
    // targets   = the genomes to compare; SynTracker names each by its file basename
    tuple val(meta), path(reference, stageAs: 'ref_in/*'), path(targets, stageAs: 'targets_in/*')

    output:
    tuple val(meta), path("syntracker/summary_output/avg_synteny_scores_*.csv"), emit: apss
    tuple val(meta), path("syntracker/summary_output/synteny_scores_per_region.csv"), emit: per_region
    tuple val(meta), path("syntracker/SynTracker_log.txt")                       , emit: log
    tuple val("${task.process}"), val('syntracker'), val(VERSION), topic: versions, emit: versions_syntracker

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SYNTRACKER_RUN only supports container profiles (docker/singularity): SynTracker is not packaged for conda."
    }
    // task.ext.args: SynTracker options, e.g. `-length 5000 --identity 97 --coverage 70`.
    def args = task.ext.args ?: ''
    // WARN: SynTracker reports no version on the CLI; keep in sync with the container tag.
    VERSION = '1.4.0'
    """
    # SynTracker takes directories of FASTA files and names each genome by its file basename.
    # Stage every file as <sample>.fasta (decompressing .gz), so the names match the PopPUNK
    # sample names that evaluate_poppunk_fastani.py's normalise_genome_id derives.
    stage_fasta() {
        local src="\$1" dest_dir="\$2" name
        name=\$(basename "\$src"); name=\${name%.gz}; name=\${name%.*}
        case "\$src" in
            *.gz) gzip -cd "\$src" > "\$dest_dir/\$name.fasta" ;;
            *)    ln -s "\$(readlink -f "\$src")" "\$dest_dir/\$name.fasta" ;;
        esac
    }
    mkdir -p ref targets
    for f in ref_in/*;     do stage_fasta "\$f" ref;     done
    for f in targets_in/*; do stage_fasta "\$f" targets; done

    # Always a fresh run (-mode new): Nextflow retries start in a clean work dir, and SynTracker's
    # continue mode preloads every finished region before forking its workers (memory-hungry).
    syntracker.py \\
        -ref ref \\
        -target targets \\
        -out syntracker \\
        -mode new \\
        -cores ${task.cpus} \\
        ${args}

    # SynTracker reports a failed R step only in its log; it still leaves header-only summary
    # tables behind. Fail the task unless the all-regions APSS table holds at least one pair.
    if [ "\$(awk 'END { print NR }' syntracker/summary_output/avg_synteny_scores_all_regions.csv)" -lt 2 ]; then
        echo "ERROR: SynTracker produced no APSS pairs; see syntracker/SynTracker_log.txt" >&2
        exit 1
    fi
    """

    stub:
    VERSION = '1.4.0'
    // Write one APSS row per target pair (all 0.95), so the clustering and evaluation steps
    // downstream run on realistic, correctly named input in stub mode.
    """
    mkdir -p syntracker/summary_output
    ls targets_in | sed -e 's/\\.gz\$//' -e 's/\\.[^.]*\$//' | sort > samples.txt
    header='"Ref_genome","Sample1","Sample2","APSS","Compared_regions"'
    echo "\$header" > syntracker/summary_output/avg_synteny_scores_all_regions.csv
    awk '{ s[NR] = \$0 } END {
            for (i = 1; i < NR; i++) for (j = i + 1; j <= NR; j++)
                printf "\\"stub_ref\\",\\"%s\\",\\"%s\\",0.95,100\\n", s[i], s[j]
         }' samples.txt >> syntracker/summary_output/avg_synteny_scores_all_regions.csv
    for n in 40 60 80 100 200; do
        # per-n tables: every pair compared on exactly n regions, as SynTracker writes them
        awk -F, -v OFS=, -v n=\$n 'NR > 1 { \$NF = n } 1' syntracker/summary_output/avg_synteny_scores_all_regions.csv > syntracker/summary_output/avg_synteny_scores_\${n}_regions.csv
    done
    echo '"Ref_genome","Sample1","Sample2","Region","Length1","Length2","Overlap","Blocks","Synteny_score"' > syntracker/summary_output/synteny_scores_per_region.csv
    touch syntracker/SynTracker_log.txt
    """
}
