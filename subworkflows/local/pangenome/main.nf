//
// PANGENOME: the gene presence/absence table of one species, over the genomes that passed SPECIESQC.
//
// 1. GGCALLER calls genes in all of the species' genomes at once on a Bifrost graph, so the calls are
//    consistent across fragmented assemblies, and writes one GFF3 per genome (renamed to the r-file
//    sample names). Its BALROG database comes from --ggcaller_db, or is downloaded once per run
//    (GGCALLER_DOWNLOADDB).
// 2. PANAROO_RUN clusters the GFFs into gene families and writes gene_presence_absence.Rtab, one
//    column per sample. Contamination is already filtered by SPECIESQC, so Panaroo runs in its
//    sensitive mode (it deletes no genes) and without refinding (--refind-mode off).
//

include { GGCALLER_DOWNLOADDB } from '../../../modules/local/ggcaller/downloaddb/main'
include { GGCALLER            } from '../../../modules/local/ggcaller/run/main'
include { PANAROO_RUN         } from '../../../modules/nf-core/panaroo/run/main'

workflow PANGENOME {

    take:
    ch_input // channel: [ val(meta), [ genome files ], path(rfile), path(labels) ]  (as POPPUNK_METHODS)

    main:

    // Only the r-file genomes: the DISCARDED ones are not staged, so they cannot invalidate -resume.
    ch_genomes = ch_input
        .map { meta, genomes, rfile, labels ->
            def qc_files = rfile.splitCsv(sep: '\t').collect { row -> file(row[1]).name } as Set
            [ meta, [ genomes ].flatten().findAll { genome -> genome.name in qc_files }, rfile ]
        }

    if ( params.ggcaller_db ) {
        ch_db = Channel.value( file( params.ggcaller_db, checkIfExists: true ) )
    } else {
        GGCALLER_DOWNLOADDB()
        ch_db = GGCALLER_DOWNLOADDB.out.db.collect()
    }

    GGCALLER( ch_genomes, ch_db )
    PANAROO_RUN( GGCALLER.out.gff )

    emit:
    rtab     = PANAROO_RUN.out.results.map { meta, files -> [ meta, files.find { it.name == 'gene_presence_absence.Rtab' } ] } // channel: [ val(meta), path(gene_presence_absence.Rtab) ]
    panaroo  = PANAROO_RUN.out.results // channel: [ val(meta), [ Panaroo result files ] ]
    gff      = GGCALLER.out.gff        // channel: [ val(meta), [ <sample>.gff ] ]
    proteins = GGCALLER.out.proteins   // channel: [ val(meta), path(*_gene_calls.faa) ]
    genes    = GGCALLER.out.genes      // channel: [ val(meta), path(*_gene_calls.ffn) ]
}
