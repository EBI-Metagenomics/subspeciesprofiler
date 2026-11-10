# Pangenome test fixture

The five FASTA files are capsular polysaccharide (CPS) loci of *Streptococcus pneumoniae* (EMBL
CR931658, CR931659, CR931660, CR931717 and CR931719, about 25 kb each). They are copied, gzipped and
renamed to their accession, from the `test/` folder of
[ggCaller](https://github.com/bacpop/ggCaller) v1.5.2 (MIT licence, Copyright (c) 2021 Sam Horsfield).

`rfile.tsv` and `labels.csv` follow the SPECIESQC output contract: sample = file name minus its
FASTA extension, path `./genomes/<file>`.
