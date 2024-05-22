This folder contains the programs necessary for the construction of a core genome and a pangenome.

TFM:

- construction_core.sh: Receives a list of NCBI identifiers (GCF) as argument. The user chooses by parameters whether to build the core or pangenome, whether to build a plasmid database and in that case, can enter the FASTA file with the sequences. The processes of downloading, annotation
downloading, annotation, plasmid removal (only for core) and obtaining alignments with Panaroo are automated. SNPs and phylogeny are also obtained.
core.

  Input:
        - accessions_list: One NCBI ID per line.
  Output:
        - Folder with downloaded genomes and their plasmid annotation and identification files-.
        - Folder with the individual gene alignments and their concatenation (core/pangenome).
        - Phylogenetic tree with the SNPs of the alignment.
  
            
DIFFERENTIAL EXPRESSION:

- pipeline.nf: It is written with Nextflow language to facilitate pipeline construction. It allows to obtain from SRA identifiers and a reference the BAM files, from which the count matrix will be built.

OTHERS:

- search_pattern_virus.py: Some code to do some exercises about pattern searching and basic work in Python.

