This folder contains the necessary scripts to build a pangenome. 

1. mejora_anotacion.sh: It receives a list of NCBI identifiers as an argument, where the user is capable of choosing whether the plasmids should be incorporated to the analysis or not. Once the genomes are downloaded, an automated annotation is performed thanks to PROKKA, which uses HMM files to detect the genes. (I used HAMAP and Pfam files). Then, the programme calls panaroo to carry out a multiple sequence alignment with all the genomes by using MAFFT aligned.

2. ejecutar_blast.sh: It creates a BLAST database from a FASTA file containing all of the UNIPROT reviewed protein sequences. Given that when panaroo is performed, some of the genes are unnamed, this is script allows us to identify which funcion/name those genes may have. For each gene_group we compare it to the created database through BLASTX and change the name of the file.

3. 
