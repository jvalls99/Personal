#!/bin/bash

# Create my own BLAST database
#makeblastdb -in uniprot_reviewed_proteins.fasta -dbtype prot -parse_seqids -input_type fasta -out proteome_blastdb

# Identify unknown gene_names by performing a BLASTx search
alignments_path="panaroo_output/aligned_gene_sequences"
blast_db="annotation_databases/proteome_blastdb"
results_blast="annotation_databases/results_blast.txt"

unnamed_files=$(ls $alignments_path | grep group)  
for group_name in $unnamed_files
do
    group_file="$alignments_path/$group_name"
    
    # Only take the first genome
    limit_line=$(cat -n $group_file | grep ">" | awk '{if (NR == 2){print $1}}')
    limit_line=$(($limit_line - 1))
    head -n $limit_line $group_file > "annotation_databases/temporary_file.txt"
    blastx -query "annotation_databases/temporary_file.txt" -db $blast_db -out $results_blast -strand "both"
    #echo -n "$group_name "
    gene_name=$(grep -F -v '...'  $results_blast | grep GN= | head -n 1 | awk -F "GN=" '{print $2}' | awk -F " " '{print $1}')
    mv $group_file "$alignments_path/$gene_name.aln.fas"

done

#rm "annotation_databases/temporary_file.txt"
