#!/bin/bash
help(){
	echo ""
	echo "Usage: $0 -d [yes/no] -a [Fasta file] -o [Output name]"
	echo "       Mandatory arguments:"
    echo "       -d: Build personal database (default = no)"
    echo "       -a: Uniprot file"
    echo "       -o:  Output panaroo"
	
	exit
}

while getopts "d:o:a:h" option; do
	case $option in
		d) database=$OPTARG;;
        o) panaroo_folder=$OPTARG;;
        a) uniprot_file=$OPTARG;;
	esac
done

#################################################
####DEFAULT VALUES IF ARGUMENT IS NOT DEFINED####
#################################################
if [[ -z $database ]];then
	database="no"
fi

mkdir -p annotation_databases
cp $uniprot_file annotation_databases/
uniprot_file=$(echo $uniprot_file | awk -F "/" '{print $NF}')
uniprot_file="annotation_databases/$uniprot_file"

# Create my own BLAST database
if [ $database == "yes" ]
then
    makeblastdb -in $uniprot_file -dbtype prot -parse_seqids -input_type fasta -out proteome_blastdb
    mv proteome_blastdb* annotation_databases
fi

# Identify unknown gene_names by performing a BLASTx search
alignments_path="$panaroo_folder/aligned_gene_sequences"
blast_db="annotation_databases/proteome_blastdb"
results_blast="annotation_databases/results_blast.txt"

# Create a folder for the alignemnts with no gaps
nogap_path="$panaroo_folder/no_gaps_aligned_gene_sequences"
mkdir $nogap_path
cp $alignments_path/* $nogap_path

# Correspondencia
corr_file="annotation_databases/corr.txt"
echo -n "" > $corr_file

# Filter
files=$(ls $nogap_path)

original_files=$(ls $panaroo_folder/aligned_gene_sequences | wc -l)
for group_name in $files
do
    fich_mod1=$(ls $nogap_path | wc -l)

    group_file="$nogap_path/$group_name"
    # Eliminate gaps and null lines
    nogap_file="$nogap_path/ng_$group_name"
    cat $group_file | tr -d "-" | tr "\n" "%" | sed 's/%\+/%/g' | sed 's/%/\n/g' >  $nogap_file
    rm $group_file

    # Only take the first genome for those with group name
    STRING="group"
    if [[ "$nogap_file" == *"$STRING"* ]]; then
        limit_line=$(cat -n $nogap_file | grep ">" | awk '{if (NR == 2){print $1}}')
        limit_line=$(($limit_line - 1))
        head -n $limit_line $nogap_file > "annotation_databases/temporary_file.txt"
        blastx -query "annotation_databases/temporary_file.txt" -db $blast_db -out $results_blast -strand "both"
        gene_name=$(grep -F -v '...'  $results_blast | grep GN= | head -n 1 | awk -F "GN=" '{print $2}' | awk -F " " '{print $1}')
        mv $nogap_file "$nogap_path/ng_$gene_name.aln.fas"
        echo "$group_name---$gene_name" >> $corr_file
    fi
done

# Renombrar
groups=$(ls $panaroo_folder/no_gaps_aligned_gene_sequences | grep group)

for group_file in $groups
do
    ori_name=$(echo $group_file | sed 's/ng_//g')
    gen_name=$(cat $corr_file | grep $ori_name | awk -F "---" '{print $2}')
    mv "$panaroo_folder/no_gaps_aligned_gene_sequences/$group_file" "$panaroo_folder/no_gaps_aligned_gene_sequences/ng_$gen_name.aln.fas"

done

