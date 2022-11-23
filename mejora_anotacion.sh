#!/bin/bash
help(){
	echo ""
	echo "Usage: $0 -g [id_genomes.txt] -p [yes/no] -f [fasta/genbank] -o [output_folder]"
	echo "       Mandatory arguments:"
	echo "       -g [File with genome identifiers. One ID per line']"
    echo "       Optional arguments:"
	echo "       -p [Plasmids contained in the genome sequence (Default: no)]"
    echo "       -f [Input format for Prokka  (Default: fasta)]"
    echo "       -o [Name of the folder to place the downloaded genomes]"
	
	exit
}


while getopts "g:p:f:o:h" option;do
	case $option in
		g) list_accesions=$OPTARG;;
        p) plasmidos=$OPTARG;;
		f) format=$OPTARG;;
        o) output=$OPTARG;;
        h) help;;
	esac
done

###########################
####MANDATORY ARGUMENTS####
###########################
if [[ -z $list_accesions ]] ;then
	help
fi

#################################################
####DEFAULT VALUES IF ARGUMENT IS NOT DEFINED####
#################################################
if [[ -z $plasmidos ]];then
	plasmidos="no"
fi

if [[ -z $format ]];then
	format="fasta"
fi

if [[ -z $output ]];then
    output="genome_folder"
fi


# Download the genomes in fasta and genebank formats.
echo "Download of genomes has started..."
ncbi-genome-download bacteria -F fasta,genbank -A $list_accesions -o $output
echo "Download has finished"


path_genomes="$output/refseq/bacteria"
genome_folders=$(ls $path_genomes)

# Unzipping genomes
for genome in $genome_folders
do
    gzip -d $path_genomes/$genome/*.gz
done

# If wished, plasmids will be eliminated from the FASTA files.


if [ "$plasmidos" == "no" ]
then
    for genome in $genome_folders
    do        
        file=$(ls $path_genomes/$genome | grep .fna)
        file_path="$path_genomes/$genome/$file"
        new_file="np_"$file
        newfile_path="$path_genomes/$genome/$new_file"
        anotaciones_path="$path_genomes/$genome/anotaciones"

        # We calculate the line where the plasmids begin to cut it off
        limit_line=$(cat -n $file_path | grep ">" | awk '{if (NR == 2){print $1}}')
        limit_line=$(($limit_line - 1))

        # Create a new file without the plasmids for each genome        
        head -n $limit_line $file_path > $newfile_path

        # Annotate the genome
        prokka --species "Klebsiella pneumoniae" --proteins "annotation_databases/uniprot_reviewed_proteins.fasta" --outdir $anotaciones_path $newfile_path 
    done
fi

if [ "$plasmidos" == "yes" ]
then
    for genome in $genome_folders
    do
        # Fasta with the cromosome and the plasmids
        file=$(ls $path/$genome | grep .fna)
        file_path="$path/$genome/$file"
        anotaciones_path="$path_genomes/$genome/anotaciones"

        # Annotate the genome
        prokka --kingdom "Bacteria" --species "Klebsiella pneumoniae" --proteins "annotation_databases/uniprot_reviewed_proteins.fasta" --outdir $anotaciones_path $file_path
        # Idenitfy the Plasmids
        plasmidfinder.py -i $file_path -o $anotaciones_path
    done
fi


# Investigate the pangenome

# Get the gff_files 
genome_folders=$(ls $path)
# Folder for gff_files
mkdir gff_folder
echo -n "" > "gff_files.txt"

for genome in $genome_folders
do
    gff_file=$(ls $path/$genome/anotaciones | grep gff )
    cp "$path/$genome/anotaciones/$gff_file" "$genome.gff" 
    echo "$genome".gff >> "gff_files.txt"
done

#Execute panaroo
if [ $plasmidos == "no" ]; then
    panaroo --input "gff_files.txt" -t 4 --clean-mode "sensitive" --aligner mafft -a core --core_threshold 0.9 -o "panaroo_output"
fi

if [ $plasmidos == "yes"]; then
    panaroo --input "gff_files.txt" -t 4 --clean-mode "sensitive" --aligner mafft -a pan --core_threshold 0.9 -o "panaroo_output"



mv *.gff gff_folder ; mv gff_files.txt gff_folder
# Encontrar los genes de los ficheros group
ejecutar_blast.sh "panaroo_output"

# Deberiamos crear un archivo de SNPss a partir del alineamiento del core (SNPsites)
cp panaroo_output/core_gene_alignment.aln .
snp-sites core_gene_alignment.aln -o snp_alignment.aln


# PlasmidFinder y Phaster.

