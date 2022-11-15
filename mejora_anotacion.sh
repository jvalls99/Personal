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


path="$output/refseq/bacteria"
genome_folders=$(ls $path)

# Unzipping genomes
for genome in $genome_folders
do
    gzip -d $path/$genome/*.gz
done

# If wished, plasmids will be eliminated from the FASTA files.


if [ "$plasmidos" == "no" ]
then
    for genome in $genome_folders
    do        
        file=$(ls $path/$genome | grep .fna)
        new_file="np_"$file

        # We calculate the line where the plasmids begin to cut it off
        limit_line=$(cat -n $path/$genome/$file | grep ">" | awk '{if (NR == 2){print $1}}')
        limit_line=$(($limit_line - 1))

        # Create a new file without the plasmids for each genome        
        head -n $limit_line $path/$genome/$file > $path/$genome/$new_file

        # Annotate the genome
        prokka --species "Klebsiella pneumoniae" --proteins "annotation_databases/uniprot_reviewed_proteins.fasta" --outdir $path/$genome/"anotaciones" $path/$genome/$new_file 
    done
fi

if [ "$plasmidos" == "yes" ]
then
    for genome in $genome_folders
    do
        # Fasta with the cromosome and the plasmids
        file=$(ls $path/$genome | grep .fna)
        # Annotate the genome
        prokka --kingdom "Bacteria" --species "Klebsiella pneumoniae" --proteins "annotation_databases/uniprot_reviewed_proteins.fasta" --outdir $path/$genome/"anotaciones" $path/$genome/$file
    done
fi


# Investigate the pangenome

# Get the gff_files 
path="genomas_descargados/refseq/bacteria"
genome_folders=$(ls $path)
# Folder for gff_files
mkdir gff_folder
gff_path="gff_folder/gff_files"
echo -n "" > $gff_path

for genome in $genome_folders
do
    gff_file=$(ls $path/$genome/anotaciones | grep gff )
    echo "$path/$genome/anotaciones/$gff_file"
    cp "$path/$genome/anotaciones/$gff_file" "$genome.gff"
    mv "$genome.gff" "gff_folder"
    echo "$genome".gff3 >> $gff_path 
done

#Execute panaroo
panaroo --input "gff_folder/gff_files" -t 4 --clean-mode "sensitive" --aligner mafft -a core --core_threshold 0.9 -o "panaroo_results"

# Renombrar groups con un blast (solo usar 1 genoma).


# Deberiamos crear un archivo de SNPss a partir del alineamiento del core (SNPsites)
cp panaroo_results/core_gene_alignment.aln .
snp-sites core_gene_alignment.aln -o snp_alignment.aln