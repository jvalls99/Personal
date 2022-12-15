#!/bin/bash
help(){
	echo ""
	echo "Usage: $0 -g [id_genomes.txt] -p [yes/no] -f [fasta/genbank] -a [uniprot_file] -o [output_folder]"
	echo "       Mandatory arguments:"
	echo "       -g [File with genome identifiers. One ID per line']"
    echo "       -a [Uniprot file]"
    echo "       Optional arguments:"
	echo "       -p [Plasmids contained in the genome sequence (Default: no)]"
    echo "       -f [Input format for Prokka  (Default: fasta)]"
    echo "       -o [Name of the folder to place the downloaded genomes]"
	
	exit
}


while getopts "g:p:f:a:o:h" option;do
	case $option in
		g) list_accesions=$OPTARG;;
        p) plasmidos=$OPTARG;;
		f) format=$OPTARG;;
        a) uniprot_file=$OPTARG;;
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
        prokka --species "Klebsiella pneumoniae"  --outdir $anotaciones_path $newfile_path --cpus 8
    done
fi

if [ "$plasmidos" == "yes" ]
then
    for genome in $genome_folders
    do
        # Fasta with the cromosome and the plasmids
        file=$(ls $path_genomes/$genome | grep .fna)
        file_path="$path_genomes/$genome/$file"
        anotaciones_path="$path_genomes/$genome/anotaciones"

        # Annotate the genome
        prokka --kingdom "Bacteria" --species "Klebsiella pneumoniae" --outdir $anotaciones_path $file_path --cpus 8
    
        # Identify Phage Sequences
        get_phaster.sh $file_path $output

        # Idenitfy the Plasmids
        plasmidfinder.py -i $file_path -o $anotaciones_path
    done
fi


# Investigate the pangenome

# Get the gff_files 
genome_folders=$(ls $path_genomes)
# Folder for gff_files
mkdir gff_folder
echo -n "" > "gff_files.txt"

for genome in $genome_folders
do
    gff_file=$(ls $path_genomes/$genome/anotaciones | grep gff )
    cp "$path_genomes/$genome/anotaciones/$gff_file" "$genome.gff" 
    echo "$genome".gff >> "gff_files.txt"
done

#Execute panaroo

echo "Building the pangenome..."

if [ $plasmidos == "no" ]; then
    a=1
    panaroo --input "gff_files.txt" -t 4 --clean-mode "sensitive" --aligner mafft -a core --core_threshold 0.95 -o "panaroo_output"
fi

if [ $plasmidos == "yes" ]; then
    panaroo --input "gff_files.txt" -t 4 --clean-mode "sensitive" --aligner mafft -a pan --core_threshold 0.95 -o "panaroo_output"
fi


mv *.gff gff_folder ; mv gff_files.txt gff_folder
# Encontrar los genes de los ficheros group
echo "Retrieving gene names from unidentified genes..."
ejecutar_blast.sh -d "yes" -o "panaroo_output" -a $uniprot_file

# Deberiamos crear un archivo de SNPss a partir del alineamiento del core (SNPsites)
cp panaroo_output/core_gene_alignment.aln .
echo "Searching for SNPs in the core alignment..."
snp-sites core_gene_alignment.aln -o snp_alignment.aln
echo "PROGRAM FINISHED SUCCESSFULLY"


