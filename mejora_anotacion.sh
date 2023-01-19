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
    # Create a Plasmid Database
    
    #esearch -db nucleotide -query "Klebsiella pneumoniae[Primary Organism] AND (refseq[filter] AND plasmid[filter]) NOT unnamed*[All Fields]" | efetch -format fasta > klebsiella_refseq_plasmids.fasta
    #cat reference_genome.fna >> klebsiella_refseq_plasmids.fasta
    #makeblastdb -in klebsiella_refseq_plasmids.fasta -input_type fasta -dbtype nucl -parse_seqids -out klebsiella_db


    # For each genome
    for genome in $genome_folders
    do  
        
        file=$(ls $path_genomes/$genome | grep .fna) # Name of the fasta file
        file_path="$path_genomes/$genome/$file"  # Path of the fasta file
        plasmid_path="$path_genomes/$genome/plasmid" # Annotation file's path
        mkdir $plasmid_path

        new_file="np_"$file   # Name of the new fasta file, once plasmids are removed
        newfile_path="$path_genomes/$genome/$new_file" # Path of the new fasta file
        

        first_line=$(cat $file_path | head -n 1)

        # If the genome is complete (Chromosome = 1 entry).
        if [[ $first_line == *"complete genome"* ]]; then
            # We calculate the line where the plasmids begin
            limit_line=$(cat -n $file_path | grep ">" | awk '{if (NR == 2){print $1}}')
            limit_line=$(($limit_line - 1))

            # Create a new file without the plasmids for each genome        
            head -n $limit_line $file_path > $newfile_path
            

        else
            blast_out="$plasmid_path/plasmid_identification.txt"
            blastn -query $file_path -db klebsiella_db -out $blast_out -max_target_seqs 1 -outfmt "6 qseqid sacc stitle"
            cat $blast_out | sort  | uniq > temp.txt
            mv temp.txt $blast_out

            #Check what entries belong to plasmids
            headers_to_remove=$(cat $blast_out | awk -F "\t" '{if (index($3, "plasmid")) print $1}')
            temp_file="temp_$file"     # Name of the temporary fasta file
            temp_file_path="$path_genomes/$genome/$temp_file"  # Path of the temporary fasta file
            # Remove plasmid entries
            python3 eliminar_plasmidos.py $file_path $temp_file_path $headers_to_remove
            mv $temp_file_path $newfile_path
            
        fi
    done

    # Annotate the genomes
    for genome in $genome_folders
    do
        no_plasmid_file=$(ls $path_genomes/$genome| grep np_)
        no_plasmid_path="$path_genomes/$genome/$no_plasmid_file"
        anotaciones_path="$path_genomes/$genome/anotaciones" # Annotation file's path
        prokka --species "Klebsiella pneumoniae"  --outdir $anotaciones_path $no_plasmid_path  --cpus 6
    done
fi


# If wished, keep the plasmids in the FASTA.
if [ "$plasmidos" == "yes" ]
then
    for genome in $genome_folders
    do
        # Fasta with the cromosome and the plasmids
        file=$(ls $path_genomes/$genome | grep .fna)
        file_path="$path_genomes/$genome/$file"
        anotaciones_path="$path_genomes/$genome/anotaciones"

        # Annotate the genome by using specificic hmm3 databases (HAMAP-Pfam)
        prokka --kingdom "Bacteria" --species "Klebsiella pneumoniae" --outdir $anotaciones_path $file_path --cpus 6
    
        # Idenitfy the Plasmids
        plasmidfinder.py -i $file_path -o $anotaciones_path
    done
fi

# Investigate the pangenome

# Get the gff_files 
get_gff.sh $output

#Execute panaroo

echo "Building the pangenome..."

if [ $plasmidos == "no" ]; then
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
echo "Searching for SNPs in the core/pan alignment..."
snp-sites core_gene_alignment.aln -o snp_alignment.aln
echo "PROGRAM FINISHED SUCCESSFULLY"


