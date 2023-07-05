#!/bin/bash
help(){
	echo ""
	echo "Usage: $0 -g [id_genomes.txt] -m [core/pangenome]  -d [plasmid_fasta] -o [output_folder] -t [Number_of_threads]"
	echo "       Mandatory arguments:"
	echo "       -g [File with genome identifiers. One ID per line']"
  	echo "       -m [Build a core genome or pangenome (Default: core)]"

    echo "       Optional arguments:"
    echo "       -f [Input format for Prokka  (Default: fasta)]"
    echo "       -o [Name of the folder to place the downloaded genomes]"
    echo "       -d [Plasmids Fasta to build local database]"
    echo "       -t [Number of threads]"
	exit
}


while getopts "g:m:d:o:t:h" option;do
	case $option in
		g) list_accesions=$OPTARG;;
        m) mode_use=$OPTARG;;
        d) plasmid_fasta=$OPTARG;;
        o) output=$OPTARG;;
        t) n_threads=$OPTARG;;
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
if [[ -z $mode_use ]];then
	mode_use="core"
fi


if [[ -z $output ]];then
    output="genome_folder"
fi

if [[ -z $n_threads ]]; then
    n_threads=$(cat /proc/cpuinfo | grep processor | wc -l)
    n_threads=$(($n_threads/2)) # Do not use all the power if not indicated
fi


# Check if mode_use is valid
if [[ $mode_use != "core" && $mode_use != "pangenome" ]]; then
    echo "Choose a valid option (core OR pangenome)"
    exit
fi

# Check if the file is prepared to create the database.

echo "Building plasmid database..."
if [[ -z $plasmid_fasta ]]; then
    # Retrieve plasmid sequences
    esearch -db nucleotide -query "Klebsiella pneumoniae[Organism] AND plasmid[Filter] NOT unnamed AND srcdb_refseq[Properties]" | efetch -format fasta > plasmids_file.fasta
            
    # Eliminate duplicates
    seqkit rmdup plasmids_file.fasta > temp.fasta
    mv temp.fasta plasmids_file.fasta
            
    # Include Reference Genome
    ncbi-genome-download bacteria -F fasta -A "GCF_000240185.1" -o referencia
    mv referencia/refseq/bacteria/GCF_000240185.1/*fna* reference_genome.fna.gz
    gzip -d reference_genome.fna.gz
    cat reference_genome.fna >> plasmids_file.fasta
            
    # Eliminate duplicates
    seqkit rmdup plasmids_file.fasta > temp.fasta
    mv temp.fasta plasmids_file.fasta
    makeblastdb -in plasmids_file.fasta -input_type fasta -dbtype nucl -parse_seqids -out klebsiella_db

else
    if [[ ! -f $plasmid_fasta ]] ; then
        echo 'File of plasmids is not there, aborting.'
        exit
    
    else
    mv $plasmid_fasta plasmids_file.fasta
    #Build database
    makeblastdb -in plasmids_file.fasta -input_type fasta -dbtype nucl -parse_seqids -out klebsiella_db
    fi
fi

echo "Done"



##############################################################################

# 1. Download the genomes in fasta and genebank formats.
echo "Download of genomes has started..."
ncbi-genome-download bacteria -F fasta,genbank -A $list_accesions -o $output
echo "Download has finished"

# 2. Unzip the genomes

path_genomes=$(readlink -f $output/refseq/bacteria/GCF*)

for genome in $path_genomes
do
    gzip -d $genome/*.gz
done

# 3A. Removal of plasmids and annotation (Optional)

if [ "$mode_use" == "core" ]
then

    for genome in $path_genomes
    do
        # Circularize assemblies
        fasta_name=$(readlink -f $genome/*fna)
        circlator fixstart $fasta_name circle  # Fix the beginning of the chromosome
        mv circle.fasta $fasta_name; rm circle* # Eliminate additional files
    
        #Sort FASTA by contig size (from bigger to smaller)
        cat $fasta_name | seqkit sort --by-length --reverse > temp.txt
        mv temp.txt $fasta_name

        # Removal of plasmids
        plasmid_path="$genome/plasmid" # Specific folder for plasmids search
        mkdir -p $plasmid_path

        first_line=$(cat $fasta_name | head -n 1)

        # 3A.1 IF GENOME IS CLOSED
        if [[ $first_line == *"complete genome"* ]]; then
            # We calculate the line where the plasmids begin
            limit_line=$(cat -n $fasta_name | grep ">" | awk '{if (NR == 2){print $1}}')
            limit_line=$(($limit_line - 1))

            # Create a new file without the plasmids for each genome        
            head -n $limit_line $fasta_name > temp.txt
            mv temp.txt $fasta_name

            #Annotate the genomes
            anotaciones_path="$genome/anotaciones"
            sudo /home/javi_vcf99/prokka/bin/prokka --species "Klebsiella pneumoniae" --outdir $anotaciones_path --cpus $n_threads $fasta_name

            # Check presence of plasmids
            plasmidfinder.py -i $fasta_name -o $plasmid_path
            
       
        # 3A.2 IF GENOME NOT CLOSED
        else
            
            blast_out="$plasmid_path/plasmid_identification.txt"
            blastn -query $fasta_name -db klebsiella_db -out $blast_out -max_target_seqs 1 -outfmt "6 qseqid sacc stitle"
            cat $blast_out | sort  | uniq > temp.txt
            mv temp.txt $blast_out

            #Check what entries belong to plasmids
            headers_to_remove=$(cat $blast_out | awk -F "\t" '{if (index($3, "plasmid")) print $1}')
            # Remove plasmid entries
            python3 /mnt/c/Users/javi_/Escritorio/MASTER_BIOINF/TFM/scripts/eliminar_plasmidos.py $fasta_name temp.txt $headers_to_remove
            mv temp.txt $fasta_name

            #Annotate the genomes
            anotaciones_path="$genome/anotaciones"
            sudo /home/javi_vcf99/prokka/bin/prokka --species "Klebsiella pneumoniae" --outdir $anotaciones_path  --cpus $n_threads $fasta_name

            # Check presence of plasmids
            plasmidfinder.py -i $fasta_name -o $anotaciones_path
        fi
    
    done
fi


# 3B If wished, keep the plasmids in the FASTA.
if [ "$$mode_use" == "pangenoma" ]
then
    for genome in $path_genomes
    do
        # Circularize assemblies
        fasta_name=$(readlink -f $genome/*fna)
        circlator fixstart $fasta_name circle  # Fix the beginning of the chromosome
        mv circle.fasta $fasta_name; rm circle* # Eliminate additional files
    
        # Sort FASTA by contig size (from bigger to smaller)
        cat $fasta_name | seqkit sort --by-length --reverse > temp.txt
        mv temp.txt $fasta_name
        
        # Annotate the genome
        anotaciones_path="$genome/anotaciones"

        # Annotate the genome by using specificic hmm3 databases (HAMAP-Pfam)
        sudo /home/javi_vcf99/prokka/bin/prokka --kingdom "Bacteria" --species "Klebsiella pneumoniae" --outdir $anotaciones_path $fasta_name --cpus $n_threads
        
        # Idenitfy the Plasmids
        plasmidfinder.py -i $fasta_name -o $anotaciones_path
    done
fi

# 4. Building the core genome/pangenome

# Get the gff_files and write the name in a file
echo -n "" > gff_files.txt
for genome in $path_genomes
do
    name_genome=$(echo $genome |awk -F "/" '{print $NF}')
	cp $genome/anotaciones/*gff $name_genome.gff
    echo "$name_genome.gff" >> gff_files.txt
done

echo "Building the pangenome..."

# If plasmids were removed, build core genome
if [ $mode_use == "core" ]; then
    panaroo --input "gff_files.txt" -t $n_threads --clean-mode "sensitive" --aligner mafft -a core --core_threshold 0.9 -o "panaroo_output" --merge_paralogs
fi

# If plasmids were not removed, build pangenome
if [ $mode_use == "pangenome" ]; then
    panaroo --input "gff_files.txt" -t $n_threads --clean-mode "sensitive" --aligner mafft -a pan --core_threshold 0.9 -o "panaroo_output" --merge_paralogs
fi

# Eliminate copies of GFF files
rm *.gff

# 5. Extract SNP alignment from core genome alignment
cp panaroo_output/core_gene_alignment.aln .
echo "Searching for SNPs in the core/pan alignment..."
snp-sites core_gene_alignment.aln -o snp_alignment.aln

# 6. Build phylogeny from SNP alignment
fconstsvar=$(snp-sites -C core_gene_alignment.aln)
iqtree -s snp_alignment.aln -m GTR+F+I+G4 -nt 6 -bb 1000 -fconst $fconstsvar -redo

# 7. End of the programme
echo "PROGRAM FINISHED"

