#/bin/bash
# Get the gff_files 
output="$1"
path_genomes="$output/refseq/bacteria"
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
