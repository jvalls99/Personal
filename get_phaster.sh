#!/bin/bash
fasta_genome="$1"
output="$2"

genome=$(echo $fasta_genome | awk -F "/" '{print $(NF -1)}')
ori_dir="$output/refseq/bacteria/$genome"

fasta_ids=$(cat $fasta_genome | grep  ">" | awk -F " " '{print $1}' | tr -d ">" )
echo $fasta_ids

for fasta_id in $fasta_ids
do
    wget "http://phaster.ca/phaster_api?acc=$fasta_id" -O output.json
    link=$(cat output.json | awk -F "," '{print $4}' | tr -d '"' | awk -F ":" '{print $2}')
    echo $link
    wget $link
done

mv *.zip $ori_dir

