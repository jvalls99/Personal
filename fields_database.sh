#!/bin/bash
db="$1"
header=$(cat $db | head -n 1)
number_fields=$(echo $header | grep -o "," | wc -l )
number_fields=$(($number_fields + 1))

for (( i=1; i<=$number_fields; i++ ))
do
    echo $header | awk -v m=$i -F "," '{print m,$m}'
done
