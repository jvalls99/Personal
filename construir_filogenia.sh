#!/bin/bash

core_genome="$1"
snp_alignment="$2"
# Porcentaje de cada base en el alineamiento
fconstsvar=$(snp-sites -C $core_genome)
iqtree2 -s $snp_alignment -m GTR+F+I+G4 -nt 6 -bb 1000 -fconst $fconstsvar -redo
