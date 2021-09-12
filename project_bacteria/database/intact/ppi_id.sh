#!/bin/bash
#SBATCH --account=def-yxia
#SBATCH --mem=12000M
#SBATCH --job-name=ppi_id
#SBATCH --output=%x-%j.out

cd ~/database/intact/

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $1,$8}' ppi_id.txt <(zcat ~/database/uniprot/idmapping_selected.tab.gz) > ppi_uniref100.tab
sed -i 's/UniRef100_//g' ppi_uniref100.tab 
cat <(cut -f 2 ppi_uniref100.tab) ppi_id.txt | sort -u > ppi_id_all.txt
mv ppi_id_all.txt ppi_id.txt && rm -f ppi_id_all.txt
