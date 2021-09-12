#!/bin/bash
#SBATCH --account=def-yxia
#SBATCH --time=3:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8000M
#SBATCH --job-name=ppi2pfam
#SBATCH --output=%x-%j.out

cd ~/database/uniprot/
grep '>' ~/database/intact/ppi.fasta | sed 's/>//g' | sort -u > ppi_uniprot
zcat idmapping.dat.gz | awk '$2 == "UniParc" {print $1,$3}' > uniprot2uniparc
awk 'BEGIN {OFS=","} NR==FNR {a[$1]; next} $1 in a {$1=$1; print $2,$1}' ppi_uniprot uniprot2uniparc > ppiwithuniparc.csv
cut -d "," -f 1 ppiwithuniparc.csv | sort -u > ppi_uniparc
awk 'BEGIN {FS=","; OFS=","} NR==FNR {a[$1]; next} $1 in a {$1=$1; print $1,$2,$4,$5,$6}' ppi_uniparc uniparc2pfam.csv > ppi_uniparc2pfam.csv
join -t , -j 1 -o 1.2,2.2,2.3,2.4,2.5 <(sort ppiwithuniparc.csv) <(sort ppi_uniparc2pfam.csv) > ppi2uniparc2pfam.csv
awk 'NR==FNR {a[$1]; next} !($1 in a)' uniprot2uniparc ppi_uniprot | xargs samtools faidx ~/database/intact/ppi.fasta > ppiwithoutuniparc.fa
interproscan.sh -appl Pfam -f TSV -i ppiwithoutuniparc.fa 
awk 'BEGIN {FS="\t"; OFS=","} {print $1,$3,$5,$7,$8}' ppiwithoutuniparc.fa.tsv > ppi2interpro2pfam.csv
echo 'seqid,length,pfamacc,pfamstart,pfamend' > ~/database/intact/ppi_pfam.csv
cat ppi2uniparc2pfam.csv ppi2interpro2pfam.csv >> ~/database/intact/ppi_pfam.csv
