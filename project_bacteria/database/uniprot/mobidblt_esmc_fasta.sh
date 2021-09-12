#!/bin/bash
#SBATCH --account=def-yxia
#SBATCH --mem=8000M
#SBATCH --job-name=mobidblt_esmc_fasta
#SBATCH --output=%x-%j.out

cd ~/database/uniprot

cut -f 1 mobidblt_esmc.tab | sort -u > mobidblt_esmc_uniprotid.txt

samtools faidx esmc.fasta 
cut -d '|' -f 3 esmc.fasta.fai > esmc.fasta.fai.temp
mv esmc.fasta.fai.temp esmc.fasta.fai

cat mobidblt_esmc_uniprotid.txt | xargs samtools faidx esmc.fasta > mobidblt_esmc.fasta
