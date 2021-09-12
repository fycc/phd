#!/bin/bash
#SBATCH --account=def-yxia
#SBATCH --time=8:00:00
#SBATCH --mem=90G
#SBATCH --job-name=getinfo_bacteria
#SBATCH --output=%x-%j.out

cd ~/database/uniprot

wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz

echo -ne "uniprotid\tgocc" > bacteria_gocc.tab
zcat uniprot_{sprot,trembl}_bacteria.dat.gz | awk 'BEGIN {ORS=""}; {if (/^ID/) {print "\n"$2"\t"} else if (/^DR/ && /GO;/ && /GO:[0-9]+; C:/) {print $3}}' >> bacteria_gocc.tab

echo -ne "uniprotid\tproteome" > bacteria_proteome.tab
zcat uniprot_{sprot,trembl}_bacteria.dat.gz | awk 'BEGIN {ORS=""}; {if (/^ID/) {print "\n"$2"\t"} else if (/^DR/ && /Proteomes/) {print $3}}' >> bacteria_proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} ($2) {print $0}' bacteria_proteome.tab > bacteria_proteome_nonempty.tab
mv bacteria_proteome_nonempty.tab bacteria_proteome.tab

tail -n +2 bacteria_proteome.tab | cut -f 1 > bacteria_uniprotid.txt
echo -e "uniprotid\tuniref100\tuniref90\tuniref50\tuniparc\ttaxid" > bacteria_uniref_uniparc_taxid.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($2 in a) {print $2,$8,$9,$10,$11,$13}' bacteria_uniprotid.txt <(zcat idmapping_selected.tab.gz) >> bacteria_uniref_uniparc_taxid.tab

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$5; next} ($1 in a) {print a[$1],$2}' bacteria_uniref_uniparc_taxid.tab bacteria_proteome.tab | tail -n +2 > bacteria_uniparc_proteome.tab
tail -n +2 bacteria_uniref_uniparc_taxid.tab | cut -f 5 > bacteria_uniparc.txt
awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$4}' bacteria_uniparc.txt uniparc2pfam.csv > bacteria_uniparc_pfam.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $2; sep[$1]=";"} END {for(k in f2) print k,f2[k]}' bacteria_uniparc_pfam.tab > bacteria_uniparc2pfam.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $2; sep[$1]=""} END {for(k in f2) print k,f2[k]}' bacteria_uniparc_proteome.tab > bacteria_uniparc2proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$2; next} ($1 in a) {print a[$1],$2}' bacteria_uniparc2pfam.tab bacteria_uniparc2proteome.tab > bacteria_pfam2proteome.tab
