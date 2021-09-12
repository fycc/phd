#!/bin/bash
#SBATCH --account=def-yxia
#SBATCH --time=4:00:00
#SBATCH --mem=30G
#SBATCH --job-name=getinfo_host
#SBATCH --output=%x-%j.out

cd ~/database/uniprot

wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_human.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_invertebrates.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_invertebrates.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_mammals.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_mammals.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_rodents.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_rodents.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_vertebrates.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_vertebrates.dat.gz

echo -ne "uniprotid\tgocc" > animal_gocc.tab
zcat uniprot_{sprot,trembl}_{human,invertebrates,mammals,rodents,vertebrates}.dat.gz | awk 'BEGIN {ORS=""}; {if (/^ID/) {print "\n"$2"\t"} else if (/^DR/ && /GO;/ && /GO:[0-9]+; C:/) {print $3}}' >> animal_gocc.tab

echo -ne "uniprotid\tproteome" > animal_proteome.tab
zcat uniprot_{sprot,trembl}_{human,invertebrates,mammals,rodents,vertebrates}.dat.gz | awk 'BEGIN {ORS=""}; {if (/^ID/) {print "\n"$2"\t"} else if (/^DR/ && /Proteomes/) {print $3}}' >> animal_proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} ($2) {print $0}' animal_proteome.tab > animal_proteome_nonempty.tab
mv animal_proteome_nonempty.tab animal_proteome.tab

tail -n +2 animal_proteome.tab | cut -f 1 > animal_uniprotid.txt
echo -e "uniprotid\tuniref100\tuniref90\tuniref50\tuniparc\ttaxid" > animal_uniref_uniparc_taxid.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($2 in a) {print $2,$8,$9,$10,$11,$13}' animal_uniprotid.txt <(zcat idmapping_selected.tab.gz) >> animal_uniref_uniparc_taxid.tab

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$5; next} ($1 in a) {print a[$1],$2}' animal_uniref_uniparc_taxid.tab animal_proteome.tab | tail -n +2 > animal_uniparc_proteome.tab
tail -n +2 animal_uniref_uniparc_taxid.tab | cut -f 5 > animal_uniparc.txt
awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$4}' animal_uniparc.txt uniparc2pfam.csv > animal_uniparc_pfam.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $2; sep[$1]=";"} END {for(k in f2) print k,f2[k]}' animal_uniparc_pfam.tab > animal_uniparc2pfam.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $2; sep[$1]=""} END {for(k in f2) print k,f2[k]}' animal_uniparc_proteome.tab > animal_uniparc2proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$2; next} ($1 in a) {print a[$1],$2}' animal_uniparc2pfam.tab animal_uniparc2proteome.tab > animal_pfam2proteome.tab

wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz
wget -nc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_plants.dat.gz

echo -ne "uniprotid\tgocc" > plant_gocc.tab
zcat uniprot_{sprot,trembl}_plants.dat.gz | awk 'BEGIN {ORS=""}; {if (/^ID/) {print "\n"$2"\t"} else if (/^DR/ && /GO;/ && /GO:[0-9]+; C:/) {print $3}}' >> plant_gocc.tab

echo -ne "uniprotid\tproteome" > plant_proteome.tab
zcat uniprot_{sprot,trembl}_plants.dat.gz | awk 'BEGIN {ORS=""}; {if (/^ID/) {print "\n"$2"\t"} else if (/^DR/ && /Proteomes/) {print $3}}' >> plant_proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} ($2) {print $0}' plant_proteome.tab > plant_proteome_nonempty.tab
mv plant_proteome_nonempty.tab plant_proteome.tab

tail -n +2 plant_proteome.tab | cut -f 1 > plant_uniprotid.txt
echo -e "uniprotid\tuniref100\tuniref90\tuniref50\tuniparc\ttaxid" > plant_uniref_uniparc_taxid.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($2 in a) {print $2,$8,$9,$10,$11,$13}' plant_uniprotid.txt <(zcat idmapping_selected.tab.gz) >> plant_uniref_uniparc_taxid.tab

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$5; next} ($1 in a) {print a[$1],$2}' plant_uniref_uniparc_taxid.tab plant_proteome.tab | tail -n +2 > plant_uniparc_proteome.tab
tail -n +2 plant_uniref_uniparc_taxid.tab | cut -f 5 > plant_uniparc.txt
awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$4}' plant_uniparc.txt uniparc2pfam.csv > plant_uniparc_pfam.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $2; sep[$1]=";"} END {for(k in f2) print k,f2[k]}' plant_uniparc_pfam.tab > plant_uniparc2pfam.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $2; sep[$1]=""} END {for(k in f2) print k,f2[k]}' plant_uniparc_proteome.tab > plant_uniparc2proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$2; next} ($1 in a) {print a[$1],$2}' plant_uniparc2pfam.tab plant_uniparc2proteome.tab > plant_pfam2proteome.tab
