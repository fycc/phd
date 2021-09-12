#!/bin/bash

wget -nc ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
wget -nc ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip.md5
md5sum -c new_taxdump.zip.md5 && unzip new_taxdump.zip
grep -wF "33208" taxidlineage.dmp | cut -f1 > animal_taxid.txt
grep -wF "2157" taxidlineage.dmp | cut -f1 > archaea_taxid.txt
grep -wF "2" taxidlineage.dmp | cut -f1 > bacteria_taxid.txt
grep -wF "4751" taxidlineage.dmp | cut -f1 > fungi_taxid.txt
grep -wF "33090" taxidlineage.dmp | cut -f1 > plant_taxid.txt
grep -wF "10239" taxidlineage.dmp | cut -f1 > virus_taxid.txt
sed -e 's/\t|\t/\t/g' -e 's/\t|//g' < names.dmp | awk 'BEGIN {FS="\t"; OFS="\t"} ($4 == "scientific name") {print $1,$2}' > names.tab
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $0}' oncoviruses.txt names.tab > oncoviruses.tab
grep -wFf oncoviruses.txt taxidlineage.dmp | cut -f1 > oncoviruses_alltaxids.txt
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $0}' oncoviruses_alltaxids.txt names.tab > oncoviruses_alltaxids.tab
sed -e 's/\t|\t/\t/g' -e 's/\t|//g' < nodes.dmp > nodes.tab
awk 'BEGIN {FS="\t"} ($3 == "species") {print $1}' nodes.tab > speciesids.txt
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $0}' speciesids.txt oncoviruses_alltaxids.tab > oncoviruses_speciesids.tab
sed -e 's/\t|\t/\t/g' -e 's/\t|//g' < rankedlineage.dmp > rankedlineage.tab
# awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$2]; next} {for(i=1; i<=NF; i++) {if($i in a) {$1=$1; print $1,$2}}}' oncoviruses.tab rankedlineage.tab > oncoviruses_alltaxids.tab
# awk 'BEGIN {FS="\t"} ($10 == "Bacteria") {print $1}' rankedlineage.tab > bacteria_taxid.txt
# awk 'BEGIN {FS="\t"} ($10 == "Viruses") {print $1}' rankedlineage.tab > virus_taxid.txt
