#!/bin/bash
#SBATCH --account=def-yxia
#SBATCH --time=2:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8000M
#SBATCH --job-name=getinfo_esmc
#SBATCH --output=%x-%j.out

cd ~/database/uniprot

# UniProt search keywords
# bacterial effector proteins by name: taxonomy:"Bacteria [2]" name:effector (name:"type 1" OR name:"type 2" OR name:"type 3" OR name:"type 4" OR name:"type 5" OR name:"type 6" OR name:"type 7" OR name:"type 8" OR name:"type 9" OR name:t*ss OR name:"secretion system")
# bacterial effector proteins by cellular location: taxonomy:"Bacteria [2]" (annotation:(type:function effector) OR locations:(note:"type 1") OR locations:(note:"type 2") OR locations:(note:"type 3") OR locations:(note:"type 4") OR locations:(note:"type 5") OR locations:(note:"type 6") OR locations:(note:"type 7") OR locations:(note:"type 8") OR locations:(note:"type 9") OR locations:(note:t*ss) OR locations:(note:"secretion system")) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 1") OR locations:(note:t1ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 2") OR locations:(note:t2ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 3") OR locations:(note:t3ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 4") OR locations:(note:t4ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 5") OR locations:(note:t5ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 6") OR locations:(note:t6ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 7") OR locations:(note:t7ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 8") OR locations:(note:t8ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# taxonomy:"Bacteria [2]" (locations:(note:"type 9") OR locations:(note:t9ss)) (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# bacterial secreted proteins: taxonomy:"Bacteria [2]" (locations:(location:"Secreted [SL-0243]") OR locations:(location:"Host [SL-0431]"))
# bacterial membrane proteins: taxonomy:"Bacteria [2]" (locations:(location:"Cell envelope [SL-0036]") OR locations:(location:"Membrane [SL-0162]"))
# bacterial cytoplasmic proteins: taxonomy:"Bacteria [2]" locations:(location:"Cytoplasm [SL-0086]")

awk 'BEGIN {FS="\t"} NR==FNR {a[$1]; next} ($6 in a) {print $2}' phibase_bacpatho_taxids.txt <(zcat uniprot_bacteria_effector_name.tab.gz uniprot_bacteria_effector_location.tab.gz uniprot_bacteria_effector_phibase.tab.gz) | sort -u > effector_uniprotid.txt
awk 'BEGIN {FS="\t"} NR==FNR {a[$1]; next} ($6 in a) {print $6}' phibase_bacpatho_taxids.txt <(zcat uniprot_bacteria_effector_name.tab.gz uniprot_bacteria_effector_location.tab.gz uniprot_bacteria_effector_phibase.tab.gz) | sort -u > effector_taxid.txt
awk 'BEGIN {FS="\t"} NR==FNR {a[$1]; next} ($6 in a) {print $2}' effector_taxid.txt <(zcat uniprot_bacteria_secreted.tab.gz) | sort -u > secreted_uniprotid_temp.txt
awk 'BEGIN {FS="\t"} NR==FNR {a[$1]; next} ($6 in a) {print $2}' effector_taxid.txt <(zcat uniprot_bacteria_membrane.tab.gz) | sort -u > membrane_uniprotid_temp.txt
awk 'BEGIN {FS="\t"} NR==FNR {a[$1]; next} ($6 in a) {print $2}' effector_taxid.txt <(zcat uniprot_bacteria_cytoplasm.tab.gz) | sort -u > cytoplasm_uniprotid_temp.txt
comm -23 secreted_uniprotid_temp.txt <(cat effector_uniprotid.txt membrane_uniprotid_temp.txt | sort -u) > secreted_uniprotid.txt
comm -23 membrane_uniprotid_temp.txt effector_uniprotid.txt > membrane_uniprotid.txt
comm -23 cytoplasm_uniprotid_temp.txt <(cat effector_uniprotid.txt secreted_uniprotid.txt membrane_uniprotid.txt | sort -u) > cytoplasm_uniprotid.txt
rm -f *_uniprotid_temp.txt
cat effector_uniprotid.txt secreted_uniprotid.txt membrane_uniprotid.txt cytoplasm_uniprotid.txt > esmc_uniprotid.txt
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($2 in a) {print $1,$2,$8,$9,$10,$11,$13}' esmc_uniprotid.txt <(zcat idmapping_selected.tab.gz) > esmc_uniref_uniparc_taxid.tab

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $0}' esmc_uniprotid.txt bacteria_proteome.tab > esmc_proteome.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {print $6,$2}' esmc_uniref_uniparc_taxid.tab > esmc_uniparc2uniprot.tab
cut -f 1 esmc_uniparc2uniprot.tab | sort -u > esmc_uniparc.txt
awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$2,$4,$5,$6}' esmc_uniparc.txt uniparc2pfam.csv > esmc_uniparc2pfam.tab
join -j 1 -o 1.2,2.2,2.3,2.4,2.5 -t $'\t' <(sort esmc_uniparc2uniprot.tab) <(sort esmc_uniparc2pfam.tab) | uniq > pfam_esmc.tab
awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $3; sep[$1]=";"} END {for(k in f2) print k,f2[k]}' pfam_esmc.tab > esmc_uniprotid2pfam.tab
awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$2,$5,$6}' esmc_uniparc.txt uniparc2mobidblt.csv > esmc_uniparc2mobidblt.tab
join -j 1 -o 1.2,2.2,2.3,2.4 -t $'\t' <(sort esmc_uniparc2uniprot.tab) <(sort esmc_uniparc2mobidblt.tab) | uniq > mobidblt_esmc.tab

# cut -f 4 esmc_uniref_uniparc_taxid.tab | sed 's/UniRef90_//g' | sort -u > esmc_uniref90_uniprotac.txt
# awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $1,$2,$8,$9,$10,$11,$13}' esmc_uniref90_uniprotac.txt <(zcat idmapping_selected.tab.gz) > esmc90_uniref_uniparc_taxid.tab
# awk 'BEGIN {FS="\t"; OFS="\t"} {print $6,$2}' esmc90_uniref_uniparc_taxid.tab > esmc90_uniparc2uniprot.tab
# cut -f 1 esmc90_uniparc2uniprot.tab | sort -u > esmc90_uniparc.txt
# awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$2,$4,$5,$6}' esmc90_uniparc.txt uniparc2pfam.csv > esmc90_uniparc2pfam.tab
# join -j 1 -o 1.2,2.2,2.3,2.4,2.5 -t $'\t' <(sort esmc90_uniparc2uniprot.tab) <(sort esmc90_uniparc2pfam.tab) | uniq > pfam_esmc90.tab
# awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $3; sep[$1]=";"} END {for(k in f2) print k,f2[k]}' pfam_esmc90.tab > esmc90_uniprotid2pfam.tab
# awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$2,$5,$6}' esmc90_uniparc.txt uniparc2mobidblt.csv > esmc90_uniparc2mobidblt.tab
# join -j 1 -o 1.2,2.2,2.3,2.4 -t $'\t' <(sort esmc90_uniparc2uniprot.tab) <(sort esmc90_uniparc2mobidblt.tab) | uniq > mobidblt_esmc90.tab

# cut -f 5 esmc_uniref_uniparc_taxid.tab | sed 's/UniRef50_//g' | sort -u > esmc_uniref50_uniprotac.txt
# awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {print $1,$2,$8,$9,$10,$11,$13}' esmc_uniref50_uniprotac.txt <(zcat idmapping_selected.tab.gz) > esmc50_uniref_uniparc_taxid.tab
# awk 'BEGIN {FS="\t"; OFS="\t"} {print $6,$2}' esmc50_uniref_uniparc_taxid.tab > esmc50_uniparc2uniprot.tab
# cut -f 1 esmc50_uniparc2uniprot.tab | sort -u > esmc50_uniparc.txt
# awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$2,$4,$5,$6}' esmc50_uniparc.txt uniparc2pfam.csv > esmc50_uniparc2pfam.tab
# join -j 1 -o 1.2,2.2,2.3,2.4,2.5 -t $'\t' <(sort esmc50_uniparc2uniprot.tab) <(sort esmc50_uniparc2pfam.tab) | uniq > pfam_esmc50.tab
# awk 'BEGIN {FS="\t"; OFS="\t"} {f2[$1]=f2[$1] sep[$1] $3; sep[$1]=";"} END {for(k in f2) print k,f2[k]}' pfam_esmc50.tab > esmc50_uniprotid2pfam.tab
# awk 'BEGIN {FS=","; OFS="\t"} NR==FNR {a[$1]; next} ($1 in a) {$1=$1; print $1,$2,$5,$6}' esmc50_uniparc.txt uniparc2mobidblt.csv > esmc50_uniparc2mobidblt.tab
# join -j 1 -o 1.2,2.2,2.3,2.4 -t $'\t' <(sort esmc50_uniparc2uniprot.tab) <(sort esmc50_uniparc2mobidblt.tab) | uniq > mobidblt_esmc50.tab
