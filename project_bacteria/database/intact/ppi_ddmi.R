#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=ppi_ddmi
#SBATCH --output=%x-%j.out

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)
library(stringr)

load("~/database/3did_elm/ddi_dmi.RData")
load("~/database/intact/ppi.RObject")
load("~/database/intact/ppi_fasta_tab.RData")

# scan protein sequences for motifs
system.time(ppi_motif <- mclapply(ppi_fasta, function(fa) all_motifs[str_detect(fa, all_motifs)], mc.cores = 4))
ppi_motif = ppi_motif[lengths(ppi_motif) > 0]

ppi_pfam = fread("~/database/intact/ppi_pfam.csv")
ppi_pfam$pfamacc = str_extract(ppi_pfam$pfamacc, "PF\\d+")
setkeyv(ppi_pfam, c("seqid", "pfamacc"))
ppi_pfam_ddmi = ppi_pfam[pfamacc %in% ddi_dmi$dmA]

#################################
### scan PPIs for DDI/DMI/MDI ###
#################################
ppi_list = split(unique(ppi[(idA %in% ppi_pfam_ddmi$seqid & idB %in% ppi_pfam_ddmi$seqid) | (idA %in% ppi_pfam_ddmi$seqid & idB %in% names(ppi_motif)) | (idB %in% ppi_pfam_ddmi$seqid & idA %in% names(ppi_motif)), .(idA, idB)]), by = c("idA", "idB"))
names(ppi_list) = sub("[.]", "_", names(ppi_list))

ddmi = function(x) {
  ddi_dmi[(dmA %in% ppi_pfam_ddmi[.(x$idA), pfamacc] & dmB %in% ppi_pfam_ddmi[.(x$idB), pfamacc]) |
            (dmA %in% ppi_pfam_ddmi[.(x$idA), pfamacc] & dmB %in% ppi_motif[[x$idB]]) |
            (dmA %in% ppi_motif[[x$idA]] & dmB %in% ppi_pfam_ddmi[.(x$idB), pfamacc])]
}

N = split(1:length(ppi_list), ceiling(seq_along(1:length(ppi_list)) / 10000))
ppi_ddmi = list()
for(i in 1:length(N)) {
  ppi_ddmi[[i]] = mclapply(ppi_list[N[[i]]], ddmi, mc.cores = 4)
}
ppi_ddmi = do.call(c, ppi_ddmi)
ppi_ddmi = ppi_ddmi[lapply(ppi_ddmi, nrow) > 0]

#####################################
### tabulate domain-resolved PPIs ###
#####################################
ddmi_dt = rbindlist(ppi_ddmi, idcol = "ID")
ddmi_dt = merge(ppi, ddmi_dt, by = "ID")
ddmi_dt$type = ddi_dmi[match(ddmi_dt[, paste(dmA, dmB, sep = "_")], paste(dmA, dmB, sep = "_")), type]

writeLines(ppi[, union(cidA, cidB)], "~/database/intact/ppi_cids.txt")
system('zgrep -wFf ~/database/intact/ppi_cids.txt ~/database/uniprot/idmapping_selected.tab.gz | cut -f 1,2 > ~/database/intact/ppi_cids_uniprotid.tab')
ppi_cids_uniprotid = fread("~/database/intact/ppi_cids_uniprotid.tab", header = FALSE)
setnames(ppi_cids_uniprotid, c("uniprotac", "uniprotid"))

ppi[, cuidA := ppi_cids_uniprotid[match(cidA, uniprotac), uniprotid]]
ppi[, cuidB := ppi_cids_uniprotid[match(cidB, uniprotac), uniprotid]]
ppi[, ID := paste(geneA, geneB, sep = "_")]
setcolorder(ppi, c("ID", "geneA", "geneB", "idA", "idB", "cidA", "cidB", "cuidA", "cuidB", "taxidA", "taxidB", "taxnameA", "taxnameB", "taxcatA", "taxcatB", "pmid"))
setorder(ppi)
ppi = unique(ppi)

ddmi_dt[, cuidA := ppi_cids_uniprotid[match(cidA, uniprotac), uniprotid]]
ddmi_dt[, cuidB := ppi_cids_uniprotid[match(cidB, uniprotac), uniprotid]]
ddmi_dt[, ID := paste(geneA, geneB, sep = "_")]
setcolorder(ddmi_dt, c("ID", "geneA", "geneB", "type", "dmA", "dmB", "idA", "idB", "cidA", "cidB", "cuidA", "cuidB", "taxidA", "taxidB", "taxnameA", "taxnameB", "taxcatA", "taxcatB", "pmid"))
setorder(ddmi_dt)
ddmi_dt = unique(ddmi_dt)

save(ppi_motif, ppi_pfam, ppi, ppi_ddmi, ddmi_dt, file = "~/database/intact/ppi_ddmi.RData")
