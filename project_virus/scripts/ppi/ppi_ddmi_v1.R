#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)
library(stringr)

load("database/3did_elm/ddi_dmi.RData")
load("database/intact/ppi.RObject")
load("database/intact/ppi_fasta_tab.RData")

# scan protein sequences for motifs
system.time(ppi_motif <- mclapply(ppi_fasta, function(fa) all_motifs[str_detect(fa, all_motifs)], mc.cores = detectCores()))

ppi_pfam = fread("database/intact/ppi_pfam.csv")
ppi_pfam$pfamacc = str_extract(ppi_pfam$pfamacc, "PF\\d{5}")
ppi_pfam = sapply(split(ppi_pfam$pfamacc, ppi_pfam$seqid), unique)

#################################
### scan PPIs for DDI/DMI/MDI ###
#################################
ddmi_scan = function(ppi, ddmi, domain, motif) {
  
  cond1 = ddmi$dmA %in% domain[[ppi$idA]] & ddmi$dmB %in% domain[[ppi$idB]]
  cond2 = ddmi$dmA %in% domain[[ppi$idA]] & ddmi$dmB %in% motif[[ppi$idB]]
  cond3 = ddmi$dmA %in% motif[[ppi$idA]] & ddmi$dmB %in% domain[[ppi$idB]]
  
  ddmi[cond1 | cond2 | cond3]
  
}

ppi_list = split(ppi[, .(idA, idB)], seq(nrow(ppi)))
names(ppi_list) = paste(ppi$idA, ppi$idB, sep = "_")
system.time(ppi_ddmi <- mclapply(ppi_list, ddmi_scan, ddmi = ddi_dmi, domain = ppi_pfam, motif = ppi_motif, mc.cores = detectCores()))
ppi_ddmi = ppi_ddmi[lapply(ppi_ddmi, nrow) > 0]

#####################################
### tabulate domain-resolved PPIs ###
#####################################
# map UniProt IDs to Gene IDs, since UniProt IDs are too specific (same sequence can have multiple UniProt IDs)
# order of preference: Gene Symbol > Gene ID (if no Gene Symbol exists) > UniProt ID (if neither Gene Symbol nor Gene ID exists)
ppi[, c("ID") := paste(idA, idB, sep = "_")]
ppi[, c("geneA", "geneB") := lapply(.SD, function(x) ppi_tab[match(x, ppi_id), `Gene names  (primary )`]), .SDcols = c("idA", "idB")]
ppi[, c("geneidA", "geneidB") := lapply(.SD, function(x) ppi_tab[match(x, ppi_id), `Cross-reference (GeneID)`]), .SDcols = c("idA", "idB")]
ppi[nchar(geneA) == 0, c("geneA") := .(geneidA)]
ppi[nchar(geneA) == 0, c("geneA") := .(cidA)]
ppi[nchar(geneB) == 0, c("geneB") := .(geneidB)]
ppi[nchar(geneB) == 0, c("geneB") := .(cidB)]
ppi[, c("geneidA", "geneidB") := NULL]
ppi = unique(ppi[geneA != geneB])

ddmi_dt = rbindlist(ppi_ddmi, idcol = "ID")
ddmi_dt = merge(ppi, ddmi_dt, by = "ID")
colsAB = c(grep("A$", names(ddmi_dt), value = TRUE), grep("B$", names(ddmi_dt), value = TRUE))
colsBA = c(grep("B$", names(ddmi_dt), value = TRUE), grep("A$", names(ddmi_dt), value = TRUE))
ddmi_dt[geneA > geneB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ddmi_dt[taxcatA %in% c("bacteria", "virus") & taxcatB == "host", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ddmi_dt$type = ddi_dmi[match(paste(ddmi_dt$dmA, ddmi_dt$dmB, sep = "_"), paste(dmA, dmB, sep = "_")), type]
ddmi_dt[, c("ID") := paste(geneA, geneB, sep = "_")]
ddmi_dt[, setdiff(names(ddmi_dt), c("ID", "type", "dmA", "dmB")) := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = c("ID", "type", "dmA", "dmB")]
cols = names(ddmi_dt)[sapply(ddmi_dt, function(x) any(grepl(";", x)))]
ddmi_dt[, c(cols) := lapply(.SD, function(x) paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")), by = c("ID", "type", "dmA", "dmB"), .SDcols = cols]
ddmi_dt[str_detect(pmid, "pubmed:\\d+"), c("pmid") := sapply(pmid, function(x) paste(sort(unique(unlist(str_extract_all(x, "pubmed:\\d+")))), collapse = ";"))]
setcolorder(ddmi_dt, c("ID", "geneA", "geneB", "type", "dmA", "dmB", "idA", "idB", "cidA", "cidB", "taxidA", "taxidB", "taxnameA", "taxnameB", "taxcatA", "taxcatB", "pmid"))
setorder(ddmi_dt)
ddmi_dt = unique(ddmi_dt)

colsAB = c(grep("A$", names(ppi), value = TRUE), grep("B$", names(ppi), value = TRUE))
colsBA = c(grep("B$", names(ppi), value = TRUE), grep("A$", names(ppi), value = TRUE))
ppi[geneA > geneB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ppi[taxcatA %in% c("bacteria", "virus") & taxcatB == "host", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ppi[, c("ID") := paste(geneA, geneB, sep = "_")]
ppi[, setdiff(names(ppi), "ID") := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = "ID"]
cols = names(ppi)[sapply(ppi, function(x) any(grepl(";", x)))]
ppi[, c(cols) := lapply(.SD, function(x) paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")), by = "ID", .SDcols = cols]
ppi[str_detect(pmid, "pubmed:\\d+"), c("pmid") := sapply(pmid, function(x) paste(sort(unique(unlist(str_extract_all(x, "pubmed:\\d+")))), collapse = ";"))]
setcolorder(ppi, c("ID", "geneA", "geneB", "idA", "idB", "cidA", "cidB", "taxidA", "taxidB", "taxnameA", "taxnameB", "taxcatA", "taxcatB", "pmid"))
setorder(ppi)
ppi = unique(ppi)

save(ppi_motif, ppi_pfam, ppi, ppi_ddmi, ddmi_dt, file = "database/intact/ppi_ddmi.RData")
