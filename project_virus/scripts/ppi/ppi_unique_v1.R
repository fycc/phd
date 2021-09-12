#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(seqinr)
library(stringr)

uniprot_regex = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\-(PRO_){0,1}\\d+){0,1}"
ppi_fasta = read.fasta("database/intact/ppi.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
names(ppi_fasta) = str_extract(names(ppi_fasta), uniprot_regex)

ppi_tab = fread("database/intact/ppi.tab", colClasses = "character", na.strings = NULL)
names(ppi_tab)[grep("yourlist", names(ppi_tab))] = "ppi_id"
names(ppi_tab)[4] = "Status"
ppi_tab = unique(ppi_tab[, list(ppi_id = unlist(strsplit(ppi_id, ","))), by = setdiff(names(ppi_tab), "ppi_id")])
genecols = c("Gene names  (primary )", "Gene names", "Cross-reference (GeneID)")
ppi_tab[, c(genecols) := lapply(.SD, str_replace, pattern = "(?=[;|,| ])(.*?)$", replacement = ""), .SDcols = genecols]
ppi_tab[is.na(`Gene names  (primary )`) | nchar(`Gene names  (primary )`) == 0, `Gene names  (primary )` := `Gene names`]
annocols = setdiff(grep("Gene ontology IDs|Cross-reference", names(ppi_tab), value = TRUE), grep("GeneID|Pfam", names(ppi_tab), value = TRUE))
ppi_tab$nPfam = nchar(ppi_tab$`Cross-reference (Pfam)`)
ppi_tab$nAnno = apply(ppi_tab[, annocols, with = FALSE], 1, function(x) sum(nchar(gsub(" |;|,", "", unlist(x))), na.rm = TRUE))
ppi_tab$nGene = apply(ppi_tab[, genecols, with = FALSE], 1, function(x) sum(nchar(gsub(" |;|,", "", unlist(x))) > 0, na.rm = TRUE))
ppi_tab = setorder(ppi_tab, -nPfam, -nAnno, -nGene, Status, Entry)
ppi_tab = unique(ppi_tab, by = "ppi_id")

load("database/intact/ppi.RObject")
ppi[grepl("-PRO_", idA), c("idA") := .(cidA)]
ppi[grepl("-PRO_", idB), c("idB") := .(cidB)]
ppi = ppi[idA %in% c(names(ppi_fasta), ppi_tab$ppi_id) & idB %in% c(names(ppi_fasta), ppi_tab$ppi_id)]
# if protein ID is not in fasta file, replace with UniProt mapped ID (in both PPI table and mapping table)
ppi[!idA %in% names(ppi_fasta), c("idA") := .(ppi_tab[match(idA, ppi_id), Entry])]
ppi[!idB %in% names(ppi_fasta), c("idB") := .(ppi_tab[match(idB, ppi_id), Entry])]
ppi_tab[!ppi_id %in% names(ppi_fasta), c("ppi_id") := .(Entry)]
stopifnot(all(c(ppi$idA, ppi$idB) %in% ppi_tab$ppi_id) & all(ppi_tab$ppi_id %in% names(ppi_fasta)))

# if multiple UniProt IDs have identical sequences (same SEGUID)
# choose a representative ID that has the most Pfam domains and annotations, and preferably reviewed
ppi_seguid = fread("database/intact/ppi_seguid.csv", colClasses = "character")
ppi_seguid$ID = str_extract(ppi_seguid$ID, uniprot_regex)
ppi_seguid = ppi_seguid[ID %in% ppi_tab$ppi_id]
ppi_seguid[, "REPID" := lapply(.SD, function(ids) ids[which.min(match(ids, ppi_tab$ppi_id))]), by = "SEGUID"]

ppi[, c("idA") := .(ppi_seguid[match(idA, ID), REPID])]
ppi[, c("idB") := .(ppi_seguid[match(idB, ID), REPID])]

# get canonical IDs (cids) from isoform IDs and remove self-interactions between different isoforms
ppi[, c("cidA", "cidB") := lapply(.SD, str_replace, pattern = "(?=-)(.*?)$", replacement = ""), .SDcols = c("idA", "idB")]
ppi = ppi[cidA != cidB]

# to remove duplicates (e.g. A_B and B_A), rearrange protein A and protein B, such that protein A will always have a smaller ID than protein B
# for host-pathogen PPIs, make sure host protein ID is in front
colsAB = c(grep("A$", names(ppi), value = TRUE), grep("B$", names(ppi), value = TRUE))
colsBA = c(grep("B$", names(ppi), value = TRUE), grep("A$", names(ppi), value = TRUE))
ppi[cidA > cidB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ppi[taxcatA %in% c("bacteria", "virus") & taxcatB == "host", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]

# for interactions involving multiple isoforms of the same pair of proteins, use canonical IDs for the pair
ppi[, setdiff(names(ppi), c("cidA", "cidB")) := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = c("cidA", "cidB")]
ppi[grepl(";", idA), c("idA") := .(cidA)]
ppi[grepl(";", idB), c("idB") := .(cidB)]
cols = names(ppi)[sapply(ppi, function(x) any(grepl(";", x)))]
ppi[, c(cols) := lapply(.SD, function(x) paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")), by = c("cidA", "cidB"), .SDcols = cols]
ppi = unique(ppi[grepl("pubmed:", pmid, fixed = TRUE) & !grepl("pubmed:20711500", pmid, fixed = TRUE)])
ppi[str_detect(pmid, "pubmed:\\d+"), c("pmid") := sapply(pmid, function(x) paste(sort(unique(unlist(str_extract_all(x, "pubmed:\\d+")))), collapse = ";"))]
setorder(ppi)
ppi = unique(ppi)

save(ppi, file = "database/intact/ppi.RObject")

ppi_fasta = ppi_fasta[unique(c(ppi$idA, ppi$cidA, ppi$idB, ppi$cidB))]
write.fasta(ppi_fasta, names(ppi_fasta), "database/intact/ppi.fasta", as.string = TRUE)
save(ppi_fasta, ppi_tab, file = "database/intact/ppi_fasta_tab.RData")
