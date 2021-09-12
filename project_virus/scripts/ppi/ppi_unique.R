#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=16000M
#SBATCH --job-name=ppi_unique
#SBATCH --output=%x-%j.out

options(stringsAsFactors = FALSE)
library(data.table)
library(seqinr)
library(stringr)

load("~/database/intact/ppi.RObject")

ppi_uniref100 = fread("~/database/intact/ppi_uniref100.tab", header = FALSE, colClasses = "character", na.strings = NULL)
setnames(ppi_uniref100, c("ppi_id", "uniref100"))

uniprot_regex = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\-(PRO_){0,1}\\d+){0,1}"
ppi_fasta = read.fasta("~/database/intact/ppi.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
names(ppi_fasta) = str_extract(names(ppi_fasta), uniprot_regex)

ppi[!idA %in% names(ppi_fasta), idA := cidA]
ppi[!idB %in% names(ppi_fasta), idB := cidB]
ppi[!idA %in% names(ppi_fasta), idA := ppi_uniref100[match(idA, ppi_id), uniref100]]
ppi[!idB %in% names(ppi_fasta), idB := ppi_uniref100[match(idB, ppi_id), uniref100]]
ppi = ppi[!is.na(idA) & !is.na(idB)]
ppi[, c("cidA", "cidB") := lapply(.SD, str_replace, pattern = "(?=-)(.*?)$", replacement = ""), .SDcols = c("idA", "idB")]
ppi[, setdiff(names(ppi), c("idA", "idB")) := lapply(.SD, function(x) paste(unique(unlist(x)), collapse = ";")), by = c("idA", "idB")]
stopifnot(all(ppi[, c(idA, idB, cidA, cidB)] %in% names(ppi_fasta)))

ppi_tab = fread("~/database/intact/ppi.tab", colClasses = "character", na.strings = NULL)
names(ppi_tab)[grep("yourlist", names(ppi_tab))] = "ppi_id"
names(ppi_tab)[4] = "Status"
ppi_tab = unique(ppi_tab[, list(ppi_id = unlist(strsplit(ppi_id, ","))), by = setdiff(names(ppi_tab), "ppi_id")])
ppi_tab[, uniref100 := ppi_uniref100[match(ppi_tab$ppi_id, ppi_id), uniref100]]
ppi_tab[, GeneSymbol := sapply(`Gene names  (primary )`, function(x) unlist(strsplit(x, ";"))[1])]
ppi_tab[is.na(GeneSymbol) | nchar(GeneSymbol) == 0, GeneSymbol := sapply(`Gene names`, function(x) unlist(strsplit(x, "[ |;]"))[1])]
entrez2symbol = fread("~/database/intact/entrez2symbol.tab", header = FALSE, colClasses = "character", na.strings = NULL)
setnames(entrez2symbol, c("EntrezID", "GeneSymbol"))
entrez2symbol[, Freq := length(EntrezID), by = "GeneSymbol"]
setorder(entrez2symbol, -Freq, GeneSymbol)
ppi_tab[is.na(GeneSymbol) | nchar(GeneSymbol) == 0, GeneSymbol := sapply(`Cross-reference (GeneID)`, function(x) entrez2symbol[EntrezID %in% unlist(strsplit(x, ";")), GeneSymbol][1])]
ppi_tab[is.na(GeneSymbol) | nchar(GeneSymbol) == 0, GeneSymbol := sapply(uniref100, function(x) ppi_tab[match(x, ppi_id), GeneSymbol])]
ppi_tab[is.na(GeneSymbol) | nchar(GeneSymbol) == 0, GeneSymbol := uniref100]
ppi_tab[is.na(GeneSymbol) | nchar(GeneSymbol) == 0, GeneSymbol := Entry]

ppi[, c("geneA", "geneB") := lapply(.SD, function(x) ppi_tab[match(x, ppi_id), GeneSymbol]), .SDcols = c("idA", "idB")]
colsAB = c(grep("A$", names(ppi), value = TRUE), grep("B$", names(ppi), value = TRUE))
colsBA = c(grep("B$", names(ppi), value = TRUE), grep("A$", names(ppi), value = TRUE))
ppi[geneA > geneB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ppi[taxcatA %in% c("archaea", "bacteria", "virus", "fungi") & taxcatB %in% c("animal", "human", "plant"), c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ppi[taxcatA != "human" & taxcatB == "human", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ppi[, pmid := sapply(.SD, function(x) paste(sort(unique(unlist(str_extract_all(unlist(strsplit(x, "[\\||;]")), "^pubmed:(.*?)[0-9]$")))), collapse = ";")), by = c("geneA", "geneB"), .SDcols = "pmid"]
ppi[, ID := paste(idA, idB, sep = "_")]
setcolorder(ppi, c("ID", "geneA", "geneB", "idA", "idB", "cidA", "cidB", "taxidA", "taxidB", "taxnameA", "taxnameB", "taxcatA", "taxcatB", "pmid"))
setorder(ppi)
ppi = unique(ppi)
save(ppi, file = "~/database/intact/ppi.RObject")

ppi_fasta = ppi_fasta[ppi[, sort(unique(c(idA, idB, cidA, cidB)))]]
write.fasta(ppi_fasta, names(ppi_fasta), "~/database/intact/ppi.fasta", as.string = TRUE)
save(ppi_fasta, ppi_tab, ppi_uniref100, file = "~/database/intact/ppi_fasta_tab.RData")
