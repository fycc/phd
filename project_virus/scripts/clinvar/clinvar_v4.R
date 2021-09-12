#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)

clinvar = fread("database/clinvar/clinvar.txt", colClasses = "character")
colnames(clinvar) = gsub("[#|()|.| ]", "", colnames(clinvar))
clinvar = clinvar[grepl("^NP_", clinvar$HGVSp), ]
clinvar$refseq = as.character(sapply(clinvar$HGVSp, function(x) unlist(strsplit(x, ":"))[1]))
clinvar$refseq = as.character(sapply(clinvar$refseq, function(x) unlist(strsplit(x, "[.]"))[1]))
writeLines(unique(clinvar$refseq), "database/clinvar/clinvar_refseq.txt")

# Download all human refseq protein sequences in fasta format
# http://www.ncbi.nlm.nih.gov/protein?term=(homo%20sapiens%5BOrganism%5D)%20AND%20srcdb_refseq%5BProperties%5D
refseq_fasta = read.fasta("database/clinvar/human_refseq.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
refseq_fasta = refseq_fasta[grep("NP_", names(refseq_fasta))]
names(refseq_fasta) = as.character(sapply(names(refseq_fasta), function(x) grep("^NP_", unlist(strsplit(x, "[|]")), value = TRUE)))
names(refseq_fasta) = as.character(sapply(names(refseq_fasta), function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_fasta = refseq_fasta[unique(clinvar$refseq)]
write.fasta(clinvar_fasta, names(clinvar_fasta), "database/clinvar/clinvar.fasta", as.string = TRUE)
# qsub -pe omp 16 database/clinvar/clinvar_pfam.sh

clinvar_r2u = read.delim("database/clinvar/clinvar_refseq2uniprot.tab")
colnames(clinvar_r2u)[grep("yourlist", colnames(clinvar_r2u))] = "refseq"
colnames(clinvar_r2u)[grep("isomap", colnames(clinvar_r2u))] = "isomap"
annocols = grepl("Gene.ontology.IDs|Cross.reference", colnames(clinvar_r2u)) & !grepl("GeneID", colnames(clinvar_r2u))

y = clinvar_r2u[duplicated(clinvar_r2u$refseq) | duplicated(clinvar_r2u$refseq, fromLast = TRUE), ]
y = y[order(y$refseq, y$Status, -y$Length, -nchar(y$Cross.reference..Pfam.), -apply(y[, annocols], 1, function(x) {
  x = unlist(x)
  x = gsub(" |;|,", "", x[!is.na(x)])
  sum(nchar(x))
})), ]
y = y[!duplicated(y$refseq), ]
clinvar_r2u = rbind.data.frame(clinvar_r2u[!clinvar_r2u$refseq %in% y$refseq, ], y)
rm(y)
clinvar_r2u$Gene.names...primary.. = as.character(sapply(clinvar_r2u$Gene.names...primary.., function(x) unlist(strsplit(x, ";| "))[1]))
clinvar_r2u$Gene.names = as.character(sapply(clinvar_r2u$Gene.names, function(x) unlist(strsplit(x, ";| "))[1]))
clinvar_r2u$Gene.names...primary..[is.na(clinvar_r2u$Gene.names...primary..) | clinvar_r2u$Gene.names...primary.. == ""] = clinvar_r2u$Gene.names[is.na(clinvar_r2u$Gene.names...primary..) | clinvar_r2u$Gene.names...primary.. == ""]
clinvar_r2u$Cross.reference..GeneID. = as.character(sapply(clinvar_r2u$Cross.reference..GeneID., function(x) unlist(strsplit(x, ";"))[1]))

clinvar_uniprot = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Entry[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$uniprot = as.character(clinvar_uniprot[clinvar$refseq])
clinvar_geneid = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Cross.reference..GeneID.[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$GeneID = as.character(clinvar_geneid[clinvar$refseq])
clinvar_genesymbol = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Gene.names...primary..[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$GeneSymbol = as.character(clinvar_genesymbol[clinvar$refseq])
clinvar$GeneSymbol[is.na(clinvar$GeneSymbol)] = as.character(sapply(clinvar$Name[is.na(clinvar$GeneSymbol)], function(x) unlist(strsplit(x, "[()]"))[2]))

clinvar = clinvar[clinvar$Assembly == "GRCh38" & grepl("OMIM:[0-9]", clinvar$PhenotypeIDs), ]
clinvar$DiseaseMIM = as.character(sapply(clinvar$PhenotypeIDs, function(x) paste(gsub("OMIM:", "", grep("^OMIM:[0-9]", unlist(strsplit(x, ";|,")), value = TRUE)), collapse = ";")))

############################
# Map SNPs to PFAM domains #
############################

clinvar = clinvar[grepl("[p.]", clinvar$HGVSp) & !grepl("[=|?|()]|ext|fs", clinvar$HGVSp), ]
clinvar$HGVSp = gsub("[*]", "Ter", clinvar$HGVSp)

library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, names(AMINO_ACID_CODE), "[+|_|=]", "del", "dup", "ins", "Ter"), collapse = "|")

hgvs2aapos = function(hgvs) {
  pos = tail(unlist(strsplit(hgvs, "p.", fixed = TRUE)), 1)
  pos = unlist(strsplit(pos, aa))
  pos = sort(as.numeric(unique(pos[nchar(pos) > 0])))
  
  if (grepl("Ter", hgvs)) paste(min(pos), nchar(clinvar_fasta[[unlist(strsplit(hgvs, "[.|:]"))[1]]]), sep = ":")
  else if (length(pos) == 1) pos
  else paste(min(pos), max(pos), sep = ":")
}

aapos = sapply(unique(clinvar$HGVSp), hgvs2aapos)
clinvar$aapos = as.character(aapos[clinvar$HGVSp])

clinvar_pfam = read.csv("database/clinvar/clinvar_pfam.csv", header = FALSE)
clinvar_pfam$V2 = as.character(sapply(clinvar_pfam$V2, function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_pfams = lapply(split(clinvar_pfam[, 2:4], clinvar_pfam$V1), unique)

clinvar = clinvar[clinvar$refseq %in% names(clinvar_pfams), ]
refseq_aapos = unique(paste(clinvar$refseq, clinvar$aapos, sep = ";"))
pfams = mclapply(refseq_aapos, function(x) {
  
  refseq = unlist(strsplit(x, ";"))[1]
  aapos = unlist(strsplit(x, ";"))[2]
  
  pfam = clinvar_pfams[[refseq]]
  i = apply(pfam, 1, function(p) any(eval(parse(text = aapos)) %in% p["V3"]:p["V4"]))
  ifelse(sum(i) == 0, NA, paste(unique(pfam$V2[i]), collapse = ";"))
  
}, mc.cores = detectCores())
names(pfams) = refseq_aapos
pfams = unlist(pfams)
clinvar$pfam = as.character(pfams[paste(clinvar$refseq, clinvar$aapos, sep = ";")])
clinvar = clinvar[!is.na(clinvar$pfam), ]

################################
# Combine humsavar and clinvar #
################################

disease_names = read.delim("database/clinvar/disease_names.txt")
colnames(disease_names)[1] = "DiseaseName"
disease_names = disease_names[!is.na(disease_names$DiseaseMIM) & disease_names$Category == "Disease", ]
disease_names = unique(disease_names[, c("DiseaseName", "DiseaseMIM")])
disease_names$DiseaseMIM = as.character(disease_names$DiseaseMIM)

load("database/clinvar/humsavar.RObject")
humsa = unique(humsavar[, list(GeneSymbol, aachange, pfam, DiseaseMIM)])
humsa$aachange = gsub("p.", "", humsa$aachange, fixed = TRUE)

clin = unique(clinvar[, list(GeneSymbol, HGVSp, pfam, DiseaseMIM)])
colnames(clin)[2] = "aachange"
clin$aachange = as.character(sapply(clin$aachange, function(x) tail(unlist(strsplit(x, "p.", fixed = TRUE)), 1)))
clin = unique(clin[, list(pfam = unlist(strsplit(pfam, ";"))), by = c(colnames(clin)[colnames(clin) != "pfam"])])
clin = unique(clin[, list(DiseaseMIM = unlist(strsplit(DiseaseMIM, ";"))), by = c(colnames(clin)[colnames(clin) != "DiseaseMIM"])])

hcvar = unique(rbind(humsa, clin))
hcvar = hcvar[hcvar$DiseaseMIM %in% disease_names$DiseaseMIM, ]
hcvar$DiseaseName = disease_names$DiseaseName[match(hcvar$DiseaseMIM, disease_names$DiseaseMIM)]
hcvar = setorder(hcvar, GeneSymbol, aachange, pfam, DiseaseMIM, DiseaseName)
rm(humsa, clin)

save.image("database/clinvar/clinvar.RData")
