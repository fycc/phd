#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

clinvar = read.delim("database/clinvar/variant_summary.txt")
colnames(clinvar)[1] = "AlleleID"
clinvar$PhenotypeIDs = gsub(",", ";", clinvar$PhenotypeIDs)
clinvar = clinvar[grepl("^NP_", clinvar$HGVS.p..), ]
clinvar$refseq = as.character(sapply(clinvar$HGVS.p.., function(x) unlist(strsplit(x, ":"))[1]))
writeLines(unique(clinvar$refseq), "database/clinvar/clinvar_refseq.txt")

clinvar_r2u = read.delim("database/clinvar/clinvar_refseq2uniprot.tab")
colnames(clinvar_r2u)[grep("yourlist", colnames(clinvar_r2u))] = "refseq"
writeLines(unique(clinvar_r2u$Entry), "database/clinvar/clinvar_uniprot.txt")

clinvar_u2u = read.delim("database/clinvar/clinvar_uniprot2uniparc.tab")
colnames(clinvar_u2u)[grep("yourlist", colnames(clinvar_u2u))] = "uniprot"

annocols = grep("Gene.ontology|Cross.reference|Involvement", colnames(clinvar_r2u))
clinvar_r2u = clinvar_r2u[order(clinvar_r2u$Status, -clinvar_r2u$Length, -nchar(clinvar_r2u$Cross.reference..Pfam.), -apply(clinvar_r2u[, annocols], 1, function(x) sum(nchar(x)))), ]
clinvar_uniprot = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Entry[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$uniprot = clinvar_uniprot[clinvar$refseq]
clinvar_uniparc = sapply(unique(clinvar$uniprot), function(x) clinvar_u2u$Entry[grep(paste("\\b", x, "\\b", sep = ""), clinvar_u2u$uniprot)[1]])
clinvar$uniparc = clinvar_uniparc[clinvar$uniprot]
clinvar = clinvar[!is.na(clinvar$uniprot) & !is.na(clinvar$uniparc), ]

load("database/clinvar/humsavar.RData")

clinvar = clinvar[clinvar$uniparc %in% ppi_pfam$uniparc, ]
clinvar = clinvar[!grepl("[=|?]", clinvar$HGVS.p..), ]

library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, "_", "del", "dup", "fs", "ins", "Ter"), collapse = "|")
positions = sapply(clinvar$HGVS.p.., function(x) {
  x = unlist(strsplit(x, "p.", fixed = TRUE))[2]
  x = unlist(strsplit(x, aa))
  paste(x[!is.na(as.numeric(x))], collapse = "_")
})

clinvar$aapos = as.character(positions)
clinvars = split(clinvar$aapos, clinvar$uniprot)
clinvars = sapply(clinvars, function(x) unique(unlist(strsplit(x, "_"))))

pfams = sapply(names(clinvars), function(x) {
  pfam = subset(ppi_pfam, V1 == x)
  z = sapply(clinvars[[x]], function(position) {
    i = apply(pfam, 1, function(p) position %in% p["V4"]:p["V5"])
    paste(pfam$V6[i], collapse = "_")
  })
  z[nchar(z) > 0]
})

pfams = pfams[lengths(pfams) > 0]
pfams = sapply(pfams, function(pfam) sapply(pfam, function(x) unique(unlist(strsplit(x, "_")))))

cvs = data.frame(unlist(pfams))
colnames(cvs) = "pfam"
cvs$uniprot = as.character(sapply(names(unlist(pfams)), function(x) unlist(strsplit(x, "[.]"))[1]))
cvs$aapos = as.character(sapply(names(unlist(pfams)), function(x) unlist(strsplit(x, "[.]"))[2]))
clinvar$pfam = apply(clinvar, 1, function(x) paste(unique(cvs$pfam[cvs$uniprot == x["uniprot"] & cvs$aapos %in% unlist(strsplit(x["aapos"], "_"))]), collapse = "_"))

# medgen_names = data.frame(CUI = as.character(sapply(readLines("database/clinvar/MedGen_names.RRF")[-1], function(x) unlist(strsplit(x, "[|]"))[1])),
#                           name = as.character(sapply(readLines("database/clinvar/MedGen_names.RRF")[-1], function(x) unlist(strsplit(x, "[|]"))[2])))
medgen_names = read.csv("database/clinvar/MedGen_names.csv")
medgen_names$CUI = paste("MedGen:", medgen_names$CUI, sep = "")
c2d = sapply(unique(clinvar$PhenotypeIDs), function(x) {
  pids = grep("^MedGen:", unlist(strsplit(x, ";")), value = TRUE)
  diseases = medgen_names$name[match(pids, medgen_names$CUI)]
  diseases = setdiff(diseases, c("not provided", "not specified"))
  diseases = paste(diseases, collapse = ";")
  diseases
})
clinvar$DiseaseNames = c2d[clinvar$PhenotypeIDs]

save(clinvar, file = "database/clinvar/clinvar.RObject")
