#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)

load("database/clinvar/humsavar.RObject")

clinvar = fread("database/clinvar/clinvar.txt", colClasses = "character")
colnames(clinvar) = gsub("[#|()|.| ]", "", colnames(clinvar))
clinvar = clinvar[grep("^NP_", clinvar$HGVSp), ]
clinvar$refseq = as.character(sapply(clinvar$HGVSp, function(x) unlist(strsplit(x, ":"))[1]))
clinvar$refseq = as.character(sapply(clinvar$refseq, function(x) unlist(strsplit(x, "[.]"))[1]))
writeLines(unique(clinvar$refseq), "database/clinvar/clinvar_refseq.txt")

# Download all human refseq protein sequences in fasta format
# http://www.ncbi.nlm.nih.gov/protein?term=(homo%20sapiens%5BOrganism%5D)%20AND%20srcdb_refseq%5BProperties%5D
human_refseq = read.fasta("database/intact/human_refseq.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
human_refseq = human_refseq[grep("NP_", names(human_refseq))]
names(human_refseq) = as.character(sapply(names(human_refseq), function(x) grep("^NP_", unlist(strsplit(x, "[|]")), value = TRUE)))
names(human_refseq) = as.character(sapply(names(human_refseq), function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_fasta = human_refseq[unique(clinvar$refseq)]
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

clinvar_uniprot = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Entry[grep(paste0("\\b", x, "\\b"), clinvar_r2u$refseq)[1]])
clinvar$uniprot = as.character(clinvar_uniprot[clinvar$refseq])
clinvar_geneid = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Cross.reference..GeneID.[grep(paste0("\\b", x, "\\b"), clinvar_r2u$refseq)[1]])
clinvar$GeneID = as.character(clinvar_geneid[clinvar$refseq])
clinvar_genesymbol = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Gene.names...primary..[grep(paste0("\\b", x, "\\b"), clinvar_r2u$refseq)[1]])
clinvar$GeneSymbol = as.character(clinvar_genesymbol[clinvar$refseq])
clinvar$GeneSymbol[is.na(clinvar$GeneSymbol)] = as.character(sapply(clinvar$Name[is.na(clinvar$GeneSymbol)], function(x) unlist(strsplit(x, "[()]"))[2]))

clinvar = clinvar[clinvar$Assembly == "GRCh38", ]

############################
# Map SNPs to PFAM domains #
############################

clinvar = clinvar[grepl("[p.]", clinvar$HGVSp) & !grepl("[=|?|()]|ext", clinvar$HGVSp), ]
clinvar = clinvar[!grepl("p.0", clinvar$HGVSp, fixed = TRUE), ]
clinvar$HGVSp = gsub("[*]", "Ter", clinvar$HGVSp)

library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, names(AMINO_ACID_CODE), "[+|_|=]", "del", "dup", "fs", "ins", "Ter"), collapse = "|")

hgvs2aapos = function(hgvs) {
  pos = tail(unlist(strsplit(hgvs, "p.", fixed = TRUE)), 1)
  pos = unlist(strsplit(pos, aa))
  pos = as.numeric(unique(pos[nchar(pos) > 0]))
  
  if (grepl("fs|Ter", hgvs)) paste(pos[1], nchar(clinvar_fasta[[unlist(strsplit(hgvs, "[.|:]"))[1]]]), sep = ":")
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

###########################################
# Combine humsavar and clinvar into hcvar #
###########################################

disease_names = fread("database/clinvar/disease_names.txt", colClasses = "character")
colnames(disease_names)[1] = "DiseaseName"
disease_names$ConceptID = paste0("UMLS:", disease_names$ConceptID)
disease_names$DiseaseMIM = paste0("OMIM:", disease_names$DiseaseMIM)
disease_names$SourceID[disease_names$SourceName == "SNOMED CT"] = paste0("SNOMED CT:", disease_names$SourceID[disease_names$SourceName == "SNOMED CT"])
disease_names$SourceID[disease_names$SourceName == "Human Phenotype Ontology"] = paste0("Human Phenotype Ontology:", disease_names$SourceID[disease_names$SourceName == "Human Phenotype Ontology"])
mesh = fread("database/clinvar/mesh.csv", colClasses = "character")
mesh$DescriptorUI = paste0("MeSH:", mesh$DescriptorUI)
orpha = fread("database/clinvar/orpha.csv", colClasses = "character")
diseases = data.table(DiseaseID = c(disease_names$ConceptID, disease_names$DiseaseMIM, disease_names$SourceID, mesh$DescriptorUI, orpha$OrphaNumber),
                      DiseaseName = c(disease_names$DiseaseName, disease_names$DiseaseName, disease_names$DiseaseName, mesh$DescriptorName, orpha$DiseaseName))
diseases = unique(diseases[grepl("^UMLS:|OMIM:|SNOMED CT:|Human Phenotype Ontology:|MeSH:|ORPHA:", diseases$DiseaseID) & !diseases$DiseaseName %in% c("not provided", "not specified"), ])

# # exclude mutations with unknown clinical significance
# nclinmut = paste(clinvar$GeneSymbol, as.character(sapply(clinvar$HGVSp, function(x) unlist(strsplit(x, "p.", fixed = TRUE))[2])), sep = "_")
# nclinmut = unique(nclinmut[clinvar$ClinicalSignificance %in% c("not provided", "other", "Uncertain significance")])
# humsa = unique(humsavar[!paste(humsavar$GeneSymbol, as.character(sapply(humsavar$aachange, function(x) unlist(strsplit(x, "p.", fixed = TRUE))[2])), sep = "_") %in% nclinmut, list(GeneSymbol, aachange, pfam, DiseaseMIM)])
# clin = unique(clinvar[!clinvar$ClinicalSignificance %in% c("not provided", "other", "Uncertain significance"), list(GeneSymbol, HGVSp, pfam, PhenotypeIDs)])

humsa = unique(humsavar[, list(GeneSymbol, aachange, pfam, DiseaseMIM)])
colnames(humsa)[colnames(humsa) == "DiseaseMIM"] = "PhenotypeIDs"
humsa$aachange = gsub("p.", "", humsa$aachange, fixed = TRUE)
humsa$DiseaseName = diseases$DiseaseName[match(humsa$PhenotypeIDs, diseases$DiseaseID)]
humsa = unique(humsa[!is.na(humsa$DiseaseName), ])

clin = unique(clinvar[, list(GeneSymbol, HGVSp, pfam, PhenotypeIDs)])
colnames(clin)[colnames(clin) == "HGVSp"] = "aachange"
clin$aachange = as.character(sapply(clin$aachange, function(x) tail(unlist(strsplit(x, "p.", fixed = TRUE)), 1)))
clin = unique(clin[, list(pfam = unlist(strsplit(pfam, ";"))), by = setdiff(colnames(clin), "pfam")])
clin = unique(clin[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";|,"))), by = setdiff(colnames(clin), "PhenotypeIDs")])
clin$PhenotypeIDs = gsub("MedGen:", "UMLS:", clin$PhenotypeIDs)
clin$PhenotypeIDs = gsub("Orphanet:ORPHA", "ORPHA:", clin$PhenotypeIDs)
clin$DiseaseName = diseases$DiseaseName[match(clin$PhenotypeIDs, diseases$DiseaseID)]
clin = unique(clin[!is.na(clin$DiseaseName), ])

hcvar = unique(rbind(humsa, clin))

# Map OMIM numbers (too specific) to OrphaNet and Disease Ontology IDs (several OMIM numbers may be combined into one term)
orpha2mim = fread("database/clinvar/orpha2mim.csv", colClasses = "character")
doid = fread("database/clinvar/doid.csv", colClasses = "character")
diseases = unique(data.table(DiseaseCode = c(orpha2mim$OrphaNumber, doid$DOID[grep("^OMIM:", doid$xref)]),
                             OMIM = c(orpha2mim$OMIM, doid$xref[grep("^OMIM:", doid$xref)])))
code2omim = sapply(split(diseases$OMIM, diseases$DiseaseCode), unique)
omim2code = sapply(split(diseases$DiseaseCode, diseases$OMIM), function(x) sort(unique(x), decreasing = TRUE))
omim2code = sapply(omim2code, function(x) x[which.max(lengths(code2omim)[x])])

hcvar = hcvar[hcvar$DiseaseMIM %in% disease_names$DiseaseMIM, ]
hcvar$DiseaseCode = hcvar$DiseaseMIM
hcvar$DiseaseCode[hcvar$DiseaseCode %in% names(omim2code)] = as.character(omim2code[hcvar$DiseaseCode[hcvar$DiseaseCode %in% names(omim2code)]])
hcvar$DiseaseName = disease_names$DiseaseName[match(hcvar$DiseaseMIM, disease_names$DiseaseMIM)]
hcvar$DiseaseName[grep("^ORPHA:", hcvar$DiseaseCode)] = orpha2mim$DiseaseName[match(grep("^ORPHA:", hcvar$DiseaseCode, value = TRUE), orpha2mim$OrphaNumber)]
hcvar$DiseaseName[grep("^DOID:", hcvar$DiseaseCode)] = doid$DiseaseName[match(grep("^DOID:", hcvar$DiseaseCode, value = TRUE), doid$DOID)]
hcvar = unique(setorder(hcvar, GeneSymbol, aachange, pfam, DiseaseMIM, DiseaseCode, DiseaseName))

rm(humsa, clin)

save.image("database/clinvar/clinvar.RData")
