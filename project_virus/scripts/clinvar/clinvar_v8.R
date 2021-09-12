#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)
library(stringr)

clinvar = fread("database/clinvar/clinvar.txt", colClasses = "character")
names(clinvar) = gsub("[^A-Z|a-z]+", "", names(clinvar))
# exclude VarIDs already included in UniProtKB humsavar.txt
clinvar = clinvar[!grepl("UniProtKB", OtherIDs, ignore.case = TRUE), ]
clinvar = unique(clinvar[, .(AlleleID, PhenotypeIDS)])
setnames(clinvar, c("AlleleID", "PhenotypeIDs"))
clinvar$PhenotypeIDs = gsub("MedGen:CN169374|MedGen:CN221809", "", clinvar$PhenotypeIDs)
clinvar$PhenotypeIDs = gsub("MedGen:", "UMLS:", clinvar$PhenotypeIDs)
clinvar$PhenotypeIDs = gsub("Orphanet:ORPHA", "ORPHA:", clinvar$PhenotypeIDs)
clinvar = unique(clinvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";"))), by = "AlleleID"])
clinvar = clinvar[nchar(PhenotypeIDs) > 0]
clinvar = unique(clinvar[, c("PhenotypeIDs") := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = "AlleleID"])

hgvs = fread("database/clinvar/hgvs4variation.txt", colClasses = "character")
names(hgvs)[grep("Symbol", names(hgvs), ignore.case = TRUE)] = "GeneSymbol"
hgvs = unique(hgvs[, .(GeneSymbol, GeneID, AlleleID, ProteinExpression)])
hgvs = hgvs[grepl("^NP_", ProteinExpression) & grepl("p.", ProteinExpression, fixed = TRUE) & !grepl("[=|?]|\\[(.*?)\\]|\\((.*?)\\)|ext", ProteinExpression) & ProteinExpression != "p.0", ]
hgvs$ProteinExpression = gsub("[*]", "Ter", hgvs$ProteinExpression)
hgvs$refseq = str_extract(hgvs$ProteinExpression, "NP_\\d+")
hgvs = merge(clinvar, hgvs, by = "AlleleID")
setorder(hgvs)
writeLines(unique(hgvs$refseq), "database/clinvar/clinvar_refseq.txt")

# Download all human refseq protein sequences in fasta format
# http://www.ncbi.nlm.nih.gov/protein?term=(homo%20sapiens%5BOrganism%5D)%20AND%20srcdb_refseq%5BProperties%5D
refseq_human = read.fasta("database/intact/refseq_human.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
refseq_human = refseq_human[grep("NP_", names(refseq_human))]
names(refseq_human) = str_extract(names(refseq_human), "NP_\\d+")
clinvar_fasta = refseq_human[unique(hgvs$refseq)]
write.fasta(clinvar_fasta, names(clinvar_fasta), "database/clinvar/clinvar.fasta", as.string = TRUE)
# qsub -pe omp 16 database/clinvar/clinvar_pfam.sh

# fill in GeneSymbol and GeneID for hgvs
clinvar_refseq = fread("database/clinvar/clinvar_refseq.tab", na.strings = NULL)
names(clinvar_refseq)[grep("yourlist", names(clinvar_refseq))] = "refseq"
names(clinvar_refseq)[grep("isomap", names(clinvar_refseq))] = "isomap"
clinvar_refseq = unique(clinvar_refseq[, list(refseq = unlist(strsplit(refseq, ","))), by = setdiff(names(clinvar_refseq), "refseq")])
genecols = c("Gene names  (primary )", "Gene names", "Cross-reference (GeneID)")
clinvar_refseq[, c(genecols) := lapply(.SD, str_replace, pattern = "(?=[;|,| ])(.*?)$", replacement = ""), .SDcols = genecols]
clinvar_refseq[is.na(`Gene names  (primary )`) | nchar(`Gene names  (primary )`) == 0, `Gene names  (primary )` := `Gene names`]
annocols = setdiff(grep("Gene ontology IDs|Cross-reference", names(clinvar_refseq), value = TRUE), grep("GeneID|Pfam", names(clinvar_refseq), value = TRUE))
clinvar_refseq$nPfam = nchar(clinvar_refseq$`Cross-reference (Pfam)`)
clinvar_refseq$nAnno = apply(clinvar_refseq[, annocols, with = FALSE], 1, function(x) sum(nchar(gsub(" |;|,", "", unlist(x))), na.rm = TRUE))
clinvar_refseq$nGene = apply(clinvar_refseq[, genecols, with = FALSE], 1, function(x) sum(nchar(gsub(" |;|,", "", unlist(x))) > 0, na.rm = TRUE))
clinvar_refseq = setorder(clinvar_refseq, -nPfam, -nAnno, -nGene, Status, Entry)
clinvar_refseq = unique(clinvar_refseq, by = "refseq")

hgvs[GeneID == "-", c("GeneID") := sapply(refseq, function(x) clinvar_refseq[match(x, refseq), `Cross-reference (GeneID)`])]
hgvs[GeneSymbol == "-", c("GeneSymbol") := sapply(refseq, function(x) clinvar_refseq[match(x, refseq), `Gene names  (primary )`])]

############################
# Map SNPs to PFAM domains #
############################

library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, names(AMINO_ACID_CODE), "[+|_|=]", "del", "dup", "fs", "ins", "Ter"), collapse = "|")

hgvsp2aapos = function(hgvsp) {
  pos = tail(unlist(strsplit(hgvsp, "p.", fixed = TRUE)), 1)
  pos = unlist(strsplit(pos, aa))
  pos = setdiff(as.numeric(pos[nchar(pos) > 0]), NA)
  
  if (length(pos) == 0) NA
  else if (grepl("fs|Ter", hgvsp)) paste(pos[1], nchar(clinvar_fasta[[str_extract(hgvsp, "NP_\\d+")]]), sep = ":")
  else if (length(pos) == 1) pos
  else paste(pos[1], pos[2], sep = ":")
}

aapos = sapply(unique(hgvs$ProteinExpression), hgvsp2aapos)
stopifnot(all(sapply(aapos, function(x) class(eval(parse(text = x)))) %in% c("integer", "numeric")))
hgvs$aapos = as.character(aapos[hgvs$ProteinExpression])

clinvar_pfam = fread("database/clinvar/clinvar_pfam.csv")
clinvar_pfam$pfamacc = str_extract(clinvar_pfam$pfamacc, "PF\\d{5}")
clinvar_pfams = lapply(split(clinvar_pfam[, .(pfamacc, pfamstart, pfamend)], clinvar_pfam$seqid), unique)

hgvs = hgvs[refseq %in% names(clinvar_pfams)]
refseq_aapos = unique(paste(hgvs$refseq, hgvs$aapos, sep = ";"))
pfams_hgvs = mclapply(refseq_aapos, function(x) {
  
  refseq = unlist(strsplit(x, ";"))[1]
  aapos = unlist(strsplit(x, ";"))[2]
  
  pfam = clinvar_pfams[[refseq]]
  i = apply(pfam, 1, function(p) any(eval(parse(text = aapos)) %in% p["pfamstart"]:p["pfamend"]))
  paste(unique(pfam[i, pfamacc]), collapse = ";")
  
}, mc.cores = detectCores())
names(pfams_hgvs) = refseq_aapos
pfams_hgvs = unlist(pfams_hgvs)
hgvs$pfam = as.character(pfams_hgvs[paste(hgvs$refseq, hgvs$aapos, sep = ";")])
hgvs = unique(hgvs[nchar(pfam) > 0])

###########################################
# Combine humsavar and clinvar into hcvar #
###########################################
load("database/clinvar/humsavar.RData")
hcvar = unique(rbind(hgvs[, .(GeneSymbol, GeneID, AlleleID, pfam, PhenotypeIDs)], humsavar[, .(GeneSymbol, GeneID, AlleleID, pfam, PhenotypeIDs)]))
hcvar = unique(hcvar[, list(pfam = unlist(strsplit(pfam, ";"))), by = setdiff(names(hcvar), "pfam")])
hcvar = unique(hcvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";|,"))), by = setdiff(names(hcvar), "PhenotypeIDs")])

# map Pfam domains to clans, to address potential biases introduced by overcounting domains with very similar sequences, structures and HMM profiles.
clan = fread("database/3did_elm/clan_membership.txt", header = FALSE)
names(clan) = c("clan", "pfam")
hcvar$clan = clan[match(hcvar$pfam, pfam), clan]
hcvar[is.na(clan), c("clan") := .(pfam)]
hcvar$genepfam = paste(hcvar$GeneSymbol, hcvar$pfam, sep = "|")
hcvar$geneclan = paste(hcvar$GeneSymbol, hcvar$clan, sep = "|")
hcvar = hcvar[grepl("\\d+$", PhenotypeIDs)]
hcvar$PhenotypeSource = str_extract(hcvar$PhenotypeIDs, "^(.*?)(?=[:])")
setcolorder(hcvar, c("GeneSymbol", "GeneID", "AlleleID", "pfam", "genepfam", "clan", "geneclan", "PhenotypeSource", "PhenotypeIDs"))

####################################################################################
# map PhenotypeIDs to Disease Ontology, MeSH, and Orphanet IDs and their top nodes #
####################################################################################
# import mapping files
mgconso = fread("database/clinvar/MGCONSO.csv", colClasses = "character")
mgconso$CUI = paste0("UMLS:", mgconso$CUI)
mgconso$SAB = gsub("HPO", "Human Phenotype Ontology", mgconso$SAB)
mgconso$SAB = gsub("MSH", "MeSH", mgconso$SAB)
mgconso$SAB = gsub("ORDO", "ORPHA", mgconso$SAB)
mgconso$SAB = gsub("SNOMEDCT_US", "SNOMED CT", mgconso$SAB)
mgconso[SAB %in% c("NCI", "SNOMED CT"), SDUI := SCUI]
mgconso$SDUI = gsub("Orphanet_", "", mgconso$SDUI)
mgconso$SDUI = paste(mgconso$SAB, mgconso$SDUI, sep = ":")
mgconso = mgconso[order(TS, STT, -ISPREF, -nchar(STR))]
setkeyv(mgconso, c("CUI", "SDUI", "SAB"))
mgconso = unique(mgconso[!CUI %in% c("UMLS:CN169374", "UMLS:CN221809") & !grepl("obsolete", STR, ignore.case = TRUE), .(CUI, SDUI, SAB, STR)])

doid2xref = fread("database/clinvar/doid2xref.csv", colClasses = "character")
doid2xref$xref = gsub("(?<=EFO)(.*?)(?=:)", "", doid2xref$xref, perl = TRUE)
doid2xref$xref = gsub("EFO:", "EFO:EFO_", doid2xref$xref)
doid2xref$xref = gsub("HP:", "Human Phenotype Ontology:HP:", doid2xref$xref)
doid2xref$xref = gsub("MSH:", "MeSH:", doid2xref$xref)
doid2xref$xref = gsub("(?<=NCI)(.*?)(?=:)", "", doid2xref$xref, perl = TRUE)
doid2xref$xref = gsub("OMM:", "OMIM:", doid2xref$xref)
doid2xref$xref = gsub("Orphanet:", "ORPHA:", doid2xref$xref)
doid2xref$xref = gsub("ORDO:", "ORPHA:", doid2xref$xref)
doid2xref$xref = gsub("(?<=SNOMEDCT)(.*?)(?=:)", "", doid2xref$xref, perl = TRUE)
doid2xref$xref = gsub("SNOMEDCT:", "SNOMED CT:", doid2xref$xref)
doid2xref$xref = gsub("UMLS_CUI:", "UMLS:", doid2xref$xref)
doid2xref$xref = gsub("^:", "", doid2xref$xref)
doid2parent = fread("database/clinvar/doid2parent.csv", colClasses = "character")

mesh = fread("database/clinvar/mesh.csv", colClasses = "character")
mesh$DescriptorUI = paste0("MeSH:", mesh$DescriptorUI)

orpha = fread("database/clinvar/orpha.csv", colClasses = "character")
orpha2xref = fread("database/clinvar/orpha2xref.csv", colClasses = "character")
orphalinear = fread("database/clinvar/orphalinear.csv", colClasses = "character")

disease_names = fread("database/clinvar/disease_names.txt", colClasses = "character")
names(disease_names)[1] = "DiseaseName"
disease_names$ConceptID = paste0("UMLS:", disease_names$ConceptID)
disease_names = disease_names[!ConceptID %in% c("UMLS:CN169374", "UMLS:CN221809")]
disease_names[nchar(DiseaseMIM) > 0, c("DiseaseMIM") := paste0("OMIM:", DiseaseMIM)]
disease_names[SourceName %in% c("Human Phenotype Ontology", "SNOMED CT"), c("SourceID") := paste(SourceName, SourceID, sep = ":")]
disease_names = data.table(DiseaseID = c(disease_names$ConceptID, disease_names[nchar(DiseaseMIM) > 0, DiseaseMIM], disease_names[SourceName %in% c("Human Phenotype Ontology", "SNOMED CT"), SourceID], mesh$DescriptorUI, orpha$OrphaNumber),
                           DiseaseName = c(disease_names$DiseaseName, disease_names[nchar(DiseaseMIM) > 0, DiseaseName], disease_names[SourceName %in% c("Human Phenotype Ontology", "SNOMED CT"), DiseaseName], mesh$DescriptorName, orpha$DiseaseName))
disease_names = unique(rbindlist(list(disease_names, mgconso[!CUI %in% disease_names$DiseaseID, .(CUI, STR)], mgconso[!SDUI %in% disease_names$DiseaseID, .(SDUI, STR)])))

hcvar$DiseaseName = disease_names[match(hcvar$PhenotypeIDs, DiseaseID), DiseaseName]
hcvar = hcvar[!is.na(DiseaseName)]

# functions to map PhenotypeIDs to lowest level (most disease-specific) nodes of Disease Ontology, MeSH and ORPHA.
setkeyv(doid2xref, "xref")
setkeyv(mgconso, "CUI")
setkeyv(orpha2xref, "ExternalID")

xref2doid = function(xid) {
  if (grepl("^UMLS:", xid)) paste(sort(union(doid2xref[.(xid), DOID], doid2xref[.(mgconso[.(xid), SDUI]), DOID])), collapse = ";")
  else paste(sort(union(doid2xref[.(xid), DOID], doid2xref[.(mgconso[SDUI == xid, CUI]), DOID])), collapse = ";")
}

xref2mesh = function(xid) {
  if (grepl("^MeSH:", xid)) xid
  else if (grepl("^UMLS:", xid)) paste(sort(unique(mgconso[SAB == "MeSH"][.(xid), SDUI])), collapse = ";")
  else paste(sort(unique(mgconso[SAB == "MeSH"][.(mgconso[SDUI == xid, CUI]), SDUI])), collapse = ";")
}

xref2orpha = function(xid) {
  if (grepl("^ORPHA:", xid)) xid
  else if (grepl("^UMLS:", xid)) paste(sort(Reduce(union, list(orpha2xref[.(xid), OrphaNumber], orpha2xref[.(mgconso[.(xid), SDUI]), OrphaNumber], mgconso[SAB == "ORPHA"][.(xid), SDUI]))), collapse = ";")
  else paste(sort(Reduce(union, list(orpha2xref[.(xid), OrphaNumber], orpha2xref[.(mgconso[SDUI == xid, CUI]), OrphaNumber], mgconso[SAB == "ORPHA"][.(mgconso[SDUI == xid, CUI]), SDUI]))), collapse = ";")
}

# functions to map PhenotypeIDs to parent terms and top levels of Disease Ontology, MeSH and ORPHA.
setkeyv(doid2parent, "DOID")
setkeyv(mesh, "DescriptorUI")
setkeyv(orphalinear, "OrphaNumber")

# some DOIDs have multiple parent nodes. some parent nodes reach the top node "DOID:4" earlier than others, which will return NA in subsequent doid2tree iterations.
doid2tree = function(x) {y = NULL; while (any(x %in% doid2parent$DOID)) {x = doid2parent[.(x), ParentID]; y = c(y, x)}; paste(setdiff(y, c("DOID:4", NA)), collapse = ";")}
mesh2tree = function(x) paste(sort(unique(mesh[.(x), TreeNumber])), collapse = ";")
orpha2top = function(x) paste(sort(unique(orphalinear[.(x), ORDOID])), collapse = ";")

# mapping PhenotypeIDs to Disease Ontology, MeSH and ORPHA IDs.
pheno2doid = sapply(unique(hcvar$PhenotypeIDs), xref2doid)
pheno2mesh = sapply(unique(hcvar$PhenotypeIDs), xref2mesh)
pheno2orpha = sapply(unique(hcvar$PhenotypeIDs), xref2orpha)

pheno2doidtree = sapply(unique(hcvar$PhenotypeIDs), function(x) {
  if (nchar(pheno2doid[[x]]) == 0) return("")
  else {
    doid = unlist(strsplit(pheno2doid[[x]], ";"))
    doidtree = as.character(sapply(doid, doid2tree))
    paste(unique(unlist(strsplit(doidtree, ";"))), collapse = ";")
  }
})

pheno2meshtree = sapply(unique(hcvar$PhenotypeIDs), function(x) {
  if (nchar(pheno2mesh[[x]]) == 0) return("")
  else {
    mesh = unlist(strsplit(pheno2mesh[[x]], ";"))
    meshtree = as.character(sapply(mesh, mesh2tree))
    paste(sort(unique(unlist(strsplit(meshtree, ";")))), collapse = ";")
  }
})

pheno2orphatop = sapply(unique(hcvar$PhenotypeIDs), function(x) {
  if (nchar(pheno2orpha[[x]]) == 0) return("")
  else {
    orpha = unlist(strsplit(pheno2orpha[[x]], ";"))
    orphatop = as.character(sapply(orpha, orpha2top))
    paste(sort(unique(unlist(strsplit(orphatop, ";")))), collapse = ";")
  }
})

hcvar$DOID = as.character(pheno2doid[hcvar$PhenotypeIDs])
hcvar$MeSH = as.character(pheno2mesh[hcvar$PhenotypeIDs])
hcvar$ORPHA = as.character(pheno2orpha[hcvar$PhenotypeIDs])

hcvar$DOIDtree = as.character(pheno2doidtree[hcvar$PhenotypeIDs])
hcvar$MeSHtree = as.character(pheno2meshtree[hcvar$PhenotypeIDs])
hcvar$ORPHAtop = as.character(pheno2orphatop[hcvar$PhenotypeIDs])

setorder(hcvar)
hcvar = unique(hcvar)

save.image("database/clinvar/clinvar.RData")
