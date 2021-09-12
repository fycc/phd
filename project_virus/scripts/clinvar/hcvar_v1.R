#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)
library(stringr)

########################################################
### map humsavar (UniProt) mutations to Pfam domains ###
########################################################
load("database/intact/ppi_fasta_tab.RData")

# function to determine column widths in fixed-width text file
colwidth = function(counterline, counter, max = 10000) {
  widths = unlist(gregexpr(paste0(" ", counter), gsub("#", "", counterline)))
  widths = c(widths[1], diff(widths), max)
  widths
}

humsavar = read.fwf(file = "database/clinvar/humsavar.txt", widths = colwidth(readLines("database/clinvar/humsavar.txt")[49], "_"),
                    skip = 49, n = length(readLines("database/clinvar/humsavar.txt")) - 54)
humsavar = as.data.table(apply(humsavar, 2, trimws))
setnames(humsavar, c("GeneSymbol", "uniprot", "AlleleID", "ProteinChange", "vartype", "MutationID", "PhenotypeIDs"))
humsavar$ProteinExpression = paste(humsavar$uniprot, humsavar$ProteinChange, sep = ":")
humsavar$GSExpression = paste(humsavar$GeneSymbol, humsavar$ProteinChange, sep = ":")
protexpr_humsavar = unique(humsavar$ProteinExpression)
gsexpr_humsavar = unique(humsavar$GSExpression)
humsavar$PhenotypeIDs = str_extract(humsavar$PhenotypeIDs, "MIM:\\d+")
humsavar = humsavar[grep("^MIM:", PhenotypeIDs)]
humsavar$PhenotypeIDs = gsub("MIM:", "OMIM:", humsavar$PhenotypeIDs)
for(j in names(humsavar)) set(humsavar, j = j, value = gsub("^-$", NA, humsavar[[j]]))
humsavar = unique(humsavar[!is.na(MutationID)])
writeLines(unique(humsavar$uniprot), "database/clinvar/humsavar_uniprot.txt")
# exclude sequences in ppi.fasta, which have already been scanned for Pfam domains
uniprot_regex = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\-(PRO_){0,1}\\d+){0,1}"
humsavar_fasta = read.fasta("database/clinvar/humsavar.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
names(humsavar_fasta) = str_extract(names(humsavar_fasta), uniprot_regex)
humsavar_fasta = humsavar_fasta[setdiff(names(humsavar_fasta), names(ppi_fasta))]
write.fasta(humsavar_fasta, names(humsavar_fasta), "database/clinvar/humsavar.fasta", as.string = TRUE)
# qsub -pe omp 16 database/clinvar/humsavar_pfam.sh
# get Gene ID from UniProt mapping file: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
uniprotmap = fread("database/HUMAN_9606_idmapping_selected.tab", colClasses = "character", header = FALSE, na.strings = NULL)
setnames(uniprotmap, c("UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90", "UniRef50", "UniParc",
                       "PIR", "NCBI-taxon", "MIM", "UniGene", "PubMed", "EMBL", "EMBL-CDS", "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Additional PubMed"))
for(j in names(uniprotmap)) set(uniprotmap, j = j, value = gsub("; ", ";", uniprotmap[[j]]))
humsavar[, c("GeneID") := uniprotmap[match(uniprot, `UniProtKB-AC`), `GeneID (EntrezGene)`]]

# map SNPs to PFAM domains
library(Biostrings)
humsavar$aapos = gsub(paste(AMINO_ACID_CODE, collapse = "|"), "", gsub("p.", "", humsavar$ProteinChange, fixed = TRUE))
stopifnot(all(sapply(humsavar$aapos, function(x) class(eval(parse(text = x)))) == "numeric"))

humsavar_pfam = fread("database/clinvar/humsavar_pfam.csv")
humsavar_pfam$seqid = str_extract(humsavar_pfam$seqid, uniprot_regex)
humsavar_pfam$pfamacc = str_extract(humsavar_pfam$pfamacc, "PF\\d{5}")
ppi_pfam = fread("database/intact/ppi_pfam.csv")
ppi_pfam$pfamacc = str_extract(ppi_pfam$pfamacc, "PF\\d{5}")
humsavar_pfam = unique(rbind(humsavar_pfam, ppi_pfam[seqid %in% humsavar$uniprot]))
humsavar_pfams = lapply(split(humsavar_pfam[, .(pfamacc, pfamstart, pfamend)], humsavar_pfam$seqid), unique)

aapos2pfam = function(protaapos, pfams) {
  prot = unlist(strsplit(protaapos, ","))[1]
  aapos = unlist(strsplit(protaapos, ","))[2]
  
  pfam = pfams[[prot]]
  i = apply(pfam, 1, function(p) any(eval(parse(text = aapos)) %in% p["pfamstart"]:p["pfamend"]))
  paste(unique(pfam[i, pfamacc]), collapse = ";")
}

humsavar = humsavar[uniprot %in% names(humsavar_pfams)]
aapos_humsavar = unique(paste(humsavar$uniprot, humsavar$aapos, sep = ","))
aapos2pfam_humsavar = mclapply(aapos_humsavar, aapos2pfam, pfams = humsavar_pfams, mc.cores = detectCores())
names(aapos2pfam_humsavar) = aapos_humsavar
aapos2pfam_humsavar = unlist(aapos2pfam_humsavar)
humsavar$pfam = as.character(aapos2pfam_humsavar[paste(humsavar$uniprot, humsavar$aapos, sep = ",")])
humsavar = humsavar[nchar(pfam) > 0]
humsavar = humsavar[, list(pfam = unlist(strsplit(pfam, ";"))), by = setdiff(names(humsavar), "pfam")]

setorder(humsavar)
humsavar = unique(humsavar)

#############################################
### map ClinVar mutations to Pfam domains ###
#############################################
clinvar = fread("database/clinvar/clinvar.txt", colClasses = "character")
names(clinvar) = gsub("[^A-Z|a-z]+", "", names(clinvar))
names(clinvar) = gsub("RSdbSNP", "MutationID", names(clinvar))
names(clinvar) = gsub("PhenotypeIDS", "PhenotypeIDs", names(clinvar))
clinvar = clinvar[MutationID != "-1" & Assembly == "GRCh38"]
clinvar$MutationID = paste0("rs", clinvar$MutationID)
clinvar[, c("ChromosomePosition") := paste(Chromosome, paste(Start, Stop, sep = "-"), sep = ":")]
clinvar = unique(clinvar[, .(ClinSigSimple, AlleleID, MutationID, PhenotypeIDs, ChromosomePosition)])
clinvar$PhenotypeIDs = gsub("MedGen:CN169374|MedGen:CN221809", "", clinvar$PhenotypeIDs)
clinvar$PhenotypeIDs = gsub("MedGen:", "UMLS:", clinvar$PhenotypeIDs)
clinvar$PhenotypeIDs = gsub("Orphanet:ORPHA", "ORPHA:", clinvar$PhenotypeIDs)
clinvar = unique(clinvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";|,"))), by = setdiff(names(clinvar), "PhenotypeIDs")])
clinvar = clinvar[nchar(PhenotypeIDs) > 1]
cid_history = fread("database/clinvar/ConceptID_history.txt")
cid_history$`Previous ConceptID` = paste0("UMLS:", cid_history$`Previous ConceptID`)
cid_history$`Current ConceptID` = paste0("UMLS:", cid_history$`Current ConceptID`)
clinvar = clinvar[!PhenotypeIDs %in% cid_history[`Current ConceptID` == "UMLS:No longer reported", `Previous ConceptID`]]
clinvar[PhenotypeIDs %in% cid_history$`Previous ConceptID`, c("PhenotypeIDs") := cid_history[match(PhenotypeIDs, `Previous ConceptID`), `Current ConceptID`]]
clinvar = unique(clinvar[, c("PhenotypeIDs") := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = setdiff(names(clinvar), "PhenotypeIDs")])
clinsig = unique(clinvar[, .(MutationID, ClinSigSimple)][order(-ClinSigSimple, na.last = TRUE)])

hgvs = fread("database/clinvar/hgvs4variation.txt", colClasses = "character")
names(hgvs)[grep("Symbol", names(hgvs), ignore.case = TRUE)] = "GeneSymbol"
hgvs = unique(hgvs[, .(GeneSymbol, GeneID, AlleleID, ProteinExpression, ProteinChange)])
hgvs$refseq = str_extract(hgvs$ProteinExpression, "[A-Z]P_\\d+")
hgvs$ProteinChange = gsub("[*]", "Ter", hgvs$ProteinChange)
hgvs$ProteinChange = gsub("Terfs|fsTer", "Ter", hgvs$ProteinChange)
hgvs$ProteinExpression = paste(hgvs$refseq, hgvs$ProteinChange, sep = ":")
hgvs$GSExpression = paste(hgvs$GeneSymbol, hgvs$ProteinChange, sep = ":")
protexpr_clinvar = unique(hgvs$ProteinExpression)
gsexpr_clinvar = unique(hgvs$GSExpression)
hgvs = hgvs[grepl("^NP_", ProteinExpression) & grepl("p.", ProteinExpression, fixed = TRUE) & !grepl("[=|?]|\\[(.*?)\\]|\\((.*?)\\)|ext", ProteinExpression) & ProteinExpression != "p.0"]

clinvar = merge(clinvar, hgvs, by = "AlleleID")
clinvar[, c("AlleleID") := NULL]
clinvar = unique(clinvar, by = c("MutationID", "ProteinExpression", "GSExpression"))

# Download all human refseq protein sequences in fasta format
# http://www.ncbi.nlm.nih.gov/protein?term=(homo%20sapiens%5BOrganism%5D)%20AND%20srcdb_refseq%5BProperties%5D
refseq_human = read.fasta("database/refseq_human.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
refseq_human = refseq_human[grep("NP_", names(refseq_human))]
names(refseq_human) = str_extract(names(refseq_human), "[A-Z]P_\\d+")
clinvar = clinvar[refseq %in% names(refseq_human)]
clinvar_fasta = refseq_human[unique(clinvar$refseq)]
write.fasta(clinvar_fasta, names(clinvar_fasta), "database/clinvar/clinvar.fasta", as.string = TRUE)
# qsub -pe omp 16 database/clinvar/clinvar_pfam.sh

# map SNPs to PFAM domains
library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, names(AMINO_ACID_CODE), "[+|_|=|>]", "del", "dup", "fs", "ins", "Ter"), collapse = "|")

hgvsp2aapos = function(hgvsp, faa) {
  prot = unlist(strsplit(hgvsp, ":p.", fixed = TRUE))[1]
  pos = unlist(strsplit(hgvsp, ":p.", fixed = TRUE))[2]
  pos = unlist(strsplit(pos, aa))
  pos = setdiff(as.numeric(pos[nchar(pos) > 0]), NA)
  
  if (length(pos) == 0) NA
  else if (grepl("fs|Ter", hgvsp)) paste(pos[1], nchar(faa[[prot]]), sep = ":")
  else if (length(pos) == 1) pos
  else paste(pos[1], pos[2], sep = ":")
}

hgvsp2aapos_clinvar = sapply(unique(clinvar$ProteinExpression), hgvsp2aapos, faa = clinvar_fasta)
stopifnot(all(sapply(hgvsp2aapos_clinvar, function(x) class(eval(parse(text = x)))) %in% c("integer", "numeric")))
clinvar$aapos = as.character(hgvsp2aapos_clinvar[clinvar$ProteinExpression])

clinvar_pfam = fread("database/clinvar/clinvar_pfam.csv")
clinvar_pfam$pfamacc = str_extract(clinvar_pfam$pfamacc, "PF\\d{5}")
clinvar_pfams = lapply(split(clinvar_pfam[, .(pfamacc, pfamstart, pfamend)], clinvar_pfam$seqid), unique)

clinvar = clinvar[refseq %in% names(clinvar_pfams)]
aapos_clinvar = unique(paste(clinvar$refseq, clinvar$aapos, sep = ","))
aapos2pfam_clinvar = mclapply(aapos_clinvar, aapos2pfam, pfams = clinvar_pfams, mc.cores = detectCores())
names(aapos2pfam_clinvar) = aapos_clinvar
aapos2pfam_clinvar = unlist(aapos2pfam_clinvar)
clinvar$pfam = as.character(aapos2pfam_clinvar[paste(clinvar$refseq, clinvar$aapos, sep = ",")])
clinvar = unique(clinvar[nchar(pfam) > 0])

###############################################
### combine humsavar and clinvar into hcvar ###
###############################################
hcvar = unique(rbind(clinvar[, .(GeneSymbol, GeneID, MutationID, ProteinExpression, GSExpression, pfam, PhenotypeIDs)], humsavar[, .(GeneSymbol, GeneID, MutationID, ProteinExpression, GSExpression, pfam, PhenotypeIDs)]))
# due to codon degeneracy, multiple MutationIDs may produce the same effect on amino acid sequence.
# for each group of MutationIDs producing the same amino acid effect, pick the MutationID
# that is clinically significant, has the highest number of distinct GeneSymbol + ProteinChange, and is alphabetically first.
hcvar$ClinSigSimple = clinsig[match(hcvar$MutationID, MutationID), ClinSigSimple]
hcvar[, c("MutationIDFreq") := uniqueN(GSExpression), by = "MutationID"]
setorder(hcvar, -ClinSigSimple, -MutationIDFreq, MutationID, na.last = TRUE)
hcvar[, c("ClinSigSimple", "MutationIDFreq") := NULL]
hcvar = unique(hcvar[, setdiff(names(hcvar), "GSExpression") := lapply(.SD, function(x) paste(unique(unlist(strsplit(x, ";"))), collapse = ";")), by = "GSExpression"])
hcvar$MutationID = str_replace(hcvar$MutationID, "(?=[;])(.*?)$", "")
hcvar = unique(hcvar[, list(pfam = unlist(strsplit(pfam, ";"))), by = setdiff(names(hcvar), "pfam")])
hcvar = unique(hcvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";"))), by = setdiff(names(hcvar), "PhenotypeIDs")])

# map Pfam domains to clans, to address potential biases introduced by overcounting domains with very similar sequences, structures and HMM profiles.
clan = fread("database/3did_elm/clan_membership.txt", header = FALSE)
names(clan) = c("clan", "pfam")
hcvar$clan = clan[match(hcvar$pfam, pfam), clan]
hcvar[is.na(clan), c("clan") := .(pfam)]
hcvar$genepfam = paste(hcvar$GeneSymbol, hcvar$pfam, sep = "|")
hcvar$geneclan = paste(hcvar$GeneSymbol, hcvar$clan, sep = "|")
hcvar = hcvar[grepl("\\d+$", PhenotypeIDs)]
hcvar$PhenotypeSource = str_extract(hcvar$PhenotypeIDs, "^(.*?)(?=[:])")
setcolorder(hcvar, c("GeneSymbol", "GeneID", "MutationID", "ProteinExpression", "GSExpression", "pfam", "clan", "genepfam", "geneclan", "PhenotypeSource", "PhenotypeIDs"))

########################################################################################
### map PhenotypeIDs to Disease Ontology, MeSH, and Orphanet IDs and their top nodes ###
########################################################################################
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
doid2xref$xref = gsub("(?<=EFO)(.*?)(?=[:])", "", doid2xref$xref, perl = TRUE)
doid2xref$xref = gsub("EFO:", "EFO:EFO_", doid2xref$xref)
doid2xref$xref = gsub("HP:", "Human Phenotype Ontology:HP:", doid2xref$xref)
doid2xref$xref = gsub("MESH:", "MeSH:", doid2xref$xref)
doid2xref$xref = gsub("(?<=NCI)(.*?)(?=[:])", "", doid2xref$xref, perl = TRUE)
doid2xref$xref = gsub("OMM:", "OMIM:", doid2xref$xref)
doid2xref$xref = gsub("Orphanet:", "ORPHA:", doid2xref$xref)
doid2xref$xref = gsub("ORDO:", "ORPHA:", doid2xref$xref)
doid2xref$xref = gsub("(?<=SNOMEDCT)(.*?)(?=[:])", "", doid2xref$xref, perl = TRUE)
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

# manually assign Disease Ontology, MeSH and ORPHA IDs to unclassified cancer terms
oncoregex = "blastoma|cancer|carcino*|glioma|leukemia|leukaemia|lymphoma|melanoma|neoplas*|sarcoma|tumour|tumor"
hcvar[grepl(oncoregex, DiseaseName, ignore.case = TRUE) & nchar(DOIDtree) == 0, c("DOIDtree") := "DOID:14566"]
hcvar[grepl(oncoregex, DiseaseName, ignore.case = TRUE) & nchar(MeSHtree) == 0, c("MeSHtree") := "C04"]
hcvar[grepl(oncoregex, DiseaseName, ignore.case = TRUE) & nchar(ORPHAtop) == 0, c("ORPHAtop") := "ORPHA:250908"]

setorder(hcvar)
hcvar_full = unique(hcvar)
hcvar = unique(hcvar[MutationID %in% clinsig[ClinSigSimple == "1", MutationID]])

save.image("database/clinvar/hcvar.RData")
