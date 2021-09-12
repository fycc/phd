#!/usr/bin/env Rscript

library(data.table)
library(parallel)
library(seqinr)
library(stringr)

clinvar = fread("database/clinvar/clinvar.txt", colClasses = "character")
names(clinvar) = gsub("[^A-Z|a-z]+", "", names(clinvar))
# exclude VarIDs already included in UniProtKB humsavar.txt
clinvar = clinvar[!grepl("UniProtKB", clinvar$OtherIDs, ignore.case = TRUE), ]

hgvs = fread("database/clinvar/hgvs4variation.txt", colClasses = "character")
names(hgvs)[grep("Symbol", names(hgvs), ignore.case = TRUE)] = "GeneSymbol"
hgvs = hgvs[hgvs$AlleleID %in% clinvar$AlleleID & grepl("^NP_", hgvs$ProteinExpression) & grepl("p.", hgvs$ProteinExpression, fixed = TRUE) & !grepl("[=|?]|\\[(.*?)\\]|\\((.*?)\\)|ext", hgvs$ProteinExpression) & hgvs$ProteinExpression != "p.0", ]
hgvs$ProteinExpression = gsub("[*]", "Ter", hgvs$ProteinExpression)
hgvs$ProteinChange = gsub("[*]", "Ter", hgvs$ProteinChange)
hgvs$refseq = as.character(sapply(hgvs$ProteinExpression, function(x) unlist(strsplit(x, "[.|:]"))[1]))
hgvs$PhenotypeIDs = clinvar$PhenotypeIDS[match(hgvs$AlleleID, clinvar$AlleleID)]
writeLines(unique(hgvs$refseq), "database/clinvar/clinvar_refseq.txt")

# Download all human refseq protein sequences in fasta format
# http://www.ncbi.nlm.nih.gov/protein?term=(homo%20sapiens%5BOrganism%5D)%20AND%20srcdb_refseq%5BProperties%5D
human_refseq = read.fasta("database/intact/human_refseq.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
human_refseq = human_refseq[grep("NP_", names(human_refseq))]
names(human_refseq) = as.character(sapply(names(human_refseq), function(x) grep("^NP_", unlist(strsplit(x, "[|]")), value = TRUE)))
names(human_refseq) = as.character(sapply(names(human_refseq), function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_fasta = human_refseq[unique(hgvs$refseq)]
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

clinvar_uniprot = sapply(unique(hgvs$refseq), function(x) clinvar_refseq$Entry[match(x, clinvar_refseq$refseq)])
clinvar_geneid = sapply(unique(hgvs$refseq), function(x) clinvar_refseq$`Cross-reference (GeneID)`[match(x, clinvar_refseq$refseq)])
clinvar_genesymbol = sapply(unique(hgvs$refseq), function(x) clinvar_refseq$`Gene names  (primary )`[match(x, clinvar_refseq$refseq)])

hgvs$uniprot = as.character(clinvar_uniprot[hgvs$refseq])
hgvs$GeneID = as.character(clinvar_geneid[hgvs$refseq])
hgvs$GeneID[is.na(hgvs$GeneID) | nchar(hgvs$GeneID) == 0] = clinvar$GeneID[match(hgvs$AlleleID[is.na(hgvs$GeneID) | nchar(hgvs$GeneID) == 0], clinvar$AlleleID)]
hgvs$GeneSymbol = as.character(clinvar_genesymbol[hgvs$refseq])
hgvs$GeneSymbol[is.na(hgvs$GeneSymbol) | nchar(hgvs$GeneSymbol) == 0] = as.character(sapply(hgvs$AlleleID[is.na(hgvs$GeneSymbol) | nchar(hgvs$GeneSymbol) == 0], function(x) {
  intersect(unlist(strsplit(clinvar$GeneSymbol[clinvar$AlleleID == x], ";")), unlist(strsplit(clinvar$Name[clinvar$AlleleID == x], "[()]"))[2])
}))

############################
# Map SNPs to PFAM domains #
############################

library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, names(AMINO_ACID_CODE), "[+|_|=]", "del", "dup", "fs", "ins", "Ter"), collapse = "|")

hgvsp2aapos = function(hgvsp) {
  pos = tail(unlist(strsplit(hgvsp, "p.", fixed = TRUE)), 1)
  pos = unlist(strsplit(pos, aa))
  pos = setdiff(as.numeric(unique(pos[nchar(pos) > 0])), NA)
  
  if (length(pos) == 0) NA
  else if (grepl("fs|Ter", hgvsp)) paste(pos[1], nchar(clinvar_fasta[[unlist(strsplit(hgvsp, "[.|:]"))[1]]]), sep = ":")
  else if (length(pos) == 1) pos
  else paste(pos[1], pos[2], sep = ":")
}

aapos = sapply(unique(hgvs$ProteinExpression), hgvsp2aapos)
stopifnot(all(sapply(aapos, function(x) class(eval(parse(text = x)))) %in% c("integer", "numeric")))
hgvs$aapos = as.character(aapos[hgvs$ProteinExpression])

clinvar_pfam = fread("database/clinvar/clinvar_pfam.csv")
clinvar_pfam$pfamacc = as.character(sapply(clinvar_pfam$pfamacc, function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_pfams = lapply(split(clinvar_pfam[, list(pfamacc, pfamstart, pfamend)], clinvar_pfam$seqid), unique)

hgvs = hgvs[hgvs$refseq %in% names(clinvar_pfams), ]
refseq_aapos = unique(paste(hgvs$refseq, hgvs$aapos, sep = ";"))
pfams_hgvs = mclapply(refseq_aapos, function(x) {
  
  refseq = unlist(strsplit(x, ";"))[1]
  aapos = unlist(strsplit(x, ";"))[2]
  
  pfam = clinvar_pfams[[refseq]]
  i = apply(pfam, 1, function(p) any(eval(parse(text = aapos)) %in% p["pfamstart"]:p["pfamend"]))
  ifelse(sum(i) == 0, NA, paste(unique(pfam$pfamacc[i]), collapse = ";"))
  
}, mc.cores = detectCores())
names(pfams_hgvs) = refseq_aapos
pfams_hgvs = unlist(pfams_hgvs)
hgvs$pfam = as.character(pfams_hgvs[paste(hgvs$refseq, hgvs$aapos, sep = ";")])
hgvs = unique(hgvs[!is.na(hgvs$pfam), ])

###########################################
# Combine humsavar and clinvar into hcvar #
###########################################

clin = unique(hgvs[, list(GeneSymbol, AlleleID, pfam, PhenotypeIDs)])
clin$PhenotypeIDs = gsub("MedGen:", "UMLS:", clin$PhenotypeIDs)
clin$PhenotypeIDs = gsub("Orphanet:ORPHA", "ORPHA:", clin$PhenotypeIDs)
clin = unique(clin[, list(pfam = unlist(strsplit(pfam, ";"))), by = setdiff(names(clin), "pfam")])
clin = unique(clin[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";|,"))), by = setdiff(names(clin), "PhenotypeIDs")])
clin = unique(clin[nchar(clin$PhenotypeIDs) > 0, ])

load("database/clinvar/humsavar.RData")
hcvar = unique(rbind(humsavar[, list(GeneSymbol, AlleleID, pfam, PhenotypeIDs)], clin))
rm(clin)

# map Pfam domains to clans, to address potential biases introduced by overcounting domains with very similar sequences, structures and HMM profiles.
clan = fread("database/3did_elm/clan_membership.txt", header = FALSE)
names(clan) = c("clan", "pfam")
hcvar$clan = clan$clan[match(hcvar$pfam, clan$pfam)]
hcvar$clan[is.na(hcvar$clan)] = hcvar$pfam[is.na(hcvar$clan)]
hcvar$genepfam = paste(hcvar$GeneSymbol, hcvar$pfam, sep = "|")
hcvar$geneclan = paste(hcvar$GeneSymbol, hcvar$clan, sep = "|")
hcvar = unique(hcvar[grepl("[0-9]$", hcvar$PhenotypeIDs) & !hcvar$PhenotypeIDs %in% c("UMLS:CN169374", "UMLS:CN221809"), list(GeneSymbol, AlleleID, pfam, clan, genepfam, geneclan, PhenotypeIDs)])
hcvar$PhenotypeSource = as.character(sapply(hcvar$PhenotypeIDs, function(x) unlist(strsplit(x, ":"))[1]))

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
mgconso$SDUI = gsub("Orphanet_", "", mgconso$SDUI)
mgconso$SDUI[mgconso$SAB == "NCI"] = mgconso$SCUI[mgconso$SAB == "NCI"]
mgconso$SDUI[mgconso$SAB == "SNOMED CT"] = mgconso$SCUI[mgconso$SAB == "SNOMED CT"]
mgconso$SDUI = paste(mgconso$SAB, mgconso$SDUI, sep = ":")
mgconso = mgconso[!mgconso$CUI %in% c("UMLS:CN169374", "UMLS:CN221809") & !grepl("obsolete", mgconso$STR, ignore.case = TRUE), list(CUI, SDUI, SAB, STR)]

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
doid2parent = fread("database/clinvar/doid2parent.csv", colClasses = "character")

mesh = fread("database/clinvar/mesh.csv", colClasses = "character")
meshtree = fread("database/clinvar/MeshTreeHierarchyWithScopeNotes.txt", colClasses = "character")

orpha2xref = fread("database/clinvar/orpha2xref.csv", colClasses = "character")
orphalinear = fread("database/clinvar/orphalinear.csv", colClasses = "character")

# functions to map PhenotypeIDs to lowest level nodes of (most disease-specific) Disease Ontology, MeSH and ORPHA.
xref2doid = function(x) {
  if (grepl("^UMLS:", x)) paste(sort(union(doid2xref$DOID[doid2xref$xref == x], doid2xref$DOID[doid2xref$xref %in% mgconso$SDUI[mgconso$CUI == x]])), collapse = ";")
  else paste(sort(union(doid2xref$DOID[doid2xref$xref == x], doid2xref$DOID[doid2xref$xref %in% mgconso$CUI[mgconso$SDUI == x]])), collapse = ";")
}

xref2mesh = function(x) {
  if (grepl("^MeSH:", x)) x
  else if (grepl("^UMLS:", x)) paste(sort(unique(mgconso$SDUI[mgconso$SAB == "MeSH" & mgconso$CUI == x])), collapse = ";")
  else paste(sort(unique(mgconso$SDUI[mgconso$SAB == "MeSH" & mgconso$CUI %in% mgconso$CUI[mgconso$SDUI == x]])), collapse = ";")
}

xref2orpha = function(x) {
  if (grepl("^ORPHA:", x)) x
  else if (grepl("^UMLS:", x)) {
    paste(sort(Reduce(union, list(orpha2xref$OrphaNumber[orpha2xref$ExternalID == x],
                                  orpha2xref$OrphaNumber[orpha2xref$ExternalID %in% mgconso$SDUI[mgconso$CUI == x]], 
                                  mgconso$SDUI[mgconso$SAB == "ORPHA" & mgconso$CUI == x]))), collapse = ";")
  }
  else {
    paste(sort(Reduce(union, list(orpha2xref$OrphaNumber[orpha2xref$ExternalID == x],
                                  orpha2xref$OrphaNumber[orpha2xref$ExternalID %in% mgconso$CUI[mgconso$SDUI == x]],
                                  mgconso$SDUI[mgconso$SAB == "ORPHA" & mgconso$CUI %in% mgconso$CUI[mgconso$SDUI == x]]))), collapse = ";")
  }
}

# functions to map PhenotypeIDs to parent terms and top levels of Disease Ontology, MeSH and ORPHA.
doid2top = function(x) {y = NULL; while (any(x %in% doid2parent$DOID)) {x = doid2parent$ParentID[doid2parent$DOID %in% x]; y = c(y, x)}; paste(setdiff(y, "DOID:4"), collapse = ";")}
mesh2top = function(x) paste(sort(unique(mesh$TreeTop[mesh$DescriptorUI == gsub("MeSH:", "", x)])), collapse = ";")
orpha2top = function(x) paste(sort(unique(orphalinear$ORDOID[orphalinear$OrphaNumber == x])), collapse = ";")

# mapping PhenotypeIDs to Disease Ontology, MeSH and ORPHA IDs.
pheno2doid = sapply(unique(hcvar$PhenotypeIDs), xref2doid)
pheno2mesh = sapply(unique(hcvar$PhenotypeIDs), xref2mesh)
pheno2orpha = sapply(unique(hcvar$PhenotypeIDs), xref2orpha)

pheno2doidtop = sapply(unique(hcvar$PhenotypeIDs), function(x) {
  if (nchar(pheno2doid[[x]]) == 0) return("")
  else {
    doid = unlist(strsplit(pheno2doid[[x]], ";"))
    doidtop = as.character(sapply(doid, doid2top))
    paste(unique(unlist(strsplit(doidtop, ";"))), collapse = ";")
  }
})

pheno2meshtop = sapply(unique(hcvar$PhenotypeIDs), function(x) {
  if (nchar(pheno2mesh[[x]]) == 0) return("")
  else {
    mesh = unlist(strsplit(pheno2mesh[[x]], ";"))
    meshtop = as.character(sapply(mesh, mesh2top))
    paste(sort(unique(unlist(strsplit(meshtop, ";")))), collapse = ";")
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

hcvar$DOIDtop = as.character(pheno2doidtop[hcvar$PhenotypeIDs])
hcvar$MeSHtop = as.character(pheno2meshtop[hcvar$PhenotypeIDs])
hcvar$ORPHAtop = as.character(pheno2orphatop[hcvar$PhenotypeIDs])

save.image("database/clinvar/clinvar.RData")
