#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(stringr)

setwd("~/database/clinvar")
load("humsavar.RObject")
load("clinvar.RObject")
load("clinvar_vep.RObject")

# combine humsavar, humsavar_rsid_vep and humsavar_hgvsp_vep
humsavar_rsid_vep = fread("humsavar_rsid_vep.tab", sep = "\t", skip = "Uploaded_variation\tLocation", colClasses = "character")
names(humsavar_rsid_vep) = sub("#Uploaded_variation", "rsID", names(humsavar_rsid_vep))
humsavar_rsid_vep[, pfam := lapply(str_extract_all(DOMAINS, "(?<=(?i)(Pfam_domain:))(PF\\d+)(?=[^0-9]|$)"), unique)]
humsavar_rsid_vep = humsavar_rsid_vep[grepl("(?i)missense", Consequence) & lengths(pfam) > 0 & USED_REF == GIVEN_REF & BIOTYPE == "protein_coding"]
humsavar_rsid_vep[, lva := paste(Location, VARIANT_CLASS, USED_REF, sep = "_")]

humsavar_hgvsp_vep = fread("humsavar_hgvsp_vep_all.tab", sep = "\t", colClasses = "character")
humsavar_hgvsp_vep[, pfam := lapply(str_extract_all(DOMAINS, "(?<=(?i)(Pfam_domain:))(PF\\d+)(?=[^0-9]|$)"), unique)]
humsavar_hgvsp_vep = humsavar_hgvsp_vep[grepl("(?i)missense", Consequence) & lengths(pfam) > 0 & USED_REF == GIVEN_REF & BIOTYPE == "protein_coding"]
humsavar_hgvsp_vep[, lva := paste(Location, VARIANT_CLASS, USED_REF, sep = "_")]
save(humsavar_rsid_vep, humsavar_hgvsp_vep, file = "humsavar_vep.RData")

hcvar_vep = rbindlist(list(clinvar_vep[, -c("VariationID"), with = FALSE], humsavar_rsid_vep[, -c("rsID"), with = FALSE], humsavar_hgvsp_vep[, -c("hgvsp"), with = FALSE]))
save(hcvar_vep, file = "hcvar_vep.RObject")

humsavar_rsid = merge(humsavar_rsid_vep[, .(rsID, pfam, lva, SOMATIC, CLIN_SIG)], humsavar[!is.na(rsID), .(rsID, GeneSymbol, PhenotypeIDs)], by = "rsID")
humsavar_hgvsp = merge(humsavar_hgvsp_vep[, .(hgvsp, pfam, lva, SOMATIC, CLIN_SIG)], humsavar[is.na(rsID), .(hgvsp, GeneSymbol, PhenotypeIDs)], by = "hgvsp")
hsvar = rbind(humsavar_rsid[, -c("rsID"), with = FALSE], humsavar_hgvsp[, -c("hgvsp"), with = FALSE])
hsvar[SOMATIC == "-", OriginSimple := "unknown"]
hsvar[SOMATIC != "-" & !grepl("1", SOMATIC), OriginSimple := "germline"]
hsvar[SOMATIC != "-" & grepl("1", SOMATIC), OriginSimple := "somatic"]
hsvar[SOMATIC != "-" & grepl("0", SOMATIC) & grepl("1", SOMATIC), OriginSimple := "germline/somatic"]
hsvar[, ClinSigSimple := ifelse(any(grepl("(?i)(pathogenic)", CLIN_SIG)), "1", "0"), by = c("lva", "PhenotypeIDs")]
hsvar = unique(hsvar[, list(pfam = unlist(pfam)), by = setdiff(names(hsvar), "pfam")])

hcvar = unique(rbind(clinvar[, .(lva, ClinSigSimple, GeneSymbol, pfam, PhenotypeIDs)], hsvar[, .(lva, ClinSigSimple, GeneSymbol, pfam, PhenotypeIDs)]))

# map Pfam domains to clans, to address potential biases introduced by overcounting domains with very similar sequences, structures and HMM profiles.
clan = fread("~/database/3did_elm/clan_membership.txt", header = FALSE)
names(clan) = c("clan", "pfam")
hcvar$clan = clan[match(hcvar$pfam, pfam), clan]
hcvar[is.na(clan), c("clan") := .(pfam)]
hcvar$genepfam = paste(hcvar$GeneSymbol, hcvar$pfam, sep = "|")
hcvar$geneclan = paste(hcvar$GeneSymbol, hcvar$clan, sep = "|")
hcvar = hcvar[grepl("\\d+$", PhenotypeIDs)]
hcvar$PhenotypeSource = str_extract(hcvar$PhenotypeIDs, "^(.*?)(?=[:])")
setcolorder(hcvar, c("lva", "ClinSigSimple", "GeneSymbol", "pfam", "clan", "genepfam", "geneclan", "PhenotypeSource", "PhenotypeIDs"))

########################################################################################
### map PhenotypeIDs to Disease Ontology, MeSH, and Orphanet IDs and their top nodes ###
########################################################################################
# import mapping files
mgconso = fread("MGCONSO.csv", colClasses = "character")
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
mgconso = unique(mgconso[!grepl("obsolete", STR, ignore.case = TRUE), .(CUI, SDUI, SAB, STR)])

doid2xref = fread("doid2xref.csv", colClasses = "character")
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
doid2parent = fread("doid2parent.csv", colClasses = "character")

mesh = fread("mesh.csv", colClasses = "character")
mesh$DescriptorUI = paste0("MeSH:", mesh$DescriptorUI)

orpha = fread("orpha.csv", colClasses = "character")
orpha2xref = fread("orpha2xref.csv", colClasses = "character")
orphalinear = fread("orphalinear.csv", colClasses = "character")

disease_names = fread("disease_names.txt", colClasses = "character")
names(disease_names)[1] = "DiseaseName"
disease_names$ConceptID = paste0("UMLS:", disease_names$ConceptID)
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
hcvar[grepl(oncoregex, DiseaseName, ignore.case = TRUE) & nchar(DOIDtree) == 0, DOIDtree := "DOID:162"]
hcvar[grepl(oncoregex, DiseaseName, ignore.case = TRUE) & nchar(MeSHtree) == 0, MeSHtree := "C04.588"]
hcvar[grepl(oncoregex, DiseaseName, ignore.case = TRUE) & nchar(ORPHAtop) == 0, ORPHAtop := "ORPHA:250908"]

setorder(hcvar)
hcvar_full = unique(hcvar)
hcvar = unique(hcvar[ClinSigSimple == "1"])

save(hcvar_full, hcvar, hsvar, file = "hcvar.RData")
