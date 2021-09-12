#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(parallel)
library(seqinr)
library(stringr)

# process clinvar_vep
clinvar_vep = fread("database/clinvar/clinvar_vep.tab", sep = "\t", skip = "Uploaded_variation\tLocation", colClasses = "character")
names(clinvar_vep) = gsub("#", "", names(clinvar_vep))
clinvar_vep[, pfam := sapply(DOMAINS, function(x) unlist(str_extract_all(x, "(?<=(?i)(Pfam_domain:))(PF\\d+)(?=[^0-9]|$)")))]
clinvar_vep = clinvar_vep[lengths(pfam) > 0 & USED_REF == GIVEN_REF]
clinvar_vep[, lva := paste(Location, VARIANT_CLASS, USED_REF, sep = "_")]
setkeyv(clinvar_vep, "Uploaded_variation")

variation_allele = fread("database/clinvar/variation_allele.txt", sep = "\t", skip = "VariationID\tType", colClasses = "character")
names(variation_allele) = gsub("#", "", names(variation_allele))
variation_allele = variation_allele[VariationID %in% clinvar_vep$Uploaded_variation]
setkeyv(variation_allele, c("AlleleID", "VariationID"))

clinvar = fread("database/clinvar/clinvar.txt", colClasses = "character")
names(clinvar) = gsub("[^A-Z|a-z]+", "", names(clinvar))
names(clinvar) = gsub("PhenotypeIDS", "PhenotypeIDs", names(clinvar))
names(clinvar) = gsub("RSdbSNP", "rsID", names(clinvar))
clinvar[rsID != "-1", rsID := paste0("rs", rsID)]
clinvar = clinvar[Assembly == "GRCh38" & AlleleID %in% variation_allele$AlleleID]
clinvar[, pfam := sapply(AlleleID, function(x) clinvar_vep[.(variation_allele[.(x), VariationID]), pfam])]
clinvar[, pfam := sapply(pfam, function(x) paste(unique(unlist(x)), collapse = ";"))]
clinvar[, lva := sapply(AlleleID, function(x) clinvar_vep[.(variation_allele[.(x), VariationID]), lva])]
clinvar[, lva := sapply(lva, function(x) paste(unique(unlist(x)), collapse = ";"))]
clinvar = unique(clinvar[, .(GeneSymbol, rsID, lva, pfam, PhenotypeIDs, ClinSigSimple)])

clinvar$PhenotypeIDs = gsub("MedGen:", "UMLS:", clinvar$PhenotypeIDs)
clinvar$PhenotypeIDs = gsub("Orphanet:ORPHA", "ORPHA:", clinvar$PhenotypeIDs)
clinvar = unique(clinvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";|,"))), by = setdiff(names(clinvar), "PhenotypeIDs")])
clinvar = clinvar[nchar(PhenotypeIDs) > 1 & PhenotypeIDs != "na"]
cid_history = fread("database/clinvar/ConceptID_history.txt")
cid_history$`Previous ConceptID` = paste0("UMLS:", cid_history$`Previous ConceptID`)
cid_history$`Current ConceptID` = paste0("UMLS:", cid_history$`Current ConceptID`)
clinvar = clinvar[!PhenotypeIDs %in% cid_history[`Current ConceptID` == "UMLS:No longer reported", `Previous ConceptID`]]
clinvar[PhenotypeIDs %in% cid_history$`Previous ConceptID`, c("PhenotypeIDs") := cid_history[match(PhenotypeIDs, `Previous ConceptID`), `Current ConceptID`]]
clinvar = unique(clinvar[, c("PhenotypeIDs") := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = setdiff(names(clinvar), "PhenotypeIDs")])
clinvar = unique(clinvar[, list(pfam = unlist(strsplit(pfam, ";"))), by = setdiff(names(clinvar), "pfam")])
clinvar = unique(clinvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";"))), by = setdiff(names(clinvar), "PhenotypeIDs")])

