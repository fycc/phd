#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(stringr)

clinvar_vep = fread("clinvar_vep.tab", sep = "\t", skip = "Uploaded_variation\tLocation", colClasses = "character")
names(clinvar_vep) = gsub("#", "", names(clinvar_vep))
names(clinvar_vep) = sub("Uploaded_variation", "VariationID", names(clinvar_vep))
clinvar_vep[, pfam := lapply(str_extract_all(DOMAINS, "(?<=(?i)(Pfam_domain:))(PF\\d+)(?=[^0-9]|$)"), unique)]
clinvar_vep = clinvar_vep[grepl("(?i)missense", Consequence) & lengths(pfam) > 0 & USED_REF == GIVEN_REF & BIOTYPE == "protein_coding"]
clinvar_vep[, lva := paste(Location, VARIANT_CLASS, USED_REF, sep = "_")]
save(clinvar_vep, file = "clinvar_vep.RObject")

variation_allele = fread("variation_allele.txt", sep = "\t", skip = "VariationID\tType", colClasses = "character")
names(variation_allele) = gsub("#", "", names(variation_allele))
variation_allele = variation_allele[VariationID %in% clinvar_vep$VariationID]
setkeyv(variation_allele, c("AlleleID", "VariationID"))

clinvar = fread("clinvar.txt", colClasses = "character")
names(clinvar) = gsub("[^A-Z|a-z]+", "", names(clinvar))
names(clinvar) = sub("PhenotypeIDS", "PhenotypeIDs", names(clinvar))
names(clinvar) = sub("RSdbSNP", "rsID", names(clinvar))
clinvar[rsID != "-1", rsID := paste0("rs", rsID)]
clinvar = unique(clinvar[Assembly == "GRCh38" & AlleleID %in% variation_allele$AlleleID, .(AlleleID, OriginSimple, ClinSigSimple, GeneSymbol, PhenotypeIDs, PhenotypeList)])
clinvar[, VariationID := sapply(AlleleID, function(x) variation_allele[.(x), VariationID])]
clinvar = unique(clinvar[, list(VariationID = unlist(VariationID)), by = setdiff(names(clinvar), "VariationID")])
clinvar[, PhenotypeIDs := gsub("MedGen:", "UMLS:", PhenotypeIDs)]
clinvar[, PhenotypeIDs := gsub("Orphanet:ORPHA", "ORPHA:", PhenotypeIDs)]
clinvar[, PhenotypeIDs := gsub("\\bna\\b|\\bUMLS:CN169374\\b|\\bUMLS:CN517202\\b", "", PhenotypeIDs)]
clinvar[, PhenotypeList := gsub("(\\bnot provided\\b|\\bnot specified\\b)([,|;]|$)", "", PhenotypeList)]
clinvar[, PhenotypeList := str_replace(PhenotypeList, "[;|,](?=$)", "")]
clinvar = unique(clinvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, "[;|,]"))), by = setdiff(names(clinvar), "PhenotypeIDs")])
clinvar = clinvar[nchar(PhenotypeIDs) > 0]
cid_history = fread("ConceptID_history.txt")
cid_history[, `Previous ConceptID` := paste0("UMLS:", `Previous ConceptID`)]
cid_history[, `Current ConceptID` := paste0("UMLS:", `Current ConceptID`)]
clinvar = clinvar[!PhenotypeIDs %in% cid_history[`Current ConceptID` == "UMLS:No longer reported", `Previous ConceptID`]]
clinvar[PhenotypeIDs %in% cid_history$`Previous ConceptID`, PhenotypeIDs := cid_history[match(PhenotypeIDs, `Previous ConceptID`), `Current ConceptID`]]
clinvar = unique(clinvar[, PhenotypeIDs := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = setdiff(names(clinvar), "PhenotypeIDs")])

clinvar = merge(clinvar_vep[, .(VariationID, SYMBOL, pfam, lva)], clinvar[, -c("AlleleID"), with = FALSE], by = "VariationID")
clinvar[grepl("[;|,]", GeneSymbol), GeneSymbol := SYMBOL]
clinvar[, SYMBOL := NULL]
clinvar = unique(clinvar[, list(pfam = unlist(pfam)), by = setdiff(names(clinvar), "pfam")])
clinvar = unique(clinvar[, list(PhenotypeIDs = unlist(strsplit(PhenotypeIDs, ";"))), by = setdiff(names(clinvar), "PhenotypeIDs")])
clinvar[, ClinSigSimple := ifelse(any(ClinSigSimple == "1"), "1", "0"), by = c("lva", "PhenotypeIDs")]
clinvar = unique(clinvar)
save(clinvar, file = "clinvar.RObject")
