#!/usr/bin/env Rscript

library(data.table)
library(GOfuncR)
library(GSEABase)
library(stringr)

# get GO annotation for Pfam domains
# from http://current.geneontology.org/ontology/external2go/pfam2go
pfam2go_go = grep("^Pfam:", readLines("~/database/3did_elm/pfam2go.txt"), value = TRUE)
pfam2go_go = unique(data.table(PFAM = str_extract(pfam2go_go, "PF\\d{5}"), GO = str_extract(pfam2go_go, "GO:\\d{7}")))

# from ftp://ftp.ebi.ac.uk/pub/databases/interpro/current/interpro.xml.gz
ipr2pf = fread("~/database/3did_elm/ipr2pf.csv")

# from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/gene_ontology.txt.gz
pfam2go_pfam = fread("~/database/3did_elm/gene_ontology.txt", header = FALSE)
setnames(pfam2go_pfam, c("PFAM", "GO", "GODescription", "GOCategory"))

pfam2go = unique(rbindlist(list(pfam2go_go, ipr2pf[xrefdb == "PFAM", .(xref, GO)], pfam2go_pfam[, .(PFAM, GO)]), use.names = FALSE))

goslim_metagenomics = getOBOCollection("~/database/3did_elm/goslim_metagenomics.obo")
goslim_generic = getOBOCollection("~/database/3did_elm/goslim_generic.obo")

parents = as.data.table(pfam2go[, get_parent_nodes(unique(GO))])
setkeyv(parents, "child_go_id")

bp = as.data.table(get_child_nodes("GO:0008150"))
pfam2goslim_metagenomics_bp = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim_metagenomics@ids, bp$child_go_id), "GO:0008150"), child_go_id]]
pfam2goslim_metagenomics_bp[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim_metagenomics@ids, bp$child_go_id), "GO:0008150")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_metagenomics_bp = unique(pfam2goslim_metagenomics_bp[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
# pfam2goslim_metagenomics_bp[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
# pfam2goslim_metagenomics_bp = unique(pfam2goslim_metagenomics_bp)
# setkeyv(pfam2goslim_metagenomics_bp, "PFAM")
pfam2goslim_metagenomics_bp[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

pfam2goslim_generic_bp = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim_generic@ids, bp$child_go_id), "GO:0008150"), child_go_id]]
pfam2goslim_generic_bp[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim_generic@ids, bp$child_go_id), "GO:0008150")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_generic_bp = unique(pfam2goslim_generic_bp[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
# pfam2goslim_generic_bp[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
# pfam2goslim_generic_bp = unique(pfam2goslim_generic_bp)
# setkeyv(pfam2goslim_generic_bp, "PFAM")
pfam2goslim_generic_bp[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

cc = as.data.table(get_child_nodes("GO:0005575"))
pfam2goslim_metagenomics_cc = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim_metagenomics@ids, cc$child_go_id), "GO:0005575"), child_go_id]]
pfam2goslim_metagenomics_cc[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim_metagenomics@ids, cc$child_go_id), "GO:0005575")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_metagenomics_cc = unique(pfam2goslim_metagenomics_cc[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
# pfam2goslim_metagenomics_cc[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
# pfam2goslim_metagenomics_cc = unique(pfam2goslim_metagenomics_cc)
# setkeyv(pfam2goslim_metagenomics_cc, "PFAM")
pfam2goslim_metagenomics_cc[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

pfam2goslim_generic_cc = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim_generic@ids, cc$child_go_id), "GO:0005575"), child_go_id]]
pfam2goslim_generic_cc[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim_generic@ids, cc$child_go_id), "GO:0005575")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_generic_cc = unique(pfam2goslim_generic_cc[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
# pfam2goslim_generic_cc[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
# pfam2goslim_generic_cc = unique(pfam2goslim_generic_cc)
# setkeyv(pfam2goslim_generic_cc, "PFAM")
pfam2goslim_generic_cc[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

mf = as.data.table(get_child_nodes("GO:0003674"))
pfam2goslim_metagenomics_mf = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim_metagenomics@ids, mf$child_go_id), "GO:0003674"), child_go_id]]
pfam2goslim_metagenomics_mf[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim_metagenomics@ids, mf$child_go_id), "GO:0003674")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_metagenomics_mf = unique(pfam2goslim_metagenomics_mf[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
# pfam2goslim_metagenomics_mf[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
# pfam2goslim_metagenomics_mf = unique(pfam2goslim_metagenomics_mf)
# setkeyv(pfam2goslim_metagenomics_mf, "PFAM")
pfam2goslim_metagenomics_mf[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

pfam2goslim_generic_mf = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim_generic@ids, mf$child_go_id), "GO:0003674"), child_go_id]]
pfam2goslim_generic_mf[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim_generic@ids, mf$child_go_id), "GO:0003674")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_generic_mf = unique(pfam2goslim_generic_mf[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
# pfam2goslim_generic_mf[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
# pfam2goslim_generic_mf = unique(pfam2goslim_generic_mf)
# setkeyv(pfam2goslim_generic_mf, "PFAM")
pfam2goslim_generic_mf[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

# PFAM domains that map to the GO term interspecies interaction (ISI)
pfam2go_isi = pfam2go[GO %in% parents[parent_go_id == "GO:0044419", child_go_id]]
pfam2go_isi[, parentGO := sapply(GO, function(x) parents[.(x), paste(parent_go_id, collapse = ";")])]
pfam2go_isi = unique(pfam2go_isi[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
pfam2go_isi[, parentGOTerm := sapply(parentGO, function(x) as.character(getGOTerm(x)))]

save(pfam2goslim_metagenomics_bp, pfam2goslim_metagenomics_cc, pfam2goslim_metagenomics_mf, pfam2goslim_generic_bp, pfam2goslim_generic_cc, pfam2goslim_generic_mf, pfam2go_isi, file = "~/database/3did_elm/pfam2goslim.RData")
