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

goslim = getOBOCollection("~/database/3did_elm/goslim_metagenomics.obo")

parents = as.data.table(pfam2go[, get_parent_nodes(unique(GO))])
setkeyv(parents, "child_go_id")

bp = as.data.table(get_child_nodes("GO:0008150"))
pfam2goslim_bp = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim@ids, bp$child_go_id), "GO:0008150"), child_go_id]]
pfam2goslim_bp[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim@ids, bp$child_go_id), "GO:0008150")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_bp = unique(pfam2goslim_bp[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
pfam2goslim_bp[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
pfam2goslim_bp = unique(pfam2goslim_bp)
setkeyv(pfam2goslim_bp, "PFAM")

cc = as.data.table(get_child_nodes("GO:0005575"))
pfam2goslim_cc = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim@ids, cc$child_go_id), "GO:0005575"), child_go_id]]
pfam2goslim_cc[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim@ids, cc$child_go_id), "GO:0005575")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_cc = unique(pfam2goslim_cc[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
pfam2goslim_cc[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
pfam2goslim_cc = unique(pfam2goslim_cc)
setkeyv(pfam2goslim_cc, "PFAM")

mf = as.data.table(get_child_nodes("GO:0003674"))
pfam2goslim_mf = pfam2go[GO %in% parents[parent_go_id %in% setdiff(intersect(goslim@ids, mf$child_go_id), "GO:0003674"), child_go_id]]
pfam2goslim_mf[, parentGO := sapply(GO, function(x) parents[parent_go_id %in% setdiff(intersect(goslim@ids, mf$child_go_id), "GO:0003674")][.(x), paste(parent_go_id, collapse = ";")])]
pfam2goslim_mf = unique(pfam2goslim_mf[, list(parentGO = unlist(strsplit(parentGO, ";"))), by = "PFAM"])
pfam2goslim_mf[, parentGO := paste(parentGO, collapse = ";"), by = "PFAM"]
pfam2goslim_mf = unique(pfam2goslim_mf)
setkeyv(pfam2goslim_mf, "PFAM")

save(pfam2goslim_bp, pfam2goslim_cc, pfam2goslim_mf, file = "~/database/3did_elm/pfam2goslim.RData")
