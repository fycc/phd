#!/usr/bin/env Rscript

library(data.table)
library(GSEABase)
library(GO.db)
library(stringr)

# get GO annotation for Pfam domains
# from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/pfam2go
pfam2go_go = grep("^Pfam:", readLines("~/database/3did_elm/pfam2go.txt"), value = TRUE)
pfam2go_go = unique(data.table(PFAM = str_extract(pfam2go_go, "PF\\d{5}"), GO = str_extract(pfam2go_go, "GO:\\d{7}")))

# from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/gene_ontology.txt.gz
pfam2go_pfam = fread("~/database/3did_elm/gene_ontology.txt", header = FALSE)
setnames(pfam2go_pfam, c("PFAM", "GO", "GODescription", "GOCategory"))

pfam2go = unique(rbind(pfam2go_go, pfam2go_pfam[, .(PFAM, GO)]))

goslim_generic = getOBOCollection("~/database/3did_elm/goslim_generic.obo")
goslim_metagenomics = getOBOCollection("~/database/3did_elm/goslim_metagenomics.obo")
goslim_pir = getOBOCollection("~/database/3did_elm/goslim_pir.obo")

GOOFFSPRING = c(as.list(GOBPOFFSPRING), as.list(GOCCOFFSPRING), as.list(GOMFOFFSPRING))

goslimgeneric2pfam = lapply(setdiff(intersect(goslim_generic@ids, names(GOOFFSPRING)), c("GO:0008150", "GO:0005575", "GO:0003674")), function(x) {
  goslim_term = GOTERM[[x]]@Term
  goslim_children = union(x, GOOFFSPRING[[x]])
  pfam = pfam2go[GO %in% goslim_children, unique(PFAM)]
  list(GOSLIM_ID = rep(x, length(pfam)), GOSLIM_TERM = rep(goslim_term, length(pfam)), PFAM = pfam)
})
goslimgeneric2pfam = rbindlist(goslimgeneric2pfam)

goslimmetagenomics2pfam = lapply(setdiff(intersect(goslim_metagenomics@ids, names(GOOFFSPRING)), c("GO:0008150", "GO:0005575", "GO:0003674")), function(x) {
  goslim_term = GOTERM[[x]]@Term
  goslim_children = union(x, GOOFFSPRING[[x]])
  pfam = pfam2go[GO %in% goslim_children, unique(PFAM)]
  list(GOSLIM_ID = rep(x, length(pfam)), GOSLIM_TERM = rep(goslim_term, length(pfam)), PFAM = pfam)
})
goslimmetagenomics2pfam = rbindlist(goslimmetagenomics2pfam)

goslimpir2pfam = lapply(setdiff(intersect(goslim_pir@ids, names(GOOFFSPRING)), c("GO:0008150", "GO:0005575", "GO:0003674")), function(x) {
  goslim_term = GOTERM[[x]]@Term
  goslim_children = union(x, GOOFFSPRING[[x]])
  pfam = pfam2go[GO %in% goslim_children, unique(PFAM)]
  list(GOSLIM_ID = rep(x, length(pfam)), GOSLIM_TERM = rep(goslim_term, length(pfam)), PFAM = pfam)
})
goslimpir2pfam = rbindlist(goslimpir2pfam)

save(pfam2go, goslimgeneric2pfam, goslimmetagenomics2pfam, goslimpir2pfam, file = "~/database/3did_elm/pfam2goslim.RData")
