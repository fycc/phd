#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(stringr)

#########################################################################
### get bacteria and virus taxonomy information from NCBI and UniProt ###
#########################################################################
taxcols = c("taxid", "kingdom", "taxname")

ncbibac = fread("ncbi_taxonomy_bacteria.csv", colClasses = "character")
ncbivir = fread("ncbi_taxonomy_virus.csv", colClasses = "character")
ncbitax = data.table(rbind(ncbibac, ncbivir), kingdom = rep(c("bacteria", "virus"), c(nrow(ncbibac), nrow(ncbivir))))
setorder(ncbitax)

uniprotbac = fread("speclist_bacteria.txt", header = FALSE, colClasses = "character")
setnames(uniprotbac, taxcols)
uniprotvir = fread("speclist_virus.txt", header = FALSE, colClasses = "character")
setnames(uniprotvir, taxcols)

bacvir = rbindlist(list(ncbitax[, taxcols, with = FALSE], uniprotbac[!taxid %in% ncbitax$taxid], uniprotvir[!taxid %in% ncbitax$taxid]), use.names = TRUE)
for(j in names(bacvir)) set(bacvir, j = j, value = gsub(";", "_", bacvir[[j]]))
setorder(bacvir)

taxmerged = fread("taxmerged.dmp", header = FALSE)
taxmerged = taxmerged[, sapply(names(taxmerged), function(x) all(is.integer(taxmerged[[x]]))), with = FALSE]
setnames(taxmerged, c("taxid_old", "taxid_new"))
for(j in names(taxmerged)) set(taxmerged, j = j, value = as.character(taxmerged[[j]]))

save(ncbitax, bacvir, taxmerged, file = "taxonomy.RData")

#############################
### parse IntAct PPI data ###
#############################
uniprot_regex = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\-(PRO_){0,1}\\d+){0,1}"

# sed -e "s/'/singlequote/g" -e 's/"/doublequote/g' < intact.txt > intact_quote_replaced.txt
intact = fread("/projectnb/xialab2/fyc/intact/intact_quote_replaced.txt", colClasses = "character")
setnames(intact, trimws(str_replace_all(names(intact), "#|\\((.*?)\\)", "")))
intact = intact[grepl("protein", `Type interactor A`) & grepl("protein", `Type interactor B`)]
intact[, setdiff(names(intact), c("ID interactor A", "ID interactor B", "Taxid interactor A", "Taxid interactor B", "Publication Identifier")) := NULL]
setcolorder(intact, c("ID interactor A", "ID interactor B", "Taxid interactor A", "Taxid interactor B", "Publication Identifier"))
setnames(intact, c("idA", "idB", "taxidA", "taxidB", "pmid"))
for(j in names(intact)) set(intact, j = j, value = gsub("singlequote", "'", intact[[j]]))
for(j in names(intact)) set(intact, j = j, value = gsub("doublequote", "", intact[[j]]))
intact[, c("idA", "idB") := lapply(.SD, toupper), .SDcols = c("idA", "idB")]
intact[, c("taxidA", "taxidB") := lapply(.SD, str_extract, pattern = "(?<=taxid:)(.*?)(?=[^0-9]+)"), .SDcols = c("taxidA", "taxidB")]
intact[taxidA %in% taxmerged$taxid_old, taxidA := taxmerged[match(taxidA, taxid_old), taxid_new]]
intact[taxidB %in% taxmerged$taxid_old, taxidB := taxmerged[match(taxidB, taxid_old), taxid_new]]
intact = intact[str_detect(idA, uniprot_regex) & str_detect(idB, uniprot_regex) & taxidA %in% c("9606", bacvir$taxid) & taxidB %in% c("9606", bacvir$taxid)]
intact[, c("idA", "idB") := lapply(.SD, str_extract, pattern = uniprot_regex), .SDcols = c("idA", "idB")]

# get canonical IDs (cids) from isoform IDs and remove self-interactions between different isoforms
intact[, c("cidA", "cidB") := lapply(.SD, str_replace, pattern = "(?=-)(.*?)$", replacement = ""), .SDcols = c("idA", "idB")]
intact = intact[cidA != cidB]

# to remove duplicates (e.g. A_B and B_A), rearrange protein A and protein B, such that protein A will always have a smaller ID than protein B
# for host-pathogen PPIs, make sure host protein ID is in front
colsAB = c(grep("A$", names(intact), value = TRUE), grep("B$", names(intact), value = TRUE))
colsBA = c(grep("B$", names(intact), value = TRUE), grep("A$", names(intact), value = TRUE))
intact[cidA > cidB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
intact[taxidA %in% bacvir$taxid & taxidB == "9606", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
intact[, c("taxnameA", "taxnameB") := lapply(.SD, function(x) ifelse(x == "9606", "Homo sapiens", bacvir[match(x, taxid), taxname])), .SDcols = c("taxidA", "taxidB")]
intact[, c("taxcatA", "taxcatB") := lapply(.SD, function(x) ifelse(x == "9606", "host", bacvir[match(x, taxid), kingdom])), .SDcols = c("taxidA", "taxidB")]
setorder(intact)
intact = unique(intact)
save(intact, file = "intact.RObject")

############################
### parse HPIDB PPI data ###
############################
hpidb = fread("hpidb2.mitab_plus.txt", colClasses = "character")
setnames(hpidb, trimws(str_replace_all(names(hpidb), "#|\\((.*?)\\)", "")))
hpidb[, setdiff(names(hpidb), c("protein_xref_1_unique", "protein_xref_2_unique", "protein_taxid_1", "protein_taxid_2", "pmid")) := NULL]
setcolorder(hpidb, c("protein_xref_1_unique", "protein_xref_2_unique", "protein_taxid_1", "protein_taxid_2", "pmid"))
setnames(hpidb, c("idA", "idB", "taxidA", "taxidB", "pmid"))
hpidb[, c("idA", "idB") := lapply(.SD, toupper), .SDcols = c("idA", "idB")]
hpidb[, c("taxidA", "taxidB") := lapply(.SD, str_extract, pattern = "(?<=taxid:)(.*?)(?=[^0-9]+)"), .SDcols = c("taxidA", "taxidB")]
hpidb[taxidA %in% taxmerged$taxid_old, taxidA := taxmerged[match(taxidA, taxid_old), taxid_new]]
hpidb[taxidB %in% taxmerged$taxid_old, taxidB := taxmerged[match(taxidB, taxid_old), taxid_new]]
hpidb = hpidb[str_detect(idA, uniprot_regex) & str_detect(idB, uniprot_regex) & taxidA %in% c("9606", bacvir$taxid) & taxidB %in% c("9606", bacvir$taxid)]
hpidb[, c("idA", "idB") := lapply(.SD, str_extract, pattern = uniprot_regex), .SDcols = c("idA", "idB")]

# get canonical IDs (cids) from isoform IDs and remove self-interactions between different isoforms
hpidb[, c("cidA", "cidB") := lapply(.SD, str_replace, pattern = "(?=-)(.*?)$", replacement = ""), .SDcols = c("idA", "idB")]
hpidb = hpidb[cidA != cidB]

# for host-pathogen PPIs, make sure host protein ID is in front
colsAB = c(grep("A$", names(hpidb), value = TRUE), grep("B$", names(hpidb), value = TRUE))
colsBA = c(grep("B$", names(hpidb), value = TRUE), grep("A$", names(hpidb), value = TRUE))
hpidb[taxidA %in% bacvir$taxid & taxidB == "9606", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
hpidb[, c("taxnameA", "taxnameB") := lapply(.SD, function(x) ifelse(x == "9606", "Homo sapiens", bacvir[match(x, taxid), taxname])), .SDcols = c("taxidA", "taxidB")]
hpidb[, c("taxcatA", "taxcatB") := lapply(.SD, function(x) ifelse(x == "9606", "host", bacvir[match(x, taxid), kingdom])), .SDcols = c("taxidA", "taxidB")]
setorder(hpidb)
hpidb = unique(hpidb)
save(hpidb, file = "hpidb.RObject")

################################
### combine IntAct and HPIDB ###
################################
ppi = rbindlist(list(intact, hpidb), use.names = TRUE)
ppi[grepl("-PRO_", idA), idA := cidA]
ppi[grepl("-PRO_", idB), idB := cidB]
ppi = ppi[str_detect(pmid, "pubmed") & !grepl("pubmed:20711500", pmid, fixed = TRUE)]
writeLines(ppi[, sort(unique(c(idA, idB, cidA, cidB)))], "ppi_id.txt")
save(ppi, file = "ppi.RObject")
