#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=12000M
#SBATCH --job-name=ppi_parse
#SBATCH --output=%x-%j.out

options(stringsAsFactors = FALSE)
library(data.table)
library(stringr)

#########################################################################
### get bacteria and virus taxonomy information from NCBI and UniProt ###
#########################################################################
archaea_taxid = readLines("~/database/taxonomy/archaea_taxid.txt")
animal_taxid = readLines("~/database/taxonomy/animal_taxid.txt")
bacteria_taxid = readLines("~/database/taxonomy/bacteria_taxid.txt")
fungi_taxid = readLines("~/database/taxonomy/fungi_taxid.txt")
plant_taxid = readLines("~/database/taxonomy/plant_taxid.txt")
virus_taxid = readLines("~/database/taxonomy/virus_taxid.txt")

taxnames = fread("~/database/taxonomy/names.tab", header = FALSE, colClasses = "character")
setnames(taxnames, c("taxid", "taxname"))
setkeyv(taxnames, "taxid")

taxmerged = fread("~/database/taxonomy/merged.dmp", header = FALSE, select = c(1, 3), colClasses = "character")
setnames(taxmerged, c("taxid_old", "taxid_new"))

#############################
### parse IntAct PPI data ###
#############################
uniprot_regex = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\-(PRO_){0,1}\\d+){0,1}"

# sed -e "s/'/singlequote/g" -e 's/"/doublequote/g' < intact.txt > intact_quote_replaced.txt
intact = fread("~/database/intact/intact_quote_replaced.txt", colClasses = "character")
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
intact = intact[str_detect(idA, uniprot_regex) & str_detect(idB, uniprot_regex) & taxidA %in% c(animal_taxid, fungi_taxid, plant_taxid, archaea_taxid, bacteria_taxid, virus_taxid) & taxidB %in% c(animal_taxid, fungi_taxid, plant_taxid, archaea_taxid, bacteria_taxid, virus_taxid)]
intact[, c("idA", "idB") := lapply(.SD, str_extract, pattern = uniprot_regex), .SDcols = c("idA", "idB")]

# get canonical IDs (cids) from isoform IDs and remove self-interactions between different isoforms
intact[, c("cidA", "cidB") := lapply(.SD, str_replace, pattern = "(?=-)(.*?)$", replacement = ""), .SDcols = c("idA", "idB")]

# to remove duplicates (e.g. A_B and B_A), rearrange protein A and protein B, such that protein A will always have a smaller ID than protein B
# for host-pathogen PPIs, make sure host protein ID is in front
colsAB = c(grep("A$", names(intact), value = TRUE), grep("B$", names(intact), value = TRUE))
colsBA = c(grep("B$", names(intact), value = TRUE), grep("A$", names(intact), value = TRUE))
intact[cidA > cidB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
intact[taxidA %in% c(archaea_taxid, bacteria_taxid, virus_taxid, fungi_taxid) & taxidB %in% c(animal_taxid, plant_taxid), c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
intact[taxidA != "9606" & taxidB == "9606", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
intact[, c("taxnameA", "taxnameB") := lapply(.SD, function(x) taxnames[.(x), taxname]), .SDcols = c("taxidA", "taxidB")]
intact[taxidA %in% archaea_taxid, taxcatA := "archaea"]
intact[taxidA %in% animal_taxid, taxcatA := "animal"]
intact[taxidA %in% bacteria_taxid, taxcatA := "bacteria"]
intact[taxidA %in% fungi_taxid, taxcatA := "fungi"]
intact[taxidA %in% plant_taxid, taxcatA := "plant"]
intact[taxidA %in% virus_taxid, taxcatA := "virus"]
intact[taxidA == "9606", taxcatA := "human"]
intact[taxidB %in% archaea_taxid, taxcatB := "archaea"]
intact[taxidB %in% animal_taxid, taxcatB := "animal"]
intact[taxidB %in% bacteria_taxid, taxcatB := "bacteria"]
intact[taxidB %in% fungi_taxid, taxcatB := "fungi"]
intact[taxidB %in% plant_taxid, taxcatB := "plant"]
intact[taxidB %in% virus_taxid, taxcatB := "virus"]
intact[taxidB == "9606", taxcatB := "human"]
setorder(intact)
intact = unique(intact)
save(intact, file = "~/database/intact/intact.RObject")

############################
### parse HPIDB PPI data ###
############################
hpidb = fread("~/database/intact/hpidb2.mitab_plus.txt", colClasses = "character")
setnames(hpidb, trimws(str_replace_all(names(hpidb), "#|\\((.*?)\\)", "")))
hpidb[, setdiff(names(hpidb), c("protein_xref_1_unique", "protein_xref_2_unique", "protein_taxid_1", "protein_taxid_2", "pmid")) := NULL]
setcolorder(hpidb, c("protein_xref_1_unique", "protein_xref_2_unique", "protein_taxid_1", "protein_taxid_2", "pmid"))
setnames(hpidb, c("idA", "idB", "taxidA", "taxidB", "pmid"))
hpidb[, c("idA", "idB") := lapply(.SD, toupper), .SDcols = c("idA", "idB")]
hpidb[, c("taxidA", "taxidB") := lapply(.SD, str_extract, pattern = "(?<=taxid:)(.*?)(?=[^0-9]+)"), .SDcols = c("taxidA", "taxidB")]
hpidb[taxidA %in% taxmerged$taxid_old, taxidA := taxmerged[match(taxidA, taxid_old), taxid_new]]
hpidb[taxidB %in% taxmerged$taxid_old, taxidB := taxmerged[match(taxidB, taxid_old), taxid_new]]
hpidb = hpidb[str_detect(idA, uniprot_regex) & str_detect(idB, uniprot_regex) & taxidA %in% c(animal_taxid, fungi_taxid, plant_taxid, archaea_taxid, bacteria_taxid, virus_taxid) & taxidB %in% c(animal_taxid, fungi_taxid, plant_taxid, archaea_taxid, bacteria_taxid, virus_taxid)]
hpidb[, c("idA", "idB") := lapply(.SD, str_extract, pattern = uniprot_regex), .SDcols = c("idA", "idB")]

# get canonical IDs (cids) from isoform IDs and remove self-interactions between different isoforms
hpidb[, c("cidA", "cidB") := lapply(.SD, str_replace, pattern = "(?=-)(.*?)$", replacement = ""), .SDcols = c("idA", "idB")]

# for host-pathogen PPIs, make sure host protein ID is in front
colsAB = c(grep("A$", names(hpidb), value = TRUE), grep("B$", names(hpidb), value = TRUE))
colsBA = c(grep("B$", names(hpidb), value = TRUE), grep("A$", names(hpidb), value = TRUE))
hpidb[cidA > cidB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
hpidb[taxidA %in% c(archaea_taxid, bacteria_taxid, virus_taxid, fungi_taxid) & taxidB %in% c(animal_taxid, plant_taxid), c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
hpidb[taxidA != "9606" & taxidB == "9606", c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
hpidb[, c("taxnameA", "taxnameB") := lapply(.SD, function(x) taxnames[.(x), taxname]), .SDcols = c("taxidA", "taxidB")]
hpidb[taxidA %in% archaea_taxid, taxcatA := "archaea"]
hpidb[taxidA %in% animal_taxid, taxcatA := "animal"]
hpidb[taxidA %in% bacteria_taxid, taxcatA := "bacteria"]
hpidb[taxidA %in% fungi_taxid, taxcatA := "fungi"]
hpidb[taxidA %in% plant_taxid, taxcatA := "plant"]
hpidb[taxidA %in% virus_taxid, taxcatA := "virus"]
hpidb[taxidA == "9606", taxcatA := "human"]
hpidb[taxidB %in% archaea_taxid, taxcatB := "archaea"]
hpidb[taxidB %in% animal_taxid, taxcatB := "animal"]
hpidb[taxidB %in% bacteria_taxid, taxcatB := "bacteria"]
hpidb[taxidB %in% fungi_taxid, taxcatB := "fungi"]
hpidb[taxidB %in% plant_taxid, taxcatB := "plant"]
hpidb[taxidB %in% virus_taxid, taxcatB := "virus"]
hpidb[taxidB == "9606", taxcatB := "human"]
setorder(hpidb)
hpidb = unique(hpidb)
save(hpidb, file = "~/database/intact/hpidb.RObject")

################################
### combine IntAct and HPIDB ###
################################
ppi = rbindlist(list(intact, hpidb), use.names = TRUE)
ppi[grepl("-PRO_", idA), idA := cidA]
ppi[grepl("-PRO_", idB), idB := cidB]
ppi = ppi[str_detect(pmid, "pubmed") & !grepl("pubmed:20711500", pmid, fixed = TRUE)]
writeLines(ppi[, sort(unique(c(idA, idB, cidA, cidB)))], "~/database/intact/ppi_id.txt")
save(ppi, file = "~/database/intact/ppi.RObject")
