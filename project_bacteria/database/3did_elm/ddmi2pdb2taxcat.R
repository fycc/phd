#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=8G
#SBATCH --job-name=ddmi2pdb2taxcat
#SBATCH --output=%x-%j.out

library(data.table)
library(stringr)

##############################
### Map DDI and DMI to PDB ###
##############################

# replace obsolete PDB IDs by superseded PDB IDs
pdb_obsolete = readLines("~/database/3did_elm/pdb_obsolete.dat")
pdb_old = sapply(pdb_obsolete[-1], function(x) {y = unlist(strsplit(x, " ")); y[nchar(y) > 0][3]})
pdb_new = sapply(pdb_obsolete[-1], function(x) {
  y = unlist(strsplit(x, " "))
  y = y[nchar(y) > 0]
  if (length(y) == 3) NA
  else paste0(y[4:length(y)], collapse = ";")
})
pdb_old2new = data.table(pdb_old = tolower(pdb_old), pdb_new = tolower(pdb_new))
pdb_old2new = unique(pdb_old2new[, list(pdb_new = unlist(strsplit(pdb_new, ";"))), by = "pdb_old"])
setkeyv(pdb_old2new, "pdb_old")

o2n = function(x) {while (x %in% pdb_old2new$pdb_old) {x = pdb_old2new[.(x), pdb_new]}; x}

# read DDI to PDB mapping file: https://3did.irbbarcelona.org/download/current/3did_flat.gz
ddi_pdb = grep("^#=ID|^#=3D", readLines("~/database/3did_elm/3did_flat"), value = TRUE)
ddi_pdb = unname(sapply(ddi_pdb, function(x) ifelse(grepl("^#=ID", x), paste0(unlist(str_extract_all(str_extract_all(x, "(?<=[(])(.*?)(?=[)])"), "PF\\d{5}")), collapse = "_"), paste0(unlist(strsplit(x, "\t|:"))[c(2, 3, 5)], collapse = "_"))))
ddi_pdb = unname(sapply(ddi_pdb, function(x) ifelse(unlist(strsplit(x, "_"))[1] %in% pdb_old2new$pdb_old, paste0(c(o2n(unlist(strsplit(x, "_"))[1]), unlist(strsplit(x, "_"))[-1]), collapse = "_"), x)))
ddi_pdb = ddi_pdb[!grepl("^NA_", ddi_pdb)]

pos_ddi = grep("PF\\d{5}_PF\\d{5}", ddi_pdb)
ddis = ddi_pdb[pos_ddi]
ddis = unname(sapply(ddis, function(x) ifelse(unlist(strsplit(x, "_"))[1] >= unlist(strsplit(x, "_"))[2], x, paste(unlist(strsplit(x, "_"))[2], unlist(strsplit(x, "_"))[1], sep = "_"))))

ddi2pdb = list()
for (i in 1:length(pos_ddi)) {
  if (i < length(pos_ddi)) ddi2pdb[[i]] = grep("PF\\d{5}", unique(ddi_pdb[seq(from = pos_ddi[i] + 1, to = pos_ddi[i + 1] - 1)]), invert = TRUE, value = TRUE)
  else ddi2pdb[[i]] = unique(ddi_pdb[seq(from = pos_ddi[i] + 1, to = length(ddi_pdb))])
}
names(ddi2pdb) = ddis
ddi2pdb = ddi2pdb[lengths(ddi2pdb) > 0]

# read DMI to PDB mapping file: https://3did.irbbarcelona.org/download/current/3did_dmi_flat.gz
dmi_pdb = grep("^#=ID|^#=3D", readLines("~/database/3did_elm/3did_dmi_flat"), value = TRUE)
dmi_pdb[grep("^#=3D", dmi_pdb)] = unname(sapply(dmi_pdb[grep("^#=3D", dmi_pdb)], function(x) paste0(unlist(strsplit(x, "\t|:"))[c(2, 3, 5)], collapse = "_")))
dmi_pdb = unname(sapply(dmi_pdb, function(x) ifelse(unlist(strsplit(x, "_"))[1] %in% pdb_old2new$pdb_old, paste0(c(o2n(unlist(strsplit(x, "_"))[1]), unlist(strsplit(x, "_"))[-1]), collapse = "_"), x)))
dmi_pdb = dmi_pdb[!grepl("^NA_", dmi_pdb)]
pos_dmi = grep("^#=ID", dmi_pdb)
dmi2pdb = list()
for (i in 1:length(pos_dmi)) {
  if (i < length(pos_dmi)) dmi2pdb[[i]] = grep("#=ID", unique(dmi_pdb[seq(from = pos_dmi[i] + 1, to = pos_dmi[i + 1] - 1)]), invert = TRUE, value = TRUE)
  else dmi2pdb[[i]] = unique(dmi_pdb[seq(from = pos_dmi[i] + 1, to = length(dmi_pdb))])
}

load("~/database/3did_elm/ddi_dmi.RData")
names(dmi2pdb) = paste(dmi_domain$PFAM_AC, dmi_motif, sep = "_")
dmi2pdb = dmi2pdb[lengths(dmi2pdb) > 0]

# check if domain assignment by 3did matches that by Pfam
pdb2pfam = fread("~/database/3did_elm/pdb_chain_pfam.tsv.gz", quote = "", colClasses = "character")
pdb_pfam_mapping = fread("~/database/3did_elm/pdb_pfam_mapping.txt")
setnames(pdb_pfam_mapping, c("PDB_ID", "CHAIN_ID", "PFAM_ACC"), c("PDB", "CHAIN", "PFAM_ID"))
pdb2pfam = rbind(pdb2pfam[, .(PDB, CHAIN, PFAM_ID)], pdb_pfam_mapping[, .(PDB, CHAIN, PFAM_ID)])
pdb2pfam[, PDB := tolower(PDB)]
pdb2pfam[, PFAM_ID := str_extract(PFAM_ID, "PF\\d{5}")]
pdb2pfam = unique(pdb2pfam)
setkeyv(pdb2pfam, c("PDB", "CHAIN", "PFAM_ID"))

pfam_check = function(ddi) {
  m = sapply(ddi2pdb[[ddi]], function(pdb_chains) {
    pdb = unlist(strsplit(pdb_chains, "_"))[1]
    chain1 = unlist(strsplit(pdb_chains, "_"))[2]
    chain2 = unlist(strsplit(pdb_chains, "_"))[3]
    all(unlist(strsplit(ddi, "_")) %in% union(pdb2pfam[.(pdb, chain1), PFAM_ID], pdb2pfam[.(pdb, chain2), PFAM_ID]))
  })
  ddi2pdb[[ddi]][m]
}

ddi2pdb = sapply(names(ddi2pdb), pfam_check)
ddi2pdb = ddi2pdb[lengths(ddi2pdb) > 0]

######################################################
### Map DDI and DMI to taxonomy IDs and categories ###
######################################################

nodes = fread("~/database/taxonomy/nodes.tab", header = FALSE, colClasses = "character")
names(nodes)[1:3] = c("taxid", "parenttaxid", "taxrank")
rankedlineage = fread("~/database/taxonomy/rankedlineage.tab", header = FALSE, colClasses = "character")
setnames(rankedlineage, c("taxid", "taxname", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"))
rankedlineage[, taxrank := nodes[match(rankedlineage$taxid, taxid), taxrank]]
rankedlineage[nchar(species) == 0 & taxrank == "species", species := taxname]
rankedlineage[nchar(genus) == 0 & taxrank == "genus", genus := taxname]
rankedlineage[nchar(family) == 0 & taxrank == "family", family := taxname]
rankedlineage[nchar(order) == 0 & taxrank == "order", order := taxname]
rankedlineage[nchar(class) == 0 & taxrank == "class", class := taxname]
rankedlineage[nchar(phylum) == 0 & taxrank == "phylum", phylum := taxname]

animal_taxid = readLines("~/database/taxonomy/animal_taxid.txt")
archaea_taxid = readLines("~/database/taxonomy/archaea_taxid.txt")
bacteria_taxid = readLines("~/database/taxonomy/bacteria_taxid.txt")
fungi_taxid = readLines("~/database/taxonomy/fungi_taxid.txt")
plant_taxid = readLines("~/database/taxonomy/plant_taxid.txt")
virus_taxid = readLines("~/database/taxonomy/virus_taxid.txt")
taxid2taxcat = data.table(taxid = c(animal_taxid, archaea_taxid, bacteria_taxid, fungi_taxid, plant_taxid, virus_taxid),
                          taxcat = rep(c("animal", "archaea", "bacteria", "fungi", "plant", "virus"), c(length(animal_taxid), length(archaea_taxid), length(bacteria_taxid), length(fungi_taxid), length(plant_taxid), length(virus_taxid))))
setkeyv(taxid2taxcat, "taxid")

# ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz
pdb2tax = fread("~/database/3did_elm/pdb_chain_taxonomy.tsv.gz", quote = "", colClasses = "character")
pdb2tax[, PDB := tolower(PDB)]
pdb2tax = unique(pdb2tax[, .(PDB, CHAIN, TAX_ID)])
pdb2tax[, taxcat := taxid2taxcat[match(TAX_ID, taxid), taxcat]]
pdb2tax[, taxcat := paste0(sort(unique(taxcat)), collapse = "_"), by = c("PDB", "CHAIN")]
pdb2tax[, TAX_ID := paste0(sort(unique(TAX_ID)), collapse = "_"), by = c("PDB", "CHAIN")]
setnames(pdb2tax, "TAX_ID", "taxid")
pdb2tax = pdb2tax[nchar(taxcat) > 0]

# ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
uniprot2tax = fread("~/database/3did_elm/pdb_chain_uniprotac.tab", quote = "", colClasses = "character")
pdb2uniprot = fread("~/database/3did_elm/pdb_chain_uniprot.tsv.gz", quote = "", colClasses = "character")
pdb2uniprot[, PDB := tolower(PDB)]
pdb2uniprot[, taxid := uniprot2tax[match(SP_PRIMARY, uniprotac), taxid]]
pdb2uniprot = unique(pdb2uniprot[!is.na(taxid), .(PDB, CHAIN, taxid)])
pdb2uniprot[, taxcat := taxid2taxcat[match(pdb2uniprot$taxid, taxid), taxcat]]
pdb2uniprot[, taxcat := paste0(sort(unique(taxcat)), collapse = "_"), by = c("PDB", "CHAIN")]
pdb2uniprot[, taxid := paste0(sort(unique(taxid)), collapse = "_"), by = c("PDB", "CHAIN")]
pdb2uniprot = pdb2uniprot[nchar(taxcat) > 0]

pdb2taxcat = rbind(pdb2tax[, .(PDB, CHAIN, taxid, taxcat)], pdb2uniprot[, .(PDB, CHAIN, taxid, taxcat)])
pdb2taxcat[, source := rep(c("pdb_chain_taxonomy", "pdb_chain_uniprot"), c(nrow(pdb2tax), nrow(pdb2uniprot)))]
setorderv(pdb2taxcat, c("PDB", "CHAIN", "source"), c(1, 1, -1))
pdb2taxcat = unique(pdb2taxcat, by = c("PDB", "CHAIN"))
setkey(pdb2taxcat)

p2t_all = data.table(pdb_chains = union(unlist(ddi2pdb), unlist(dmi2pdb)))
p2t_all[, pdb := unlist(strsplit(pdb_chains, "_"))[1], by = "pdb_chains"]
p2t_all[, chain1 := unlist(strsplit(pdb_chains, "_"))[2], by = "pdb_chains"]
p2t_all[, chain2 := unlist(strsplit(pdb_chains, "_"))[3], by = "pdb_chains"]
p2t_all[, taxid1 := pdb2taxcat[.(pdb, chain1), taxid]]
p2t_all[, taxid2 := pdb2taxcat[.(pdb, chain2), taxid]]
p2t_all[, taxcat1 := pdb2taxcat[.(pdb, chain1), taxcat]]
p2t_all[, taxcat2 := pdb2taxcat[.(pdb, chain2), taxcat]]
p2t_all[taxcat1 == taxcat2, taxcat := taxcat1]
p2t_all[taxcat1 != taxcat2, taxcat := paste0(sort(unique(unlist(strsplit(c(taxcat1, taxcat2), "_")))), collapse = "_"), by = "pdb_chains"]
p2t_all[grepl("animal|fungi|plant", taxcat1) & grepl("animal|fungi|plant", taxcat2), taxcat := "eukaryote"]
# for DDIs between host proteins and bacterial effector proteins, replace "animal/plant_bacteria" with "host_effector"
effector = rbindlist(list(fread("~/database/uniprot/uniprot_bacteria_effector_name.tab.gz", select = c("Entry", "Status", "Length", "Annotation")),
                          fread("~/database/uniprot/uniprot_bacteria_effector_location.tab.gz", select = c("Entry", "Status", "Length", "Annotation")),
                          fread("~/database/uniprot/uniprot_bacteria_effector_phibase.tab.gz", select = c("Entry", "Status", "Length", "Annotation"))))
pdb2uniprot = fread("~/database/3did_elm/pdb_chain_uniprot.tsv.gz", quote = "", colClasses = "character")
pdb2uniprot[, PDB := tolower(PDB)]
setkeyv(pdb2uniprot, c("PDB", "CHAIN"))
pathobac_pdb_chains = union(p2t_all[grepl("animal|plant", taxcat1) & grepl("bacteria", taxcat2), pdb2uniprot[.(pdb, chain2), any(SP_PRIMARY %in% effector$Entry)], by = "pdb_chains"][V1 == TRUE, pdb_chains],
                            p2t_all[grepl("animal|plant", taxcat2) & grepl("bacteria", taxcat1), pdb2uniprot[.(pdb, chain1), any(SP_PRIMARY %in% effector$Entry)], by = "pdb_chains"][V1 == TRUE, pdb_chains])
p2t_all[pdb_chains %in% pathobac_pdb_chains, taxcat := "host_effector"]
setkeyv(p2t_all, "pdb_chains")
rm(nodes, rankedlineage)
gc(TRUE)

ddi2pdb_all = ddi2pdb
ddi2taxcat_all = lapply(ddi2pdb_all, function(x) p2t_all[.(x), sort(unique(taxcat))])
ddi2taxcat_all = ddi2taxcat_all[lengths(ddi2taxcat_all) > 0]

dmi2pdb_all = dmi2pdb
dmi2taxcat_all = lapply(dmi2pdb_all, function(x) p2t_all[.(x), sort(unique(taxcat))])
dmi2taxcat_all = dmi2taxcat_all[lengths(dmi2taxcat_all) > 0]

################################################################
### Keep DDIs between different chains of different proteins ###
################################################################

pdb_nprot = fread("~/database/3did_elm/pdb_nprot.csv", fill = TRUE)
while(length(ind <- which(pdb_nprot$`Entry ID` == "" | is.na(pdb_nprot$`Entry ID`))) > 0){
  pdb_nprot$`Entry ID`[ind] <- pdb_nprot$`Entry ID`[ind - 1]
}
pdb_nprot[, `Entry ID` := tolower(`Entry ID`)]
pdb_nprot[, `Chain ID` := as.character(`Chain ID`)]
pdb_nprot[, `Entity ID` := as.character(`Entity ID`)]
pdb_nprot = pdb_nprot[`Entity Polymer Type` == "Protein", .(`Entry ID`, `Entity Polymer Type`, `Chain ID`, `Entity ID`)]

p2t = p2t_all[pdb %in% pdb_nprot$`Entry ID` & chain1 != chain2]
p2t[, entity1 := pdb_nprot[`Entry ID` == pdb][grepl(paste0("\\b", chain1, "\\b"), `Chain ID`), unique(`Entity ID`)], by = c("pdb", "chain1")]
p2t[, entity2 := pdb_nprot[`Entry ID` == pdb][grepl(paste0("\\b", chain2, "\\b"), `Chain ID`), unique(`Entity ID`)], by = c("pdb", "chain2")]
p2t = p2t[!is.na(entity1) & !is.na(entity2) & entity1 != entity2]
setkeyv(p2t, "pdb_chains")

ddi2pdb = lapply(ddi2pdb_all, intersect, p2t$pdb_chains)
ddi2pdb = ddi2pdb[lengths(ddi2pdb) > 0]
ddi2taxcat = lapply(ddi2pdb, function(x) p2t[.(x), sort(unique(taxcat))])
ddi2taxcat = ddi2taxcat[lengths(ddi2taxcat) > 0]

dmi2pdb = lapply(dmi2pdb_all, intersect, p2t$pdb_chains)
dmi2pdb = dmi2pdb[lengths(dmi2pdb) > 0]
dmi2taxcat = lapply(dmi2pdb, function(x) p2t[.(x), sort(unique(taxcat))])
dmi2taxcat = dmi2taxcat[lengths(dmi2taxcat) > 0]

# for each host-effector DDI derived from PDB structure, retrieve host domain and effector domain
ddi_host_effector = names(ddi2taxcat)[sapply(ddi2taxcat, function(x) any(x %in% "host_effector"))]
pdb_chain2effector_pfam = function(pdb_chain) pdb2pfam[.(unlist(strsplit(pdb_chain, "_"))[1], p2t[pdb_chains == pdb_chain, ifelse(grepl("bacteria", taxcat1), chain1, chain2)]), PFAM_ID]
pdb_chain2host_pfam = function(pdb_chain) pdb2pfam[.(unlist(strsplit(pdb_chain, "_"))[1], p2t[pdb_chains == pdb_chain, ifelse(grepl("animal|plant", taxcat1), chain1, chain2)]), PFAM_ID]
pfam_pdb_bh_effector = sapply(ddi_host_effector, function(ddi) {
  intersect(unlist(sapply(intersect(ddi2pdb[[ddi]], p2t[taxcat == "host_effector", pdb_chains]), pdb_chain2effector_pfam)), unlist(strsplit(ddi, "_")))
})
pfam_pdb_bh_host = sapply(ddi_host_effector, function(ddi) {
  intersect(unlist(sapply(intersect(ddi2pdb[[ddi]], p2t[taxcat == "host_effector", pdb_chains]), pdb_chain2host_pfam)), unlist(strsplit(ddi, "_")))
})

save(pdb2pfam, pdb2taxcat, pdb_nprot, pathobac_pdb_chains, p2t_all, ddi2pdb_all, ddi2taxcat_all, dmi2pdb_all, dmi2taxcat_all,
     p2t, ddi2pdb, ddi2taxcat, dmi2pdb, dmi2taxcat, ddi_host_effector, pfam_pdb_bh_effector, pfam_pdb_bh_host, file = "~/database/3did_elm/ddmi2pdb2taxcat.RData")
