#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --time=3:00:00
#SBATCH --mem=60G
#SBATCH --job-name=pfam_host_other
#SBATCH --output=%x-%j.out

library(data.table)
library(RcppGreedySetCover)
library(stringr)

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

refnrprot = fread("~/database/uniprot/uniprot_refnr_proteomes.tab")
ppmembership = fread("~/database/uniprot/PPMembership.txt")

pfam2proteome = function(DT, proteomes_organism) {
  unique(rbind(DT[!grepl(";", pfam), list(proteomes = unlist(strsplit(proteomes, ";"))), by = "pfam"],
               DT[grepl(";", pfam), list(pfam = unlist(strsplit(pfam, ";"))), by = "proteomes"], use.names = TRUE)[, list(proteomes = paste0(sort(intersect(proteomes_organism, unlist(strsplit(proteomes, ";")))), collapse = ";")), by = "pfam"][nchar(proteomes) > 0])
}

# To minimize spurious domain matches due to cross-species contamination of genome sequences,
# retain domains that are found in at least three host (bacterial) proteomes, one of which must be a reference proteome or belong to a pan proteome.

proteomes_animal = fread("~/database/uniprot/proteomes_animal.tab", colClasses = "character")
proteomes_plant = fread("~/database/uniprot/proteomes_plant.tab", colClasses = "character")
proteomes_fungi = fread("~/database/uniprot/proteomes_fungi.tab", colClasses = "character")
proteomes_host = rbindlist(list(proteomes_animal, proteomes_plant, proteomes_fungi))
proteomes_host[, speciesname := rankedlineage[match(`Organism ID`, taxid), species]]
proteomes_host[, speciesid := rankedlineage[match(speciesname, species), taxid]]
proteomes_host[, speciesname := NULL]
proteomes_host[, genusname := rankedlineage[match(`Organism ID`, taxid), genus]]
proteomes_host[, genusid := rankedlineage[match(genusname, genus), taxid]]
proteomes_host[, genusname := NULL]
proteomes_host[, refnr_proteome := ifelse(`Proteome ID` %in% refnrprot$`Proteome ID`, `Proteome ID`, "")]
proteomes_host = proteomes_host[!is.na(speciesid)]
phibase_host_speciesids = proteomes_host[, intersect(speciesid, readLines("~/database/uniprot/phibase_host_speciesids.txt"))]
phibase_host_proteomes = proteomes_host[speciesid %in% phibase_host_speciesids, sort(unique(`Proteome ID`))]

animal_pfam2proteome = fread("~/database/uniprot/animal_pfam2proteome.tab", header = FALSE)
setnames(animal_pfam2proteome, c("pfam", "proteomes"))
plant_pfam2proteome = fread("~/database/uniprot/plant_pfam2proteome.tab", header = FALSE)
setnames(plant_pfam2proteome, c("pfam", "proteomes"))
fungi_pfam2proteome = fread("~/database/uniprot/fungi_pfam2proteome.tab", header = FALSE)
setnames(fungi_pfam2proteome, c("pfam", "proteomes"))
host_pfam2proteome = rbindlist(list(animal_pfam2proteome, plant_pfam2proteome, fungi_pfam2proteome))
host_pfam2proteome = unique(host_pfam2proteome[nchar(proteomes) > 0])
pfam2proteome_host = pfam2proteome(host_pfam2proteome, proteomes_host$`Proteome ID`)
setnames(pfam2proteome_host, "proteomes", "host_proteomes")
pfam2proteome_host[, host_refnr_proteomes := paste0(intersect(unlist(strsplit(host_proteomes, ";")), refnrprot$`Proteome ID`), collapse = ";"), by = "pfam"]
pfam2proteome_host[, host_speciesids := proteomes_host[`Proteome ID` %in% unlist(strsplit(host_proteomes, ";")), paste0(unique(speciesid), collapse = ";")], by = "pfam"]
pfam2proteome_host[, host_proteomes_phibase := paste0(intersect(unlist(strsplit(host_proteomes, ";")), phibase_host_proteomes), collapse = ";"), by = "pfam"]
pfam2proteome_host[, host_refnr_proteomes_phibase := paste0(intersect(unlist(strsplit(host_refnr_proteomes, ";")), phibase_host_proteomes), collapse = ";"), by = "pfam"]
pfam2proteome_host[, host_speciesids_phibase := paste0(intersect(unlist(strsplit(host_speciesids, ";")), phibase_host_speciesids), collapse = ";"), by = "pfam"]
pfam2proteome_host[, N_host_proteomes := ifelse(nchar(host_proteomes) == 0, 0, str_count(host_proteomes, ";") + 1)]
pfam2proteome_host[, N_host_refnr_proteomes := ifelse(nchar(host_refnr_proteomes) == 0, 0, str_count(host_refnr_proteomes, ";") + 1)]
pfam2proteome_host[, N_host_species := ifelse(nchar(host_speciesids) == 0, 0, str_count(host_speciesids, ";") + 1)]
pfam2proteome_host[, N_host_proteomes_phibase := ifelse(nchar(host_proteomes_phibase) == 0, 0, str_count(host_proteomes_phibase, ";") + 1)]
pfam2proteome_host[, N_host_refnr_proteomes_phibase := ifelse(nchar(host_refnr_proteomes_phibase) == 0, 0, str_count(host_refnr_proteomes_phibase, ";") + 1)]
pfam2proteome_host[, N_host_species_phibase := ifelse(nchar(host_speciesids_phibase) == 0, 0, str_count(host_speciesids_phibase, ";") + 1)]
pfam2proteome_host = pfam2proteome_host[N_host_proteomes > 2 & N_host_refnr_proteomes > 0]
setkey(pfam2proteome_host, pfam)

proteomes_bacteria = fread("~/database/uniprot/proteomes_bacteria.tab", colClasses = "character")
proteomes_bacteria[, speciesname := rankedlineage[match(`Organism ID`, taxid), species]]
proteomes_bacteria[, speciesid := rankedlineage[match(speciesname, species), taxid]]
proteomes_bacteria[, speciesname := NULL]
proteomes_bacteria[, genusname := rankedlineage[match(`Organism ID`, taxid), genus]]
proteomes_bacteria[, genusid := rankedlineage[match(genusname, genus), taxid]]
proteomes_bacteria[, genusname := NULL]
proteomes_bacteria[, refnr_proteome := ifelse(`Proteome ID` %in% refnrprot$`Proteome ID`, `Proteome ID`, "")]
proteomes_bacteria[, pan_proteome := ppmembership[match(`Proteome ID`, PPMember), PP]]
proteomes_bacteria[is.na(pan_proteome), pan_proteome := ""]
proteomes_bacteria = proteomes_bacteria[!is.na(speciesid)]

# pathogenic species = PHI-base pathogen species
phibase_bacpatho_speciesids = proteomes_bacteria[!is.na(genusid), intersect(speciesid, readLines("~/database/uniprot/phibase_bacpatho_speciesids.txt"))]
phibase_bacpatho_genera = rankedlineage[taxid %in% phibase_bacpatho_speciesids, unique(genus)]
phibase_bacpatho_proteomes = proteomes_bacteria[speciesid %in% phibase_bacpatho_speciesids, sort(unique(`Proteome ID`))]
# non-pathogenic species = all other species in the same genera as PHI-base pathogen species
phibase_bacnonpatho_speciesids = proteomes_bacteria[!is.na(genusid), intersect(speciesid, readLines("~/database/uniprot/phibase_bacnonpatho_speciesids.txt"))]
phibase_bacnonpatho_proteomes = proteomes_bacteria[speciesid %in% phibase_bacnonpatho_speciesids, sort(unique(`Proteome ID`))]

# All pathogenic species = PHI-base pathogen species + effector-encoding species in the same genera
effector = rbindlist(list(fread("~/database/uniprot/uniprot_bacteria_effector_name.tab.gz", colClasses = "character"),
                          fread("~/database/uniprot/uniprot_bacteria_effector_location.tab.gz", colClasses = "character"),
                          fread("~/database/uniprot/uniprot_bacteria_effector_phibase.tab.gz", colClasses = "character")))
effector[, speciesname := rankedlineage[match(`Organism ID`, taxid), species]]
effector[, speciesid := rankedlineage[match(speciesname, species), taxid]]
bacpatho_speciesids = rankedlineage[genus %in% phibase_bacpatho_genera & taxrank == "species", intersect(taxid, union(effector$speciesid, phibase_bacpatho_speciesids))]
bacpatho_speciesids = proteomes_bacteria[, intersect(speciesid, bacpatho_speciesids)]
bacpatho_proteomes = proteomes_bacteria[speciesid %in% bacpatho_speciesids, sort(unique(`Proteome ID`))]
# All non-pathogenic species = effector-less species in the same genera as PHI-base pathogen species
bacnonpatho_speciesids = rankedlineage[genus %in% phibase_bacpatho_genera & taxrank == "species", setdiff(taxid, bacpatho_speciesids)]
bacnonpatho_speciesids = proteomes_bacteria[, intersect(speciesid, bacnonpatho_speciesids)]
bacnonpatho_proteomes = proteomes_bacteria[speciesid %in% bacnonpatho_speciesids, sort(unique(`Proteome ID`))]

bacteria_pfam2proteome = fread("~/database/uniprot/bacteria_pfam2proteome.tab", header = FALSE)
setnames(bacteria_pfam2proteome, c("pfam", "proteomes"))
bacteria_pfam2proteome = unique(bacteria_pfam2proteome[nchar(proteomes) > 0])
pfam2proteome_bacteria = pfam2proteome(bacteria_pfam2proteome, proteomes_bacteria$`Proteome ID`)
setnames(pfam2proteome_bacteria, "proteomes", "bacteria_proteomes")
pfam2proteome_bacteria[, bacteria_refnr_proteomes := paste0(intersect(unlist(strsplit(bacteria_proteomes, ";")), refnrprot$`Proteome ID`), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_pan_proteomes := ppmembership[PPMember %in% unlist(strsplit(bacteria_proteomes, ";")), paste0(sort(unique(PP)), collapse = ";")], by = "pfam"]
pfam2proteome_bacteria[, bacteria_speciesids := proteomes_bacteria[`Proteome ID` %in% unlist(strsplit(bacteria_proteomes, ";")), paste0(unique(speciesid), collapse = ";")], by = "pfam"]
pfam2proteome_bacteria[, bacteria_proteomes_patho := paste0(intersect(unlist(strsplit(bacteria_proteomes, ";")), phibase_bacpatho_proteomes), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_refnr_proteomes_patho := paste0(intersect(unlist(strsplit(bacteria_refnr_proteomes, ";")), phibase_bacpatho_proteomes), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_pan_proteomes_patho := ppmembership[PPMember %in% unlist(strsplit(bacteria_proteomes_patho, ";")), paste0(sort(unique(PP)), collapse = ";")], by = "pfam"]
pfam2proteome_bacteria[, bacteria_speciesids_patho := paste0(intersect(unlist(strsplit(bacteria_speciesids, ";")), phibase_bacpatho_speciesids), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_genera_patho := rankedlineage[taxid %in% unlist(strsplit(bacteria_speciesids_patho, ";")), paste0(unique(genus), collapse = ";")], by = "pfam"]
pfam2proteome_bacteria[, bacteria_proteomes_nonpatho := paste0(intersect(unlist(strsplit(bacteria_proteomes, ";")), bacnonpatho_proteomes), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_refnr_proteomes_nonpatho := paste0(intersect(unlist(strsplit(bacteria_refnr_proteomes, ";")), bacnonpatho_proteomes), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_pan_proteomes_nonpatho := ppmembership[PPMember %in% unlist(strsplit(bacteria_proteomes_nonpatho, ";")), paste0(sort(unique(PP)), collapse = ";")], by = "pfam"]
pfam2proteome_bacteria[, bacteria_speciesids_nonpatho := paste0(intersect(unlist(strsplit(bacteria_speciesids, ";")), bacnonpatho_speciesids), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, bacteria_genera_nonpatho := rankedlineage[taxid %in% unlist(strsplit(bacteria_speciesids_nonpatho, ";")), paste0(unique(genus), collapse = ";")], by = "pfam"]
pfam2proteome_bacteria[, bacteria_genera_pnp := paste0(intersect(unlist(strsplit(bacteria_genera_patho, ";")), unlist(strsplit(bacteria_genera_nonpatho, ";"))), collapse = ";"), by = "pfam"]
pfam2proteome_bacteria[, N_bacteria_proteomes := ifelse(nchar(bacteria_proteomes) == 0, 0, str_count(bacteria_proteomes, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_refnr_proteomes := ifelse(nchar(bacteria_refnr_proteomes) == 0, 0, str_count(bacteria_refnr_proteomes, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_pan_proteomes := ifelse(nchar(bacteria_pan_proteomes) == 0, 0, str_count(bacteria_pan_proteomes, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_species := ifelse(nchar(bacteria_speciesids) == 0, 0, str_count(bacteria_speciesids, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_proteomes_patho := ifelse(nchar(bacteria_proteomes_patho) == 0, 0, str_count(bacteria_proteomes_patho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_refnr_proteomes_patho := ifelse(nchar(bacteria_refnr_proteomes_patho) == 0, 0, str_count(bacteria_refnr_proteomes_patho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_pan_proteomes_patho := ifelse(nchar(bacteria_pan_proteomes_patho) == 0, 0, str_count(bacteria_pan_proteomes_patho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_species_patho := ifelse(nchar(bacteria_speciesids_patho) == 0, 0, str_count(bacteria_speciesids_patho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_proteomes_nonpatho := ifelse(nchar(bacteria_proteomes_nonpatho) == 0, 0, str_count(bacteria_proteomes_nonpatho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_refnr_proteomes_nonpatho := ifelse(nchar(bacteria_refnr_proteomes_nonpatho) == 0, 0, str_count(bacteria_refnr_proteomes_nonpatho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_pan_proteomes_nonpatho := ifelse(nchar(bacteria_pan_proteomes_nonpatho) == 0, 0, str_count(bacteria_pan_proteomes_nonpatho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_species_nonpatho := ifelse(nchar(bacteria_speciesids_nonpatho) == 0, 0, str_count(bacteria_speciesids_nonpatho, ";") + 1)]
pfam2proteome_bacteria[, N_bacteria_genera_pnp := ifelse(nchar(bacteria_genera_pnp) == 0, 0, str_count(bacteria_genera_pnp, ";") + 1)]
pfam2proteome_bacteria = pfam2proteome_bacteria[N_bacteria_proteomes > 2 & (N_bacteria_refnr_proteomes > 0 | N_bacteria_pan_proteomes > 0)]
setkey(pfam2proteome_bacteria, pfam)

# number of pathogenic and non-pathogenic species within the same genus
genus_patho = rankedlineage[taxid %in% pfam2proteome_bacteria[, unique(unlist(strsplit(bacteria_speciesids_patho, ";")))], genus]
genus_nonpatho = rankedlineage[taxid %in% pfam2proteome_bacteria[, unique(unlist(strsplit(bacteria_speciesids_nonpatho, ";")))], genus]
genus_tab = t(rbind(table(genus_patho)[intersect(genus_patho, genus_nonpatho)], table(genus_nonpatho)[intersect(genus_patho, genus_nonpatho)]))
genus_tab = cbind.data.frame(intersect(genus_patho, genus_nonpatho), genus_tab)
rownames(genus_tab) = NULL
setnames(genus_tab, c("genus", "pathogenic species", "non-pathogenic species"))
setorder(genus_tab, genus)

save(pfam2proteome_host, pfam2proteome_bacteria, proteomes_host, proteomes_bacteria, genus_tab, file = "~/database/uniprot/pfam_host_other.RData")
