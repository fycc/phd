#!/usr/bin/env Rscript

library(data.table)

phibase = fread("~/database/uniprot/phi-base_current.csv", sep = ",", colClasses = "character")

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

phibase_host_speciesids = rankedlineage[kingdom %in% c("Metazoa", "Viridiplantae") & taxrank == "species", intersect(taxid, unname(unlist(fread("~/database/uniprot/host_species_list.csv", select = 1, colClasses = "character"))))]
writeLines(phibase_host_speciesids, "~/database/uniprot/phibase_host_speciesids.txt")
phibase_host_speciesnames = rankedlineage[match(phibase_host_speciesids, taxid), unique(taxname)]
writeLines(phibase_host_speciesnames, "~/database/uniprot/phibase_host_speciesnames.txt")

phibase_bacpatho_speciesids = rankedlineage[superkingdom == "Bacteria" & taxrank == "species", intersect(taxid, unname(unlist(fread("~/database/uniprot/pathogen_species_list.csv", select = 1, colClasses = "character"))))]
writeLines(phibase_bacpatho_speciesids, "~/database/uniprot/phibase_bacpatho_speciesids.txt")
phibase_bacpatho_speciesnames = rankedlineage[match(phibase_bacpatho_speciesids, taxid), unique(taxname)]
writeLines(phibase_bacpatho_speciesnames, "~/database/uniprot/phibase_bacpatho_speciesnames.txt")
phibase_bacpatho_taxids = rankedlineage[species %in% phibase_bacpatho_speciesnames, taxid]
writeLines(phibase_bacpatho_taxids, "~/database/uniprot/phibase_bacpatho_taxids.txt")
phibase_bacpatho_genera = rankedlineage[species %in% phibase_bacpatho_speciesnames, sort(unique(genus))]
writeLines(phibase_bacpatho_genera, "~/database/uniprot/phibase_bacpatho_genera.txt")
phibase_bacnonpatho_speciesids = rankedlineage[genus %in% phibase_bacpatho_genera & taxrank == "species", setdiff(taxid, phibase_bacpatho_speciesids)]
writeLines(phibase_bacnonpatho_speciesids, "~/database/uniprot/phibase_bacnonpatho_speciesids.txt")
phibase_bacnonpatho_speciesnames = rankedlineage[match(phibase_bacnonpatho_speciesids, taxid), taxname]
writeLines(phibase_bacnonpatho_speciesnames, "~/database/uniprot/phibase_bacnonpatho_speciesnames.txt")
phibase_effector_uniprotac = phibase[grepl("(?i)uniprot", `Protein ID source`) & grepl("(?i)effector", paste(`Gene Function`, `Mutant Phenotype`)), unique(`Protein ID`)]
writeLines(phibase_effector_uniprotac, "~/database/uniprot/phibase_effector_uniprotac.txt")
