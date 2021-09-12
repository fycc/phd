#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=8G
#SBATCH --job-name=host_bacteria_pfam_ddi2species
#SBATCH --output=%x-%j.out

library(data.table)
library(stringr)

# function to compute log odds ratio and weight (Mantel-Haenszel method)
or_weight = function(a11, a1, a21, a2) {
  a12 = a1 - a11
  a22 = a2 - a21
  if (prod(a11, a12, a21, a22) == 0) return( list(or = (a11 + 0.5) * (a22 + 0.5) / (a12 + 0.5) / (a21 + 0.5), weight = (a12 + 0.5) * (a21 + 0.5) / (a11 + a22 + a12 + a21 + 2)) )
  else return( list(or = a11 * a22 / a12 / a21, weight = a12 * a21 / (a11 + a22 + a12 + a21)) )
}

############################################################
### Enrichment of domain in host vs. bacterial proteomes ###
############################################################

# Log odds ratio of domain frequency between host and bacterial species
load("~/database/uniprot/pfam_host_other.RData")

pfam2proteome_host_bacteria = merge(pfam2proteome_host[, .(pfam, host_speciesids, N_host_species)], pfam2proteome_bacteria[, .(pfam, bacteria_speciesids, N_bacteria_species)], all = TRUE, by = "pfam")
pfam2proteome_host_bacteria[is.na(host_speciesids), host_speciesids := ""]
pfam2proteome_host_bacteria[is.na(N_host_species), N_host_species := 0L]
pfam2proteome_host_bacteria[is.na(bacteria_speciesids), bacteria_speciesids := ""]
pfam2proteome_host_bacteria[is.na(N_bacteria_species), N_bacteria_species := 0L]
N_host_species_total = pfam2proteome_host_bacteria[, uniqueN(unlist(strsplit(host_speciesids, ";")))]
N_bacteria_species_total = pfam2proteome_host_bacteria[, uniqueN(unlist(strsplit(bacteria_speciesids, ";")))]
pfam2proteome_host_bacteria[, c("host2bacteria_species_or", "host2bacteria_species_weight") := or_weight(N_host_species, N_host_species_total, N_bacteria_species, N_bacteria_species_total), by = "pfam"]
setkeyv(pfam2proteome_host_bacteria, "pfam")

###########################################################################################
### Enrichment of domain in pathogenic vs. non-pathogenic species within the same genus ###
###########################################################################################

# Log odds ratio of domain frequency between pathogenic and non-pathogenic species within the same genus
N_bacteria_species_patho_total = pfam2proteome_bacteria[N_bacteria_genera_pnp > 0, uniqueN(unlist(strsplit(bacteria_speciesids_patho, ";")))]
N_bacteria_species_nonpatho_total = pfam2proteome_bacteria[N_bacteria_genera_pnp > 0, uniqueN(unlist(strsplit(bacteria_speciesids_nonpatho, ";")))]
pfam2proteome_bacteria[N_bacteria_genera_pnp > 0, c("patho2nonpatho_species_or", "patho2nonpatho_species_weight") := or_weight(N_bacteria_species_patho, N_bacteria_species_patho_total, N_bacteria_species_nonpatho, N_bacteria_species_nonpatho_total), by = "pfam"]
setkeyv(pfam2proteome_bacteria, "pfam")

###########################################################################################
### Domain's propensity for mediating eukaryote-endogenous vs. bacteria-endogenous DDIs ###
###########################################################################################

### Domains with experimental evidence of mediating PPIs ###

load("~/database/3did_elm/ddmi2pdb2taxcat.RData")
load("~/database/intact/ppi_ddmi.RData")

# when resolving PPIs to DDIs, give highest confidence to interchain DDIs, derived from PDB structures consisting of at least two distinct protein entities.
ddi_dt = ddmi_dt[type == "DDI" & dmA %in% pfam2proteome_host_bacteria$pfam & dmB %in% pfam2proteome_host_bacteria$pfam & taxcatA %in% c("human", "animal", "plant", "fungi", "bacteria") & taxcatB %in% c("human", "animal", "plant", "fungi", "bacteria")]
ddi_dt_intertemp = ddi_dt[, all(paste(dmA, dmB, sep = "_") %in% names(ddi2pdb) | paste(dmB, dmA, sep = "_") %in% names(ddi2pdb)), by = c("cidA", "cidB")][V1 == TRUE, paste(cidA, cidB, sep = "_")]
ddi_dt_intratemp = ddi_dt[, !any(paste(dmA, dmB, sep = "_") %in% names(ddi2pdb) | paste(dmB, dmA, sep = "_") %in% names(ddi2pdb)), by = c("cidA", "cidB")][V1 == TRUE, paste(cidA, cidB, sep = "_")]
ddi_dt_mixedtemp = ddi_dt[, any(paste(dmA, dmB, sep = "_") %in% names(ddi2pdb) | paste(dmB, dmA, sep = "_") %in% names(ddi2pdb)) & any(!paste(dmA, dmB, sep = "_") %in% names(ddi2pdb) & !paste(dmB, dmA, sep = "_") %in% names(ddi2pdb)), by = c("cidA", "cidB")][V1 == TRUE, paste(cidA, cidB, sep = "_")]
ddi_dt = ddi_dt[paste(cidA, cidB, sep = "_") %in% ddi_dt_intertemp | (paste(cidA, cidB, sep = "_") %in% ddi_dt_intratemp & cidA != cidB) |
                  (paste(cidA, cidB, sep = "_") %in% ddi_dt_mixedtemp & (paste(dmA, dmB, sep = "_") %in% names(ddi2pdb) | paste(dmB, dmA, sep = "_") %in% names(ddi2pdb)))]
ppi2ddi = ddi_dt[, sort(unique(ifelse(dmA >= dmB, paste(dmA, dmB, sep = "_"), paste(dmB, dmA, sep = "_"))))]
rm(ppi, ppi_ddmi, ppi_motif, ppi_pfam, ddmi_dt)
gc(TRUE)

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

ddi_dt = unique(ddi_dt[, .(cidA, cidB, dmA, dmB, taxidA, taxidB, taxcatA, taxcatB, pmid)])
ddi_dt = ddi_dt[, list(taxidA = unlist(strsplit(taxidA, ";"))), by = setdiff(names(ddi_dt), "taxidA")]
ddi_dt = ddi_dt[, list(taxidB = unlist(strsplit(taxidB, ";"))), by = setdiff(names(ddi_dt), "taxidB")]
ddi_dt[, speciesA := rankedlineage[match(taxidA, taxid), species]]
ddi_dt[, speciesB := rankedlineage[match(taxidB, taxid), species]]
ddi_dt[, speciesidA := rankedlineage[match(speciesA, species), taxid]]
ddi_dt[, speciesidB := rankedlineage[match(speciesB, species), taxid]]
rm(nodes, rankedlineage)
gc(TRUE)

ddi_dt = unique(ddi_dt[!is.na(speciesA) & !is.na(speciesB), .(cidA, cidB, dmA, dmB, speciesA, speciesB, speciesidA, speciesidB, taxcatA, taxcatB, pmid)])
colsAB = c(grep("A$", names(ddi_dt), value = TRUE), grep("B$", names(ddi_dt), value = TRUE))
colsBA = c(grep("B$", names(ddi_dt), value = TRUE), grep("A$", names(ddi_dt), value = TRUE))
ddi_dt[dmA < dmB, c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ddi_dt[taxcatA == "bacteria" & taxcatB %in% c("human", "animal", "plant"), c(colsAB) := .SD[, colsBA, with = FALSE], .SDcols = colsAB]
ddi_dt = unique(ddi_dt)
ddi_hh = ddi_dt[taxcatA %in% c("human", "animal", "plant", "fungi") & taxcatB %in% c("human", "animal", "plant", "fungi")]
ddi_bb = ddi_dt[taxcatA == "bacteria" & taxcatB == "bacteria"]
ddi_bh = ddi_dt[taxcatA %in% c("human", "animal", "plant") & taxcatB == "bacteria"]
nonpatho_species = c("Bacillus amyloliquefaciens", "Bacillus licheniformis", "Bacillus pumilus", "Bacillus subtilis", "Escherichia coli", "Nostoc sp. PCC 7119", "Saccharopolyspora erythraea", "Streptomyces albogriseolus", "Streptomyces griseus", "Streptomyces tendae", "Synechococcus elongatus", "Synechocystis sp. PCC 6803", "Thermoactinomyces vulgaris", "Thermotoga maritima")
effector = rbindlist(list(fread("~/database/uniprot/uniprot_bacteria_effector_name.tab.gz", select = c("Entry", "Status", "Length", "Annotation")),
                          fread("~/database/uniprot/uniprot_bacteria_effector_location.tab.gz", select = c("Entry", "Status", "Length", "Annotation")),
                          fread("~/database/uniprot/uniprot_bacteria_effector_phibase.tab.gz", select = c("Entry", "Status", "Length", "Annotation"))))
ddi_bh = ddi_bh[cidB %in% effector$Entry | !speciesB %in% nonpatho_species]

# eukaryotic domains that mediate eukaryote-endogenous PPIs
pfam_ppi_hh = ddi_hh[, intersect(pfam2proteome_host$pfam, c(dmA, dmB))]
# bacterial domains that mediate bacteria-endogenous PPIs
pfam_ppi_bb = ddi_bb[, intersect(pfam2proteome_bacteria$pfam, c(dmA, dmB))]
# bacterial domains that mediate bacteria-host exogenous PPIs
pfam_ppi_bh = ddi_bh[, intersect(pfam2proteome_bacteria$pfam, dmB)]

### Domains with structural evidence of mediating DDIs ###

# eukaryotic domains with structural evidence of mediating eukaryote-endogenous DDIs
pfam_pdb_hh = intersect(pfam2proteome_host$pfam, unlist(strsplit(names(ddi2taxcat)[sapply(ddi2taxcat, function(x) any(x %in% "eukaryote"))], "_")))
pfam_pdb_hh_all = intersect(pfam2proteome_host$pfam, unlist(strsplit(names(ddi2taxcat_all)[sapply(ddi2taxcat_all, function(x) any(x %in% "eukaryote"))], "_")))
# bacterial domains with structural evidence of mediating bacteria-endogenous DDIs
pfam_pdb_bb = intersect(pfam2proteome_bacteria$pfam, unlist(strsplit(names(ddi2taxcat)[sapply(ddi2taxcat, function(x) any(x %in% "bacteria"))], "_")))
pfam_pdb_bb_all = intersect(pfam2proteome_bacteria$pfam, unlist(strsplit(names(ddi2taxcat_all)[sapply(ddi2taxcat_all, function(x) any(x %in% "bacteria"))], "_")))

load("~/database/3did_elm/ddi_dmi.RData")

host_bacteria_ddi2species = ddi_dmi[type == "DDI" & dmA >= dmB]
host_bacteria_ddi2species[, domain_pair := paste(dmA, dmB, sep = "_")]
host_bacteria_ddi2species[domain_pair %in% names(ddi2pdb_all), type := "intraprotein"]
host_bacteria_ddi2species[domain_pair %in% names(ddi2pdb), type := "interprotein"]
host_bacteria_ddi2species[domain_pair %in% ppi2ddi, type := "ppi"]
host_bacteria_ddi2species[, N_host_species_pfams_both := length(intersect(pfam2proteome_host[pfam == dmA, unlist(strsplit(host_speciesids, ";"))],
                                                                          pfam2proteome_host[pfam == dmB, unlist(strsplit(host_speciesids, ";"))])), by = "domain_pair"]
host_bacteria_ddi2species[, N_host_species_pfams_either := length(union(pfam2proteome_host[pfam == dmA, unlist(strsplit(host_speciesids, ";"))],
                                                                        pfam2proteome_host[pfam == dmB, unlist(strsplit(host_speciesids, ";"))])), by = "domain_pair"]
host_bacteria_ddi2species[, N_bacteria_species_pfams_both := length(intersect(pfam2proteome_bacteria[pfam == dmA, unlist(strsplit(bacteria_speciesids, ";"))],
                                                                              pfam2proteome_bacteria[pfam == dmB, unlist(strsplit(bacteria_speciesids, ";"))])), by = "domain_pair"]
host_bacteria_ddi2species[, N_bacteria_species_pfams_either := length(union(pfam2proteome_bacteria[pfam == dmA, unlist(strsplit(bacteria_speciesids, ";"))],
                                                                            pfam2proteome_bacteria[pfam == dmB, unlist(strsplit(bacteria_speciesids, ";"))])), by = "domain_pair"]
host_bacteria_ddi2species[, c("host2bacteria_ddi_or", "host2bacteria_ddi_weight") := or_weight(N_host_species_pfams_both, N_host_species_pfams_either, N_bacteria_species_pfams_both, N_bacteria_species_pfams_either), by = "domain_pair"]
setkeyv(host_bacteria_ddi2species, "domain_pair")

##########################################################################################
### DDIs in host-bacteria exogenous PPI network tend to arise via convergent evolution ###
##########################################################################################

ddi_hh_nonself = unique(ddi_hh[cidA != cidB, .(cidA, dmA, cidB, dmB)])

# Are similar domains used by two host proteins in binding to the same domain of another host protein?
same_domain_ddi_hh = function(x) {
  # proteins and domains targeted by protein A
  i = union(ddi_hh_nonself[cidA == x, paste(cidB, dmB, sep = "_")], ddi_hh_nonself[cidB == x, paste(cidA, dmA, sep = "_")])
  # domains used by protein A to target other proteins and domains
  j = union(ddi_hh_nonself[cidA == x, paste(dmA, cidB, dmB, sep = "_")], ddi_hh_nonself[cidB == x, paste(dmB, cidA, dmA, sep = "_")])
  # do other proteins target the same proteins and domains as protein A?
  ka = ddi_hh_nonself[cidA != x & paste(cidB, dmB, sep = "_") %in% i]
  kb = ddi_hh_nonself[cidB != x & paste(cidA, dmA, sep = "_") %in% i]
  # do other proteins use the same domains to target the same proteins and domains as protein A?
  if (nrow(ka) > 0) list(same_domain = TRUE, same_ddi = ka[, any(paste(dmA, cidB, dmB, sep = "_") %in% j)])
  else if (nrow(kb) > 0) list(same_domain = TRUE, same_ddi = kb[, any(paste(dmB, cidA, dmA, sep = "_") %in% j)])
  else list(same_domain = FALSE, same_ddi = NA)
}
sdd_hh_list = lapply(ddi_hh_nonself[, union(cidA, cidB)], same_domain_ddi_hh)
names(sdd_hh_list) = ddi_hh_nonself[, union(cidA, cidB)]
sd_hh = names(sdd_hh_list)[sapply(sdd_hh_list, function(x) x[["same_domain"]] == TRUE)]
sdd_hh = names(sdd_hh_list)[sapply(sdd_hh_list, function(x) x[["same_domain"]] == TRUE & x[["same_ddi"]] == TRUE)]

# Are similar domains used by a bacterial protein and a host protein in binding to the same domain of another host protein?
same_domain_ddi_bh = function(x) {
  # host proteins and domains targeted by bacterial protein x
  i = ddi_bh[cidB == x, unique(paste(cidA, dmA, sep = "_"))]
  # domains used by bacterial protein x to target host proteins and domains
  j = ddi_bh[cidB == x, unique(paste(dmB, cidA, dmA, sep = "_"))]
  # do other host proteins target the same proteins and domains as bacterial protein x?
  ka = ddi_hh_nonself[cidA != x & paste(cidB, dmB, sep = "_") %in% i]
  kb = ddi_hh_nonself[cidB != x & paste(cidA, dmA, sep = "_") %in% i]
  # do other host proteins use the same domains to target the same proteins and domains as bacterial protein x?
  if (nrow(ka) > 0) list(same_domain = TRUE, same_ddi = ka[, any(paste(dmA, cidB, dmB, sep = "_") %in% j)])
  else if (nrow(kb) > 0) list(same_domain = TRUE, same_ddi = kb[, any(paste(dmB, cidA, dmA, sep = "_") %in% j)])
  else list(same_domain = FALSE, same_ddi = NA)
}

sdd_bh_list = lapply(ddi_bh[, unique(cidB)], same_domain_ddi_bh)
names(sdd_bh_list) = ddi_bh[, unique(cidB)]
sd_bh = names(sdd_bh_list)[sapply(sdd_bh_list, function(x) x[["same_domain"]] == TRUE)]
sdd_bh = names(sdd_bh_list)[sapply(sdd_bh_list, function(x) x[["same_domain"]] == TRUE & x[["same_ddi"]] == TRUE)]

save(pfam2proteome_host, pfam2proteome_bacteria, pfam2proteome_host_bacteria, ddi_dt, ddi_hh, ddi_bb, ddi_bh, sd_hh, sdd_hh, sd_bh, sdd_bh,
     pfam_ppi_hh, pfam_ppi_bb, pfam_ppi_bh, pfam_pdb_hh, pfam_pdb_hh_all, pfam_pdb_bb, pfam_pdb_bb_all, pfam_pdb_bh_effector, pfam_pdb_bh_host,
     host_bacteria_ddi2species, file = "~/database/uniprot/host_bacteria_pfam_ddi2species.RData")
