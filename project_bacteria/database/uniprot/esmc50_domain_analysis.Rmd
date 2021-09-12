library(data.table)
library(ggplot2)
library(knitr)
library(stringr)

setwd("~/project_bacteria/")

load("database/uniprot/esmc50_domain.RData")
load("database/uniprot/host_bacteria_pfam_ddi2species.RData")

####################################################################################################
### Experimental evidence that pathogens target host domains involved in host-specific processes ###
####################################################################################################

div_evol = data.table(Domain = ddi_bh[dmB %in% pfam2proteome_host$pfam, sort(intersect(dmB, pfam_ppi_hh))])
div_evol[, ddi_hh_exc := paste0(sort(setdiff(ddi_hh[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")], ddi_bb[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")])), collapse = "; "), by = "Domain"]
div_evol[, ddi_hh_bb := paste0(sort(intersect(ddi_hh[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")], ddi_bb[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")])), collapse = "; "), by = "Domain"]
div_evol[, num_ddi_hh_exc := ifelse(nchar(ddi_hh_exc) == 0, 0, str_count(ddi_hh_exc, "; ") + 1), by = "Domain"]
div_evol[, num_ddi_hh_bb := ifelse(nchar(ddi_hh_bb) == 0, 0, str_count(ddi_hh_bb, "; ") + 1), by = "Domain"]
div_evol[, eukexc := num_ddi_hh_bb == 0]
div_evol[, euk2bac := (num_ddi_hh_exc + 0.001) / (num_ddi_hh_bb + 0.001)]
setorder(div_evol, -eukexc, -euk2bac)
div_evol = div_evol[, .(Domain, ddi_hh_exc, ddi_hh_bb, num_ddi_hh_exc, num_ddi_hh_bb)]
div_evol[, baceuk_ppi := ddi_bh[match(Domain, dmB), paste0(paste(cidA, cidB, sep = "_"), " (", pmid, ")")]]
div_evol[, baceuk_ppi := gsub("pubmed:", "", baceuk_ppi)]
setnames(div_evol, c("Horizontally-acquired host domain", "DDIs mapped to eukaryote-endogenous PPIs", "DDIs mapped to eukaryote-endogenous and bacteria-endogenous PPIs", "# DDIs mapped to eukaryote-endogenous PPIs", "# DDIs mapped to eukaryote-endogenous and bacteria-endogenous PPIs", "Host-pathogen PPIs (PMIDs)"))
write.csv(div_evol, "figs_tabs/Supplementary Table 1.csv", quote = FALSE, row.names = FALSE)

conv_evol = data.table(Domain = ddi_bh[!dmB %in% pfam2proteome_host$pfam, sort(intersect(dmA, pfam_ppi_hh))])
conv_evol[, ddi_hh_exc := paste0(sort(setdiff(ddi_hh[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")], ddi_bb[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")])), collapse = "; "), by = "Domain"]
conv_evol[, ddi_hh_bb := paste0(sort(intersect(ddi_hh[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")], ddi_bb[dmA %in% Domain | dmB %in% Domain, paste(dmA, dmB, sep = "_")])), collapse = "; "), by = "Domain"]
conv_evol[, num_ddi_hh_exc := ifelse(nchar(ddi_hh_exc) == 0, 0, str_count(ddi_hh_exc, "; ") + 1), by = "Domain"]
conv_evol[, num_ddi_hh_bb := ifelse(nchar(ddi_hh_bb) == 0, 0, str_count(ddi_hh_bb, "; ") + 1), by = "Domain"]
conv_evol[, eukexc := num_ddi_hh_bb == 0]
conv_evol[, euk2bac := (num_ddi_hh_exc + 0.001) / (num_ddi_hh_bb + 0.001)]
setorder(conv_evol, -eukexc, -euk2bac)
conv_evol = conv_evol[, .(Domain, ddi_hh_exc, ddi_hh_bb, num_ddi_hh_exc, num_ddi_hh_bb)]
conv_evol[, baceuk_ppi := ddi_bh[match(Domain, dmA), paste0(paste(cidA, cidB, sep = "_"), " (", pmid, ")")]]
conv_evol[, baceuk_ppi := gsub("pubmed:", "", baceuk_ppi)]
setnames(conv_evol, c("Convergently-targeted host domain", "DDIs mapped to eukaryote-endogenous PPIs", "DDIs mapped to eukaryote-endogenous and bacteria-endogenous PPIs", "# DDIs mapped to eukaryote-endogenous PPIs", "# DDIs mapped to eukaryote-endogenous and bacteria-endogenous PPIs", "Host-pathogen PPIs (PMIDs)"))
write.csv(conv_evol, "figs_tabs/Supplementary Table 2.csv", quote = FALSE, row.names = FALSE)

######################################
### Interacting eukaryotic domains ###
######################################

# bacterial protein contains eukaryotic domain
esmc50_domain[, N_hom := sum(unlist(strsplit(pfam, ";")) %in% pfam2proteome_host$pfam), by = "uniprotid"]
# eukaryotic domain mediates DDIs in eukaryotes
esmc50_domain[, N_hom_ddi_hh := sum(unlist(strsplit(pfam, ";")) %in% pfam_ppi_hh), by = "uniprotid"]
# eukaryotic domain mediates DDIs in eukaryotes only
pfam_ppi_hh_exc = setdiff(pfam_ppi_hh, c(pfam_ppi_bb, pfam_pdb_bb_all))
esmc50_domain[, N_hom_ddi_hh_exc := sum(unlist(strsplit(pfam, ";")) %in% pfam_ppi_hh_exc), by = "uniprotid"]
esmc50_domain[, pfam_hom_ddi_hh_exc := paste0(intersect(unlist(strsplit(pfam, ";")), pfam_ppi_hh_exc), collapse = ";"), by = "uniprotid"]

esmc50_domain[, table(N_hom == 0), by = "type"]
esmc50_domain[N_hom > 0, table(N_hom_ddi_hh == 0), by = "type"]
esmc50_domain[N_hom_ddi_hh > 0, table(N_hom_ddi_hh_exc == 0), by = "type"]
fisher.test(matrix(esmc50_domain[N_hom_ddi_hh > 0, table(N_hom_ddi_hh_exc == 0), by = "type"][, V1], ncol = 2))

ddi_hh_exc_tab = esmc50_domain[N_hom_ddi_hh > 0, as.integer(table(N_hom_ddi_hh_exc == 0)), by = "type"]
setnames(ddi_hh_exc_tab, "V1", "prop")
ddi_hh_exc_tab[, "contain eukaryotic-like domains that" := rep(c("mediate PPIs exclusively in eukaryotes", "also mediate DDIs in bacteria"), 2)]
ddi_hh_exc_tab[, lab_pos := cumsum(prop) - 0.5 * prop, by = "type"]

ggplot(ddi_hh_exc_tab[type == "Effector"], aes(x = 2, y = prop, fill = `contain eukaryotic-like domains that`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab_pos, label = prop), color = "white", size = 10) +
  scale_fill_manual(values = c("cyan", "blue")) +
  theme_void() +
  xlim(0.5, 2.5)

ggplot(ddi_hh_exc_tab[type == "Non-effector"], aes(x = 2, y = prop, fill = `contain eukaryotic-like domains that`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab_pos, label = prop), color = "white", size = 10) +
  scale_fill_manual(values = c("cyan", "blue")) +
  theme_void() +
  xlim(0.5, 2.5)

nodes = fread("database/taxonomy/nodes.tab", header = FALSE, colClasses = "character")
names(nodes)[1:3] = c("taxid", "parenttaxid", "taxrank")
rankedlineage = fread("database/taxonomy/rankedlineage.tab", header = FALSE, colClasses = "character")
setnames(rankedlineage, c("taxid", "taxname", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"))
rankedlineage[, taxrank := nodes[match(rankedlineage$taxid, taxid), taxrank]]
rankedlineage[nchar(species) == 0 & taxrank == "species", species := taxname]
rankedlineage[nchar(genus) == 0 & taxrank == "genus", genus := taxname]
rankedlineage[nchar(family) == 0 & taxrank == "family", family := taxname]
rankedlineage[nchar(order) == 0 & taxrank == "order", order := taxname]
rankedlineage[nchar(class) == 0 & taxrank == "class", class := taxname]
rankedlineage[nchar(phylum) == 0 & taxrank == "phylum", phylum := taxname]

esmc50_domain_speciesids = rankedlineage[species %in% esmc50_domain$species & taxrank == "species", unique(taxid)]

effector_hom_ddi_hh_exc = esmc50_domain[type == "Effector" & N_hom_ddi_hh_exc > 0, .(uniprotac, species, pfam, pfam_hom_ddi_hh_exc)]
# effector domains mediating PPIs exclusively in eukaryotes in host-endogenous PPIs, also identified from experimental host-pathogen PPIs
esmc50_domain[type == "Effector" & N_hom_ddi_hh_exc > 0, intersect(unlist(strsplit(pfam_hom_ddi_hh_exc, ";")), div_evol$`Horizontally-acquired host domain`)]
# effector domains mediating PPIs exclusively in eukaryotes in host-endogenous PPIs, but not identified from experimental host-pathogen PPIs
esmc50_domain[type == "Effector" & N_hom_ddi_hh_exc > 0, setdiff(unlist(strsplit(pfam_hom_ddi_hh_exc, ";")), div_evol$`Horizontally-acquired host domain`)]
effector_hom_ddi_hh_exc[, pfam := gsub(";", "; ", pfam)]
effector_hom_ddi_hh_exc[, pfam_speciesids_patho := {x = unlist(strsplit(pfam, "; ")); paste0(Reduce(intersect, lapply(x, function(y) pfam2proteome_bacteria[pfam == y, unlist(strsplit(bacteria_speciesids_patho, ";"))])), collapse = "; ")}, by = "pfam"]
effector_hom_ddi_hh_exc[, pfam_speciesids_patho := paste0(intersect(esmc50_domain_speciesids, unlist(strsplit(pfam_speciesids_patho, "; "))), collapse = "; "), by = "pfam"]
effector_hom_ddi_hh_exc[, N_pfam_speciesids_patho := str_count(pfam_speciesids_patho, ";") + 1]
effector_hom_ddi_hh_exc[, pfam_speciesids_patho := NULL]
effector_hom_ddi_hh_exc[, pfam_hom_ddi_hh_exc := gsub(";", "; ", pfam_hom_ddi_hh_exc)]
setorder(effector_hom_ddi_hh_exc, pfam_hom_ddi_hh_exc)
setcolorder(effector_hom_ddi_hh_exc, c("uniprotac", "species", "pfam", "N_pfam_speciesids_patho", "pfam_hom_ddi_hh_exc"))
setnames(effector_hom_ddi_hh_exc, c("UniProt Accession", "Representative species", "Domains", "# Pathogenic spp. encoding proteins with domain(s)", "Domains mediating PPIs exclusively in eukaryotes"))
write.csv(effector_hom_ddi_hh_exc, "figs_tabs/Table 1.csv", quote = FALSE, row.names = FALSE)
kable(effector_hom_ddi_hh_exc, "html")

# divergent evolution in host-endogenous PPI network vs. convergent evolution in host-bacteria PPI network
domain_usage = ddi_hh_exc_tab
domain_usage[, type := rep(c("host-endogenous PPI network", "host-bacteria PPI network"), each = 2)]
domain_usage[, prop := c(length(sdd_hh), length(sd_hh) - length(sdd_hh), length(sdd_bh), length(sd_bh) - length(sdd_bh))]
names(domain_usage)[3] = "domains used by two proteins to bind\nto the same domain of a common target"
domain_usage[, "domains used by two proteins to bind\nto the same domain of a common target" := rep(c("same type", "different types"), length.out = nrow(domain_usage))]
domain_usage[, lab_pos := cumsum(prop) - 0.5 * prop, by = "type"]

ggplot(domain_usage[type == "host-endogenous PPI network"], aes(x = 2, y = prop, fill = `domains used by two proteins to bind\nto the same domain of a common target`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab_pos, label = prop), color = "white", size = 10) +
  scale_fill_manual(values = c("gray", "black")) +
  theme_void() +
  xlim(0.5, 2.5)

ggplot(domain_usage[type == "host-bacteria PPI network"], aes(x = 2, y = prop, fill = `domains used by two proteins to bind\nto the same domain of a common target`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab_pos, label = prop), color = "white", size = 10) +
  scale_fill_manual(values = c("gray", "black")) +
  theme_void() +
  xlim(0.5, 2.5)

esmc_mean_h2b_ddi_hh = esmc50_domain[N_hom_ddi_hh > N_hom_ddi_hh_exc, {
  i = host_bacteria_ddi2species[dmA != dmB & type != "DDI", unique(unlist(lapply(intersect(unlist(strsplit(pfam, ";")), setdiff(pfam_ppi_hh, pfam_ppi_hh_exc)), grep, x = domain_pair)))]
  j = host_bacteria_ddi2species[dmA == dmB & type != "DDI", Reduce(intersect, list(unlist(strsplit(pfam, ";")), dmA, setdiff(pfam_ppi_hh, pfam_ppi_hh_exc)))]
  h2b = rbind(host_bacteria_ddi2species[dmA != dmB & type != "DDI"][i, .(host2bacteria_ddi_or, host2bacteria_ddi_weight)],
              pfam2proteome_host_bacteria[pfam %in% j, .(host2bacteria_species_or, host2bacteria_species_weight)], use.names = FALSE)
  wmlor = h2b[, log(sum(host2bacteria_ddi_or * host2bacteria_ddi_weight) / sum(host2bacteria_ddi_weight))]
  list(wmlor = wmlor)}, by = "uniprotid"]
esmc_mean_h2b_ddi_hh[, wmlor := scale(wmlor, scale = FALSE)]
esmc50_domain[, mean_h2b_ddi_hh := esmc_mean_h2b_ddi_hh[match(esmc50_domain$uniprotid, uniprotid), wmlor]]

e1 = unique(esmc50_domain[N_hom_ddi_hh > N_hom_ddi_hh_exc & !is.na(mean_h2b_ddi_hh)], by = "mean_h2b_ddi_hh")
e1[, table(type)]
ggplot(e1, aes(x = type, y = mean_h2b_ddi_hh, fill = type)) + geom_violin(trim = TRUE) + geom_boxplot(width = 0.1) + labs(title = "Bacterial proteins containing domains that can mediate DDIs\nin both eukaryotes and bacteria", x = "Type of bacterial protein") + ylab("Mean log-odds ratio of component domain(s) co-occurring with interacting domains\nin eukaryotes vs. in bacteria")
e1[, mean(mean_h2b_ddi_hh), by = "type"]
e1[, wilcox.test(mean_h2b_ddi_hh ~ type)]

effector_hom_ddi_hh_pri = e1[type == "Effector"]
effector_hom_ddi_hh_pri[, pfam_hom_ddi_hh_bb := paste0(intersect(unlist(strsplit(pfam, ";")), setdiff(pfam_ppi_hh, pfam_ppi_hh_exc)), collapse = "; "), by = "uniprotid"]
effector_hom_ddi_hh_pri = effector_hom_ddi_hh_pri[, .(uniprotac, species, pfam, pfam_hom_ddi_hh_bb, mean_h2b_ddi_hh)]
effector_hom_ddi_hh_pri[, pfam := gsub(";", "; ", pfam)]
effector_hom_ddi_hh_pri[, pfam_speciesids_patho := {x = unlist(strsplit(pfam, "; ")); paste0(Reduce(intersect, lapply(x, function(y) pfam2proteome_bacteria[pfam == y, unlist(strsplit(bacteria_speciesids_patho, ";"))])), collapse = "; ")}, by = "pfam"]
effector_hom_ddi_hh_pri[, pfam_speciesids_patho := paste0(intersect(esmc50_domain_speciesids, unlist(strsplit(pfam_speciesids_patho, "; "))), collapse = "; "), by = "pfam"]
effector_hom_ddi_hh_pri[, N_pfam_speciesids_patho := str_count(pfam_speciesids_patho, ";") + 1]
effector_hom_ddi_hh_pri[, pfam_speciesids_patho := NULL]
setorder(effector_hom_ddi_hh_pri, -mean_h2b_ddi_hh)
effector_hom_ddi_hh_pri = effector_hom_ddi_hh_pri[1:10]
effector_hom_ddi_hh_pri[, mean_h2b_ddi_hh := round(mean_h2b_ddi_hh, 1)]
setcolorder(effector_hom_ddi_hh_pri, c("uniprotac", "species", "pfam", "pfam_hom_ddi_hh_bb", "N_pfam_speciesids_patho", "mean_h2b_ddi_hh"))
setnames(effector_hom_ddi_hh_pri, c("UniProt Accession", "Species", "Domains", "Domains with DDI partners in both eukaryotes and bacteria", "# Pathogenic spp. encoding proteins with domain(s)", "Log odds ratio of domains co-occurring with DDI partners in eukaryotes vs. in bacteria"))
write.csv(effector_hom_ddi_hh_pri, "figs_tabs/Table 3.csv", quote = FALSE, row.names = FALSE)
kable(effector_hom_ddi_hh_pri, "html")

##########################################
### Non-interacting eukaryotic domains ###
##########################################

pfam_ddi_hh_all = host_bacteria_ddi2species[N_host_species_pfams_both > 3, union(dmA, dmB)]
esmc50_domain[N_hom > N_hom_ddi_hh, mean_h2b_spp_noddi_hh := pfam2proteome_host_bacteria[pfam %in% .SD[, intersect(unlist(strsplit(pfam, ";")), setdiff(pfam2proteome_host$pfam, pfam_ddi_hh_all))], log(sum(host2bacteria_species_or * host2bacteria_species_weight) / sum(host2bacteria_species_weight))], by = "uniprotid"]
esmc50_domain[N_hom > N_hom_ddi_hh, mean_h2b_spp_noddi_hh := scale(mean_h2b_spp_noddi_hh, scale = FALSE)]
esmc50_domain[N_hom > N_hom_ddi_hh & !is.na(mean_h2b_spp_noddi_hh), table(type)]

ggplot(esmc50_domain[N_hom > N_hom_ddi_hh & !is.na(mean_h2b_spp_noddi_hh)], aes(x = type, y = mean_h2b_spp_noddi_hh, fill = type)) + geom_violin(trim = TRUE) + geom_boxplot(width = 0.1) + labs(title = "Bacterial proteins containing eukaryotic-like domains\nthat are unlikely to mediate DDIs in eukaryotes", x = "Type of bacterial protein") + ylab("Mean log-odds ratio of finding domain in host vs. bacteria")
esmc50_domain[N_hom > N_hom_ddi_hh & !is.na(mean_h2b_spp_noddi_hh), mean(mean_h2b_spp_noddi_hh), by = "type"]
esmc50_domain[N_hom > N_hom_ddi_hh & !is.na(mean_h2b_spp_noddi_hh), wilcox.test(mean_h2b_spp_noddi_hh ~ type)]

##############################################################################################################
### Bacteria-exclusive domain hijacks eukaryotic domains otherwise exclusively involved in eukaryotic DDIs ###
##############################################################################################################

pfam_bacexc = setdiff(pfam2proteome_bacteria$pfam, pfam2proteome_host$pfam)
pfam_bacexc_ppi_bh = ddi_bh[dmA %in% pfam_ppi_hh, intersect(pfam_bacexc, dmB)]
pfam_bacexc_ppi_bh_hh_exc = ddi_bh[dmA %in% pfam_ppi_hh_exc, intersect(pfam_bacexc, dmB)]

# all prokaryotic domains with the potential to target eukaryotic domains
pfam_bacexc_ddi_bh = union(host_bacteria_ddi2species[type %in% c("ppi", "interprotein") & dmA %in% pfam_bacexc & dmB %in% pfam_ppi_hh, dmA],
                           host_bacteria_ddi2species[type %in% c("ppi", "interprotein") & dmB %in% pfam_bacexc & dmA %in% pfam_ppi_hh, dmB])

pfam_bacexc_ddi_bh_hh_exc = union(host_bacteria_ddi2species[type %in% c("ppi", "interprotein") & dmA %in% pfam_bacexc & dmB %in% pfam_ppi_hh_exc, dmA],
                                  host_bacteria_ddi2species[type %in% c("ppi", "interprotein") & dmB %in% pfam_bacexc & dmA %in% pfam_ppi_hh_exc, dmB])

esmc50_domain[N_hom == 0, N_bacexc_ddi_bh := sum(unlist(strsplit(pfam, ";")) %in% c(pfam_bacexc_ppi_bh, pfam_bacexc_ddi_bh)), by = "uniprotid"]
esmc50_domain[N_hom == 0 & N_bacexc_ddi_bh > 0, N_bacexc_ddi_bh_hh_exc := sum(unlist(strsplit(pfam, ";")) %in% c(pfam_bacexc_ppi_bh_hh_exc, pfam_bacexc_ddi_bh_hh_exc)), by = "uniprotid"]

esmc50_domain[N_hom == 0, table(N_bacexc_ddi_bh == 0), by = "type"]
esmc50_domain[N_hom == 0 & N_bacexc_ddi_bh > 0, table(N_bacexc_ddi_bh_hh_exc == 0), by = "type"]

effector_conv_ddi_hh_exc = esmc50_domain[N_hom == 0 & N_bacexc_ddi_bh_hh_exc > 0 & type == "Effector", .(uniprotid, species, pfam)]
effector_conv_ddi_hh_exc[, pb := paste0(sort(intersect(unlist(strsplit(pfam, ";")), c(pfam_bacexc_ppi_bh_hh_exc, pfam_bacexc_ddi_bh_hh_exc))), collapse = "; "), by = "uniprotid"]
effector_conv_ddi_hh_exc[, ph := paste0(sort(union(host_bacteria_ddi2species[type %in% c("ppi", "interprotein") & dmA %in% unlist(strsplit(pb, "; ")) & dmB %in% pfam_ppi_hh_exc, dmB],
                                                   host_bacteria_ddi2species[type %in% c("ppi", "interprotein") & dmB %in% unlist(strsplit(pb, "; ")) & dmA %in% pfam_ppi_hh_exc, dmA])), collapse = "; "), by = "uniprotid"]
effector_conv_ddi_hh_exc[, pfam := gsub(";", "; ", pfam)]
setorder(effector_conv_ddi_hh_exc, ph, pb)
setnames(effector_conv_ddi_hh_exc, c("UniProt Accession", "Species", "Domains", "Pathogen domains convergently targeting host domains", "Targeted host domains involved in eukaryote-specific PPIs"))
write.csv(effector_conv_ddi_hh_exc, "figs_tabs/Supplementary Table 3.csv", quote = FALSE, row.names = FALSE)

pfam_bacexc_ddi_bh_tab = esmc50_domain[N_hom == 0 & N_bacexc_ddi_bh > 0, as.integer(table(N_bacexc_ddi_bh_hh_exc == 0)), by = "type"]
setnames(pfam_bacexc_ddi_bh_tab, "V1", "prop")
pfam_bacexc_ddi_bh_tab[, "contain bacteria-exclusive domains\ntargeting eukaryotic domains that" := rep(c("mediate PPIs exclusively in eukaryotes", "also mediate DDIs in bacteria"), 2)]
pfam_bacexc_ddi_bh_tab[, lab_pos := cumsum(prop) - 0.5 * prop, by = "type"]

ggplot(pfam_bacexc_ddi_bh_tab[type == "Effector"], aes(x = 2, y = prop, fill = `contain bacteria-exclusive domains\ntargeting eukaryotic domains that`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab_pos, label = prop), color = "white", size = 10) +
  scale_fill_manual(values = c("pink", "red")) +
  theme_void() +
  xlim(0.5, 2.5)

ggplot(pfam_bacexc_ddi_bh_tab[type == "Non-effector"], aes(x = 2, y = prop, fill = `contain bacteria-exclusive domains\ntargeting eukaryotic domains that`)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab_pos, label = prop), color = "white", size = 10) +
  scale_fill_manual(values = c("pink", "red")) +
  theme_void() +
  xlim(0.5, 2.5)
