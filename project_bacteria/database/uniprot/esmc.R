#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=16G
#SBATCH --job-name=esmc
#SBATCH --output=%x-%j.out

library(data.table)

effector = rbindlist(list(fread("~/database/uniprot/uniprot_bacteria_effector_name.tab.gz", select = c("Entry", "Status", "Length", "Annotation")),
                          fread("~/database/uniprot/uniprot_bacteria_effector_location.tab.gz", select = c("Entry", "Status", "Length", "Annotation")),
                          fread("~/database/uniprot/uniprot_bacteria_effector_phibase.tab.gz", select = c("Entry", "Status", "Length", "Annotation"))))
secreted = fread("~/database/uniprot/uniprot_bacteria_secreted.tab.gz", select = c("Entry", "Status", "Length", "Annotation"))
membrane = fread("~/database/uniprot/uniprot_bacteria_membrane.tab.gz", select = c("Entry", "Status", "Length", "Annotation"))
cytoplasm = fread("~/database/uniprot/uniprot_bacteria_cytoplasm.tab.gz", select = c("Entry", "Status", "Length", "Annotation"))

esmc = fread("~/database/uniprot/esmc_uniref_uniparc_taxid.tab", header = FALSE, colClasses = "character")
setnames(esmc, c("uniprotac", "uniprotid", "uniref100", "uniref90", "uniref50", "uniparc", "taxid"))
esmc[uniprotac %in% cytoplasm$Entry, type := "Cytoplasm"]
esmc[uniprotac %in% secreted$Entry, type := "Secreted non-effector"]
esmc[uniprotac %in% membrane$Entry, type := "Membrane"]
esmc[uniprotac %in% effector$Entry, type := "Effector"]
esmc = esmc[!is.na(type)]
esmc[, uniprotac_eq_uniref50 := (paste0("UniRef50_", uniprotac) == uniref50)]
esmc[, uniprotac_eq_uniref90 := (paste0("UniRef90_", uniprotac) == uniref90)]
esmc[type == "Cytoplasm", status := cytoplasm[match(uniprotac, Entry), Status]]
esmc[type == "Secreted non-effector", status := secreted[match(uniprotac, Entry), Status]]
esmc[type == "Membrane", status := membrane[match(uniprotac, Entry), Status]]
esmc[type == "Effector", status := effector[match(uniprotac, Entry), Status]]
esmc[type == "Cytoplasm", annotation := cytoplasm[match(uniprotac, Entry), Annotation]]
esmc[type == "Secreted non-effector", annotation := secreted[match(uniprotac, Entry), Annotation]]
esmc[type == "Membrane", annotation := membrane[match(uniprotac, Entry), Annotation]]
esmc[type == "Effector", annotation := effector[match(uniprotac, Entry), Annotation]]
esmc[type == "Cytoplasm", length := cytoplasm[match(uniprotac, Entry), Length]]
esmc[type == "Secreted non-effector", length := secreted[match(uniprotac, Entry), Length]]
esmc[type == "Membrane", length := membrane[match(uniprotac, Entry), Length]]
esmc[type == "Effector", length := effector[match(uniprotac, Entry), Length]]
esmc[type != "Effector", type := "Non-effector"]
esmc[, type := factor(type, levels = c("Non-effector", "Effector"))]

refnrprot = fread("~/database/uniprot/uniprot_refnr_proteomes.tab")
ppmembership = fread("~/database/uniprot/PPMembership.txt")
esmc_proteome = fread("~/database/uniprot/esmc_proteome.tab", header = FALSE)
setnames(esmc_proteome, c("uniprotid", "proteome"))
esmc_proteome[, refnr_proteome := paste0(intersect(refnrprot$`Proteome ID`, unlist(strsplit(proteome, ";"))), collapse = ";"), by = "proteome"]
esmc_proteome[, pan_proteome := ppmembership[PPMember %in% unlist(strsplit(proteome, ";")), paste0(unique(PP), collapse = ";")], by = "proteome"]
esmc = merge(esmc, esmc_proteome, by = "uniprotid", all = TRUE)
esmc[nchar(refnr_proteome) > 0 & nchar(pan_proteome) > 0, proteome_quality := 4L]
esmc[nchar(refnr_proteome) > 0 & nchar(pan_proteome) == 0, proteome_quality := 3L]
esmc[nchar(refnr_proteome) == 0 & nchar(pan_proteome) > 0, proteome_quality := 2L]
esmc[nchar(refnr_proteome) == 0 & nchar(pan_proteome) == 0, proteome_quality := 1L]
esmc[is.na(proteome), proteome_quality := 0L]

rm(cytoplasm, effector, secreted, membrane, refnrprot, ppmembership)
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

esmc[, taxname := rankedlineage[match(esmc$taxid, taxid), taxname]]
esmc[, species := rankedlineage[match(esmc$taxid, taxid), species]]
esmc[, genus := rankedlineage[match(esmc$taxid, taxid), genus]]

rm(nodes, rankedlineage)
gc(TRUE)

load("~/database/uniprot/pfam_host_other.RData")

esmc_pfam = fread("~/database/uniprot/esmc_uniprotid2pfam.tab", header = FALSE)
setnames(esmc_pfam, c("uniprotid", "pfam"))
esmc_pfam[, pfam := paste0(sort(unique(unlist(strsplit(pfam, ";")))), collapse = ";"), by = "pfam"]
esmc_pfam = esmc_pfam[!pfam %in% esmc_pfam[, any(!unlist(strsplit(pfam, ";")) %in% pfam2proteome_bacteria$pfam), by = "pfam"][V1 == TRUE, pfam]]

esmc[, pfam := esmc_pfam[match(esmc$uniprotid, uniprotid), pfam]]
esmc[!is.na(pfam), N_pfams := uniqueN(unlist(strsplit(pfam, ";"))), by = "uniprotid"]
esmc[is.na(pfam), N_pfams := 0]
setorder(esmc, -type, -uniprotac_eq_uniref50, -uniprotac_eq_uniref90, -annotation, status, -proteome_quality, -N_pfams, -length, na.last = TRUE)
save(esmc, file = "~/database/uniprot/esmc.RObject")
