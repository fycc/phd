#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --time=3:00:00
#SBATCH --mem=16G
#SBATCH --job-name=esmc_elm
#SBATCH --output=%x-%j.out

library(data.table)
library(stringr)
library(stringdist)

load("~/database/uniprot/esmc.RObject")
load("~/database/3did_elm/ddi_dmi.RData")

esmc[uniref50 %in% esmc[, any(type == "Effector"), by = "uniref50"][V1 == TRUE, uniref50], type := "Effector"]
setorder(esmc, -type, -uniprotac_eq_uniref50, -uniprotac_eq_uniref90, -annotation, status, -proteome_quality, -N_pfams, -length, na.last = TRUE)
esmc = unique(esmc, by = "uniparc")

esmc_elm = fread("~/database/uniprot/esmc_motif.csv", sep = ",", header = FALSE)
setnames(esmc_elm, c("uniprotid", "length", "motifs"))
esmc_elm = esmc_elm[uniprotid %in% esmc$uniprotid]
esmc_elm[, pos_in_esmc := esmc[, match(esmc_elm$uniprotid, uniprotid)]]
setorder(esmc_elm, pos_in_esmc)
esmc_elm[, pos_in_esmc := NULL]
esmc_elm[, c("type", "species", "pfam", "uniref50", "uniref90", "uniparc", "uniprotac_eq_uniref50", "uniprotac_eq_uniref90", "annotation", "status", "proteome_quality") := esmc[match(esmc_elm$uniprotid, uniprotid), .(type, species, pfam, uniref50, uniref90, uniparc, uniprotac_eq_uniref50, uniprotac_eq_uniref90, annotation, status, proteome_quality)]]
rm(esmc)
gc(TRUE)

# ELMs merged if targeting the same Pfam domains
elm_domains = lapply(elm_motifs, function(x) ddi_dmi[type == "DMI" & dmB == x, gsub("PF", "", sort(unique(dmA)))])
names(elm_domains) = elm_motifs
elm_domains_dist = seq_distmatrix(elm_domains, elm_domains, method = "jaccard")
elm_clusters = cutree(hclust(as.dist(elm_domains_dist), method = "single"), k = uniqueN(unlist(lapply(elm_motifs, function(x) ddi_dmi[type == "DMI" & dmB == x, paste0(sort(unique(dmA)), collapse = "_")]))))
elm_clusters = sapply(elm_clusters, function(x) sprintf("c%05d", as.integer(x)))
names(elm_clusters) = elm_motifs

esmc_elm = esmc_elm[species %in% esmc_elm[, any(type == "Effector") & any(type == "Non-effector"), by = "species"][V1 == TRUE, species]]
esmc_elm[, motifs := paste0(unlist(strsplit(motifs, ";"))[match(elm_motifs, all_motifs)], collapse = ";"), by = "motifs"]
esmc_elm[, N_motifs := sum(as.integer(unlist(strsplit(motifs, ";")))), by = "motifs"]
esmc_elm[, density_motifs := N_motifs / length]
esmc_elm[N_motifs == 0, types_motifs := ""]
esmc_elm[N_motifs > 0, types_motifs := paste0(mapply(sprintf, "m%05d", which(unlist(strsplit(motifs, ";")) != "0")), collapse = ";"), by = "motifs"]
esmc_elm[, N_types_motifs := uniqueN(unlist(strsplit(types_motifs, ";"))), by = "types_motifs"]
esmc_elm[, N_motifs_merged := {
  i = as.integer(unlist(strsplit(motifs, ";")))
  i = i[i > 0]
  j = elm_clusters[unlist(strsplit(motifs, ";")) != "0"]
  round(sum(by(i, j, median)))
}, by = "motifs"]
esmc_elm[, density_motifs_merged := N_motifs_merged / length]
esmc_elm[, types_motifs_merged := paste0(sort(unique(elm_clusters[unlist(strsplit(motifs, ";")) != "0"])), collapse = ";"), by = "motifs"]
esmc_elm[, N_types_motifs_merged := uniqueN(unlist(strsplit(types_motifs_merged, ";"))), by = "types_motifs_merged"]

save(esmc_elm, file = "~/database/uniprot/esmc_elm.RObject")
