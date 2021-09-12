#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=8G
#SBATCH --job-name=esmc_mobidblt_elm
#SBATCH --output=%x-%j.out

library(data.table)
library(stringr)
library(stringdist)

load("~/database/uniprot/esmc.RObject")
load("~/database/uniprot/esmc_mobidblt_motif.RObject")
load("~/database/3did_elm/ddi_dmi.RData")

esmc[uniref50 %in% esmc[, any(type == "Effector"), by = "uniref50"][V1 == TRUE, uniref50], type := "Effector"]
setorder(esmc, -type, -uniprotac_eq_uniref50, -uniprotac_eq_uniref90, -annotation, status, -proteome_quality, -N_pfams, -length, na.last = TRUE)
esmc = unique(esmc, by = "uniparc")

esmc_mobidblt_elm = esmc_mobidblt_motif[uniprotid %in% esmc$uniprotid]
esmc_mobidblt_elm[, pos_in_esmc := esmc[, match(esmc_mobidblt_elm$uniprotid, uniprotid)]]
setorder(esmc_mobidblt_elm, pos_in_esmc)
esmc_mobidblt_elm[, pos_in_esmc := NULL]
esmc_mobidblt_elm[, c("type", "species", "pfam", "uniref50", "uniref90", "uniparc", "uniprotac_eq_uniref50", "uniprotac_eq_uniref90", "annotation", "status", "proteome_quality") := esmc[match(esmc_mobidblt_elm$uniprotid, uniprotid), .(type, species, pfam, uniref50, uniref90, uniparc, uniprotac_eq_uniref50, uniprotac_eq_uniref90, annotation, status, proteome_quality)]]
rm(esmc)
gc(TRUE)

# ELMs merged if targeting the same Pfam domains
elm_domains = lapply(elm_motifs, function(x) ddi_dmi[type == "DMI" & dmB == x, gsub("PF", "", sort(unique(dmA)))])
names(elm_domains) = elm_motifs
elm_domains_dist = seq_distmatrix(elm_domains, elm_domains, method = "jaccard")
elm_clusters = cutree(hclust(as.dist(elm_domains_dist), method = "single"), k = uniqueN(unlist(lapply(elm_motifs, function(x) ddi_dmi[type == "DMI" & dmB == x, paste0(sort(unique(dmA)), collapse = "_")]))))
elm_clusters = sapply(elm_clusters, function(x) sprintf("c%05d", as.integer(x)))
names(elm_clusters) = elm_motifs

esmc_mobidblt_elm = esmc_mobidblt_elm[species %in% esmc_mobidblt_elm[, any(type == "Effector") & any(type == "Non-effector"), by = "species"][V1 == TRUE, species]]
esmc_mobidblt_elm[, mobidblt_motifs := paste0(unlist(strsplit(mobidblt_motifs, ";"))[match(elm_motifs, all_motifs)], collapse = ";"), by = "mobidblt_motifs"]
esmc_mobidblt_elm[, N_mobidblt_motifs := sum(as.double(unlist(strsplit(mobidblt_motifs, ";")))), by = "mobidblt_motifs"]
esmc_mobidblt_elm[, density_mobidblt_motifs := N_mobidblt_motifs / mobidblt_residues]
esmc_mobidblt_elm[N_mobidblt_motifs == 0, types_mobidblt_motifs := ""]
esmc_mobidblt_elm[N_mobidblt_motifs > 0, types_mobidblt_motifs := paste0(mapply(sprintf, "m%05d", which(unlist(strsplit(mobidblt_motifs, ";")) != "0")), collapse = ";"), by = "mobidblt_motifs"]
esmc_mobidblt_elm[, N_types_mobidblt_motifs := uniqueN(unlist(strsplit(types_mobidblt_motifs, ";"))), by = "types_mobidblt_motifs"]
esmc_mobidblt_elm[, N_mobidblt_motifs_merged := {
  i = as.double(unlist(strsplit(mobidblt_motifs, ";")))
  i = i[i > 0]
  j = elm_clusters[unlist(strsplit(mobidblt_motifs, ";")) != "0"]
  round(sum(by(i, j, median)))
}, by = "mobidblt_motifs"]
esmc_mobidblt_elm[, density_mobidblt_motifs_merged := N_mobidblt_motifs_merged / mobidblt_residues]
esmc_mobidblt_elm[, types_mobidblt_motifs_merged := paste0(sort(unique(elm_clusters[unlist(strsplit(mobidblt_motifs, ";")) != "0"])), collapse = ";"), by = "mobidblt_motifs"]
esmc_mobidblt_elm[, N_types_mobidblt_motifs_merged := uniqueN(unlist(strsplit(types_mobidblt_motifs_merged, ";"))), by = "types_mobidblt_motifs_merged"]

save(esmc_mobidblt_elm, file = "~/database/uniprot/esmc_mobidblt_elm.RObject")
