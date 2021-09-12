#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --mem=8G
#SBATCH --job-name=esmc50_domain
#SBATCH --output=%x-%j.out

library(data.table)
library(stringr)

load("~/database/uniprot/esmc.RObject")
load("~/database/uniprot/pfam_host_other.RData")

# for domain signatures that are common to effectors and non-effectors
# find out if they are enriched in effectors or non-effectors
esmc[uniref50 %in% esmc[, any(type == "Effector"), by = "uniref50"][V1 == TRUE, uniref50], type := "Effector"]
esmc = esmc[!is.na(pfam)]
esmc = esmc[!pfam %in% esmc[, any(!unlist(strsplit(pfam, ";")) %in% pfam2proteome_bacteria$pfam), by = "pfam"][V1 == TRUE, pfam]]
esmc = esmc[species %in% esmc[, any(type == "Effector") & any(type == "Non-effector"), by = "species"][V1 == TRUE, species]]
esmc = unique(esmc, by = c("uniref50", "pfam"))
pfam_common = esmc[, any(type == "Effector") & any(type == "Non-effector"), by = "pfam"][V1 == TRUE, pfam]
N_eff_total = esmc[, table(type)[["Effector"]]]
N_noneff_total = esmc[, table(type)[["Non-effector"]]]
pfam_common_freq = as.data.table(t(sapply(pfam_common, function(x) {
  if (grepl(";", x)) esmc[Reduce(intersect, lapply(unlist(strsplit(x, ";")), function(y) esmc[, grep(y, pfam, fixed = TRUE)])), table(type)]
  else esmc[grep(x, pfam, fixed = TRUE), table(type)]
})))
pfam_common_freq[, "Domain composition" := pfam_common]
pfam_common_freq[, c("Odds ratio", "p-value") := fisher.test(matrix(c(Effector, N_eff_total - Effector, `Non-effector`, N_noneff_total - `Non-effector`), ncol = 2))[c("estimate", "p.value")], by = "pfam_common"]
pfam_common_freq[, "q-value" := p.adjust(`p-value`, method = "BH")]
pfam_common_freq[`Odds ratio` > 1 & `q-value` < 0.1, Enrichment := "Enriched in effectors"]
pfam_common_freq[`Odds ratio` < 1 & `q-value` < 0.1, Enrichment := "Depleted in effectors"]
pfam_common_freq[`q-value` > 0.1, Enrichment := "Statistically insignificant"]
pfam_common_freq[, Enrichment := factor(Enrichment, levels = c("Statistically insignificant", "Depleted in effectors", "Enriched in effectors"))]
setorderv(pfam_common_freq, c("Enrichment", "Odds ratio"), order = c(-1, -1), na.last = TRUE)
pfam_common_freq[, c("Odds ratio", "p-value", "q-value") := lapply(.SD, signif, 2), .SDcols = c("Odds ratio", "p-value", "q-value")]
write.table(pfam_common_freq, file = "~/database/uniprot/esmc50_pfam_common_freq.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

esmc = esmc[!pfam %in% pfam_common_freq[`q-value` > 0.1, `Domain composition`]]
esmc = esmc[!(pfam %in% pfam_common_freq[`Odds ratio` > 1 & `q-value` < 0.1, `Domain composition`] & type == "Non-effector")]
esmc = esmc[!(pfam %in% pfam_common_freq[`Odds ratio` < 1 & `q-value` < 0.1, `Domain composition`] & type == "Effector")]
esmc50_domain = unique(esmc, by = "pfam")
esmc50_domain = unique(esmc50_domain, by = "uniref50")

# GO slim functional analysis of domain signatures
load("~/database/3did_elm/pfam2goslim.RData")

esmc50_domain[, goslim_pir := goslimpir2pfam[PFAM %in% unlist(strsplit(pfam, ";")), paste0(unique(GOSLIM_ID), collapse = ";")], by = "uniprotid"]

N_domsigs_eff = esmc50_domain[nchar(goslim_pir) > 0, table(type)[["Effector"]]]
N_domsigs_noneff = esmc50_domain[nchar(goslim_pir) > 0, table(type)[["Non-effector"]]]
esmc50_domain_goslimpir = esmc50_domain[, unique(unlist(strsplit(goslim_pir, ";")))]

esmc50_domain_goslimpir_tab = as.data.table(t(sapply(esmc50_domain_goslimpir, function(x) esmc50_domain[grepl(x, goslim_pir, fixed = TRUE), table(type)])))
esmc50_domain_goslimpir_tab[, "GO ID" := esmc50_domain_goslimpir]
esmc50_domain_goslimpir_tab[, "GO Term" := str_to_sentence(goslimpir2pfam[match(`GO ID`, GOSLIM_ID), GOSLIM_TERM])]
goslim_pir_eff_exc = esmc50_domain_goslimpir_tab[`Non-effector` == 0 & Effector > 0, `GO ID`]
goslim_pir_noneff_exc = esmc50_domain_goslimpir_tab[`Non-effector` > 0 & Effector == 0, `GO ID`]
esmc50_domain_goslimpir_tab[Effector > 0 & `Non-effector` > 0, c("Odds ratio", "p-value") := fisher.test(matrix(c(Effector, N_domsigs_eff - Effector, `Non-effector`, N_domsigs_noneff - `Non-effector`), ncol = 2))[c("estimate", "p.value")], by = "GO ID"]
esmc50_domain_goslimpir_tab[Effector > 0 & `Non-effector` > 0, "q-value" := p.adjust(`p-value`, method = "BH")]
esmc50_domain_goslimpir_tab[`GO ID` %in% goslim_pir_eff_exc, Enrichment := "Exclusive to effectors"]
esmc50_domain_goslimpir_tab[`GO ID` %in% goslim_pir_noneff_exc, Enrichment := "Exclusive to non-effectors"]
esmc50_domain_goslimpir_tab[`Odds ratio` > 1 & `q-value` < 0.1, Enrichment := "Enriched in effectors"]
esmc50_domain_goslimpir_tab[`Odds ratio` < 1 & `q-value` < 0.1, Enrichment := "Depleted in effectors"]
esmc50_domain_goslimpir_tab[`q-value` > 0.1, Enrichment := "Statistically insignificant"]
esmc50_domain_goslimpir_tab[, Enrichment := factor(Enrichment, levels = c("Statistically insignificant", "Exclusive to non-effectors", "Depleted in effectors", "Enriched in effectors", "Exclusive to effectors"))]
setorderv(esmc50_domain_goslimpir_tab, c("Enrichment", "Odds ratio", "GO Term"), order = c(-1, -1, 1), na.last = TRUE)
esmc50_domain_goslimpir_tab[, c("Odds ratio", "p-value", "q-value") := lapply(.SD, signif, 2), .SDcols = c("Odds ratio", "p-value", "q-value")]
esmc50_domain_goslimpir_tab[, names(esmc50_domain_goslimpir_tab)[1:2] := NULL]
write.table(esmc50_domain_goslimpir_tab, file = "~/database/uniprot/esmc50_effector_goslimpir.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(esmc50_domain, pfam_common_freq, esmc50_domain_goslimpir_tab, file = "~/database/uniprot/esmc50_domain.RData")
