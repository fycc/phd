#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(stringr)

# InterPro file mapping Pfam to external IDs, also giving total and taxonomy-specific protein counts
ipr2pf = fread("~/database/3did_elm/ipr2pf.csv")
setkeyv(ipr2pf, "xref")

# Pfam domain accession, ID and description
domain_anno = data.table(PFAM_AC = grep("^#=GF AC", readLines("~/database/3did_elm/Pfam-A.hmm.dat"), value = TRUE),
                         PFAM_ID = grep("^#=GF ID", readLines("~/database/3did_elm/Pfam-A.hmm.dat"), value = TRUE),
                         PFAM_DE = grep("^#=GF DE", readLines("~/database/3did_elm/Pfam-A.hmm.dat"), value = TRUE))
domain_anno$PFAM_AC = str_replace(trimws(sub("#=GF AC", "", domain_anno$PFAM_AC, fixed = TRUE)), "(?=[.])(.*?)$", "")
domain_anno$PFAM_ID = trimws(sub("#=GF ID", "", domain_anno$PFAM_ID, fixed = TRUE))
domain_anno$PFAM_DE = trimws(sub("#=GF DE", "", domain_anno$PFAM_DE, fixed = TRUE))

# get domain-domain interactions from 3did and Pfam
ddi_3did = grep("Pfam", readLines("~/database/3did_elm/3did_flat"), value = TRUE)
ddi_3did = cbind(setDT(transpose(str_extract_all(str_extract_all(ddi_3did, "(?<=[(])(.*?)(?=[)])"), "PF\\d{5}"))),
                 setDT(transpose(lapply(ddi_3did, function(x) unlist(strsplit(x, "\t"))[2:3]))))
setnames(ddi_3did, c("PFAM_AC_1", "PFAM_AC_2", "PFAM_ID_1", "PFAM_ID_2"))
ddi_pfam = fread("~/database/3did_elm/pfamA_interactions.txt", header = FALSE)
# make sure all columns contain valid Pfam Accession Numbers (e.g. PF12345)
stopifnot(all(str_detect(unlist(ddi_3did[, .(PFAM_AC_1, PFAM_AC_2)], ddi_pfam), "PF\\d{5}$")))
ddi = unique(rbindlist(list(ddi_3did[, .(PFAM_AC_1, PFAM_AC_2)], ddi_pfam), use.names = FALSE))
setnames(ddi, c("dmA", "dmB"))
ddi = unique(rbindlist(list(ddi[, .(dmA, dmB)], ddi[, .(dmB, dmA)]), use.names = FALSE))
ddi$type = "DDI"

# get domain-motif interactions from 3did
dmi_domain = data.table(PFAM_ID = unname(sapply(grep("Pfam", readLines("~/database/3did_elm/3did_dmi_flat"), value = TRUE), function(x) unlist(strsplit(x, "\t"))[2])))
dmi_domain[, PFAM_AC := domain_anno[match(dmi_domain$PFAM_ID, PFAM_ID), PFAM_AC]]
dmi_domain[is.na(PFAM_AC), PFAM_AC := ddi_3did[match(PFAM_ID, PFAM_ID_1), PFAM_AC_1]]
dmi_domain[is.na(PFAM_AC), PFAM_AC := ddi_3did[match(PFAM_ID, PFAM_ID_2), PFAM_AC_2]]
dmi_motif = grep("^#=PT", readLines("~/database/3did_elm/3did_dmi_flat"), value = TRUE)
dmi_motif = trimws(str_extract(dmi_motif, "(?<=\t)(.*?)(?=[(])"))
stopifnot(all(str_detect(dmi_domain$PFAM_AC, "PF\\d{5}$")) & length(dmi_domain$PFAM_AC) == length(dmi_motif))
dmi_3did = unique(data.table("dmA" = c(dmi_domain$PFAM_AC, dmi_motif), "dmB" = c(dmi_motif, dmi_domain$PFAM_AC), "type" = rep(c("DMI", "MDI"), each = length(dmi_motif))))

# get domain-motif interactions from ELM
elm_classes = as.data.table(read.delim("~/database/3did_elm/elm_classes.tsv", skip = 5))
elm = fread("~/database/3did_elm/elm_interaction_domains.tsv")
setnames(elm, "ELM identifier", "ELMIdentifier")
# convert non-Pfam IDs to Pfam IDs, based on InterPro mapping
xref2pfam = function(x) {
  if (grepl("^IPR", x)) ipr2pf[interpro %in% x & xrefdb == "PFAM", paste0(sort(unique(xref)), collapse = ";")]
  else ipr2pf[interpro %in% ipr2pf[xref %in% x, interpro] & xrefdb == "PFAM", paste0(sort(unique(xref)), collapse = ";")]
}
elm[!str_detect(`Interaction Domain Id`, "PF\\d{5}$"), `Interaction Domain Id` := lapply(`Interaction Domain Id`, xref2pfam), by = "Interaction Domain Id"]
elm = elm[, list(`Interaction Domain Id` = unlist(strsplit(`Interaction Domain Id`, ";"))), by = setdiff(names(elm), "Interaction Domain Id")]
elm[, Regex := elm_classes[match(elm$ELMIdentifier, ELMIdentifier), Regex]]
elm[, ELMType := sapply(ELMIdentifier, function(x) unlist(strsplit(x, "_"))[1])]
elm = elm[!is.na(Regex)]
setkey(elm)
stopifnot(all(str_detect(elm$`Interaction Domain Id`, "PF\\d{5}$")))
dmi_elm = unique(data.table("dmA" = c(elm$`Interaction Domain Id`, elm$Regex), "dmB" = c(elm$Regex, elm$`Interaction Domain Id`), "type" = rep(c("DMI", "MDI"), each = length(elm$Regex))))

# combine ddi and dmi, replace obsolete Pfam IDs with current (merged) Pfam IDs
pfam_dead = readLines("~/database/3did_elm/pfam_dead.txt")
oldpfam = trimws(sub("#=GF AC", "", grep("#=GF AC", pfam_dead, value = TRUE)))
newpfam = trimws(sub("#=GF FW", "", grep("#=GF FW", pfam_dead, value = TRUE)))
pfam_dead = data.table(oldpfam, newpfam)
pfam_dead = unique(pfam_dead[grepl("PF\\d{5}$", oldpfam) & grepl("PF\\d{5}$", newpfam)])
setkeyv(pfam_dead, "oldpfam")

ddi_dmi = rbindlist(list(ddi, dmi_3did, dmi_elm), use.names = TRUE)
ddi_dmi[type %in% c("DDI", "DMI") & !dmA %in% domain_anno$PFAM_AC, dmA := lapply(dmA, function(x) pfam_dead[.(x), paste0(sort(unique(newpfam)), collapse = ";")]), by = "dmA"]
ddi_dmi[type %in% c("DDI", "MDI") & !dmB %in% domain_anno$PFAM_AC, dmB := lapply(dmB, function(x) pfam_dead[.(x), paste0(sort(unique(newpfam)), collapse = ";")]), by = "dmB"]
ddi_dmi = ddi_dmi[, list(dmA = unlist(strsplit(dmA, ";"))), by = setdiff(names(ddi_dmi), "dmA")]
ddi_dmi = ddi_dmi[, list(dmB = unlist(strsplit(dmB, ";"))), by = setdiff(names(ddi_dmi), "dmB")]
ddi_dmi = unique(ddi_dmi)
setkey(ddi_dmi)
stopifnot(ddi_dmi[type %in% c("DDI", "DMI"), all(dmA %in% domain_anno$PFAM_AC)] & ddi_dmi[type %in% c("DDI", "MDI"), all(dmB %in% domain_anno$PFAM_AC)])

elm_motifs = elm[, sort(unique(Regex))]
writeLines(elm_motifs, "~/database/3did_elm/elm_motifs.txt")
all_motifs = sort(union(ddi_dmi[type == "MDI", dmA], ddi_dmi[type == "DMI", dmB]))
writeLines(all_motifs, "~/database/3did_elm/all_motifs.txt")

save(all_motifs, ddi_dmi, dmi_domain, dmi_motif, domain_anno, elm, elm_motifs, file = "~/database/3did_elm/ddi_dmi.RData")
