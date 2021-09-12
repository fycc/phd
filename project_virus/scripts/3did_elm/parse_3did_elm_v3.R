#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(stringr)

# InterPro file mapping Pfam to external IDs, also giving total and taxonomy-specific protein counts
ipr2pf = fread("ipr2pf.csv")

# Pfam domain accession, ID and description
domain_anno = data.table(PFAM_AC = grep("^#=GF AC", readLines("Pfam-A.hmm.dat"), value = TRUE),
                         PFAM_ID = grep("^#=GF ID", readLines("Pfam-A.hmm.dat"), value = TRUE),
                         PFAM_DE = grep("^#=GF DE", readLines("Pfam-A.hmm.dat"), value = TRUE))
domain_anno$PFAM_AC = str_replace(trimws(sub("#=GF AC", "", domain_anno$PFAM_AC, fixed = TRUE)), "(?=[.])(.*?)$", "")
domain_anno$PFAM_ID = trimws(sub("#=GF ID", "", domain_anno$PFAM_ID, fixed = TRUE))
domain_anno$PFAM_DE = trimws(sub("#=GF DE", "", domain_anno$PFAM_DE, fixed = TRUE))
# # get GO annotation for Pfam domains
# # from http://geneontology.org/external2go/pfam2go
# pfam2go_go = grep("^Pfam:", readLines("pfam2go.txt"), value = TRUE)
# pfam2go_go = unique(data.table(PFAM = str_extract(pfam2go_go, "PF\\d{5}"), GO = str_extract(pfam2go_go, "GO:\\d{7}")))
# 
# # from ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz
# ipr2pf = fread("ipr2pf.csv")
# 
# # from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/database_files/gene_ontology.txt.gz
# pfam2go_pfam = fread("gene_ontology.txt", header = FALSE)
# setnames(pfam2go_pfam, c("PFAM", "GO", "GODescription", "GOCategory"))
# 
# pfam2go = unique(rbindlist(list(pfam2go_go, unique(ipr2pf[xrefdb == "PFAM", .(xref, GO)]), pfam2go_pfam[, .(PFAM, GO)])))
# setkey(pfam2go)
# pfam2go = unique(pfam2go[, c("GO") := lapply(.SD, function(x) paste(sort(unique(x)), collapse = ";")), by = "PFAM"])
# 
# domain_anno$GO = pfam2go[match(domain_anno$PFAM_AC, PFAM), GO]
setkeyv(ipr2pf, "xref")
domain_anno[PFAM_AC %in% ipr2pf[xrefdb == "PFAM", xref], c("GO") := lapply(.SD, function(x) paste(sort(unique(ipr2pf[.(x), GO])), collapse = ";")), by = setdiff(names(domain_anno), "PFAM_AC")]
domain_anno[PFAM_AC %in% ipr2pf[xrefdb == "PFAM", xref], c("GOCategory") := lapply(.SD, function(x) paste(sort(unique(ipr2pf[.(x), GOCategory])), collapse = ";")), by = setdiff(names(domain_anno), "PFAM_AC")]
domain_anno[PFAM_AC %in% ipr2pf[xrefdb == "PFAM", xref], c("GODescription") := lapply(.SD, function(x) paste(sort(unique(ipr2pf[.(x), GODescription])), collapse = ";")), by = setdiff(names(domain_anno), "PFAM_AC")]
setkey(domain_anno)

# get domain-domain interactions from 3did, iPfam and Pfam
ddi_3did = grep("Pfam", readLines("3did_flat"), value = TRUE)
ddi_3did = setDT(transpose(str_extract_all(str_extract_all(ddi_3did, "(?<=[(])(.*?)(?=[)])"), "PF\\d{5}")))
ddi_ipfam_het = fread("heterodomain_interaction.txt", header = FALSE, skip = 1, select = c(1, 3))
ddi_ipfam_hom = fread("homodomain_interaction.txt", header = FALSE, skip = 1, select = 1L)
ddi_ipfam_hom$V2 = ddi_ipfam_hom$V1
ddi_pfam = fread("pfamA_interactions.txt", header = FALSE)
# ddi_domine = fread("domine-tables-2.0/INTERACTION.txt", header = FALSE, sep = "|", select = c(1, 2))
# make sure all columns contain valid Pfam Accession Numbers (e.g. PF12345)
stopifnot(all(sapply(list(ddi_3did, ddi_ipfam_het, ddi_ipfam_hom, ddi_pfam), function(DT) all(str_detect(unlist(DT), "PF\\d{5}$")))))
ddi = unique(rbindlist(list(ddi_3did, ddi_ipfam_het, ddi_ipfam_hom, ddi_pfam), use.names = FALSE))
setnames(ddi, c("dmA", "dmB"))
ddi = unique(rbindlist(list(ddi[, .(dmA, dmB)], ddi[, .(dmB, dmA)]), use.names = FALSE))
ddi$type = "DDI"

# get domain-motif interactions from 3did
dmi_domain = grep("Pfam", readLines("3did_dmi_flat"), value = TRUE)
dmi_domain = str_extract(dmi_domain, "(?<=\t)(.*?)(?=\t)")
dmi_domain = gsub("TFIIH_BTF_p62_N", "PH_TFIIH", dmi_domain, fixed = TRUE)
dmi_domain = gsub("Ribosomal_L18ae", "Ribosomal_L18A", dmi_domain, fixed = TRUE)
dmi_domain = domain_anno[match(dmi_domain, PFAM_ID), PFAM_AC]
dmi_motif = grep("^#=PT", readLines("3did_dmi_flat"), value = TRUE)
dmi_motif = trimws(str_extract(dmi_motif, "(?<=\t)(.*?)(?=[(])"))
stopifnot(all(str_detect(dmi_domain, "PF\\d{5}$")) & length(dmi_domain) == length(dmi_motif))
dmi_3did = unique(data.table("dmA" = c(dmi_domain, dmi_motif), "dmB" = c(dmi_motif, dmi_domain), "type" = rep(c("DMI", "MDI"), each = length(dmi_domain))))

# get domain-motif interactions from ELM
elm_classes = as.data.table(read.delim("elm_classes.tsv", skip = 5))
elm_instances = as.data.table(read.delim("elm_instances.tsv", skip = 5))
elm_instances = elm_instances[InstanceLogic == "true positive"]
elm = fread("elm_interaction_domains.tsv")
# convert non-Pfam IDs to Pfam IDs, based on InterPro mapping
xref2pf = function(x) {
  if(grepl("^IPR", x)) paste(sort(unique(ipr2pf[interpro %in% x & xrefdb == "PFAM", xref])), collapse = ";")
  else paste(sort(unique(ipr2pf[interpro %in% ipr2pf[xref %in% x, interpro] & xrefdb == "PFAM", xref])), collapse = ";")
}
elm[!str_detect(`Interaction Domain Id`, "PF\\d{5}"), c("Interaction Domain Id") := lapply(.SD, xref2pf), by = setdiff(names(elm), "Interaction Domain Id")]
elm = elm[, list(`Interaction Domain Id` = unlist(strsplit(`Interaction Domain Id`, ";"))), by = setdiff(names(elm), "Interaction Domain Id")]
elm[, c("Regex") := .(elm_classes[match(`ELM identifier`, ELMIdentifier), Regex])]
elm[, c("ELMType") := .(elm_instances[match(`ELM identifier`, ELMIdentifier), ELMType])]
elm[, c("InstanceLogic") := .(elm_instances[match(`ELM identifier`, ELMIdentifier), InstanceLogic])]
elm = elm[!is.na(Regex)]
setkey(elm)
stopifnot(all(str_detect(elm$`Interaction Domain Id`, "PF\\d{5}$")))
dmi_elm = unique(data.table("dmA" = c(elm$`Interaction Domain Id`, elm$Regex), "dmB" = c(elm$Regex, elm$`Interaction Domain Id`), "type" = rep(c("DMI", "MDI"), each = length(elm$Regex))))

# combine ddi and dmi, replace obsolete Pfam IDs with current (merged) Pfam IDs
pfam_dead = readLines("pfam_dead.txt")
oldpfam = trimws(sub("#=GF AC", "", grep("#=GF AC", pfam_dead, value = TRUE)))
newpfam = trimws(sub("#=GF FW", "", grep("#=GF FW", pfam_dead, value = TRUE)))
pfam_dead = data.table(oldpfam, newpfam)
pfam_dead = unique(pfam_dead[grepl("^PF", oldpfam) & grepl("^PF", newpfam)])
setkey(pfam_dead)

ddi_dmi = rbindlist(list(ddi, dmi_3did, dmi_elm), use.names = TRUE)
ddi_dmi[type %in% c("DDI", "DMI") & !dmA %in% domain_anno$PFAM_AC, c("dmA") := lapply(.SD, function(x) paste(sort(unique(pfam_dead[.(x), newpfam])), collapse = ";")), by = setdiff(names(ddi_dmi), "dmA")]
ddi_dmi[type %in% c("DDI", "MDI") & !dmB %in% domain_anno$PFAM_AC, c("dmB") := lapply(.SD, function(x) paste(sort(unique(pfam_dead[.(x), newpfam])), collapse = ";")), by = setdiff(names(ddi_dmi), "dmB")]
ddi_dmi = ddi_dmi[, list(dmA = unlist(strsplit(dmA, ";"))), by = setdiff(names(ddi_dmi), "dmA")]
ddi_dmi = ddi_dmi[, list(dmB = unlist(strsplit(dmB, ";"))), by = setdiff(names(ddi_dmi), "dmB")]
ddi_dmi = unique(ddi_dmi)
setkey(ddi_dmi)

all_motifs = sort(union(ddi_dmi[type == "MDI", dmA], ddi_dmi[type == "DMI", dmB]))

save(all_motifs, ddi_dmi, domain_anno, elm, file = "ddi_dmi.RData")
