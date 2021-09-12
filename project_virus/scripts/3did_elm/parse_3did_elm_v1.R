#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# read in Pfam domain accession, ID and description
domain_anno = data.frame(PFAM_AC = grep("^#=GF AC", readLines("database/3did_elm/Pfam-A.hmm.dat"), value = TRUE),
                         PFAM_ID = grep("^#=GF ID", readLines("database/3did_elm/Pfam-A.hmm.dat"), value = TRUE),
                         PFAM_DE = grep("^#=GF DE", readLines("database/3did_elm/Pfam-A.hmm.dat"), value = TRUE))
domain_anno$PFAM_AC = sapply(domain_anno$PFAM_AC, function(x) unlist(strsplit(trimws(gsub("#=GF AC", "", x, fixed = TRUE)), "[.]"))[1])
domain_anno$PFAM_ID = sapply(domain_anno$PFAM_ID, function(x) trimws(gsub("#=GF ID", "", x, fixed = TRUE)))
domain_anno$PFAM_DE = sapply(domain_anno$PFAM_DE, function(x) trimws(gsub("#=GF DE", "", x, fixed = TRUE)))

# get domain-domain interactions from 3did and iPfam
ddi_3did = readLines("database/3did_elm/3did_flat")
ddi_3did = ddi_3did[grep("^#=ID", ddi_3did)]
ddi_3did = matrix(sapply(ddi_3did, function(x) unlist(strsplit(gsub("@Pfam|[()]| |[.]\\d+", "", x), "\t"))[4:5]), byrow = TRUE, ncol = 2)
ddi_3did_pairs = union(paste(ddi_3did[, 1], ddi_3did[, 2], sep = "_"), paste(ddi_3did[, 2], ddi_3did[, 1], sep = "_"))

ddi_ipfam_het = read.table("database/3did_elm/heterodomain_interaction.txt", header = TRUE)
ddi_ipfam_het = ddi_ipfam_het[, c("family", "family.1")]
ddi_ipfam_hom = read.table("database/3did_elm/homodomain_interaction.txt", header = TRUE)
ddi_ipfam_hom = ddi_ipfam_hom$family
ddi_ipfam_pairs = unique(c(paste(ddi_ipfam_het$family, ddi_ipfam_het$family.1, sep = "_"), paste(ddi_ipfam_het$family.1, ddi_ipfam_het$family, sep = "_"), paste(ddi_ipfam_hom, ddi_ipfam_hom, sep = "_")))

ddi_pfam = read.table("database/3did_elm/pfamA_interactions.txt")
ddi_pfam_pairs = union(paste(ddi_pfam[, 1], ddi_pfam[, 2], sep = "_"), paste(ddi_pfam[, 2], ddi_pfam[, 1], sep = "_"))

# ddi_domine = readLines("database/3did_elm/domine-tables-2.0/INTERACTION.txt")
# ddi_domine = as.character(sapply(ddi_domine, function(x) unlist(strsplit(x, "[|]"))[1:2]))
# ddi_domine = matrix(ddi_domine, ncol = 2, byrow = TRUE)
# ddi_domine_pairs = union(paste(ddi_domine[, 1], ddi_domine[, 2], sep = "_"),
#                          paste(ddi_domine[, 2], ddi_domine[, 1], sep = "_"))

# ddi_pairs = sort(unique(c(ddi_3did_pairs, ddi_ipfam_pairs, ddi_domine_pairs)))
ddi_pairs = sort(unique(c(ddi_3did_pairs, ddi_ipfam_pairs, ddi_pfam_pairs)))
ddi_mat = matrix(unlist(strsplit(ddi_pairs, "_")), ncol = 2, byrow = TRUE)
ddi_mat = cbind(ddi_mat, rep("DDI", nrow(ddi_mat)))
ddi_mat = ddi_mat[ddi_mat[, 1] %in% domain_anno$PFAM_AC & ddi_mat[, 2] %in% domain_anno$PFAM_AC, ]
ddi_mat = unique(ddi_mat)

# get domain-motif interactions from 3did and ELM
dmi = readLines("database/3did_elm/3did_dmi_flat")
dmi_domain = dmi[grep("^#=ID", dmi)]
dmi_domain = as.character(sapply(dmi_domain, function(x) unlist(strsplit(x, "\t"))[2]))
dmi_domain[dmi_domain == "TFIIH_BTF_p62_N"] = "PH_TFIIH"
dmi_domain[dmi_domain == "Ribosomal_L18ae"] = "Ribosomal_L18A"
dmi_domain = domain_anno$PFAM_AC[match(dmi_domain, domain_anno$PFAM_ID)]
dmi_motif = dmi[grep("^#=PT", dmi)]
dmi_motif = as.character(sapply(dmi_motif, function(x) unlist(strsplit(x, "\t"))[2]))
dmi_motif = as.character(sapply(dmi_motif, function(x) unlist(strsplit(x, " "))[1]))
stopifnot(length(dmi_domain) == length(dmi_motif))
dmi_mat = matrix(c(dmi_domain, dmi_motif), ncol = 2)
dmi_mat = dmi_mat[!is.na(dmi_mat[, 1]) & !is.na(dmi_mat[, 2]), ]
dmi_mat = rbind(dmi_mat, dmi_mat[, 2:1])
dmi_mat = cbind(dmi_mat, rep(c("DMI", "MDI"), each = nrow(dmi_mat)/2))
dmi_mat = unique(dmi_mat)

ddi_dmi = unique(rbind(ddi_mat, dmi_mat))

# require(XML)
# elm_classes = readHTMLTable("http://elm.eu.org/elms/browse_elms.html")
# elm_classes = as.data.frame(elm_classes, stringsAsFactors = FALSE)
# colnames(elm_classes) = gsub("NULL.", "", colnames(elm_classes), fixed = TRUE)
# elm_classes$ELM.Identifier = gsub("[updated]", "", elm_classes$ELM.Identifier, fixed = TRUE)
# elm_classes$ELM.Identifier = gsub("[new]", "", elm_classes$ELM.Identifier, fixed = TRUE)

elm_classes = read.delim("database/3did_elm/elm_classes.tsv", skip = 5)
elm_instances = read.delim("database/3did_elm/elm_instances.tsv", skip = 5)
# elm_neg = elm_instances$ELMIdentifier[elm_instances$InstanceLogic != "true positive"]
# elm_neg = elm_classes$Regex[elm_classes$ELMIdentifier %in% elm_neg]
elm_instances = subset(elm_instances, InstanceLogic == "true positive")
elm_interaction_domains = read.delim("database/3did_elm/elm_interaction_domains.tsv")

elm = data.frame(elm_interaction_domains,
                 Regex = elm_classes$Regex[match(elm_interaction_domains$ELM.identifier, elm_classes$ELMIdentifier)],
                 ELMType = elm_instances$ELMType[match(elm_interaction_domains$ELM.identifier, elm_instances$ELMIdentifier)])
elm = subset(elm, !is.na(ELMType) & Interaction.Domain.Id %in% domain_anno$PFAM_AC)
# Look at the frequency each motif is a "true positive", "false positive", "true negative", and "unknown".
# elm_instances = read.delim("database/3did_elm/elm_instances.tsv", skip = 5)
# p = elm_instances[elm_instances$ELMIdentifier %in% elm$ELM.identifier, ]
# tab = tapply(p$InstanceLogic, p$ELMIdentifier, table)

dmi_elm = rbind(as.matrix(elm[, c("Interaction.Domain.Id", "Regex")]), as.matrix(elm[, c("Regex", "Interaction.Domain.Id")]))
dmi_elm = cbind(dmi_elm, rep(c("DMI", "MDI"), each = nrow(dmi_elm)/2))

# ddi_dmi = rbind(ddi_mat, dmi_elm)
ddi_dmi = rbind(ddi_dmi, dmi_elm)
ddi_dmi = ddi_dmi[ddi_dmi[, 1][ddi_dmi[, 3] == "DMI"] %in% domain_anno$PFAM_AC & ddi_dmi[, 2][ddi_dmi[, 3] == "MDI"] %in% domain_anno$PFAM_AC, ]
# ddi_dmi = ddi_dmi[!ddi_dmi[, 1] %in% elm_neg & !ddi_dmi[, 2] %in% elm_neg, ]
ddi_dmi = unique(ddi_dmi)
dimnames(ddi_dmi) = NULL

# domain_anno = rbind(as.matrix(domain_anno), as.matrix(elm_interaction_domains[, c("Interaction.Domain.Id", "Interaction.Domain.Description", "Interaction.Domain.Name")]))
# domain_anno = as.data.frame(domain_anno[!duplicated(domain_anno[, 1:2]), ])
# rownames(domain_anno) = NULL
# colnames(domain_anno) = c("domain_acc", "domain_name", "domain_desc")

all_motifs = union(ddi_dmi[, 1][ddi_dmi[, 3] == "MDI"], ddi_dmi[, 2][ddi_dmi[, 3] == "DMI"])

save(ddi_dmi, elm, all_motifs, domain_anno, file = "database/3did_elm/ddi_dmi.RData")
