#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(parallel)
library(seqinr)
library(stringr)

# function to expand data frame where the entries of one column may be collapsed by a symbol (e.g. "_" or ",")
expand.df = function(dat, colname, collapse) {
  u = grep(collapse, dat[, colname])
  if(length(u) == 0) dat
  else {
    for(i in 1:length(u)) {
      y = dat[u[i], ]
      require(stringr)
      counts = str_count(y[, colname], collapse)
      for(j in 1:(counts + 1)) {
        z = y
        z[, colname] = unlist(strsplit(z[, colname], collapse))[j]
        dat = rbind(dat, z)
      }
    }
    dat = dat[-grep(collapse, dat[, colname]), ]
    dat
  }
}

clinvar = read.delim("database/clinvar/variant_summary.txt")
clinvar = as.data.frame(apply(clinvar, 2, trimws))
colnames(clinvar)[1] = "AlleleID"
clinvar$PhenotypeIDs = gsub(",", ";", clinvar$PhenotypeIDs)
clinvar = clinvar[grepl("^NP_", clinvar$HGVS.p..), ]
clinvar$refseq = as.character(sapply(clinvar$HGVS.p.., function(x) unlist(strsplit(x, ":"))[1]))
clinvar$refseq = as.character(sapply(clinvar$refseq, function(x) unlist(strsplit(x, "[.]"))[1]))
writeLines(unique(clinvar$refseq), "database/clinvar/clinvar_refseq.txt")

# Download all human refseq protein sequences in fasta format
# http://www.ncbi.nlm.nih.gov/protein?term=(homo%20sapiens%5BOrganism%5D)%20AND%20srcdb_refseq%5BProperties%5D
refseq_fasta = read.fasta("database/clinvar/human_refseq.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
refseq_fasta = refseq_fasta[grep("NP_", names(refseq_fasta))]
names(refseq_fasta) = as.character(sapply(names(refseq_fasta), function(x) grep("^NP_", unlist(strsplit(x, "[|]")), value = TRUE)))
names(refseq_fasta) = as.character(sapply(names(refseq_fasta), function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_fasta = refseq_fasta[unique(clinvar$refseq)]
write.fasta(clinvar_fasta, names(clinvar_fasta), "database/clinvar/clinvar.fasta", as.string = TRUE)
# qsub -pe omp 16 database/clinvar/clinvar_pfam.sh
# awk 'NR > 28 {print $1,$2,$3,$4,$5,$6,$7}' < database/clinvar/clinvar.pfam > database/clinvar/clinvar_pfam.tab

clinvar_r2u = read.delim("database/clinvar/clinvar_refseq2uniprot.tab")
colnames(clinvar_r2u)[grep("yourlist", colnames(clinvar_r2u))] = "refseq"
colnames(clinvar_r2u)[grep("isomap", colnames(clinvar_r2u))] = "isomap"
annocols = grepl("Gene.ontology.IDs|Cross.reference", colnames(clinvar_r2u)) & !grepl("GeneID", colnames(clinvar_r2u))

y = clinvar_r2u[duplicated(clinvar_r2u$refseq) | duplicated(clinvar_r2u$refseq, fromLast = TRUE), ]
y = y[order(y$refseq, y$Status, -y$Length, -nchar(y$Cross.reference..Pfam.), -apply(y[, annocols], 1, function(x) {
  x = unlist(x)
  x = gsub(" |;|,", "", x[!is.na(x)])
  sum(nchar(x))
})), ]
y = y[!duplicated(y$refseq), ]
clinvar_r2u = rbind.data.frame(clinvar_r2u[!clinvar_r2u$refseq %in% y$refseq, ], y)
rm(y)
clinvar_r2u$Gene.names...primary.. = as.character(sapply(clinvar_r2u$Gene.names...primary.., function(x) unlist(strsplit(x, ";"))[1]))
clinvar_r2u$Cross.reference..GeneID. = as.character(sapply(clinvar_r2u$Cross.reference..GeneID., function(x) unlist(strsplit(x, ";"))[1]))

clinvar_uniprot = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Entry[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$uniprot = as.character(clinvar_uniprot[clinvar$refseq])
clinvar_geneid = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Cross.reference..GeneID.[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$GeneID = as.character(clinvar_geneid[clinvar$refseq])
clinvar_genesymbol = sapply(unique(clinvar$refseq), function(x) clinvar_r2u$Gene.names...primary..[grep(paste("\\b", x, "\\b", sep = ""), clinvar_r2u$refseq)[1]])
clinvar$GeneSymbol = as.character(clinvar_genesymbol[clinvar$refseq])
clinvar$GeneSymbol[is.na(clinvar$GeneSymbol)] = as.character(sapply(clinvar$Name[is.na(clinvar$GeneSymbol)], function(x) unlist(strsplit(x, "[()]"))[2]))

load("database/clinvar/humsavar.RData")
clinvar = clinvar[grep("MedGen:", clinvar$PhenotypeIDs), ]
pid2cid = sapply(unique(clinvar$PhenotypeIDs), function(x) {
  p2c = gsub("MedGen:", "", grep("^MedGen:", unlist(strsplit(x, ";")), value = TRUE))
  paste(intersect(p2c, cid_disease$ConceptID), collapse = ";")
})
pid2cid = pid2cid[lengths(pid2cid) > 0 & nchar(pid2cid) > 0]
clinvar$ConceptID = as.character(pid2cid[clinvar$PhenotypeIDs])
clinvar = clinvar[!is.na(clinvar$ConceptID), ]

cid2disease = sapply(unique(clinvar$ConceptID), function(x) {
  cids = unlist(strsplit(x, ";"))
  paste(unique(as.character(cid_disease$DiseaseName[match(cids, cid_disease$ConceptID)])), collapse = ";")
})
clinvar$disease = as.character(cid2disease[clinvar$ConceptID])

############################
# Map SNPs to PFAM domains #
############################

clinvar = clinvar[grepl("p.", clinvar$HGVS.p.., fixed = TRUE), ]
# clinvar = clinvar[!grepl("[=|?]", clinvar$HGVS.p..), ]
clinvar = clinvar[!grepl("[?]", clinvar$HGVS.p..), ]
clinvar = clinvar[!grepl("fs", clinvar$HGVS.p..), ]
clinvar = clinvar[!grepl("ext", clinvar$HGVS.p..), ]
clinvar = clinvar[!grepl("[()]", clinvar$HGVS.p..), ]
clinvar$HGVS.p.. = gsub("[*]", "Ter", clinvar$HGVS.p..)

library(Biostrings)
aa = paste(c(AMINO_ACID_CODE, names(AMINO_ACID_CODE), "[+|_|=]", "del", "dup", "ins", "Ter"), collapse = "|")

hgvs2aapos = function(hgvs) {
  pos = tail(unlist(strsplit(hgvs, "p.", fixed = TRUE)), 1)
  pos = unlist(strsplit(pos, aa))
  pos = sort(as.numeric(unique(pos[nchar(pos) > 0])))
  
  if (grepl("Ter", hgvs)) paste(min(pos), nchar(clinvar_fasta[[unlist(strsplit(hgvs, "[.|:]"))[1]]]), sep = ":")
  else if (length(pos) == 1) pos
  else paste(min(pos), max(pos), sep = ":")
}

aapos = sapply(unique(clinvar$HGVS.p..), hgvs2aapos)
clinvar$aapos = as.character(aapos[clinvar$HGVS.p..])

clinvar_pfam = read.table("database/clinvar/clinvar_pfam.tab")
clinvar_pfam$V6 = as.character(sapply(clinvar_pfam$V6, function(x) unlist(strsplit(x, "[.]"))[1]))
clinvar_pfams = split(clinvar_pfam[, 4:6], clinvar_pfam$V1)

clinvar = clinvar[clinvar$refseq %in% names(clinvar_pfams), ]
refseq_aapos = unique(paste(clinvar$refseq, clinvar$aapos, sep = ";"))
pfams = mclapply(refseq_aapos, function(x) {
  
  refseq = unlist(strsplit(x, ";"))[1]
  aapos = unlist(strsplit(x, ";"))[2]
  
  pfam = clinvar_pfams[[refseq]]
  i = apply(pfam, 1, function(p) any(eval(parse(text = aapos)) %in% p["V4"]:p["V5"]))
  ifelse(sum(i) == 0, NA, paste(unique(pfam$V6[i]), collapse = ";"))
  
}, mc.cores = detectCores())
names(pfams) = refseq_aapos
pfams = unlist(pfams)
clinvar$pfam = as.character(pfams[paste(clinvar$refseq, clinvar$aapos, sep = ";")])
clinvar = clinvar[!is.na(clinvar$pfam), ]
clinvar = expand.df(dat = clinvar, colname = "pfam", collapse = ";")

snpeff = read.delim("database/clinvar/Galaxy2-[SnpEff_on_clinvar_vcf].vcf", skip = 71)
colnames(snpeff)[1] = "CHROM"
snpeff$ID = gsub("rs", "", snpeff$ID)
nmd = snpeff[grep("NMD=", snpeff$INFO), ]
clinvar$nmd = clinvar$RS...dbSNP. %in% nmd$ID
clinvar = unique(clinvar)

################################
# Combine humsavar and clinvar #
################################

oneRowPerId = function(x, ids) {
  stopifnot(all(x[,1] %in% ids))
  x = x[ apply(x[,-1,drop=FALSE], 1, function(y) any(y!="")), ]
  d = lapply(2:ncol(x), function(i) {
    r  = character(length(ids))
    v  = sapply(split(x[,i], x[,1]), unique)
    v  = sapply(v, paste, collapse=";")
    # mt = match(names(v), ids)
    # r[mt] = v
    mt = match(ids, names(v))
    r[!is.na(mt)] = v[mt[!is.na(mt)]]
    r[r==""] = NA
    return(I(r))
  })
  names(d) = colnames(x)[2:ncol(x)]
  do.call(data.frame, d)
}

humsa = unique(humsavar[, c("gene", "uniprot", "aapos", "pfam", "ConceptID")])
clin = unique(clinvar[, c("GeneSymbol", "uniprot", "aapos", "pfam", "ConceptID")])
colnames(clin)[1] = "gene"
hcvar = unique(rbind.data.frame(humsa, clin))
hcvar = data.frame(IDs = paste(hcvar$gene, hcvar$pfam, sep = "_"), hcvar)
# hcvar = data.frame(IDs = paste(hcvar$uniprot, hcvar$pfam, sep = "_"), hcvar)
hcvar = oneRowPerId(hcvar, unique(hcvar$IDs))
hcvar$ConceptID = as.character(sapply(hcvar$ConceptID, function(x) paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")))
hcvar$aapos_unique = str_count(hcvar$aapos, ";") + 1

find_contig = function(contig) {
  contig = unlist(strsplit(contig, ";"))
  contig = unique(sort(unname(unlist(sapply(contig, function(x) eval(parse(text = x)))))))
  start = c(1, which(diff(contig) != 1 & diff(contig) != 0) + 1)
  end = c(start - 1, length(contig))
  paste(mapply(function(x, y) ifelse(x == y, x, paste(x, y, sep = ":")), contig[start], contig[end]), collapse = ";")
}
hcvar$aapos_final = as.character(sapply(hcvar$aapos, find_contig))
hcvar = hcvar[order(hcvar$gene, hcvar$pfam), c("gene", "uniprot", "aapos", "aapos_unique", "aapos_final", "pfam", "ConceptID")]
rm(humsa, clin)

hcvar$DiseaseNames = as.character(sapply(hcvar$ConceptID, function(x) {
  cids = unlist(strsplit(x, ";"))
  paste(as.character(cid_disease$DiseaseName[match(cids, cid_disease$ConceptID)]), collapse = ";")
}))

hcvar$meshtree = as.character(sapply(hcvar$ConceptID, function(x) {
  cids = unlist(strsplit(x, ";"))
  trees = unique(cid_disease$meshtree[match(cids, cid_disease$ConceptID)])
  ifelse(all(is.na(trees)), NA, paste(trees[!is.na(trees)], collapse = ";"))
}))

hcvar$meshtreetop = as.character(sapply(hcvar$ConceptID, function(x) {
  cids = unlist(strsplit(x, ";"))
  treetops = unique(cid_disease$meshtreetop[match(cids, cid_disease$ConceptID)])
  ifelse(all(is.na(treetops)), NA, paste(treetops[!is.na(treetops)], collapse = ";"))
}))

###########################################################
# Intersect disease SNPs with virally implicated diseases #
###########################################################

vds_cids = sapply(vds, function(vd) {
  sapply(vd, function(cid) grep(paste("\\b", cid, "\\b", sep = ""), hcvar$ConceptID), simplify = FALSE)
})
vds_cids = sapply(vds_cids, function(x) x[lengths(x) > 0])
vds_cids = vds_cids[lengths(vds_cids) > 0]

Virus = rep(names(vds_cids), sapply(vds_cids, function(x) length(unlist(x))))
ConceptID = rep(unname(unlist(sapply(vds_cids, names))), rapply(vds_cids, length))
DiseaseName = as.character(cid_disease$DiseaseName[match(ConceptID, cid_disease$ConceptID)])
hcvar_vid = unique(data.frame(Virus, ConceptID, DiseaseName, hcvar[unlist(vds_cids), c("gene", "uniprot", "pfam")]))
hcvar_vid$meshtreetop = cid_disease$meshtreetop[match(hcvar_vid$ConceptID, cid_disease$ConceptID)]

save.image("database/clinvar/clinvar.RData")
