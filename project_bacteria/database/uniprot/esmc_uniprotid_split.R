#!/usr/bin/env Rscript

esmc_uniprotid = readLines("~/database/uniprot/esmc_uniprotid.txt")
l = split(esmc_uniprotid, ceiling(seq_along(esmc_uniprotid)/100000))
for(i in 1:length(l)) {writeLines(l[[i]], paste0("~/database/uniprot/esmc_uniprotid_", i, ".txt"))}
