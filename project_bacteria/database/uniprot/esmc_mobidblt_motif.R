#!/usr/bin/env Rscript
#SBATCH --account=def-yxia
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8000M
#SBATCH --job-name=esmc_mobidblt_motif
#SBATCH --output=%x-%j.out

library(data.table)
library(parallel)
library(seqinr)
library(stringr)

mobidblt_esmc_fasta = read.fasta("~/database/uniprot/mobidblt_esmc.fasta", seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
esmc_mobidblt_motif = fread("~/database/uniprot/mobidblt_esmc.tab", header = FALSE)
setnames(esmc_mobidblt_motif, c("uniprotid", "length", "mobidblt_start", "mobidblt_end"))
esmc_mobidblt_motif = esmc_mobidblt_motif[uniprotid %in% names(mobidblt_esmc_fasta)]
setkeyv(esmc_mobidblt_motif, "uniprotid")

# Get fraction of disordered residues and number of motifs in disordered regions (MobiDB-lite) of bacterial proteins
all_motifs = readLines("~/database/3did_elm/all_motifs.txt")
uids = esmc_mobidblt_motif[, unique(uniprotid)]
uids_list = split(uids, ceiling(seq_along(uids) / 100000))
residues_motifs = list()
for (i in 1:length(uids_list)) {
  residues_motifs[[i]] = mclapply(uids_list[[i]], function(u) {
    disordered = esmc_mobidblt_motif[.(u), sort(unique(unlist(lapply(paste(mobidblt_start, mobidblt_end, sep = ":"), function(x) eval(parse(text = x))))))]
    residues = length(disordered)
    disordered = split(disordered, cumsum(c(0, diff(disordered) > 1)))
    # get coordinates of all motif matches in full protein
    s = str_locate_all(mobidblt_esmc_fasta[[u]], all_motifs)
    # find coordinates that are inside disordered regions
    L = sapply(s[lengths(s) > 0], function(x) {
      sapply(disordered, function(y) apply(x, 1, function(j) all(eval(parse(text = paste(j["start"], j["end"], sep = ":"))) %in% y)))
    })
    motifs = lengths(s) / 2
    motifs[motifs > 0] = sapply(L, sum)
    list(residues = residues, motifs = paste0(motifs, collapse = ";"))
  }, mc.cores = 4)
}
residues_motifs = unlist(residues_motifs, recursive = FALSE)
names(residues_motifs) = uids
esmc_mobidblt_motif = unique(esmc_mobidblt_motif[, .(uniprotid, length)])
esmc_mobidblt_motif[, mobidblt_residues := sapply(residues_motifs, "[[", "residues")]
esmc_mobidblt_motif[, frac_mobidblt := mobidblt_residues / length]
esmc_mobidblt_motif[, mobidblt_motifs := sapply(residues_motifs, "[[", "motifs")]

save(residues_motifs, file = "~/database/uniprot/esmc_mobidblt_residues_motifs.RObject")
save(esmc_mobidblt_motif, file = "~/database/uniprot/esmc_mobidblt_motif.RObject")
