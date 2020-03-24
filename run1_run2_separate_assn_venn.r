library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
require(tidyverse)
require(VennDiargram)

get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

setwd("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320")

seqtab.nochim.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/dada2_seq_table_no_chim.rds")
seqtab.nochim.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/dada2_seq_table_no_chim.rds")

#set sample names to match metadata
rownames(seqtab.nochim.r1) = unname(sapply(rownames(seqtab.nochim.r1), get.sample.name))
rownames(seqtab.nochim.r1) = paste0("runOne", rownames(seqtab.nochim.r1))
rownames(seqtab.nochim.r2) = unname(sapply(rownames(seqtab.nochim.r2), get.sample.name))
rownames(seqtab.nochim.r2) = paste0("runTwo", rownames(seqtab.nochim.r2))

seqtab.nochim.r1.t = t(seqtab.nochim.r1)
seqtab.nochim.r2.t = t(seqtab.nochim.r2)

cbind(data.frame(run1 = rownames(seqtab.nochim.r1.t)), run2 = rownames(seqtab.nochim.r2.t)))

data.frame(rowSums(seqtab.nochim.r1.t)
