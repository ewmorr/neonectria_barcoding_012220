library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
require(tidyverse)

get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

setwd("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320")
system("mkdir separate_denoising_per_run")

seqtab.nochim.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/taxa_w_bootstraps.rds")

seqtab.nochim.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/taxa_w_bootstraps.rds")

#set sample names to match metadata
rownames(seqtab.nochim.r1) = unname(sapply(rownames(seqtab.nochim.r1), get.sample.name))
rownames(seqtab.nochim.r1) = paste0("runOne", rownames(seqtab.nochim.r1))
rownames(seqtab.nochim.r2) = unname(sapply(rownames(seqtab.nochim.r2), get.sample.name))
rownames(seqtab.nochim.r2) = paste0("runTwo", rownames(seqtab.nochim.r2))

#join tables
#full
seqtab.nochim.join_by_ASVseq = full_join(
    data.frame(asv.seq = rownames(t(seqtab.nochim.r1)), t(seqtab.nochim.r1)),
    data.frame(asv.seq = rownames(t(seqtab.nochim.r2)), t(seqtab.nochim.r2)),
    by = "asv.seq"
)
rownames(seqtab.nochim.join_by_ASVseq) = seqtab.nochim.join_by_ASVseq$asv.seq
seqtab.nochim.join_by_ASVseq$asv.seq = NULL
seqtab.nochim.join_by_ASVseq[is.na(seqtab.nochim.join_by_ASVseq)] = 0

#shared
seqtab.nochim.join_by_ASVseq.sharedASVs = inner_join(
data.frame(asv.seq = rownames(t(seqtab.nochim.r1)), t(seqtab.nochim.r1)),
data.frame(asv.seq = rownames(t(seqtab.nochim.r2)), t(seqtab.nochim.r2)),
by = "asv.seq"
)
rownames(seqtab.nochim.join_by_ASVseq.sharedASVs) = seqtab.nochim.join_by_ASVseq.sharedASVs$asv.seq
seqtab.nochim.join_by_ASVseq.sharedASVs$asv.seq = NULL

#run1 uniques
seqtab.nochim.join_by_ASVseq.run1uniqASVs = anti_join(
data.frame(asv.seq = rownames(t(seqtab.nochim.r1)), t(seqtab.nochim.r1)),
data.frame(asv.seq = rownames(t(seqtab.nochim.r2)), t(seqtab.nochim.r2)),
by = "asv.seq"
)
rownames(seqtab.nochim.join_by_ASVseq.run1uniqASVs) = seqtab.nochim.join_by_ASVseq.run1uniqASVs$asv.seq
seqtab.nochim.join_by_ASVseq.run1uniqASVs$asv.seq = NULL

#run2 uniques
seqtab.nochim.join_by_ASVseq.run2uniqASVs = anti_join(
data.frame(asv.seq = rownames(t(seqtab.nochim.r2)), t(seqtab.nochim.r2)),
data.frame(asv.seq = rownames(t(seqtab.nochim.r1)), t(seqtab.nochim.r1)),
by = "asv.seq"
)
rownames(seqtab.nochim.join_by_ASVseq.run2uniqASVs) = seqtab.nochim.join_by_ASVseq.run2uniqASVs$asv.seq
seqtab.nochim.join_by_ASVseq.run2uniqASVs$asv.seq = NULL


#counts
total_ASVs = seqtab.nochim.join_by_ASVseq[rowSums(seqtab.nochim.join_by_ASVseq) > 0,]  %>% rownames %>% length
shared_ASVs = seqtab.nochim.join_by_ASVseq.sharedASVs[rowSums(seqtab.nochim.join_by_ASVseq.sharedASVs) > 0,]  %>% rownames %>% length

run1_ASVs = t(seqtab.nochim.r1)[rowSums(t(seqtab.nochim.r1)) > 0,] %>% rownames %>% length
run2_ASVs = t(seqtab.nochim.r2)[rowSums(t(seqtab.nochim.r2)) > 0,] %>% rownames %>% length

run1_uniq = seqtab.nochim.join_by_ASVseq.run1uniqASVs[rowSums(seqtab.nochim.join_by_ASVseq.run1uniqASVs) > 0,]  %>% rownames %>% length
run2_uniq = seqtab.nochim.join_by_ASVseq.run2uniqASVs[rowSums(seqtab.nochim.join_by_ASVseq.run2uniqASVs) > 0,]  %>% rownames %>% length

#Venn Diagram
require(VennDiagram)
venn.diagram(list(run1 = colnames(seqtab.nochim.r1), run2 = colnames(seqtab.nochim.r2)), filename = "shared_ASVs_Venn.tiff")

#alternate way to get names
rownames(seqtab.nochim.join_by_ASVseq[
    rowSums(seqtab.nochim.join_by_ASVseq[ ,1:length(rownames(seqtab.nochim.r1))
    ]) > 0
,])
rownames( seqtab.nochim.join_by_ASVseq[
    rowSums(seqtab.nochim.join_by_ASVseq[ ,length(rownames(seqtab.nochim.r1)) + 1:length(rownames(seqtab.nochim.r2))
    ]) > 0
,])
