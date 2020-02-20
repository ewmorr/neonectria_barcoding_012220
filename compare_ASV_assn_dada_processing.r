library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
require(tidyverse)

get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

setwd("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320")
system("mkdir dada2_processing_compare")

seqtab.nochim.pool_sep.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.pool_sep.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/taxa_w_bootstraps.rds")

seqtab.nochim.pool_sep.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.pool_sep.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/taxa_w_bootstraps.rds")

#set sample names to match metadata
rownames(seqtab.nochim.pool_sep.r1) = unname(sapply(rownames(seqtab.nochim.pool_sep.r1), get.sample.name))
rownames(seqtab.nochim.pool_sep.r1) = paste0("runOne", rownames(seqtab.nochim.pool_sep.r1))
rownames(seqtab.nochim.pool_sep.r2) = unname(sapply(rownames(seqtab.nochim.pool_sep.r2), get.sample.name))
rownames(seqtab.nochim.pool_sep.r2) = paste0("runTwo", rownames(seqtab.nochim.pool_sep.r2))

#join tables
#full
seqtab.nochim.pool_sep = full_join(
    data.frame(asv.seq = rownames(t(seqtab.nochim.pool_sep.r1)), t(seqtab.nochim.pool_sep.r1)),
    data.frame(asv.seq = rownames(t(seqtab.nochim.pool_sep.r2)), t(seqtab.nochim.pool_sep.r2)),
    by = "asv.seq"
)
rownames(seqtab.nochim.pool_sep) = seqtab.nochim.pool_sep$asv.seq
seqtab.nochim.pool_sep$asv.seq = NULL
seqtab.nochim.pool_sep[is.na(seqtab.nochim.pool_sep)] = 0

seq_tab.nochim.files_cat = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.files_cat = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat/intermediate_RDS/taxa_w_bootstraps.rds")

seq_tab.nochim.files_sep = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.files_sep = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep/intermediate_RDS/taxa_w_bootstraps.rds")


#Venn Diagram
require(VennDiagram)
venn.diagram(list(sep_run_pool = colnames(t(seqtab.nochim.pool_sep)), global_pool = colnames(seq_tab.nochim.files_sep), cat_seq_files = colnames(seq_tab.nochim.files_cat)), filename = "compare_dada_processing_figs/shared_ASVs_Venn.tiff")

#ASV richness by method

asv_richness_method = data.frame(
    method = c("sep_run_pool", "global_pool", "cat_seq_files"),
    asv_richness = c(
        length(colnames(t(seqtab.nochim.pool_sep))),
        length(colnames(seq_tab.nochim.files_sep)),
        length(colnames(seq_tab.nochim.files_cat))
    )
)

ggplot(asv_richness_method)


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


#alternate way to get names
rownames(seqtab.nochim.join_by_ASVseq[
    rowSums(seqtab.nochim.join_by_ASVseq[ ,1:length(rownames(seqtab.nochim.r1))
    ]) > 0
,])
rownames( seqtab.nochim.join_by_ASVseq[
    rowSums(seqtab.nochim.join_by_ASVseq[ ,length(rownames(seqtab.nochim.r1)) + 1:length(rownames(seqtab.nochim.r2))
    ]) > 0
,])
