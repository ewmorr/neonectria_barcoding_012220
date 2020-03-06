library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
require(tidyverse)


get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]
fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    #reformat zeros
    l <- gsub("0e\\+00","0",l)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
}

setwd("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320")
system("mkdir dada2_processing_compare")

sample_pairs = read.table("sample_data/run1_run2_sample_pairs.txt")
colnames(sample_pairs) = c("metadata.label", "sample.runOne", "sample.runTwo")


#########################
#Separate pooling by run#
seqtab.nochim.pool_sep.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.pool_sep.r1 = readRDS("~/GARNAS_neonectria_barcoding_091819/intermediate_RDS/taxa_w_bootstraps.rds")

seqtab.nochim.pool_sep.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.pool_sep.r2 = readRDS("~/GARNAS_neonectria_barcoding_012220/intermediate_RDS/taxa_w_bootstraps.rds")

#Join by ASV seq
seqtab.nochim.pool_sep = full_join(
data.frame(asv.seq = rownames(t(seqtab.nochim.pool_sep.r1)), t(seqtab.nochim.pool_sep.r1)),
data.frame(asv.seq = rownames(t(seqtab.nochim.pool_sep.r2)), t(seqtab.nochim.pool_sep.r2)),
by = "asv.seq"
)
rownames(seqtab.nochim.pool_sep) = seqtab.nochim.pool_sep$asv.seq
seqtab.nochim.pool_sep$asv.seq = NULL
seqtab.nochim.pool_sep[is.na(seqtab.nochim.pool_sep)] = 0

#Sum ASV counts by sample pairs
#define empty DF
seqtab.nochim.pool_sep.sample_sum = data.frame(
    matrix(
        data = NA,
        nrow = length(rownames(seqtab.nochim.pool_sep)),
        ncol = length(sample_pairs$metadata.label)
    )
)
colnames(seqtab.nochim.pool_sep.sample_sum) = sample_pairs$metadata.label
rownames(seqtab.nochim.pool_sep.sample_sum) = rownames(seqtab.nochim.pool_sep)


for(i in 1:length(sample_pairs$metadata.label)){
    #check for colanme exists
    if(sample_pairs[i,2] %in% colnames(seqtab.nochim.pool_sep) & sample_pairs[i,3] %in% colnames(seqtab.nochim.pool_sep))
        for(u in 1:length(rownames(seqtab.nochim.pool_sep))){
            seqtab.nochim.pool_sep.sample_sum[u,i] = sum(
                seqtab.nochim.pool_sep[u, colnames(seqtab.nochim.pool_sep) == sample_pairs[i,2] ],
                seqtab.nochim.pool_sep[u, colnames(seqtab.nochim.pool_sep) == sample_pairs[i,3] ]
                )
        }else if(!sample_pairs[i,2] %in% colnames(seqtab.nochim.pool_sep) & !sample_pairs[i,3] %in% colnames(seqtab.nochim.pool_sep)){
            next
        }else if(sample_pairs[i,2] %in% colnames(seqtab.nochim.pool_sep) & !sample_pairs[i,3] %in% colnames(seqtab.nochim.pool_sep)){
            for(u in 1:length(rownames(seqtab.nochim.pool_sep))){
                seqtab.nochim.pool_sep.sample_sum[u,i] = seqtab.nochim.pool_sep[u, colnames(seqtab.nochim.pool_sep) == sample_pairs[i,2] ]
            }
        }else if(!sample_pairs[i,2] %in% colnames(seqtab.nochim.pool_sep) & sample_pairs[i,3] %in% colnames(seqtab.nochim.pool_sep)){
            for(u in 1:length(rownames(seqtab.nochim.pool_sep))){
                seqtab.nochim.pool_sep.sample_sum[u,i] = seqtab.nochim.pool_sep[u, colnames(seqtab.nochim.pool_sep) == sample_pairs[i,3] ]
            }
        }
}

seqtab.nochim.pool_sep.sample_sum[is.na(seqtab.nochim.pool_sep.sample_sum)] = 0
seqtab.nochim.pool_sep.sample_sum = seqtab.nochim.pool_sep.sample_sum[,colSums(seqtab.nochim.pool_sep.sample_sum) > 0]

saveRDS(seqtab.nochim.pool_sep.sample_sum, "run1_run2_dada_compare/sep_run_pool/intermediate_RDS/dada2_seq_table_no_chim_run1run2_pool_sep_samples_summed.rds")

################
#Global pooling#
seqtab.nochim.files_sep = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.files_sep = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep/intermediate_RDS/taxa_w_bootstraps.rds")

seqtab.nochim.files_sep.t = t(seqtab.nochim.files_sep)
#Sum ASV counts by sample pairs
#define empty DF
seqtab.nochim.files_sep.sample_sum = data.frame(
matrix(
data = NA,
nrow = length(rownames(seqtab.nochim.files_sep.t)),
ncol = length(sample_pairs$metadata.label)
)
)
colnames(seqtab.nochim.files_sep.sample_sum) = sample_pairs$metadata.label
rownames(seqtab.nochim.files_sep.sample_sum) = rownames(seqtab.nochim.files_sep.t)

for(i in 1:length(sample_pairs$metadata.label)){
    #check for colanme exists
    if(sample_pairs[i,2] %in% colnames(seqtab.nochim.files_sep.t) & sample_pairs[i,3] %in% colnames(seqtab.nochim.files_sep.t))
    for(u in 1:length(rownames(seqtab.nochim.files_sep.t))){
        seqtab.nochim.files_sep.sample_sum[u,i] = sum(
        seqtab.nochim.files_sep.t[u, colnames(seqtab.nochim.files_sep.t) == sample_pairs[i,2] ][1], #added this last index to deal with the repeat sampling from runOne, this  only looks at the first sequencing event
        seqtab.nochim.files_sep.t[u, colnames(seqtab.nochim.files_sep.t) == sample_pairs[i,3] ][1]
        )
    }else if(!sample_pairs[i,2] %in% colnames(seqtab.nochim.files_sep.t) & !sample_pairs[i,3] %in% colnames(seqtab.nochim.files_sep.t)){
        next
    }else if(sample_pairs[i,2] %in% colnames(seqtab.nochim.files_sep.t) & !sample_pairs[i,3] %in% colnames(seqtab.nochim.files_sep.t)){
        for(u in 1:length(rownames(seqtab.nochim.files_sep.t))){
            seqtab.nochim.files_sep.sample_sum[u,i] = seqtab.nochim.files_sep.t[u, colnames(seqtab.nochim.files_sep.t) == sample_pairs[i,2] ][1]
        }
    }else if(!sample_pairs[i,2] %in% colnames(seqtab.nochim.files_sep.t) & sample_pairs[i,3] %in% colnames(seqtab.nochim.files_sep.t)){
        for(u in 1:length(rownames(seqtab.nochim.files_sep.t))){
            seqtab.nochim.files_sep.sample_sum[u,i] = seqtab.nochim.files_sep.t[u, colnames(seqtab.nochim.files_sep.t) == sample_pairs[i,3] ][1]
        }
    }
}

seqtab.nochim.files_sep.sample_sum[is.na(seqtab.nochim.files_sep.sample_sum)] = 0
seqtab.nochim.files_sep.sample_sum = seqtab.nochim.files_sep.sample_sum[,colSums(seqtab.nochim.files_sep.sample_sum) > 0]

saveRDS(seqtab.nochim.files_sep.sample_sum, "run1_run2_dada_compare/files_sep/intermediate_RDS/dada2_seq_table_no_chim_global_pool_samples_summed.rds")

####################
#Files concatenated#
seq_tab.nochim.files_cat = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.files_cat = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat/intermediate_RDS/taxa_w_bootstraps.rds")

