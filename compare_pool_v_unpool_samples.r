library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
require(tidyverse)

source("~/ggplot_theme.txt")

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

get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

system("mkdir pool_v_unpool_figs")

seq_tab.nochim.files_sep.pool_T = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep/intermediate_RDS/dada2_seq_table_no_chim.rds")
seq_tab.nochim.files_sep.pool_F = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep_pool_F/intermediate_RDS/dada2_seq_table_no_chim.rds")
seq_tab.nochim.files_cat.pool_T = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat/intermediate_RDS/dada2_seq_table_no_chim.rds")
seq_tab.nochim.files_cat.pool_F = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat_pool_F/intermediate_RDS/dada2_seq_table_no_chim.rds")
taxa.w_bootstraps.files_sep.pool_T = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep/intermediate_RDS/taxa_w_bootstraps.rds")
taxa.w_bootstraps.files_sep.pool_F = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_sep_pool_F/intermediate_RDS/taxa_w_bootstraps.rds")
taxa.w_bootstraps.files_cat.pool_T = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat/intermediate_RDS/taxa_w_bootstraps.rds")
taxa.w_bootstraps.files_cat.pool_F = readRDS("~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_dada_compare/files_cat_pool_F/intermediate_RDS/taxa_w_bootstraps.rds")

####################################
#Venn of shared ASVs across pooling#
require(VennDiagram)
#files_sep
venn.diagram(list(pool_T = colnames(seq_tab.nochim.files_sep.pool_T), pool_F = colnames(seq_tab.nochim.files_sep.pool_F)), filename = "pool_v_unpool_figs/files_sep_pooling_shared_ASVs_Venn.tiff")
#files_cat
venn.diagram(list(pool_T = colnames(seq_tab.nochim.files_cat.pool_T), pool_F = colnames(seq_tab.nochim.files_cat.pool_F)), filename = "pool_v_unpool_figs/files_cat_pooling_shared_ASVs_Venn.tiff")

#########################################################
#ASV frequency (i.e., count of ASVs by number of samples#

files_sep.pool_T.ASV_freq = colSums(seq_tab.nochim.files_sep.pool_T > 0) %>% unname
files_sep.pool_F.ASV_freq = colSums(seq_tab.nochim.files_sep.pool_F > 0) %>% unname
files_cat.pool_T.ASV_freq = colSums(seq_tab.nochim.files_cat.pool_T > 0) %>% unname
files_cat.pool_F.ASV_freq = colSums(seq_tab.nochim.files_cat.pool_F > 0) %>% unname

files_sep.asv_freq = rbind(
    data.frame(pooling = "pool_T", sample_count = files_sep.pool_T.ASV_freq),
    data.frame(pooling = "pool_F", sample_count = files_sep.pool_F.ASV_freq)
)
files_cat.asv_freq = rbind(
data.frame(pooling = "pool_T", sample_count = files_cat.pool_T.ASV_freq),
data.frame(pooling = "pool_F", sample_count = files_cat.pool_F.ASV_freq)
)

p = ggplot(files_sep.asv_freq, aes(sample_count)) +
geom_histogram(stat = "bin", binwidth = 1) +
facet_wrap(~pooling) +
scale_y_sqrt(breaks = c(10,100,250,500,750)) +
scale_x_continuous(breaks = c(1,100,200,300)) +
labs(x = "ASV frequency (no. samples observed)", y = "ASV count") +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_ASV_sample_freq.pdf", width = 8, height = 4)
p
dev.off()

p = ggplot(files_cat.asv_freq, aes(sample_count)) +
geom_histogram(stat = "bin", binwidth = 1) +
facet_wrap(~pooling) +
scale_y_sqrt(breaks = c(10,100,250,500,750,1000)) +
scale_x_continuous(breaks = c(1,50,100,150)) +
labs(x = "ASV frequency (no. samples observed)", y = "ASV count") +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_ASV_sample_freq.pdf", width = 8, height = 4)
p
dev.off()

#######################
#Rank abundance curves#

files_sep.pool_T.rank_abd = colSums(seq_tab.nochim.files_sep.pool_T) %>% unname %>% sort(decreasing = T)
files_sep.pool_F.rank_abd = colSums(seq_tab.nochim.files_sep.pool_F) %>% unname %>% sort(decreasing = T)
files_cat.pool_T.rank_abd = colSums(seq_tab.nochim.files_cat.pool_T) %>% unname %>% sort(decreasing = T)
files_cat.pool_F.rank_abd = colSums(seq_tab.nochim.files_cat.pool_F) %>% unname %>% sort(decreasing = T)

files_sep.rank_abd = rbind(
data.frame(pooling = "pool_T", count = files_sep.pool_T.rank_abd, rank = seq(1,length(files_sep.pool_T.rank_abd), by = 1)),
data.frame(pooling = "pool_F", count = files_sep.pool_F.rank_abd, rank = seq(1,length(files_sep.pool_F.rank_abd), by = 1))
)
files_cat.rank_abd = rbind(
data.frame(pooling = "pool_T", count = files_cat.pool_T.rank_abd, rank = seq(1,length(files_cat.pool_T.rank_abd), by = 1)),
data.frame(pooling = "pool_F", count = files_cat.pool_F.rank_abd, rank = seq(1,length(files_cat.pool_F.rank_abd), by = 1))
)

p = ggplot(files_sep.rank_abd, aes(rank, count +1)) +
geom_col(width = 1) +
facet_wrap(~pooling) +
my_gg_theme +
labs(y = "sequence count + 1") +
scale_y_log10(labels = fancy_scientific)

pdf("pool_v_unpool_figs/files_sep_rank_abd.pdf", width = 8, height = 4)
p
dev.off()

p = ggplot(files_cat.rank_abd, aes(rank, count +1)) +
geom_col(width = 1) +
facet_wrap(~pooling) +
my_gg_theme +
labs(y = "sequence count + 1") +
scale_y_log10(labels = fancy_scientific)

pdf("pool_v_unpool_figs/files_cat_rank_abd.pdf", width = 8, height = 4)
p
dev.off()

####################################
#set sample names to match metadata#

rownames(seq_tab.nochim.files_sep.pool_T) = paste0("poolT", rownames(seq_tab.nochim.files_sep.pool_T))

rownames(seq_tab.nochim.files_sep.pool_F) = paste0("poolF", rownames(seq_tab.nochim.files_sep.pool_F))

rownames(seq_tab.nochim.files_cat.pool_T) = paste0("poolT", rownames(seq_tab.nochim.files_cat.pool_T))

rownames(seq_tab.nochim.files_cat.pool_F) = paste0("poolF", rownames(seq_tab.nochim.files_cat.pool_F))

#Read metadata
id_bench_map.files_sep = read.table("sample_data/sample_mapping_pooling_files_sep.txt", header = T)
id_bench_map.files_cat = read.table("sample_data/sample_mapping_pooling_files_cat.txt", header = T)
id_bench_map.files_sep$Pool = as.character(id_bench_map.files_sep$Pool)
id_bench_map.files_cat$Pool = as.character(id_bench_map.files_cat$Pool)

metadata_map = read.table("sample_data/metadata.txt", header = T)
survey_dat = read.table("sample_data/trees_site_survey_data.txt", header = T, sep = "\t")
neo_cov = read.table("sample_data/plug_neonectria_coverage.txt", header = T)
site_info = read.csv("sample_data/site_info.csv")

#colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
metadata_ordered.files_sep = full_join(metadata_map, id_bench_map.files_sep)
metadata_ordered.files_cat = full_join(metadata_map, id_bench_map.files_cat)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>% left_join(., site_info, by = "Site")

full_metadata.files_sep = full_join(metadata_ordered.files_sep, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))
full_metadata.files_cat = full_join(metadata_ordered.files_cat, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))


#############
#join tables#

#full count tab FILES_SEP
seqtab.nochim.files_sep.join_by_ASVseq = full_join(
    data.frame(asv.seq = rownames(t(seq_tab.nochim.files_sep.pool_T)), t(seq_tab.nochim.files_sep.pool_T)),
    data.frame(asv.seq = rownames(t(seq_tab.nochim.files_sep.pool_F)), t(seq_tab.nochim.files_sep.pool_F)),
    by = "asv.seq"
)
rownames(seqtab.nochim.files_sep.join_by_ASVseq) = seqtab.nochim.files_sep.join_by_ASVseq$asv.seq
seqtab.nochim.files_sep.join_by_ASVseq$asv.seq = NULL
seqtab.nochim.files_sep.join_by_ASVseq[is.na(seqtab.nochim.files_sep.join_by_ASVseq)] = 0

#full count tab FILES_CAT
seqtab.nochim.files_cat.join_by_ASVseq = full_join(
data.frame(asv.seq = rownames(t(seq_tab.nochim.files_cat.pool_T)), t(seq_tab.nochim.files_cat.pool_T)),
data.frame(asv.seq = rownames(t(seq_tab.nochim.files_cat.pool_F)), t(seq_tab.nochim.files_cat.pool_F)),
by = "asv.seq"
)
rownames(seqtab.nochim.files_cat.join_by_ASVseq) = seqtab.nochim.files_cat.join_by_ASVseq$asv.seq
seqtab.nochim.files_cat.join_by_ASVseq$asv.seq = NULL
seqtab.nochim.files_cat.join_by_ASVseq[is.na(seqtab.nochim.files_cat.join_by_ASVseq)] = 0

#full taxa table
#Its seems after doing some joins that some taxonomic names are not matching between different seqs, and this may be attribued at the Kingdom level (based on error from full_join). Will need to explore this

#FILES_SEP
asv_tax.files_sep.join_by_ASVseq = rbind(
    data.frame(asv.seq = rownames(taxa.w_bootstraps.files_sep.pool_T$tax), taxa.w_bootstraps.files_sep.pool_T$tax),
    data.frame(asv.seq = rownames(taxa.w_bootstraps.files_sep.pool_F$tax), taxa.w_bootstraps.files_sep.pool_F$tax)
) %>% unique

#convert to character before replacing incorrect names
asv_tax.files_sep.join_by_ASVseq = asv_tax.files_sep.join_by_ASVseq %>%
mutate_all(as.character)

#convert Cylindrocapron faginatum to Neonectria faginata
for(i in 1:length(rownames(asv_tax.files_sep.join_by_ASVseq))){
    if(is.na(asv_tax.files_sep.join_by_ASVseq[i,7]) == T | is.na(asv_tax.files_sep.join_by_ASVseq[i,8]) == T){
        next
    }else
    if(asv_tax.files_sep.join_by_ASVseq[i,7] == "g__Cylindrocarpon" & asv_tax.files_sep.join_by_ASVseq[i,8] == "s__faginatum"){
        print("T")
        asv_tax.files_sep.join_by_ASVseq[i,7] = "g__Neonectria"
        asv_tax.files_sep.join_by_ASVseq[i,8] = "s__faginata"
    }
}

seqtab.nochim.files_sep.join_by_ASVseq.w_tax = full_join(
    data.frame(asv.seq = rownames(seqtab.nochim.files_sep.join_by_ASVseq), seqtab.nochim.files_sep.join_by_ASVseq),
    data.frame(asv.seq = rownames(asv_tax.files_sep.join_by_ASVseq), asv_tax.files_sep.join_by_ASVseq),
    by = "asv.seq"
)

#FILES_CAT
asv_tax.files_cat.join_by_ASVseq = rbind(
data.frame(asv.seq = rownames(taxa.w_bootstraps.files_cat.pool_T$tax), taxa.w_bootstraps.files_cat.pool_T$tax),
data.frame(asv.seq = rownames(taxa.w_bootstraps.files_cat.pool_F$tax), taxa.w_bootstraps.files_cat.pool_F$tax)
) %>% unique

#convert to character before replacing incorrect names
asv_tax.files_cat.join_by_ASVseq = asv_tax.files_cat.join_by_ASVseq %>%
mutate_all(as.character)

#convert Cylindrocapron faginatum to Neonectria faginata
for(i in 1:length(rownames(asv_tax.files_cat.join_by_ASVseq))){
    if(is.na(asv_tax.files_cat.join_by_ASVseq[i,7]) == T | is.na(asv_tax.files_cat.join_by_ASVseq[i,8]) == T){
        next
    }else
    if(asv_tax.files_cat.join_by_ASVseq[i,7] == "g__Cylindrocarpon" & asv_tax.files_cat.join_by_ASVseq[i,8] == "s__faginatum"){
        print("T")
        asv_tax.files_cat.join_by_ASVseq[i,7] = "g__Neonectria"
        asv_tax.files_cat.join_by_ASVseq[i,8] = "s__faginata"
    }
}

seqtab.nochim.files_cat.join_by_ASVseq.w_tax = full_join(
data.frame(asv.seq = rownames(seqtab.nochim.files_cat.join_by_ASVseq), seqtab.nochim.files_cat.join_by_ASVseq),
data.frame(asv.seq = rownames(asv_tax.files_cat.join_by_ASVseq), asv_tax.files_cat.join_by_ASVseq),
by = "asv.seq"
)

###################
#Neonectria ASVs

#FILES_SEP
asv_tax.files_sep.join_by_ASVseq.Nf = filter(asv_tax.files_sep.join_by_ASVseq, Genus == "g__Neonectria" & Species == "s__faginata")
asv_tax.files_sep.join_by_ASVseq.Nf$asv.seq %>% length
asv_tax.files_sep.join_by_ASVseq.Nf$asv.seq %>% unique %>% length

asv_tax.files_sep.join_by_ASVseq.Nd = filter(asv_tax.files_sep.join_by_ASVseq, Genus == "g__Neonectria" & Species == "s__ditissima")
asv_tax.files_sep.join_by_ASVseq.Nd$asv.seq %>% length
asv_tax.files_sep.join_by_ASVseq.Nd$asv.seq %>% unique %>% length


#FILES_CAT
asv_tax.files_cat.join_by_ASVseq.Nf = filter(asv_tax.files_cat.join_by_ASVseq, Genus == "g__Neonectria" & Species == "s__faginata")
asv_tax.files_cat.join_by_ASVseq.Nf$asv.seq %>% length
asv_tax.files_cat.join_by_ASVseq.Nf$asv.seq %>% unique %>% length

asv_tax.files_cat.join_by_ASVseq.Nd = filter(asv_tax.files_cat.join_by_ASVseq, Genus == "g__Neonectria" & Species == "s__ditissima")
asv_tax.files_cat.join_by_ASVseq.Nd$asv.seq %>% length
asv_tax.files_cat.join_by_ASVseq.Nd$asv.seq %>% unique %>% length


########
#Run through comps with files_sep before modify for Files_cat

#FILES_SEP

######################################################################
#tables of pool_F vs. pool_T values by sample (i.e. metadata.label)#

#total sequence counts after full dada2
sample_counts = data.frame(sample = colnames(seqtab.nochim.files_sep.join_by_ASVseq), count = colSums(seqtab.nochim.files_sep.join_by_ASVseq))
sample_counts.meta = left_join(
    full_metadata.files_sep %>%
        filter(seq.rep == "n" & locus == "ITS2" & resequenced == "n" & !is.na(metadata.label)) %>%
        select(c("metadata.label", "sample", "run.seq", "Pool")) ,
    sample_counts,
    by = "sample"
)
sample_counts.meta[is.na(sample_counts.meta)] = 0

sample_counts.meta.wide_run = pivot_wider(sample_counts.meta, id_cols = c(run.seq, metadata.label), names_from = Pool, values_from = count)

qqnorm(residuals(lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run)))
plot(residuals(lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run)))
summary(lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run))

fit = lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run)

p = ggplot(sample_counts.meta.wide_run, aes(x = pool_F, y = pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_log10(labels = fancy_scientific) +
#scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F sequence count", y = "pool_T sequence count",
title = paste(
"Sequence count by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
" Slope = ",signif(fit$coef[[2]], 2),
" P =",signif(summary(fit)$coef[2,4], 3)
)
) +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_pooling_seq_count.pdf", width = 8, height = 6)
print(p)
dev.off()

#############
#Nf rel.abd.#
seqtab.nochim.join_by_ASVseq.Nf = data.frame(asv.seq = rownames(seqtab.nochim.files_sep.join_by_ASVseq), seqtab.nochim.files_sep.join_by_ASVseq) %>%
    filter(asv.seq %in% asv_tax.files_sep.join_by_ASVseq.Nf$asv.seq)
#seqtab.nochim.join_by_ASVseq.Nf[is.na(seqtab.nochim.join_by_ASVseq.Nf)] = 0
seqtab.nochim.join_by_ASVseq.Nf$asv.seq = NULL

Nf_count = data.frame(
Nf = colSums(seqtab.nochim.join_by_ASVseq.Nf),
sample = colnames(seqtab.nochim.files_sep.join_by_ASVseq),
richness = colSums(seqtab.nochim.files_sep.join_by_ASVseq > 0),
count = colSums(seqtab.nochim.files_sep.join_by_ASVseq)
)

Nf_count.meta = left_join(
full_metadata.files_sep %>%
filter(seq.rep == "n" & locus == "ITS2" & resequenced == "n" & !is.na(metadata.label)) %>%
select(c("metadata.label", "sample", "run.seq", "Pool")) ,
Nf_count,
by = "sample"
)
Nf_count.meta[is.na(Nf_count.meta)] = 0

Nf_count.meta.wide = pivot_wider(Nf_count.meta, id_cols = c(metadata.label, run.seq), names_from = Pool, values_from = c(count, richness, Nf))

p = ggplot(Nf_count.meta.wide %>% filter(count_pool_F != 0 & count_pool_T != 0), aes(Nf_pool_F/count_pool_F, Nf_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_log10(labels = fancy_scientific) +
#scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F N. faginata rel. abd.", y = "pool_T N. faginata rel. abd.") +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_pooling_Nf_RA.pdf", width = 8, height = 6)
print(p)
dev.off()


#############
#Nd rel.abd.#
seqtab.nochim.join_by_ASVseq.Nd = data.frame(asv.seq = rownames(seqtab.nochim.files_sep.join_by_ASVseq), seqtab.nochim.files_sep.join_by_ASVseq) %>%
filter(asv.seq %in% asv_tax.files_sep.join_by_ASVseq.Nd$asv.seq)
seqtab.nochim.join_by_ASVseq.Nd[is.na(seqtab.nochim.join_by_ASVseq.Nd)] = 0
seqtab.nochim.join_by_ASVseq.Nd$asv.seq = NULL

Nd_count = data.frame(
Nd = colSums(seqtab.nochim.join_by_ASVseq.Nd),
sample = colnames(seqtab.nochim.files_sep.join_by_ASVseq),
richness = colSums(seqtab.nochim.files_sep.join_by_ASVseq > 0),
count = colSums(seqtab.nochim.files_sep.join_by_ASVseq)
)

Nd_count.meta = left_join(
full_metadata.files_sep %>%
filter(seq.rep == "n" & locus == "ITS2" & resequenced == "n" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "Pool", "run.seq")) ,
Nd_count,
by = "sample"
)
Nd_count.meta[is.na(Nd_count.meta)] = 0

Nd_count.meta.wide = pivot_wider(Nd_count.meta, id_cols = c(run.seq,metadata.label), names_from = Pool, values_from = c(count, richness, Nd))

p = ggplot(Nd_count.meta.wide %>% filter(count_pool_T != 0 & count_pool_F != 0), aes(Nd_pool_F/count_pool_F, Nd_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_sqrt(labels = fancy_scientific) +
#scale_y_sqrt(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F N. ditissima rel. abd.", y = "pool_T N. ditissima rel. abd.") +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_Nd_RA.pdf", width = 8, height = 6)
print(p)
dev.off()


########################################################
#Pairwise asv_tabs for rarefying samples by run minimum#
########################################################
require(vegan)

#table for selecting sample pairs
sample_map.wide_run = full_metadata.files_sep %>% filter(!is.na(metadata.label) & resequenced == "n" & locus == "ITS2") %>% select(c(sample, run.seq, metadata.label, Pool))  %>%
    pivot_wider(., id_cols = c(run.seq,metadata.label), names_from = Pool, values_from = sample)

sample_map.wide_run = sample_map.wide_run %>% filter(!is.na(pool_F) & !is.na(pool_T) )

seqtab.nochim.join_by_ASVseq.t = t(seqtab.nochim.files_sep.join_by_ASVseq)
seqtab.nochim.join_by_ASVseq.t.sample = data.frame(sample = rownames(seqtab.nochim.join_by_ASVseq.t), seqtab.nochim.join_by_ASVseq.t)

#empty matrix with rows as asvs and cols as samples
seqtab.nochim.sample_rare_by_run_min = matrix(
    nrow = length(rownames(seqtab.nochim.join_by_ASVseq.t)),
    ncol = length(colnames(seqtab.nochim.join_by_ASVseq.t)),
    dimnames = list(rownames(seqtab.nochim.join_by_ASVseq.t), colnames(seqtab.nochim.join_by_ASVseq.t))
)

#simplest method is to call row by var name (bc dplyr filter NSE variable interpolation)
#This will break if the rownames in the empty matrix are not properly names, so need to run the above call if the matrix is modified
for(i in 1:length(sample_map.wide_run$metadata.label)){
    #skip sample pairs where one is missing
    if(any(rownames(seqtab.nochim.join_by_ASVseq.t) == sample_map.wide_run$pool_F[i]) == FALSE |
    any(rownames(seqtab.nochim.join_by_ASVseq.t) == sample_map.wide_run$pool_T[i]) == FALSE){
        next
    }
    
    asv_tab.temp.pool_F = seqtab.nochim.join_by_ASVseq.t[as.character(sample_map.wide_run$pool_F[i]), ]
    asv_tab.temp.pool_T = seqtab.nochim.join_by_ASVseq.t[as.character(sample_map.wide_run$pool_T[i]), ]
    min_seqs = min(sum(asv_tab.temp.pool_F), sum(asv_tab.temp.pool_T))
    
    asv_tab.temp = rbind(asv_tab.temp.pool_F, asv_tab.temp.pool_T)
    rownames(asv_tab.temp) = c(as.character(sample_map.wide_run$pool_F[i]), as.character(sample_map.wide_run$pool_T[i]))
    asv_tab.temp.rare = rrarefy(asv_tab.temp, min_seqs)
    seqtab.nochim.sample_rare_by_run_min[c(as.character(sample_map.wide_run$pool_F[i]), as.character(sample_map.wide_run$pool_T[i])),] =
        asv_tab.temp.rare[c(as.character(sample_map.wide_run$pool_F[i]), as.character(sample_map.wide_run$pool_T[i])),]
}

seqtab.nochim.sample_rare_by_run_min.rmNA = drop_na(as.data.frame(seqtab.nochim.sample_rare_by_run_min))


#####################################################
#Calculate vals and create pairwise DFs for rarefied#

######################
#rarefied Simoson's D#
seqtab.nochim.sample_rare_by_run_min.rmNA.simp = diversity(seqtab.nochim.sample_rare_by_run_min.rmNA, "simpson")


asv_simpson = data.frame(
sample = names(seqtab.nochim.sample_rare_by_run_min.rmNA.simp),
simp.D = unname(1 - seqtab.nochim.sample_rare_by_run_min.rmNA.simp),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_simpson.meta = left_join(asv_simpson,
full_metadata.files_sep %>%
filter(resequenced == "n" & locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "run.seq", "Pool")),
by = "sample"
)

asv_simpson.meta.wide = pivot_wider(asv_simpson.meta, id_cols = c(metadata.label, run.seq), names_from = Pool, values_from = c(simp.D, count))

qqnorm(residuals(lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100))))
plot(residuals(lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100))))
summary(lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100)))

fit = lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100) )

p = ggplot(asv_simpson.meta.wide  %>% filter(count_pool_F >= 100), aes(simp.D_pool_F, simp.D_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "pool_F Simpson D", y = "pool_T Simpson D",
    title = paste(
        "Simpson D by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
        " Slope = ",signif(fit$coef[[2]], 2),
        " P =",signif(summary(fit)$coef[2,4], 3)
    )
)

pdf("pool_v_unpool_figs/files_sep_Simpson_D.pdf", width = 8, height = 6)
print(p)
dev.off()

######################
#rarefied Shannon H#
seqtab.nochim.sample_rare_by_run_min.rmNA.shan = diversity(seqtab.nochim.sample_rare_by_run_min.rmNA, "shannon")


asv_shannon = data.frame(
sample = names(seqtab.nochim.sample_rare_by_run_min.rmNA.shan),
shannon = unname(seqtab.nochim.sample_rare_by_run_min.rmNA.shan),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_shannon.meta = left_join(asv_shannon,
full_metadata.files_sep %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "run.seq", "Pool")),
by = "sample"
)

asv_shannon.meta.wide = pivot_wider(asv_shannon.meta, id_cols = c(run.seq, metadata.label), names_from = Pool, values_from = c(shannon, count))

qqnorm(residuals(lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100))))
plot(residuals(lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100))))
summary(lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100)))

fit = lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100) )

p = ggplot(asv_shannon.meta.wide %>% filter(count_pool_F >= 100), aes(shannon_pool_F, shannon_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "pool_F Shannon H", y = "pool_T Shannon H",
title = paste(
    "Shannon H by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
    " Slope = ",signif(fit$coef[[2]], 2),
    " P =",signif(summary(fit)$coef[2,4], 3)
    )
)


pdf("pool_v_unpool_figs/files_sep_shannon.pdf", width = 8, height = 6)
print(p)
dev.off()


##############
#asv richness#
asv_richness = data.frame(
sample = rownames(seqtab.nochim.sample_rare_by_run_min.rmNA),
richness = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA > 0),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_richness.meta = left_join(asv_richness,
full_metadata.files_sep %>%
filter(seq.rep == "n" & locus == "ITS2") %>%
select(c("metadata.label", "sample", "run.seq", "Pool")) ,
by = "sample"
)
asv_richness.meta[is.na(asv_richness.meta)] = 0

asv_richness.meta.wide = pivot_wider(asv_richness.meta, id_cols = c(run.seq, metadata.label), names_from = Pool, values_from = c(count, richness))

qqnorm(residuals(lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
plot(residuals(lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
summary(lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

fit = lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))

require(gridExtra)

p1 = ggplot(asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(richness_pool_F, richness_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "pool_F ASV richness", y = "pool_T ASV richness",
    title = paste(
        "ASV richness by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
        " Slope = ",signif(fit$coef[[2]], 2),
        " P =",signif(summary(fit)$coef[2,4], 3)
    )
)

p2 = ggplot(asv_richness.meta.wide %>% filter(count_pool_F != 0 & count_pool_T != 0), aes(richness_pool_F, richness_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F ASV richness log10+1", y = "pool_T ASV richness log10+1", title = "ASV richness by sample log10") +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_richness.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#############
#Nf rel.abd.#
seqtab.nochim.join_by_ASVseq.Nf = data.frame(asv.seq = rownames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)), t(seqtab.nochim.sample_rare_by_run_min.rmNA)) %>%
filter(asv.seq %in% asv_tax.files_sep.join_by_ASVseq.Nf$asv.seq)
seqtab.nochim.join_by_ASVseq.Nf[is.na(seqtab.nochim.join_by_ASVseq.Nf)] = 0
seqtab.nochim.join_by_ASVseq.Nf$asv.seq = NULL

Nf_count = data.frame(
Nf = colSums(seqtab.nochim.join_by_ASVseq.Nf),
sample = colnames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)),
richness = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA) > 0),
count = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA))
)

Nf_count.meta = left_join(
Nf_count,
full_metadata.files_sep %>%
filter(resequenced == "n" & locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "run.seq", "Pool")) ,
by = "sample"
)
Nf_count.meta[is.na(Nf_count.meta)] = 0

Nf_count.meta.wide = pivot_wider(Nf_count.meta, id_cols = c(run.seq,metadata.label), names_from = Pool, values_from = c(count, richness, Nf))

qqnorm(residuals(lm(log10(Nf_pool_T+1) ~ log10(Nf_pool_F+1), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
plot(residuals(lm(log10(Nf_pool_T+1) ~ log10(Nf_pool_F+1), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
summary(lm(log10(Nf_pool_T+1) ~ log10(Nf_pool_F+1), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

qqnorm(residuals(lm(Nf_pool_T/count_pool_T ~ I(Nf_pool_F/count_pool_F), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
summary(lm(Nf_pool_T/count_pool_T ~ I(Nf_pool_F/count_pool_F), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

#proportional counts
Nf.binomial_glm = glm(Nf_pool_T/count_pool_T ~ I(Nf_pool_F/count_pool_F*100), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), family = binomial)
1-(Nf.binomial_glm$deviance/Nf.binomial_glm$null.deviance
)
1 - pchisq(summary(Nf.binomial_glm)$deviance,
summary(Nf.binomial_glm)$df.residual
)

p1 = ggplot(Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nf_pool_F/count_pool_F, Nf_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
my_gg_theme +
labs(x = "pool_F N. faginata rel. abd.", y = "pool_T N. faginata rel. abd.",
title = paste(
"N. faginata rel abd by sample\nR2 = ", signif(1-(Nf.binomial_glm$deviance/Nf.binomial_glm$null.deviance), 2),
" P =",signif(summary(Nf.binomial_glm)$coef[2,4], 3)
)
)

p2 = ggplot(Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nf_pool_F+1, Nf_pool_T+1)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "pool_F N. faginata count log10+1", y = "pool_T N. faginata count log10+1", title = "N. faginata counts") +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_Nf_RA_rarefied.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#############
#Nd rel.abd.#
seqtab.nochim.join_by_ASVseq.Nd = data.frame(asv.seq = rownames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)), t(seqtab.nochim.sample_rare_by_run_min.rmNA)) %>%
filter(asv.seq %in% asv_tax.files_sep.join_by_ASVseq.Nd$asv.seq)
seqtab.nochim.join_by_ASVseq.Nd[is.na(seqtab.nochim.join_by_ASVseq.Nd)] = 0
seqtab.nochim.join_by_ASVseq.Nd$asv.seq = NULL

Nd_count = data.frame(
Nd = colSums(seqtab.nochim.join_by_ASVseq.Nd),
sample = colnames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)),
richness = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA) > 0),
count = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA))
)

Nd_count.meta = left_join(
Nd_count,
full_metadata.files_sep %>%
filter(resequenced == "n" & locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "run.seq", "Pool")) ,
by = "sample"
)
Nd_count.meta[is.na(Nd_count.meta)] = 0

Nd_count.meta.wide = pivot_wider(Nd_count.meta, id_cols = c(metadata.label,run.seq), names_from = Pool, values_from = c(count, richness, Nd))

qqnorm(residuals(lm(log10(Nd_pool_T+1) ~ log10(Nd_pool_F+1), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100 & Nd_pool_F+Nd_pool_T != 0))))
plot(residuals(lm(log10(Nd_pool_T+1) ~ log10(Nd_pool_F+1), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100 & Nd_pool_F+Nd_pool_T != 0))))
summary(lm(log10(Nd_pool_T+1) ~ log10(Nd_pool_F+1), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100 & Nd_pool_F+Nd_pool_T != 0)))
summary(lm(Nd_pool_T/count_pool_T ~ I(Nd_pool_F/count_pool_F), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

#proportional counts
Nd.binomial_glm = glm(Nd_pool_T/count_pool_T ~ I(Nd_pool_F/count_pool_F), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), family = binomial)
1-(Nd.binomial_glm$deviance/Nd.binomial_glm$null.deviance
)
1 - pchisq(summary(Nd.binomial_glm)$deviance,
summary(Nd.binomial_glm)$df.residual
)

p1 = ggplot(Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nd_pool_F/count_pool_F, Nd_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
#scale_x_sqrt(labels = fancy_scientific) +
#scale_y_sqrt(labels = fancy_scientific) +
my_gg_theme +
labs(x = "pool_F N. ditissima rel. abd.", y = "pool_T N. ditissima rel. abd.",
title = paste(
"N. ditissima rel abd by sample\nR2 = ", signif(1-(Nd.binomial_glm$deviance/Nd.binomial_glm$null.deviance), 2),
" P =",signif(summary(Nd.binomial_glm)$coef[2,4], 1)
)
)

p2 = ggplot(Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nd_pool_F+1, Nd_pool_T+1)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "pool_F N. ditissima count log10+1", y = "pool_T N. ditissima count log10+1", title = "N. ditissima counts") +
my_gg_theme

pdf("pool_v_unpool_figs/files_sep_Nd_RA_rarefied.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()





###############################################
#Calculate Nf and Nd detection agreement between pooling#

#Nf

Nf_pooling_agree = data.frame(
    pooling = c("pool T and F \nNf present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNf absent"),
    count = vector(length = 4, mode = "numeric")
)
Nf_pooling_agree$pooling = factor(Nf_pooling_agree$pooling, levels = c("pool T and F \nNf present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNf absent"))

agree_or_not = function(x){
    temp.mat = c(0,0,0,0)
    if(x[1] > 0 & x[2] > 0){
        temp.mat[1] = temp.mat[1]+1
    }else if(x[1] == 0 & x[2] == 0){
        temp.mat[4] = temp.mat[4]+1
    }else if(x[1] > 0 & x[2] == 0){
        temp.mat[2] = temp.mat[2]+1
    }else if(x[1] == 0 & x[2] > 0){
        temp.mat[3] = temp.mat[3]+1
    }
    return(temp.mat)
}

Nf_pooling_agree[,2] = rowSums(apply(data.frame(Nf_count.meta.wide$Nf_pool_T, Nf_count.meta.wide$Nf_pool_F), 1, agree_or_not) )

p1 = ggplot(Nf_pooling_agree, aes(x = pooling, y = count)) +
geom_col() +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1)
) +
labs(x = "", y = "sample count")

pdf("pool_v_unpool_figs/Nf_detection_by_pooling_files_sep.pdf")
p1
dev.off()

#Nd

Nd_pooling_agree = data.frame(
pooling = c("pool T and F \nNd present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNd absent"),
count = vector(length = 4, mode = "numeric")
)
Nd_pooling_agree$pooling = factor(Nd_pooling_agree$pooling, levels = c("pool T and F \nNd present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNd absent"))

agree_or_not = function(x){
    temp.mat = c(0,0,0,0)
    if(x[1] > 0 & x[2] > 0){
        temp.mat[1] = temp.mat[1]+1
    }else if(x[1] == 0 & x[2] == 0){
        temp.mat[4] = temp.mat[4]+1
    }else if(x[1] > 0 & x[2] == 0){
        temp.mat[2] = temp.mat[2]+1
    }else if(x[1] == 0 & x[2] > 0){
        temp.mat[3] = temp.mat[3]+1
    }
    return(temp.mat)
}

Nd_pooling_agree[,2] = rowSums(apply(data.frame(Nd_count.meta.wide$Nd_pool_T, Nd_count.meta.wide$Nd_pool_F), 1, agree_or_not) )

p1 = ggplot(Nd_pooling_agree, aes(x = pooling, y = count)) +
geom_col() +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1)
) +
labs(x = "", y = "sample count")

pdf("pool_v_unpool_figs/Nd_detection_by_pooling_files_sep.pdf")
p1
dev.off()







#################################################
#################################################
#################################################
#CODE FROM ABOVE MODIFIED FOR FILES_CAT
#################################################
#################################################
#################################################

#FILES_CAT_START

######################################################################
#tables of pool_F vs. pool_T values by sample (i.e. metadata.label)#

#total sequence counts after full dada2
sample_counts = data.frame(sample = colnames(seqtab.nochim.files_cat.join_by_ASVseq), count = colSums(seqtab.nochim.files_cat.join_by_ASVseq))
sample_counts.meta = left_join(
full_metadata.files_cat %>%
filter(locus == "ITS2", !is.na(metadata.label)) %>%
select(c("metadata.label", "sample", "Pool")) ,
sample_counts,
by = "sample"
)
sample_counts.meta[is.na(sample_counts.meta)] = 0

sample_counts.meta.wide_run = pivot_wider(sample_counts.meta, id_cols = metadata.label, names_from = Pool, values_from = count)

qqnorm(residuals(lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run)))
plot(residuals(lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run)))
summary(lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run))

fit = lm(pool_T ~ pool_F, data = sample_counts.meta.wide_run)

p = ggplot(sample_counts.meta.wide_run, aes(x = pool_F, y = pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_log10(labels = fancy_scientific) +
#scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F sequence count", y = "pool_T sequence count",
title = paste(
"Sequence count by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
" Slope = ",signif(fit$coef[[2]], 2),
" P =",signif(summary(fit)$coef[2,4], 3)
)
) +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_pooling_seq_count.pdf", width = 8, height = 6)
print(p)
dev.off()

#############
#Nf rel.abd.#
seqtab.nochim.join_by_ASVseq.Nf = data.frame(asv.seq = rownames(seqtab.nochim.files_cat.join_by_ASVseq), seqtab.nochim.files_cat.join_by_ASVseq) %>%
filter(asv.seq %in% asv_tax.files_cat.join_by_ASVseq.Nf$asv.seq)
#seqtab.nochim.join_by_ASVseq.Nf[is.na(seqtab.nochim.join_by_ASVseq.Nf)] = 0
seqtab.nochim.join_by_ASVseq.Nf$asv.seq = NULL

Nf_count = data.frame(
Nf = colSums(seqtab.nochim.join_by_ASVseq.Nf),
sample = colnames(seqtab.nochim.files_cat.join_by_ASVseq),
richness = colSums(seqtab.nochim.files_cat.join_by_ASVseq > 0),
count = colSums(seqtab.nochim.files_cat.join_by_ASVseq)
)

Nf_count.meta = left_join(
full_metadata.files_cat %>%
filter(locus == "ITS2" & !is.na(metadata.label)) %>%
select(c("metadata.label", "sample", "Pool")) ,
Nf_count,
by = "sample"
)
Nf_count.meta[is.na(Nf_count.meta)] = 0

Nf_count.meta.wide = pivot_wider(Nf_count.meta, id_cols = metadata.label, names_from = Pool, values_from = c(count, richness, Nf))

p = ggplot(Nf_count.meta.wide %>% filter(count_pool_F != 0 & count_pool_T != 0), aes(Nf_pool_F/count_pool_F, Nf_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_log10(labels = fancy_scientific) +
#scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F N. faginata rel. abd.", y = "pool_T N. faginata rel. abd.") +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_pooling_Nf_RA.pdf", width = 8, height = 6)
print(p)
dev.off()


#############
#Nd rel.abd.#
seqtab.nochim.join_by_ASVseq.Nd = data.frame(asv.seq = rownames(seqtab.nochim.files_cat.join_by_ASVseq), seqtab.nochim.files_cat.join_by_ASVseq) %>%
filter(asv.seq %in% asv_tax.files_cat.join_by_ASVseq.Nd$asv.seq)
seqtab.nochim.join_by_ASVseq.Nd[is.na(seqtab.nochim.join_by_ASVseq.Nd)] = 0
seqtab.nochim.join_by_ASVseq.Nd$asv.seq = NULL

Nd_count = data.frame(
Nd = colSums(seqtab.nochim.join_by_ASVseq.Nd),
sample = colnames(seqtab.nochim.files_cat.join_by_ASVseq),
richness = colSums(seqtab.nochim.files_cat.join_by_ASVseq > 0),
count = colSums(seqtab.nochim.files_cat.join_by_ASVseq)
)

Nd_count.meta = left_join(
full_metadata.files_cat %>%
filter(locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "Pool")) ,
Nd_count,
by = "sample"
)
Nd_count.meta[is.na(Nd_count.meta)] = 0

Nd_count.meta.wide = pivot_wider(Nd_count.meta, id_cols = metadata.label, names_from = Pool, values_from = c(count, richness, Nd))

p = ggplot(Nd_count.meta.wide %>% filter(count_pool_T != 0 & count_pool_F != 0), aes(Nd_pool_F/count_pool_F, Nd_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
#scale_x_sqrt(labels = fancy_scientific) +
#scale_y_sqrt(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F N. ditissima rel. abd.", y = "pool_T N. ditissima rel. abd.") +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_Nd_RA.pdf", width = 8, height = 6)
print(p)
dev.off()


########################################################
#Pairwise asv_tabs for rarefying samples by run minimum#
########################################################
require(vegan)

#table for selecting sample pairs
sample_map.wide_run = full_metadata.files_cat %>% filter(!is.na(metadata.label) & locus == "ITS2") %>% select(c(sample, metadata.label, Pool))  %>%
pivot_wider(., id_cols = c(metadata.label), names_from = Pool, values_from = sample)

sample_map.wide_run = sample_map.wide_run %>% filter(!is.na(pool_F) & !is.na(pool_T) )

seqtab.nochim.join_by_ASVseq.t = t(seqtab.nochim.files_cat.join_by_ASVseq)
seqtab.nochim.join_by_ASVseq.t.sample = data.frame(sample = rownames(seqtab.nochim.join_by_ASVseq.t), seqtab.nochim.join_by_ASVseq.t)

#empty matrix with rows as asvs and cols as samples
seqtab.nochim.sample_rare_by_run_min = matrix(
nrow = length(rownames(seqtab.nochim.join_by_ASVseq.t)),
ncol = length(colnames(seqtab.nochim.join_by_ASVseq.t)),
dimnames = list(rownames(seqtab.nochim.join_by_ASVseq.t), colnames(seqtab.nochim.join_by_ASVseq.t))
)

#simplest method is to call row by var name (bc dplyr filter NSE variable interpolation)
#This will break if the rownames in the empty matrix are not properly names, so need to run the above call if the matrix is modified
for(i in 1:length(sample_map.wide_run$metadata.label)){
    #skip sample pairs where one is missing
    if(any(rownames(seqtab.nochim.join_by_ASVseq.t) == sample_map.wide_run$pool_F[i]) == FALSE |
    any(rownames(seqtab.nochim.join_by_ASVseq.t) == sample_map.wide_run$pool_T[i]) == FALSE){
        next
    }
    
    asv_tab.temp.pool_F = seqtab.nochim.join_by_ASVseq.t[as.character(sample_map.wide_run$pool_F[i]), ]
    asv_tab.temp.pool_T = seqtab.nochim.join_by_ASVseq.t[as.character(sample_map.wide_run$pool_T[i]), ]
    min_seqs = min(sum(asv_tab.temp.pool_F), sum(asv_tab.temp.pool_T))
    
    asv_tab.temp = rbind(asv_tab.temp.pool_F, asv_tab.temp.pool_T)
    rownames(asv_tab.temp) = c(as.character(sample_map.wide_run$pool_F[i]), as.character(sample_map.wide_run$pool_T[i]))
    asv_tab.temp.rare = rrarefy(asv_tab.temp, min_seqs)
    seqtab.nochim.sample_rare_by_run_min[c(as.character(sample_map.wide_run$pool_F[i]), as.character(sample_map.wide_run$pool_T[i])),] =
    asv_tab.temp.rare[c(as.character(sample_map.wide_run$pool_F[i]), as.character(sample_map.wide_run$pool_T[i])),]
}

seqtab.nochim.sample_rare_by_run_min.rmNA = drop_na(as.data.frame(seqtab.nochim.sample_rare_by_run_min))


#####################################################
#Calculate vals and create pairwise DFs for rarefied#

######################
#rarefied Simoson's D#
seqtab.nochim.sample_rare_by_run_min.rmNA.simp = diversity(seqtab.nochim.sample_rare_by_run_min.rmNA, "simpson")


asv_simpson = data.frame(
sample = names(seqtab.nochim.sample_rare_by_run_min.rmNA.simp),
simp.D = unname(1 - seqtab.nochim.sample_rare_by_run_min.rmNA.simp),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_simpson.meta = left_join(asv_simpson,
full_metadata.files_cat %>%
filter(locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "Pool")),
by = "sample"
)

asv_simpson.meta.wide = pivot_wider(asv_simpson.meta, id_cols = metadata.label, names_from = Pool, values_from = c(simp.D, count))

qqnorm(residuals(lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100))))
plot(residuals(lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100))))
summary(lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100)))

fit = lm(simp.D_pool_T ~ simp.D_pool_F, data = asv_simpson.meta.wide %>% filter(count_pool_F >= 100) )

p = ggplot(asv_simpson.meta.wide  %>% filter(count_pool_F >= 100), aes(simp.D_pool_F, simp.D_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "pool_F Simpson D", y = "pool_T Simpson D",
title = paste(
"Simpson D by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
" Slope = ",signif(fit$coef[[2]], 2),
" P =",signif(summary(fit)$coef[2,4], 3)
)
)

pdf("pool_v_unpool_figs/files_cat_Simpson_D.pdf", width = 8, height = 6)
print(p)
dev.off()

######################
#rarefied Shannon H#
seqtab.nochim.sample_rare_by_run_min.rmNA.shan = diversity(seqtab.nochim.sample_rare_by_run_min.rmNA, "shannon")


asv_shannon = data.frame(
sample = names(seqtab.nochim.sample_rare_by_run_min.rmNA.shan),
shannon = unname(seqtab.nochim.sample_rare_by_run_min.rmNA.shan),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_shannon.meta = left_join(asv_shannon,
full_metadata.files_cat %>%
filter(!is.na(metadata.label) & locus == "ITS2") %>%
dplyr::select(c("metadata.label", "sample", "Pool")),
by = "sample"
)

asv_shannon.meta.wide = pivot_wider(asv_shannon.meta, id_cols = metadata.label, names_from = Pool, values_from = c(shannon, count))

qqnorm(residuals(lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100))))
plot(residuals(lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100))))
summary(lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100)))

fit = lm(shannon_pool_T ~ shannon_pool_F, data = asv_shannon.meta.wide %>% filter(count_pool_F >= 100) )

p = ggplot(asv_shannon.meta.wide %>% filter(count_pool_F >= 100), aes(shannon_pool_F, shannon_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "pool_F Shannon H", y = "pool_T Shannon H",
title = paste(
"Shannon H by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
" Slope = ",signif(fit$coef[[2]], 2),
" P =",signif(summary(fit)$coef[2,4], 3)
)
)


pdf("pool_v_unpool_figs/files_cat_shannon.pdf", width = 8, height = 6)
print(p)
dev.off()


##############
#asv richness#
asv_richness = data.frame(
sample = rownames(seqtab.nochim.sample_rare_by_run_min.rmNA),
richness = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA > 0),
count = rowSums(seqtab.nochim.sample_rare_by_run_min.rmNA)
)

asv_richness.meta = left_join(asv_richness,
full_metadata.files_cat %>%
filter(!is.na(metadata.label) & locus == "ITS2") %>%
select(c("metadata.label", "sample", "Pool")) ,
by = "sample"
)
asv_richness.meta[is.na(asv_richness.meta)] = 0

asv_richness.meta.wide = pivot_wider(asv_richness.meta, id_cols = metadata.label, names_from = Pool, values_from = c(count, richness))

qqnorm(residuals(lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
plot(residuals(lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
summary(lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

fit = lm(richness_pool_T ~ richness_pool_F, data = asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))

require(gridExtra)

p1 = ggplot(asv_richness.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(richness_pool_F, richness_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
geom_smooth(method = "lm") +
my_gg_theme +
labs(x = "pool_F ASV richness", y = "pool_T ASV richness",
title = paste(
"ASV richness by sample\nR2 = ", signif(summary(fit)$adj.r.squared, 2),
" Slope = ",signif(fit$coef[[2]], 2),
" P =",signif(summary(fit)$coef[2,4], 3)
)
)

p2 = ggplot(asv_richness.meta.wide %>% filter(count_pool_F != 0 & count_pool_T != 0), aes(richness_pool_F, richness_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
geom_smooth(method = "lm") +
labs(x = "pool_F ASV richness log10+1", y = "pool_T ASV richness log10+1", title = "ASV richness by sample log10") +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_richness.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#############
#Nf rel.abd.#
seqtab.nochim.join_by_ASVseq.Nf = data.frame(asv.seq = rownames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)), t(seqtab.nochim.sample_rare_by_run_min.rmNA)) %>%
filter(asv.seq %in% asv_tax.files_cat.join_by_ASVseq.Nf$asv.seq)
seqtab.nochim.join_by_ASVseq.Nf[is.na(seqtab.nochim.join_by_ASVseq.Nf)] = 0
seqtab.nochim.join_by_ASVseq.Nf$asv.seq = NULL

Nf_count = data.frame(
Nf = colSums(seqtab.nochim.join_by_ASVseq.Nf),
sample = colnames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)),
richness = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA) > 0),
count = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA))
)

Nf_count.meta = left_join(
Nf_count,
full_metadata.files_cat %>%
filter(locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "Pool")) ,
by = "sample"
)
Nf_count.meta[is.na(Nf_count.meta)] = 0

Nf_count.meta.wide = pivot_wider(Nf_count.meta, id_cols = metadata.label, names_from = Pool, values_from = c(count, richness, Nf))

qqnorm(residuals(lm(log10(Nf_pool_T+1) ~ log10(Nf_pool_F+1), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
plot(residuals(lm(log10(Nf_pool_T+1) ~ log10(Nf_pool_F+1), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
summary(lm(log10(Nf_pool_T+1) ~ log10(Nf_pool_F+1), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

qqnorm(residuals(lm(Nf_pool_T/count_pool_T ~ I(Nf_pool_F/count_pool_F), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100))))
summary(lm(Nf_pool_T/count_pool_T ~ I(Nf_pool_F/count_pool_F), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

#proportional counts
Nf.binomial_glm = glm(Nf_pool_T/count_pool_T ~ I(Nf_pool_F/count_pool_F*100), data = Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), family = binomial)
1-(Nf.binomial_glm$deviance/Nf.binomial_glm$null.deviance
)
1 - pchisq(summary(Nf.binomial_glm)$deviance,
summary(Nf.binomial_glm)$df.residual
)

p1 = ggplot(Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nf_pool_F/count_pool_F, Nf_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
my_gg_theme +
labs(x = "pool_F N. faginata rel. abd.", y = "pool_T N. faginata rel. abd.",
title = paste(
"N. faginata rel abd by sample\nR2 = ", signif(1-(Nf.binomial_glm$deviance/Nf.binomial_glm$null.deviance), 2),
" P =",signif(summary(Nf.binomial_glm)$coef[2,4], 3)
)
)

p2 = ggplot(Nf_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nf_pool_F+1, Nf_pool_T+1)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "pool_F N. faginata count log10+1", y = "pool_T N. faginata count log10+1", title = "N. faginata counts") +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_Nf_RA_rarefied.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#############
#Nd rel.abd.#
seqtab.nochim.join_by_ASVseq.Nd = data.frame(asv.seq = rownames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)), t(seqtab.nochim.sample_rare_by_run_min.rmNA)) %>%
filter(asv.seq %in% asv_tax.files_cat.join_by_ASVseq.Nd$asv.seq)
seqtab.nochim.join_by_ASVseq.Nd[is.na(seqtab.nochim.join_by_ASVseq.Nd)] = 0
seqtab.nochim.join_by_ASVseq.Nd$asv.seq = NULL

Nd_count = data.frame(
Nd = colSums(seqtab.nochim.join_by_ASVseq.Nd),
sample = colnames(t(seqtab.nochim.sample_rare_by_run_min.rmNA)),
richness = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA) > 0),
count = colSums(t(seqtab.nochim.sample_rare_by_run_min.rmNA))
)

Nd_count.meta = left_join(
Nd_count,
full_metadata.files_cat %>%
filter(locus == "ITS2" & !is.na(metadata.label)) %>%
dplyr::select(c("metadata.label", "sample", "Pool")) ,
by = "sample"
)
Nd_count.meta[is.na(Nd_count.meta)] = 0

Nd_count.meta.wide = pivot_wider(Nd_count.meta, id_cols = metadata.label, names_from = Pool, values_from = c(count, richness, Nd))

qqnorm(residuals(lm(log10(Nd_pool_T+1) ~ log10(Nd_pool_F+1), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100 & Nd_pool_F+Nd_pool_T != 0))))
plot(residuals(lm(log10(Nd_pool_T+1) ~ log10(Nd_pool_F+1), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100 & Nd_pool_F+Nd_pool_T != 0))))
summary(lm(log10(Nd_pool_T+1) ~ log10(Nd_pool_F+1), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100 & Nd_pool_F+Nd_pool_T != 0)))
summary(lm(Nd_pool_T/count_pool_T ~ I(Nd_pool_F/count_pool_F), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100)))

#proportional counts
Nd.binomial_glm = glm(Nd_pool_T/count_pool_T ~ I(Nd_pool_F/count_pool_F), data = Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), family = binomial)
1-(Nd.binomial_glm$deviance/Nd.binomial_glm$null.deviance
)
1 - pchisq(summary(Nd.binomial_glm)$deviance,
summary(Nd.binomial_glm)$df.residual
)

p1 = ggplot(Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nd_pool_F/count_pool_F, Nd_pool_T/count_pool_T)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
#scale_x_sqrt(labels = fancy_scientific) +
#scale_y_sqrt(labels = fancy_scientific) +
my_gg_theme +
labs(x = "pool_F N. ditissima rel. abd.", y = "pool_T N. ditissima rel. abd.",
title = paste(
"N. ditissima rel abd by sample\nR2 = ", signif(1-(Nd.binomial_glm$deviance/Nd.binomial_glm$null.deviance), 2),
" P =",signif(summary(Nd.binomial_glm)$coef[2,4], 1)
)
)

p2 = ggplot(Nd_count.meta.wide %>% filter(count_pool_F >= 100 & count_pool_T >= 100), aes(Nd_pool_F+1, Nd_pool_T+1)) +
geom_abline(slope = 1, intercept = 0) +
geom_smooth(method = "lm") +
geom_point(alpha = 0.5) +
scale_x_log10(labels = fancy_scientific) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "pool_F N. ditissima count log10+1", y = "pool_T N. ditissima count log10+1", title = "N. ditissima counts") +
my_gg_theme

pdf("pool_v_unpool_figs/files_cat_Nd_RA_rarefied.pdf", width = 16, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()





###############################################
#Calculate Nf and Nd detection agreement between pooling#

#Nf

Nf_pooling_agree = data.frame(
pooling = c("pool T and F \nNf present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNf absent"),
count = vector(length = 4, mode = "numeric")
)
Nf_pooling_agree$pooling = factor(Nf_pooling_agree$pooling, levels = c("pool T and F \nNf present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNf absent"))

agree_or_not = function(x){
    temp.mat = c(0,0,0,0)
    if(x[1] > 0 & x[2] > 0){
        temp.mat[1] = temp.mat[1]+1
    }else if(x[1] == 0 & x[2] == 0){
        temp.mat[4] = temp.mat[4]+1
    }else if(x[1] > 0 & x[2] == 0){
        temp.mat[2] = temp.mat[2]+1
    }else if(x[1] == 0 & x[2] > 0){
        temp.mat[3] = temp.mat[3]+1
    }
    return(temp.mat)
}

Nf_pooling_agree[,2] = rowSums(apply(data.frame(Nf_count.meta.wide$Nf_pool_T, Nf_count.meta.wide$Nf_pool_F), 1, agree_or_not) )

p1 = ggplot(Nf_pooling_agree, aes(x = pooling, y = count)) +
geom_col() +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1)
) +
labs(x = "", y = "sample count")

pdf("pool_v_unpool_figs/Nf_detection_by_pooling_files_cat.pdf")
p1
dev.off()

#Nd

Nd_pooling_agree = data.frame(
pooling = c("pool T and F \nNd present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNd absent"),
count = vector(length = 4, mode = "numeric")
)
Nd_pooling_agree$pooling = factor(Nd_pooling_agree$pooling, levels = c("pool T and F \nNd present", "pool T present\npool F absent", "pool T absent\npool F present", "pool T and F\nNd absent"))

agree_or_not = function(x){
    temp.mat = c(0,0,0,0)
    if(x[1] > 0 & x[2] > 0){
        temp.mat[1] = temp.mat[1]+1
    }else if(x[1] == 0 & x[2] == 0){
        temp.mat[4] = temp.mat[4]+1
    }else if(x[1] > 0 & x[2] == 0){
        temp.mat[2] = temp.mat[2]+1
    }else if(x[1] == 0 & x[2] > 0){
        temp.mat[3] = temp.mat[3]+1
    }
    return(temp.mat)
}

Nd_pooling_agree[,2] = rowSums(apply(data.frame(Nd_count.meta.wide$Nd_pool_T, Nd_count.meta.wide$Nd_pool_F), 1, agree_or_not) )

p1 = ggplot(Nd_pooling_agree, aes(x = pooling, y = count)) +
geom_col() +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1)
) +
labs(x = "", y = "sample count")

pdf("pool_v_unpool_figs/Nd_detection_by_pooling_files_cat.pdf")
p1
dev.off()


