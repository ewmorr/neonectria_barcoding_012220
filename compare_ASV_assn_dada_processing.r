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


#seq_tab.nochim.files_cat[rowSums(seq_tab.nochim.files_cat) > 1000,] %>% rownames

###############
#Venn Diagram#
require(VennDiagram)
venn.diagram(list(sep_run_pool = colnames(t(seqtab.nochim.pool_sep)), global_pool = colnames(seq_tab.nochim.files_sep), cat_seq_files = colnames(seq_tab.nochim.files_cat)), filename = "compare_dada_processing_figs/shared_ASVs_Venn.tiff")



#Code for ASV frequency plots -- this is done on tables that have been summed at the level of sample pairs so that comparisons are 1:1

seqtab.nochim.pool_sep.sample_sum = readRDS("run1_run2_dada_compare/sep_run_pool/intermediate_RDS/dada2_seq_table_no_chim_run1run2_pool_sep_samples_summed.rds") %>% t

seqtab.nochim.files_sep.sample_sum = readRDS("run1_run2_dada_compare/files_sep/intermediate_RDS/dada2_seq_table_no_chim_global_pool_samples_summed.rds") %>% t

#filter samples based on shared across methods

seqtab.nochim.pool_sep.sample_sum = seqtab.nochim.pool_sep.sample_sum[rownames(seqtab.nochim.pool_sep.sample_sum) %in% rownames(seq_tab.nochim.files_cat),]
seqtab.nochim.files_sep.sample_sum = seqtab.nochim.files_sep.sample_sum[rownames(seqtab.nochim.files_sep.sample_sum) %in% rownames(seq_tab.nochim.files_cat),]

#filter out zero count ASVS
seqtab.nochim.pool_sep.sample_sum = seqtab.nochim.pool_sep.sample_sum[,colSums(seqtab.nochim.pool_sep.sample_sum) > 0]
seqtab.nochim.files_sep.sample_sum = seqtab.nochim.files_sep.sample_sum[,colSums(seqtab.nochim.files_sep.sample_sum) > 0]

#files_cat already loaded

pool_sep.ASV_freq = colSums(seqtab.nochim.pool_sep.sample_sum > 0) %>% unname
files_sep.ASV_freq = colSums(seqtab.nochim.files_sep.sample_sum > 0) %>% unname
files_cat.ASV_freq = colSums(seq_tab.nochim.files_cat > 0) %>% unname

sum(pool_sep.ASV_freq == 1)
sum(files_sep.ASV_freq == 1)
sum(files_cat.ASV_freq == 1)

dada_processing.asv_freq = rbind(
data.frame(processing = "sep run pool", sample_count = pool_sep.ASV_freq),
data.frame(processing = "global pool", sample_count = files_sep.ASV_freq),
data.frame(processing = "cat seq files", sample_count = files_cat.ASV_freq)
)


p = ggplot(dada_processing.asv_freq, aes(sample_count)) +
geom_histogram(stat = "bin", binwidth = 1) +
facet_wrap(~processing) +
scale_y_sqrt(breaks = c(10,100,200,300,400,500)) +
scale_x_continuous(breaks = c(1,50,100,150)) +
labs(x = "ASV frequency (no. samples observed)", y = "ASV count") +
my_gg_theme

pdf("compare_dada_processing_figs/ASV_sample_freq.pdf", width = 12, height = 4)
p
dev.off()


#Rank abundance (Even though this is done below for tables that have not been summed by sample pairs)
pool_sep.rank_abd = colSums(seqtab.nochim.pool_sep.sample_sum) %>% unname %>% sort(decreasing = T)
files_sep.rank_abd = colSums(seqtab.nochim.files_sep.sample_sum) %>% unname %>% sort(decreasing = T)
files_cat.rank_abd = colSums(seq_tab.nochim.files_cat) %>% unname %>% sort(decreasing = T)

dada_processing.rank_abd = rbind(
data.frame(processing = "sep run pool", count = pool_sep.rank_abd, rank = seq(1,length(pool_sep.rank_abd), by = 1)),
data.frame(processing = "global pool", count = files_sep.rank_abd, rank = seq(1,length(files_sep.rank_abd), by = 1)),
data.frame(processing = "cat seq files", count = files_cat.rank_abd, rank = seq(1,length(files_cat.rank_abd), by = 1))
)

p = ggplot(dada_processing.rank_abd, aes(rank, count +1)) +
geom_col(width = 1) +
facet_wrap(~processing) +
my_gg_theme +
labs(y = "sequence count + 1") +
scale_y_log10(labels = fancy_scientific)

pdf("compare_dada_processing_figs/ASV_rank_abundance_sample_pairs_summed.pdf", width = 12, height = 4)
p
dev.off()

#############
#Venn diagram exclude singletons
venn.diagram(list(sep_run_pool = colnames(seqtab.nochim.pool_sep.sample_sum)[colSums(seqtab.nochim.files_sep.sample_sum > 0) > 1], global_pool = colnames(seqtab.nochim.files_sep.sample_sum)[colSums(seqtab.nochim.files_sep.sample_sum > 0) > 1], cat_seq_files = colnames(seq_tab.nochim.files_cat)[colSums(seq_tab.nochim.files_cat > 0) > 1]), filename = "compare_dada_processing_figs/shared_ASVs_Venn_no_singletons.tiff")


########################
#ASV richness by method#

#taking numbers from Venn
asv_richness_method = data.frame(
    method = c(rep("sep run\npool", 4), rep("global\npool", 4), rep("cat seq\nfiles", 4)),
    shared_with = c("all", "global pool", "cat seq file", "unique",
            "all", "sep run pool", "cat seq file", "unique",
            "all", "sep run pool", "global pool", "unique"
    ),
    asv_richness = c(
        1291, 41, 10, 144,
        1291, 41, 28, 17,
        1291, 10, 28, 10
    )
)

asv_richness_method$shared_with = factor(asv_richness_method$shared_with, levels = rev(c("all", "cat seq file", "global pool", "sep run pool", "unique")))

require(RColorBrewer)
require(gridExtra)

cols = list("all" = "light grey", "cat seq file" = brewer.pal(8, "Dark2")[2], "global pool" = brewer.pal(8, "Dark2")[1], "sep run pool" = brewer.pal(8, "Dark2")[6], "unique" = brewer.pal(8, "Dark2")[4])

p1 = ggplot(asv_richness_method, aes(method, asv_richness, fill = shared_with)) +
geom_col() +
scale_fill_manual(values = cols, guide = F) +
my_gg_theme +
labs(x = "", y = "number ASVs") +
theme(
legend.title = element_text(size = 20)
)

p2 = ggplot(asv_richness_method, aes(method, asv_richness, fill = shared_with)) +
geom_col(position = "fill") +
scale_fill_manual(values = cols) +
my_gg_theme +
labs(x = "", y = "number ASVs") +
theme(
    legend.title = element_text(size = 20)
)

pdf("compare_dada_processing_figs/asv_richness_shared_col.pdf", width = 12, height = 4)
grid.arrange(p1,p2, ncol = 2, widths = c(.4,.6))
dev.off()

################################################
#rank abundance curves with indicator of method#
################################################

pool_sep.rank_abd = data.frame(
#    method = "sep run pool" ,
    asv.seq = colnames(t(seqtab.nochim.pool_sep)),
    count.pool_sep = colSums(t(seqtab.nochim.pool_sep)),
    pool_sep = 1
)

global_pool.rank_abd = data.frame(
#    method = "global pool",
    asv.seq = colnames(seq_tab.nochim.files_sep),
    count.global_pool = colSums(seq_tab.nochim.files_sep),
    global_pool = 1
)

cat_seq_files.rank_abd = data.frame(
#    method = "cat seq files",
    asv.seq = colnames(seq_tab.nochim.files_cat),
    count.cat_seq_files = colSums(seq_tab.nochim.files_cat),
    cat_seq_files = 1
)

rank_abd.method_join = full_join(pool_sep.rank_abd, global_pool.rank_abd, by = "asv.seq") %>%
    full_join(., cat_seq_files.rank_abd, by = "asv.seq")

rank_abd.method_join[is.na(rank_abd.method_join)] = 0
rank_abd.method_join$shared_with = vector(mode = "character", length = length(rank_abd.method_join$asv.seq))


rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 3, colnames(rank_abd.method_join) == "shared_with"] = "all"
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 1, colnames(rank_abd.method_join) == "shared_with"] = "unique"



#Rerun these with dif flags before each filter to get the appropriate labeling
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 2 & rank_abd.method_join$global_pool == 1, colnames(rank_abd.method_join) == "shared_with"] = "global pool"
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 2 & rank_abd.method_join$global_pool == 0, colnames(rank_abd.method_join) == "shared_with"] = "sep run pool"

temp.rank_abd.1 = rank_abd.method_join %>% filter(cat_seq_files == 1) %>% arrange(desc(count.cat_seq_files))
temp.rank_abd.1$rank = as.numeric(rownames(temp.rank_abd.1))

#Rerun these with dif flags before each filter to get the appropriate labeling
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 2 & rank_abd.method_join$cat_seq_files == 1, colnames(rank_abd.method_join) == "shared_with"] = "cat seq file"
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 2 & rank_abd.method_join$cat_seq_files == 0, colnames(rank_abd.method_join) == "shared_with"] = "sep run pool"

temp.rank_abd.2 = rank_abd.method_join %>% filter(global_pool == 1) %>% arrange(desc(count.global_pool))
temp.rank_abd.2$rank = as.numeric(rownames(temp.rank_abd.2))

#Rerun these with dif flags before each filter to get the appropriate labeling
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 2 & rank_abd.method_join$cat_seq_files == 1, colnames(rank_abd.method_join) == "shared_with"] = "cat seq file"
rank_abd.method_join[(rank_abd.method_join$pool_sep + rank_abd.method_join$global_pool + rank_abd.method_join$cat_seq_files) == 2 & rank_abd.method_join$cat_seq_files == 0, colnames(rank_abd.method_join) == "shared_with"] = "global pool"

temp.rank_abd.3 = rank_abd.method_join %>% filter(pool_sep == 1) %>% arrange(desc(count.pool_sep))
temp.rank_abd.3$rank = as.numeric(rownames(temp.rank_abd.3))

temp.rank_abd = rbind(
    data.frame(method = "cat seq file", rank = temp.rank_abd.1 %>% select(rank), count = temp.rank_abd.1$count.cat_seq_files, shared_with = temp.rank_abd.1$shared_with),
    data.frame(method = "global pool", rank = temp.rank_abd.2 %>% select(rank), count = temp.rank_abd.2$count.global_pool, shared_with = temp.rank_abd.2$shared_with),
    data.frame(method = "sep run pool", rank = temp.rank_abd.3 %>% select(rank), count = temp.rank_abd.3$count.pool_sep, shared_with = temp.rank_abd.3$shared_with)
)

temp.rank_abd$shared_with = factor(temp.rank_abd$shared_with, levels = c("all", "cat seq file", "global pool", "sep run pool", "unique"))


cols = list("all" = rgb(211/255,211/255,211/255,alpha = 0.5), "cat seq file" = brewer.pal(8, "Dark2")[2], "global pool" = brewer.pal(8, "Dark2")[1], "sep run pool" = brewer.pal(8, "Dark2")[6], "unique" = brewer.pal(8, "Dark2")[4])

p1 = ggplot(temp.rank_abd, aes(rank, count, fill = shared_with)) +
geom_col(width = 1) +
facet_wrap(~method, ncol = 1) +
scale_fill_manual(values = cols) +
scale_y_log10(labels = fancy_scientific) +
labs(y = "ASV sequence count") +
my_gg_theme +
theme(legend.title = element_text(size = 20))


p1 = ggplot() +
geom_col(data = temp.rank_abd %>% filter(shared_with == "all"), aes(rank, count, fill = shared_with), width = 1) +
geom_col(data = temp.rank_abd %>% filter(shared_with != "all"), aes(rank, count, fill = shared_with), width = 1) +
facet_wrap(~method, ncol = 1) +
scale_fill_manual(values = cols) +
scale_y_log10(labels = fancy_scientific) +
scale_x_continuous(expand = c(0,20), breaks = c(1,500,1000,1500)) +
labs(y = "ASV sequence count") +
my_gg_theme +
theme(legend.title = element_text(size = 20))


pdf("compare_dada_processing_figs/rank_abundance.pdf", width = 16, height = 6)
p1
dev.off()







