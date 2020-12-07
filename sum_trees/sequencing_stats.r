require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)

source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
#this pulls in objects:
#asv_tab
#asv_tax
#id_bench_map

#joins metadata files to get metadata object:
#full_metadata

#creates negatives and controls only asv_tab (long format):
#asv_tab.negatives.long

#creates new object with lowest informative taxonomic level and a character asv_tax with unknown instead of NA
#asv_informative_taxa
#asv_tax.char

#and calc duration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

#######
#PLOTS#
######
#NMDS#
asv_tab.gt1k = asv_tab[,colSums(asv_tab) > 999]
sum(rowSums(asv_tab.gt1k) > 0)
colSums(asv_tab.gt1k > 0) %>% sort
colSums(asv_tab.gt1k > 0) %>% median
##########################
#1000 seqs per sample min#
##########################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
asv_tab.gt1K.rare = asv_tab.gt1K.rare[,colSums(asv_tab.gt1K.rare) > 0]

richness_rarefied = rowSums(asv_tab.gt1K.rare > 0)
richness_rarefied = data.frame("sample" = names(richness_rarefied), richness_rarefied)

full_meatadata.sorted = left_join(
    data.frame(sample = rownames(asv_tab.gt1K.rare)),
    full_metadata
) %>% left_join(., Nf_v_Nd.bin)

full_meatadata.sorted = left_join(full_meatadata.sorted, richness_rarefied)
summarized_richness = full_meatadata.sorted %>% group_by(Site) %>% summarize(mean = mean(richness_rarefied), sd = sd(richness_rarefied))
full_meatadata.sorted$summarized_richness %>% median


sample_rarefied_richness = rarefy(t(asv_tab[,colSums(asv_tab) > 1000]), sample = 1000)
sample_rarefied_richness.metadata = left_join(
data.frame(sample = names(sample_rarefied_richness), richness = sample_rarefied_richness),
full_metadata,
by = "sample"
)

