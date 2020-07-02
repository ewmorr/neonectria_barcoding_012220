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

##########################
#1000 seqs per sample min#
##########################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")

full_meatadata.sorted = left_join(
    data.frame(sample = rownames(asv_tab.gt1K.rare)),
    full_metadata
) %>% left_join(., Nf_v_Nd.bin)

#site
adonis(log10(asv_tab.gt1K.rare+1) ~ full_meatadata.sorted$Site)
adonis(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T) ~ full_meatadata.sorted$Site)

#Neo occurence
adonis(log10(asv_tab.gt1K.rare+1) ~ as.factor(Nf)*as.factor(Nd), data = full_meatadata.sorted)
adonis(log10(asv_tab.gt1K.rare+1) ~ as.factor(Nd)*as.factor(Nf), data = full_meatadata.sorted)
adonis(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T) ~ as.factor(Nf)*as.factor(Nd), data = full_meatadata.sorted)
adonis(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T) ~ as.factor(Nd)*as.factor(Nf), data = full_meatadata.sorted)

#both
adonis(log10(asv_tab.gt1K.rare+1) ~ (Site*as.factor(Nf)*as.factor(Nd))^2, data = full_meatadata.sorted)
adonis(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T) ~ (Site*as.factor(Nf)*as.factor(Nd))^2, data = full_meatadata.sorted)

