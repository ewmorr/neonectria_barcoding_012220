require(tidyverse)
require(vegan)
require(indicspecies)

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
asv_tab.gt1K.rare.bin = asv_tab.gt1K.rare
asv_tab.gt1K.rare.bin[asv_tab.gt1K.rare.bin > 0] = 1

full_metadata.sorted = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare)),
full_metadata
) %>% left_join(., Nf_v_Nd.bin)

#for ISA only consider single factor level (not combinations) for disease severity vars
duleg = F

levels_of = c(0, 1, 2, 3, 4, 5)
for(i in levels_of){
    
}

restcomb = c(1,2,3,4,5,6,7,)

#vars to use from full_metadata.sorted
#duration_infection NeoFruiting RaisedCanker  Wax Xylococcus TreeCond quartile occurence

neo_fruit = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$NeoFruiting, duleg = F)
summary(neo_fruit)
summary(neo_fruit, indvalcom = T)


#In order to set restcomb to an ecologically meaningful grouping can first extract the different combinations of levels from $comb of the multipatt object
#The indies of this vector will then be the indices of the groups to look at (pick as desired)
possible_groupings = colnames(neo_fruit$comb)

RaisedCanker = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$RaisedCanker, duleg = F)
summary(RaisedCanker)

Wax = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Wax, duleg = F)
summary(Wax)

Xylococcus = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Xylococcus, duleg = F)
summary(Xylococcus)

TreeCond = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$TreeCond, duleg = F)
summary(TreeCond)

TreeCond = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$TreeCond, duleg = F)
summary(TreeCond)

quartile = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$quartile, duleg = F)
summary(quartile)

occurrence = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$occurence, duleg = T)
summary(occurrence)

Nf = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Nf, duleg = T)
summary(Nf)

Nd = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Nd, duleg = T)
summary(Nd)

Site = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Site, duleg = F, max.order = 3)
summary(Site)

asv_tax[rownames(asv_tax) == "ASV_29", ]

















