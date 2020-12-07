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

#vars to use from full_metadata.sorted
#duration_infection NeoFruiting RaisedCanker  Wax Xylococcus TreeCond quartile occurence

neo_fruit = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$NeoFruiting, duleg = F)
summary(neo_fruit)
summary(neo_fruit, indvalcom = T)


#In order to set restcomb to an ecologically meaningful grouping can first extract the different combinations of levels from $comb of the multipatt object
#The indices of this vector will then be the indices of the groups to look at (pick as desired)
#With ordinal variable sensible combinations include consecutive groupings
#0, 1, 2, 3, 4, 5, 0+1, 0+1+2, 2+3+4+5, 3+4+5, 4+5
possible_groupings = colnames(neo_fruit$comb)

restcomb = c(1,2,3,4,5,6,7, 22, 56, 41, 21)

neo_fruit = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$NeoFruiting, duleg = F, restcomb = restcomb)
summary(neo_fruit)


RaisedCanker = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$RaisedCanker, duleg = F)
possible_groupings = colnames(RaisedCanker$comb)
#0, 1, 2, 3,0+1, 1+2+3,2+3
restcomb = c(1,2,3,4,5, 14,10)

RaisedCanker = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$RaisedCanker, duleg = F, restcomb = restcomb)
summary(RaisedCanker)

Wax = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Wax, duleg = F)
possible_groupings = colnames(Wax$comb)
#0, 1, 2, 3, 4, 5, 0+1, 0+1+2, 2+3+4+5, 3+4+5, 4+5
restcomb = c(1,2,3,4,5,6,7, 22, 56, 41, 21)
Wax = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Wax, duleg = F, restcomb = restcomb)
summary(Wax)

Xylococcus = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Xylococcus, duleg = F)
possible_groupings = colnames(Xylococcus$comb)
#0,1,2,0+1,1+2
restcomb = c(1,2,3,4,6)
Xylococcus = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Xylococcus, duleg = F, restcomb = restcomb)
summary(Xylococcus)

TreeCond = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$TreeCond, duleg = F)
possible_groupings = colnames(TreeCond$comb)
#0,1,2,3, 0+1, 1+2+3, 2+3
restcomb = c(1,2,3,4,5,14, 10)
TreeCond = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$TreeCond, duleg = F, restcomb = restcomb)
summary(TreeCond)

quartile = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$quartile, duleg = F)
possible_groupings = colnames(quartile$comb)
#lower, inter, upper, inter+lower, inter+upper
restcomb = c(1,2,3,4,5)
quartile = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$quartile, duleg = F, restcomb = restcomb)
summary(quartile)

occurrence = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$occurence, duleg = T)
possible_groupings = colnames(occurrence$comb)
#all
summary(occurrence)

Nf = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Nf, duleg = T)
summary(Nf)

Nd = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Nd, duleg = T)
summary(Nd)

Site = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$Site, duleg = F, max.order = 3)
summary(Site)

asv_tax[rownames(asv_tax) == "ASV_29", ]



######################################################
#Write sig tables to parse sig spp associations for multiple category ASVs

write.table(neo_fruit$sign, "ISA_spp_tables/neo_fruit.txt", quote = F, sep = "\t", row.names = T)
write.table(RaisedCanker$sign, "ISA_spp_tables/RaisedCanker.txt", quote = F, sep = "\t", row.names = T)
write.table(Wax$sign, "ISA_spp_tables/Wax.txt", quote = F, sep = "\t", row.names = T)
write.table(TreeCond$sign, "ISA_spp_tables/TreeCond.txt", quote = F, sep = "\t", row.names = T)
write.table(occurrence$sign, "ISA_spp_tables/occurrence.txt", quote = F, sep = "\t", row.names = T)

