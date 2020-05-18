require(vegan)

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.site_sum.r")
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

#and calcuduration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

min_depth = min(colSums(asv_tab))

#min depth is 168227
######
#NMDS#


##########################
#min seqs per site min#
asv_tab.gtMin = asv_tab[colSums(asv_tab) >= min_depth]
asv_tab.gtMin = asv_tab.gtMin[rowSums(asv_tab.gtMin) >=1, ]
asv_tab.gtMin.rare = rrarefy(t(asv_tab.gtMin), min_depth)

length(colnames(asv_tab.gtMin.rare))

saveRDS(asv_tab.gtMin.rare, file = "intermediate_RDS/asv_tab.gtMin.rare.site_sum.rds")


asv_tab.gtMin.rare.mds = metaMDS(log10(asv_tab.gtMin.rare+1), autotransform = F, k = 2)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
#Add metadata

saveRDS(asv_tab.gtMin.rare.mds, file = "intermediate_RDS/asv_tab.gtMin.rare.mds.site_sum.rds")

##########################
#min seqs per site min#

asv_tab.gtMin.rare.mds.bin = metaMDS(vegdist(asv_tab.gtMin.rare, method = "bray", binary = T), autotransform = F, k = 2)

saveRDS(asv_tab.gtMin.rare.mds.bin, file = "intermediate_RDS/asv_tab.gtMin.rare.mds.bin.site_sum.rds")

