require(vegan)

#read data
source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")
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

#Sum asv counts at the tree level
source("~/repo/neonectria_barcoding_012220/sum_trees/sum_ASV_counts_by_tree.r")

asv_tab = asv_tab.tree_sum

######
#NMDS#

#asv_tab = asv_tab[rowSums(asv_tab > 0) > 2,]
#asv_tab = asv_tab[rowSums(asv_tab > 0) > 3,]
#asv_tab = asv_tab[rowSums(asv_tab > 0) > n_samples*0.025,]

#858 taxa originally (this is after LULU filtering and removing taxa that occur in only one sample
#736 taxa remaining at filtering based on > 1 frequence
#610 taxa remaining at filtering based on > 3 frequence
#500 taxa remaining at filtering based on > 2.5% frequencey
#262 taxa remaining at filtering based on > 5% frequencey
#122 taxa remaining at filtering based on > 10% frequencey
#See plot in folder rarefaction_ana

#stricter filtering does not seem to improve convergence. Stick with filtering out those spp that only occur in one sample
#may revisit this

##########################
#1000 seqs per sample min#
asv_tab.gt1K = asv_tab[colSums(asv_tab) >= 1000]
asv_tab.gt1K = asv_tab.gt1K[rowSums(asv_tab.gt1K) >=1, ]
asv_tab.gt1K.rare = rrarefy(t(asv_tab.gt1K), 1000)

length(colnames(asv_tab.gt1K.rare))

saveRDS(asv_tab.gt1K.rare, file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")


asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = 3)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
while(asv_tab.gt1K.rare.mds$converged != TRUE){
    asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = 3, previous.best = asv_tab.gt1K.rare.mds)
}
#Add metadata

saveRDS(asv_tab.gt1K.rare.mds, file = "intermediate_RDS/asv_tab.gt1K.rare.mds.tree_sum.rds")

##########################
#5000 seqs per sample min#
asv_tab.gt5K = asv_tab[colSums(asv_tab) >= 5000]
asv_tab.gt5K = asv_tab.gt5K[rowSums(asv_tab.gt5K) >=1, ]
asv_tab.gt5K.rare = rrarefy(t(asv_tab.gt5K), 5000)

saveRDS(asv_tab.gt5K.rare, file = "intermediate_RDS/asv_tab.gt5K.rare.tree_sum.rds")

asv_tab.gt5K.rare.mds = metaMDS(log10(asv_tab.gt5K.rare+1), autotransform = F, k = 3)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
while(asv_tab.gt5K.rare.mds$converged != TRUE){
    asv_tab.gt5K.rare.mds = metaMDS(log10(asv_tab.gt5K.rare+1), autotransform = F, k = 3, previous.best = asv_tab.gt5K.rare.mds)
}

saveRDS(asv_tab.gt5K.rare.mds, file = "intermediate_RDS/asv_tab.gt5K.rare.mds.tree_sum.rds")

########
#binary#

##########################
#1000 seqs per sample min#

asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = 3)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
while(asv_tab.gt1K.rare.mds.bin$converged != TRUE){
    asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = 3, previous.best = asv_tab.gt1K.rare.mds.bin)
}
#Add metadata

saveRDS(asv_tab.gt1K.rare.mds.bin, file = "intermediate_RDS/asv_tab.gt1K.rare.mds.bin.tree_sum.rds")

##########################
#5000 seqs per sample min#

asv_tab.gt5K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt5K.rare, method = "bray", binary = F), autotransform = F, k = 3)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
while(asv_tab.gt5K.rare.mds.bin$converged != TRUE){
    asv_tab.gt5K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt5K.rare, method = "bray", binary = F), autotransform = F, k = 3, previous.best = asv_tab.gt5K.rare.mds.bin)
}

saveRDS(asv_tab.gt5K.rare.mds.bin, file = "intermediate_RDS/asv_tab.gt5K.rare.mds.bin.tree_sum.rds")
