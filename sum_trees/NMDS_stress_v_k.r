require(vegan)
require(ggplot2)
source("~/ggplot_theme.txt")

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

stress_v_k = data.frame(
    k = vector(mode="numeric", length = 6),
    stress = vector(mode="numeric", length = 6),
    run_count =  vector(mode="numeric", length = 6)
)

for(i in 2:7){
    print(paste("running k =", i))
    run_count = 1
    asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = i)
    #No convergence reached with k = 2, solution reached with k = 3
    #also tried binary without filtering to 5k seqs per sample, no convergence
    while(asv_tab.gt1K.rare.mds$converged != TRUE & run_count <= 25){
        run_count = run_count+1
        asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = i, previous.best = asv_tab.gt1K.rare.mds)
    }
    if(asv_tab.gt1K.rare.mds$converged != TRUE){
        stress_v_k$stress[i-1] = 0
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
    else{
        stress_v_k$stress[i-1] = asv_tab.gt1K.rare.mds$stress
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
}

p1 = ggplot(stress_v_k, aes(k, stress)) +
geom_point() +
labs(title = "Abundance weighted (log10+1 sequence count)\nStress = 0 -- no convergence after 26 runs") +
my_gg_theme +
theme(plot.title = element_text(size = 15))



stress_v_k = data.frame(
k = vector(mode="numeric", length = 6),
stress = vector(mode="numeric", length = 6),
run_count =  vector(mode="numeric", length = 6)
)

for(i in 2:7){
    print(paste("running k =", i))
    run_count = 1
    asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = i)
    #No convergence reached with k = 2, solution reached with k = 3
    #also tried binary without filtering to 5k seqs per sample, no convergence
    while(asv_tab.gt1K.rare.mds.bin$converged != TRUE & run_count <= 25){
        run_count = run_count+1
        asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = i, previous.best = asv_tab.gt1K.rare.mds.bin)
    }
    if(asv_tab.gt1K.rare.mds.bin$converged != TRUE){
        stress_v_k$stress[i-1] = 0
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
    else{
        stress_v_k$stress[i-1] = asv_tab.gt1K.rare.mds.bin$stress
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
}

require(ggplot2)

p2 = ggplot(stress_v_k, aes(k, stress)) +
geom_point() +
labs(title = "Presence-absence\nStress = 0 -- no convergence after 26 runs") +
my_gg_theme +
theme(plot.title = element_text(size = 15))


require(gridExtra)
pdf("NMDS_fits/k_v_stress_1Kseqs_noSingletons.sum_trees.pdf", width = 12, height = 4)
grid.arrange(p1,p2,ncol = 2)
dev.off()


##############################################################
#filter out doubletons

asv_tab = asv_tab[rowSums(asv_tab > 0) > 2,]

asv_tab.gt1K = asv_tab[colSums(asv_tab) >= 1000]
asv_tab.gt1K = asv_tab.gt1K[rowSums(asv_tab.gt1K) >=1, ]
asv_tab.gt1K.rare = rrarefy(t(asv_tab.gt1K), 1000)

stress_v_k = data.frame(
k = vector(mode="numeric", length = 6),
stress = vector(mode="numeric", length = 6),
run_count =  vector(mode="numeric", length = 6)
)

for(i in 2:7){
    print(paste("running k =", i))
    run_count = 1
    asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = i)
    #No convergence reached with k = 2, solution reached with k = 3
    #also tried binary without filtering to 5k seqs per sample, no convergence
    while(asv_tab.gt1K.rare.mds$converged != TRUE & run_count <= 25){
        run_count = run_count+1
        asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = i, previous.best = asv_tab.gt1K.rare.mds)
    }
    if(asv_tab.gt1K.rare.mds$converged != TRUE){
        stress_v_k$stress[i-1] = 0
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
    else{
        stress_v_k$stress[i-1] = asv_tab.gt1K.rare.mds$stress
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
}


p1 = ggplot(stress_v_k, aes(k, stress)) +
geom_point() +
labs(title = "Abundance weighted (log10+1 sequence count)\nStress = 0 -- no convergence after 26 runs") +
my_gg_theme +
theme(plot.title = element_text(size = 15))



stress_v_k = data.frame(
k = vector(mode="numeric", length = 6),
stress = vector(mode="numeric", length = 6),
run_count =  vector(mode="numeric", length = 6)
)

for(i in 2:7){
    print(paste("running k =", i))
    run_count = 1
    asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = i)
    #No convergence reached with k = 2, solution reached with k = 3
    #also tried binary without filtering to 5k seqs per sample, no convergence
    while(asv_tab.gt1K.rare.mds.bin$converged != TRUE & run_count <= 25){
        run_count = run_count+1
        asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = i, previous.best = asv_tab.gt1K.rare.mds.bin)
    }
    if(asv_tab.gt1K.rare.mds.bin$converged != TRUE){
        stress_v_k$stress[i-1] = 0
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
    else{
        stress_v_k$stress[i-1] = asv_tab.gt1K.rare.mds.bin$stress
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
}

p2 = ggplot(stress_v_k, aes(k, stress)) +
geom_point() +
labs(title = "Presence-absence\nStress = 0 -- no convergence after 26 runs") +
my_gg_theme +
theme(plot.title = element_text(size = 15))

pdf("NMDS_fits/k_v_stress_1Kseqs_noDoubletons.sum_trees.pdf", width = 12, height = 4)
grid.arrange(p1,p2,ncol = 2)
dev.off()


##############################################################
#filter out tripletons

asv_tab = asv_tab[rowSums(asv_tab > 0) > 3,]

asv_tab.gt1K = asv_tab[colSums(asv_tab) >= 1000]
asv_tab.gt1K = asv_tab.gt1K[rowSums(asv_tab.gt1K) >=1, ]
asv_tab.gt1K.rare = rrarefy(t(asv_tab.gt1K), 1000)

stress_v_k = data.frame(
k = vector(mode="numeric", length = 6),
stress = vector(mode="numeric", length = 6),
run_count =  vector(mode="numeric", length = 6)
)

for(i in 2:7){
    print(paste("running k =", i))
    run_count = 1
    asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = i)
    #No convergence reached with k = 2, solution reached with k = 3
    #also tried binary without filtering to 5k seqs per sample, no convergence
    while(asv_tab.gt1K.rare.mds$converged != TRUE & run_count <= 25){
        run_count = run_count+1
        asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = i, previous.best = asv_tab.gt1K.rare.mds)
    }
    if(asv_tab.gt1K.rare.mds$converged != TRUE){
        stress_v_k$stress[i-1] = 0
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
    else{
        stress_v_k$stress[i-1] = asv_tab.gt1K.rare.mds$stress
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
}


p1 = ggplot(stress_v_k, aes(k, stress)) +
geom_point() +
labs(title = "Abundance weighted (log10+1 sequence count)\nStress = 0 -- no convergence after 26 runs") +
my_gg_theme +
theme(plot.title = element_text(size = 15))



stress_v_k = data.frame(
k = vector(mode="numeric", length = 6),
stress = vector(mode="numeric", length = 6),
run_count =  vector(mode="numeric", length = 6)
)

for(i in 2:7){
    print(paste("running k =", i))
    run_count = 1
    asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = i)
    #No convergence reached with k = 2, solution reached with k = 3
    #also tried binary without filtering to 5k seqs per sample, no convergence
    while(asv_tab.gt1K.rare.mds.bin$converged != TRUE & run_count <= 25){
        run_count = run_count+1
        asv_tab.gt1K.rare.mds.bin = metaMDS(vegdist(asv_tab.gt1K.rare, method = "bray", binary = T), autotransform = F, k = i, previous.best = asv_tab.gt1K.rare.mds.bin)
    }
    if(asv_tab.gt1K.rare.mds.bin$converged != TRUE){
        stress_v_k$stress[i-1] = 0
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
    else{
        stress_v_k$stress[i-1] = asv_tab.gt1K.rare.mds.bin$stress
        stress_v_k$k[i-1] = i
        stress_v_k$run_count[i-1] = run_count
    }
}

p2 = ggplot(stress_v_k, aes(k, stress)) +
geom_point() +
labs(title = "Presence-absence\nStress = 0 -- no convergence after 26 runs") +
my_gg_theme +
theme(plot.title = element_text(size = 15))

pdf("NMDS_fits/k_v_stress_1Kseqs_noTripletons.sum_trees.pdf", width = 12, height = 4)
grid.arrange(p1,p2,ncol = 2)
dev.off()
