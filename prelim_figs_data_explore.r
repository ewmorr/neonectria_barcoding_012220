require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)

my_gg_theme = theme(
    panel.background = element_rect(fill='white', colour='black'),
    panel.grid.major=element_blank(),
    panel.grid.minor= element_blank(),
    text=element_text(family="sans"),
    axis.text=element_text(size=15, color="black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust=0, size=20),
    axis.title = element_text(size=17),
    legend.title = element_blank(),
    legend.text = element_text(size = 19),
    strip.text = element_text(size = 15),
    axis.title.x = element_text(margin = margin(t= 10)),
    axis.title.y = element_text(margin = margin(r=10))
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#read data
source("~/repo/neonectria_barcoding_012220/read_ASV_dat.r")
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

#and calculates neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

system("mkdir prelim_figs")

#######
#PLOTS#

#negatives/control by type with taxonomy
asv.negatives.plot_dat = left_join(asv_tab.negatives.long, full_metadata, by = "sample") %>%
    left_join(.,
        data.frame(ASV = names(asv_informative_taxa),
            taxon = paste(names(asv_informative_taxa), asv_informative_taxa, sep = " ")
        ),
    by = "ASV"
    )

#Get color palette
palette_len = asv.negatives.plot_dat %>% filter(count > 0) %>% select(taxon) %>% unique %>% rownames %>% length
getPalette = colorRampPalette(brewer.pal(n = 12, name = "Paired"))

p = ggplot(asv.negatives.plot_dat, aes(x = sample, y = count, fill = taxon)) +
geom_col(color = "black") +
facet_wrap(~bench.control, scales = "free") +
#scale_fill_brewer(palette = "Paired") +
scale_fill_manual(values = getPalette(palette_len)) +
my_gg_theme +
theme(
legend.text = element_text(size = 14),
axis.text.x = element_text(angle = 55, hjust = 1)
)

pdf("prelim_figs/controls_stacked_taxa.pdf", width = 14, height = 8)
print(p)
dev.off()

#Neo occurence across sites, trees, plugs
pdf("prelim_figs/Neo_occurence_min_100_seqs_per_sample.pdf", width = 10)
ggplot(Nf_v_Nd.bin.metadata %>% filter(total_seqs > 100 & bench.control == "n"), aes(as.factor(Tree), ..count.., fill = occurence, color = as.factor(Plug))) +
geom_histogram(stat = "count", width = 0.9) +
facet_wrap(lat~State, scales = "free_x", ncol = 5) +
scale_fill_manual(values = rev(cbPalette[1:4])) +
scale_color_manual(values = rep("black", length((Nf_v_Nd.bin.metadata %>% filter(total_seqs > 100 & bench.control == "n"))$Plug) ), guide = F) +
labs(x = "Trees", y = "count (plugs)", fill = "Neonectria\noccurence") +
my_gg_theme +
theme(axis.text.x = element_blank(), legend.title = element_text(size = 18))
dev.off()

#100 seq minimum NF/ND frequency by site
sites = Nf_v_Nd.bin.metadata %>% filter(bench.control == "n") %>% select(Site) %>% unique

Nf_Nd_site_freq = data.frame(Site = vector(mode = "character", length = length(sites$Site)),
Nf = vector(mode = "numeric", length = length(sites$Site)),
Nd = vector(mode = "numeric", length = length(sites$Site)),
total = vector(mode = "numeric", length = length(sites$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(sites$Site)){
    temp_tab = Nf_v_Nd.bin.metadata %>% filter(Site == sites$Site[i] & total_seqs >= 100)
    Nf_Nd_site_freq$Nf[i] = (filter(temp_tab, Nf > 0)) %>% nrow
    Nf_Nd_site_freq$Nd[i] = (filter(temp_tab, Nd > 0)) %>% nrow
    Nf_Nd_site_freq$total[i] = nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(sites$Site[i])
}

Nf_Nd_site_freq.site_info = left_join(Nf_Nd_site_freq, site_info)

p1 = ggplot(Nf_Nd_site_freq.site_info, aes(lat, Nd/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "Latitude") +
scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nf/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Latitude") +
scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_lat_100_min.pdf", width = 10, height = 4)
grid.arrange(p1,p2,ncol= 2, widths = c(0.45,0.55))
dev.off()

######
#NMDS#


##########################
#1000 seqs per sample min#
asv_tab.gt1K = asv_tab[colSums(asv_tab) >= 1000]
asv_tab.gt1K = asv_tab.gt1K[rowSums(asv_tab.gt1K) >=1, ]
asv_tab.gt1K.rare = rrarefy(t(asv_tab.gt1K), 1000)

asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = 3)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
while(asv_tab.gt1K.rare.mds$converged != TRUE){
    asv_tab.gt1K.rare.mds = metaMDS(log10(asv_tab.gt1K.rare+1), autotransform = F, k = 3, previous.best = asv_tab.gt1K.rare.mds)
}
#Add metadata

asv_tab.gt1K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare.mds$points), asv_tab.gt1K.rare.mds$points),
full_metadata, by = "sample")

asv_tab.gt1K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt1K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

#PLOT

p = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = lat, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "Latitude")

pdf("prelim_figs/NMDS_neo_occurence_by_lat_1Kmin.pdf", width = 8, height = 6)
p
dev.off()

##############################
#Neonectria frequency by site#

sites_nmds = asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n") %>% select(Site) %>% unique

Nf_Nd_site_freq = data.frame(Site = vector(mode = "character", length = length(sites_nmds$Site)),
Nf = vector(mode = "numeric", length = length(sites_nmds$Site)),
Nd = vector(mode = "numeric", length = length(sites_nmds$Site)),
total = vector(mode = "numeric", length = length(sites_nmds$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(sites_nmds$Site)){
    temp_tab = asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(Site == sites_nmds$Site[i])
    Nf_Nd_site_freq$Nf[i] = (filter(temp_tab, Nf > 0)) %>% nrow
    Nf_Nd_site_freq$Nd[i] = (filter(temp_tab, Nd > 0)) %>% nrow
    Nf_Nd_site_freq$total[i] = nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(sites_nmds$Site[i])
}

#This line probably not necessary because already added full metadata to NMDS, but haven't tested in this script
Nf_Nd_site_freq.site_info = left_join(Nf_Nd_site_freq, site_info)

p1 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nd/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "Latitude") +
scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nf/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Latitude") +
scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_lat_1K_min.pdf", width = 10, height = 4)
grid.arrange(p1,p2,ncol= 2, widths = c(0.45,0.55))
dev.off()


#####################
#Richness by lat etc#

sample_richness = data.frame(sample = names(apply(asv_tab.gt1K.rare,1,function(x) sum(x > 0))), richness =  apply(asv_tab.gt1K.rare,1,function(x) sum(x > 0)))
sample_richness.metadata = left_join(sample_richness, id_bench_map) %>%
left_join(., metadata_map)


sample_richness.metadata.site_info = left_join(sample_richness.metadata, Nf_v_Nd.bin) %>% left_join(., site_info, by = "Site")

pdf("prelim_figs/ASV_richness_by_lat_1K_min.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(lat, richness, group = lat, color = state_prov)) +
geom_boxplot(width = .25) +
labs(y = "ASV richness\nper 1K sequences", x = "Latitude") +
scale_color_brewer(palette = "Dark2") +
my_gg_theme
dev.off()

pdf("prelim_figs/ASV_richness_by_neo_occurence_1K_min.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(occurence, richness)) +
geom_boxplot(width = .25) +
labs(y = "ASV richness\nper 1K sequences", x = "Neonectria occurence") +
scale_color_manual(values = cbPalette, guide = F) +
scale_x_discrete(labels = c("both", "N. ditissima", "N. faginata", "none")) +
my_gg_theme
dev.off()

##########################
#5000 seqs per sample min#
asv_tab.gt5K = asv_tab[colSums(asv_tab) >= 5000]
asv_tab.gt5K = asv_tab.gt5K[rowSums(asv_tab.gt5K) >=1, ]
asv_tab.gt5K.rare = rrarefy(t(asv_tab.gt5K), 5000)

asv_tab.gt5K.rare.mds = metaMDS(log10(asv_tab.gt5K.rare+1), autotransform = F, k = 3)
#No convergence reached with k = 2, solution reached with k = 3
#also tried binary without filtering to 5k seqs per sample, no convergence
while(asv_tab.gt5K.rare.mds$converged != TRUE){
    asv_tab.gt5K.rare.mds = metaMDS(log10(asv_tab.gt5K.rare+1), autotransform = F, k = 3, previous.best = asv_tab.gt5K.rare.mds)
}
#Add metadata

asv_tab.gt5K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt5K.rare.mds$points), asv_tab.gt5K.rare.mds$points),
full_metadata, by = "sample")

asv_tab.gt5K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt5K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

#PLOT

p = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = lat, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "Latitude")

pdf("prelim_figs/NMDS_neo_occurence_by_lat_5Kmin.pdf", width = 8, height = 6)
p
dev.off()


##############################
#Neonectria frequency by site#

sites_nmds = asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n") %>% select(Site) %>% unique

Nf_Nd_site_freq = data.frame(Site = vector(mode = "character", length = length(sites_nmds$Site)),
Nf = vector(mode = "numeric", length = length(sites_nmds$Site)),
Nd = vector(mode = "numeric", length = length(sites_nmds$Site)),
total = vector(mode = "numeric", length = length(sites_nmds$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(sites_nmds$Site)){
    temp_tab = asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(Site == sites_nmds$Site[i])
    Nf_Nd_site_freq$Nf[i] = (filter(temp_tab, Nf > 0)) %>% nrow
    Nf_Nd_site_freq$Nd[i] = (filter(temp_tab, Nd > 0)) %>% nrow
    Nf_Nd_site_freq$total[i] = nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(sites_nmds$Site[i])
}

Nf_Nd_site_freq.site_info = left_join(Nf_Nd_site_freq, site_info)

p1 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nd/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "Latitude") +
scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nf/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Latitude") +
scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_lat_5K_min.pdf", width = 10, height = 4)
grid.arrange(p1,p2,ncol= 2, widths = c(0.45,0.55))
dev.off()


#####################
#Richness by lat etc#

sample_richness = data.frame(sample = names(apply(asv_tab.gt5K.rare,1,function(x) sum(x > 0))), richness =  apply(asv_tab.gt5K.rare,1,function(x) sum(x > 0)))
sample_richness.metadata = left_join(sample_richness, id_bench_map) %>%
left_join(., metadata_map)


sample_richness.metadata.site_info = left_join(sample_richness.metadata, Nf_v_Nd.bin) %>% left_join(., site_info, by = "Site")

pdf("prelim_figs/ASV_richness_by_lat_5K_min.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(lat, richness, group = lat, color = state_prov)) +
geom_boxplot(width = .25) +
labs(y = "ASV richness\nper 5K sequences", x = "Latitude") +
scale_color_brewer(palette = "Dark2") +
my_gg_theme
dev.off()

pdf("prelim_figs/ASV_richness_by_neo_occurence_5K_min.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(occurence, richness)) +
geom_boxplot(width = .25) +
labs(y = "ASV richness\nper 5K sequences", x = "Neonectria occurence") +
scale_color_manual(values = cbPalette) +
scale_x_discrete(labels = c("both", "N. ditissima", "N. faginata", "none")) +
my_gg_theme
dev.off()


