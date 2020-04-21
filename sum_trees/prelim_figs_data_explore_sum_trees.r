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

#and calcuduration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

#######
#PLOTS#

#Neo occurence across sites, trees, plugs
pdf("prelim_figs/Neo_occurence_min_100_seqs_per_sample_dur_inf_MAT.pdf", width = 16)
ggplot(Nf_v_Nd.bin.metadata %>% filter(total_seqs > 100 & bench.control == "n"), aes(as.factor(Tree), ..count.., fill = occurence)) +
geom_histogram(stat = "count", width = 0.9, color = "black") +
facet_wrap(duration_infection~signif(MAT, 2), scales = "free_x", ncol = 5) +
scale_fill_manual(values = rev(cbPalette[1:4])) +
labs(x = "Trees", y = "count (plugs)") +
my_gg_theme +
theme(axis.text.x = element_blank())
dev.off()


######
#NMDS#

##########################
#1000 seqs per sample min#

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.rds")
#NMDS
asv_tab.gt1K.rare.mds = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.mds.rds")

#Add metadata

asv_tab.gt1K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare.mds$points), asv_tab.gt1K.rare.mds$points),
full_metadata, by = "sample")

asv_tab.gt1K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt1K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

#PLOT
asv_tab.gt1K.rare.mds.metadata.neoOcurence$Site = factor(asv_tab.gt1K.rare.mds.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

#climate/continental scale vars

p1 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = duration_infection, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Infection\nduration") +
labs(title = "Infection duration")

p2 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, shape = occurence, color = Site)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
)+
labs(title = "Site", x = "")

p3 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = tmin, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Min.\ntemp.")+
labs(title = "Average minimum temperature", y = "", x = "")

p4 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = MAT, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "MAT")+
labs(title = "MAT", y = "")

p5 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = tmax, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Max.\ntemp.")+
labs(title = "Average maximum temperature", y = "", x = "")

p6 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = ppt, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Annual\nprecip.")+
labs(title = "Average annual precipitation", y = "")


pdf("NMDS_fits/NMDS_neo_occurence_by_climate_1Kmin.pdf", width = 22, height = 12)
grid.arrange(p2,p3,p5,p1,p4,p6,nrow = 2)
dev.off()

#Disease severity vars

p1 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, shape = occurence, color = Site_mean.dbh)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "DBH")+
labs(title = "DBH (cm)", x = "")

p2 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.NeoFruiting, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Neonectria\nfruting")+
labs(title = "Mean neonectria fruiting (0-5)", y = "", x = "")

p3 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.Wax, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Wax")+
labs(title = "Mean wax load (0-5)", y = "")

p4 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.RaisedCanker, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Cankers")+
labs(title = "Mean cankers", y = "", x = "")

p5 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.TreeCond, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Tree\ncondition")+
labs(title = "Mean tree condition (0-5)", y = "")

p6 = ggplot(asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.Xylococcus, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Xylococcus")+
labs(title = "Mean Xylococcus", y = "")


pdf("NMDS_fits/NMDS_neo_occurence_by_disease_severity_1Kmin.pdf", width = 22, height = 12)
grid.arrange(p1,p2,p3,p4,p5,p6,nrow = 2)
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

Nf_Nd_site_freq.site_info = left_join(Nf_Nd_site_freq, site_info) %>%
left_join(., site_climate) %>%
left_join(., site_means)

p1 = ggplot(Nf_Nd_site_freq.site_info , aes(duration_infection, Nd/total)) +
geom_point(size = 2) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(duration_infection, Nf/total)) +
geom_point(size = 2) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Infection duration (yrs)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p3 = ggplot(Nf_Nd_site_freq.site_info , aes(MAT, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p4 = ggplot(Nf_Nd_site_freq.site_info , aes(MAT, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "MAT") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p5 = ggplot(Nf_Nd_site_freq.site_info , aes(ppt, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p6 = ggplot(Nf_Nd_site_freq.site_info , aes(ppt, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Annual precip (cm)") +
expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p7 = ggplot(Nf_Nd_site_freq.site_info , aes(tmax, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p8 = ggplot(Nf_Nd_site_freq.site_info , aes(tmax, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Average max. temp.") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

p9 = ggplot(Nf_Nd_site_freq.site_info , aes(tmin, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p10 = ggplot(Nf_Nd_site_freq.site_info , aes(tmin, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Average min. temp.") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_duration_infection_and_climate_1K_min.pdf", width = 25, height = 8)
grid.arrange(p1,p7,p3,p9,p5,p2,p8,p4,p10,p6,nrow= 2)
dev.off()

p1 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.dbh, Nd/total)) +
geom_point(size = 2) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.dbh, Nf/total)) +
geom_point(size = 2) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Mean dbh (cm)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p3 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.NeoFruiting, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p4 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.NeoFruiting, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean neonectria fruiting (0-5)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p5 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.Wax, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p6 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.Wax, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean wax load (0-5)") +
#expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p7 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.RaisedCanker, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p8 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.RaisedCanker, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean cankers") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

p9 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.TreeCond, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p10 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.TreeCond, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean tree condition (0-5)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_disease_severity_1K_min.pdf", width = 25, height = 8)
grid.arrange(p1,p3,p5,p7,p9,p2,p4,p6,p8,p10,nrow= 2)
dev.off()

#####################
#Richness by duration_infection etc#

sample_richness = data.frame(sample = names(apply(asv_tab.gt1K.rare,1,function(x) sum(x > 0))), richness =  apply(asv_tab.gt1K.rare,1,function(x) sum(x > 0)))
sample_richness.metadata = left_join(sample_richness, id_bench_map) %>%
left_join(., metadata_map)


sample_richness.metadata.site_info = left_join(sample_richness.metadata, Nf_v_Nd.bin) %>%
    left_join(., site_info, by = "Site") %>%
left_join(., site_climate) %>%
left_join(., site_means, by = "Site")

sample_richness.metadata.site_info$Site = factor(sample_richness.metadata.site_info$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

#continental scale vars

p1 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 1K sequences", x = "Site") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme +
scale_x_discrete(labels = c(rep("", 10)))

p2 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(duration_infection, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 1K sequences", x = "Infection duration") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p3 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(tmin, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Average minimum temperature") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p4 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(MAT, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "MAT") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p5 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(tmax, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Average maximum temperature") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

p6 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(ppt, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Average annual precipitation") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

pdf("prelim_figs/ASV_richness_by_climate_1K_min.pdf", width = 24, height = 8)
grid.arrange(p1,p3,p5,p2,p4,p6,nrow = 2, widths = c(0.31,0.31,0.38))
dev.off()

#site/diseaxe severity vars


p1 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.dbh, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 1K sequences", x = "DBH") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme +
scale_x_discrete(labels = c(rep("", 10)))

p2 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.NeoFruiting, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 1K sequences", x = "Mean neonectria fruiting (0-5)") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p3 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.Wax, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Mean wax load (0-5)") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p4 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.RaisedCanker, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Mean cankers") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p5 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.TreeCond, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Meantree condition (0-5)") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

p6 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.Xylococcus, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Mean Xylococcus") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

pdf("prelim_figs/ASV_richness_by_disease_severity_1K_min.pdf", width = 24, height = 8)
grid.arrange(p1,p3,p5,p2,p4,p6,nrow = 2, widths = c(0.31,0.31,0.38))
dev.off()

#Neo occurence

pdf("prelim_figs/ASV_richness_by_neo_occurence_1K_min.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(occurence, richness)) +
geom_boxplot(width = .25) +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 1K sequences", x = "Neonectria occurence") +
scale_color_manual(values = cbPalette, guide = F) +
scale_x_discrete(labels = c("both", "N. ditissima", "N. faginata", "none")) +
my_gg_theme
dev.off()

##########################
#5000 seqs per sample min#

#rarefied table
asv_tab.gt5K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt5K.rare.rds")
#NMDS
asv_tab.gt5K.rare.mds = readRDS(file = "intermediate_RDS/asv_tab.gt5K.rare.mds.rds")

#Add metadata

asv_tab.gt5K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt5K.rare.mds$points), asv_tab.gt5K.rare.mds$points),
full_metadata, by = "sample")

asv_tab.gt5K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt5K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

#PLOT
asv_tab.gt5K.rare.mds.metadata.neoOcurence$Site = factor(asv_tab.gt5K.rare.mds.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

#climate/continental scale vars

p1 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = duration_infection, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Infection\nduration") +
labs(title = "Infection duration")

p2 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, shape = occurence, color = Site)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
)+
labs(title = "Site", x = "")

p3 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = tmin, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Min.\ntemp.")+
labs(title = "Average minimum temperature", y = "", x = "")

p4 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = MAT, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "MAT")+
labs(title = "MAT", y = "")

p5 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = tmax, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Max.\ntemp.")+
labs(title = "Average maximum temperature", y = "", x = "")

p6 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = ppt, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Annual\nprecip.")+
labs(title = "Average annual precipitation", y = "")


pdf("NMDS_fits/NMDS_neo_occurence_by_climate_5Kmin.pdf", width = 22, height = 12)
grid.arrange(p2,p3,p5,p1,p4,p6,nrow = 2)
dev.off()

#Disease severity vars

p1 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, shape = occurence, color = Site_mean.dbh)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "DBH")+
labs(title = "DBH (cm)", x = "")

p2 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.NeoFruiting, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Neonectria\nfruting")+
labs(title = "Mean neonectria fruiting (0-5)", y = "", x = "")

p3 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.Wax, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Wax")+
labs(title = "Mean wax load (0-5)", y = "")

p4 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.RaisedCanker, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Cankers")+
labs(title = "Mean cankers", y = "", x = "")

p5 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.TreeCond, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Tree\ncondition")+
labs(title = "Mean tree condition (0-5)", y = "")

p6 = ggplot(asv_tab.gt5K.rare.mds.metadata.neoOcurence %>% filter(bench.control == "n"), aes(MDS1, MDS2, color = Site_mean.Xylococcus, shape = occurence)) +
geom_point(size = 5) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(16,15,17,3), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "Xylococcus")+
labs(title = "Mean Xylococcus", y = "")


pdf("NMDS_fits/NMDS_neo_occurence_by_disease_severity_5Kmin.pdf", width = 22, height = 12)
grid.arrange(p1,p2,p3,p4,p5,p6,nrow = 2)
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

Nf_Nd_site_freq.site_info = left_join(Nf_Nd_site_freq, site_info) %>%
left_join(., site_climate) %>%
left_join(., site_means)

p1 = ggplot(Nf_Nd_site_freq.site_info , aes(duration_infection, Nd/total)) +
geom_point(size = 2) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(duration_infection, Nf/total)) +
geom_point(size = 2) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Infection duration (yrs)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p3 = ggplot(Nf_Nd_site_freq.site_info , aes(MAT, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p4 = ggplot(Nf_Nd_site_freq.site_info , aes(MAT, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "MAT") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p5 = ggplot(Nf_Nd_site_freq.site_info , aes(ppt, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p6 = ggplot(Nf_Nd_site_freq.site_info , aes(ppt, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Annual precip (cm)") +
expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p7 = ggplot(Nf_Nd_site_freq.site_info , aes(tmax, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p8 = ggplot(Nf_Nd_site_freq.site_info , aes(tmax, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Average max. temp.") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

p9 = ggplot(Nf_Nd_site_freq.site_info , aes(tmin, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p10 = ggplot(Nf_Nd_site_freq.site_info , aes(tmin, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Average min. temp.") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_duration_infection_and_climate_5K_min.pdf", width = 25, height = 8)
grid.arrange(p1,p7,p3,p9,p5,p2,p8,p4,p10,p6,nrow= 2)
dev.off()

p1 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.dbh, Nd/total)) +
geom_point(size = 2) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.dbh, Nf/total)) +
geom_point(size = 2) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Mean dbh (cm)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p3 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.NeoFruiting, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p4 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.NeoFruiting, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean neonectria fruiting (0-5)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p5 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.Wax, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p6 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.Wax, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean wax load (0-5)") +
#expand_limits(x=c(750,2050)) +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme


p7 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.RaisedCanker, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p8 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.RaisedCanker, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean cankers") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

p9 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.TreeCond, Nd/total)) +
geom_point(size = 2) +
labs(y = "", x = "") +
#scale_color_brewer(palette = "Dark2", guide = F) +
my_gg_theme

p10 = ggplot(Nf_Nd_site_freq.site_info , aes(Site_mean.TreeCond, Nf/total)) +
geom_point(size = 2) +
labs(y = "", x = "Mean tree condition (0-5)") +
#scale_color_brewer(palette = "Dark2") +
my_gg_theme

pdf("prelim_figs/Neonectria_frequency_by_disease_severity_5K_min.pdf", width = 25, height = 8)
grid.arrange(p1,p3,p5,p7,p9,p2,p4,p6,p8,p10,nrow= 2)
dev.off()

#####################
#Richness by duration_infection etc#

sample_richness = data.frame(sample = names(apply(asv_tab.gt5K.rare,1,function(x) sum(x > 0))), richness =  apply(asv_tab.gt5K.rare,1,function(x) sum(x > 0)))
sample_richness.metadata = left_join(sample_richness, id_bench_map) %>%
left_join(., metadata_map)


sample_richness.metadata.site_info = left_join(sample_richness.metadata, Nf_v_Nd.bin) %>%
left_join(., site_info, by = "Site") %>%
left_join(., site_climate) %>%
left_join(., site_means, by = "Site")

sample_richness.metadata.site_info$Site = factor(sample_richness.metadata.site_info$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

#continental scale vars

p1 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 5K sequences", x = "Site") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme +
scale_x_discrete(labels = c(rep("", 10)))

p2 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(duration_infection, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 5K sequences", x = "Infection duration") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p3 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(tmin, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Average minimum temperature") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p4 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(MAT, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "MAT") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p5 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(tmax, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Average maximum temperature") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

p6 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(ppt, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Average annual precipitation") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

pdf("prelim_figs/ASV_richness_by_climate_5K_min.pdf", width = 24, height = 8)
grid.arrange(p1,p3,p5,p2,p4,p6,nrow = 2, widths = c(0.31,0.31,0.38))
dev.off()

#site/diseaxe severity vars


p1 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.dbh, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 5K sequences", x = "DBH") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme +
scale_x_discrete(labels = c(rep("", 10)))

p2 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.NeoFruiting, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 5K sequences", x = "Mean neonectria fruiting (0-5)") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p3 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.Wax, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Mean wax load (0-5)") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p4 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.RaisedCanker, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Mean cankers") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"),
guide = F
) +
my_gg_theme

p5 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.TreeCond, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Meantree condition (0-5)") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

p6 = ggplot(sample_richness.metadata.site_info %>% filter(bench.control == "n"), aes(Site_mean.Xylococcus, richness, group = Site, color = Site)) +
geom_boxplot() +
#facet_wrap(~run.seq) +
labs(y = "", x = "Mean Xylococcus") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")
) +
my_gg_theme

pdf("prelim_figs/ASV_richness_by_disease_severity_5K_min.pdf", width = 24, height = 8)
grid.arrange(p1,p3,p5,p2,p4,p6,nrow = 2, widths = c(0.31,0.31,0.38))
dev.off()

#Neo occurence

pdf("prelim_figs/ASV_richness_by_neo_occurence_5K_min.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(occurence, richness)) +
geom_boxplot(width = .25) +
#facet_wrap(~run.seq) +
labs(y = "ASV richness\nper 5K sequences", x = "Neonectria occurence") +
scale_color_manual(values = cbPalette, guide = F) +
scale_x_discrete(labels = c("both", "N. ditissima", "N. faginata", "none")) +
my_gg_theme
dev.off()
