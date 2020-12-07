
require(tidyverse)
require(vegan)
require(indicspecies)
require(Ternary)

source("~/ggplot_theme.txt")
#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

#filter by 1K min seqs and set binary
lt_1K_samps = full_metadata %>%
filter(total_seqs < 1000) %>%
dplyr::select("sample")

asv_tab.gt1K = t(asv_tab[,!colnames(asv_tab) %in% lt_1K_samps$sample])
asv_tab.gt1K.ten_perc_ASVs = asv_tab.gt1K[,colSums(asv_tab.gt1K > 0) > 10]
#convert to binary
asv_tab.gt1K.ten_perc_ASVs[asv_tab.gt1K.ten_perc_ASVs > 0] = 1
#180 ASVs

full_metadata.sorted = left_join(
    data.frame(sample = rownames(asv_tab.gt1K.ten_perc_ASVs)),
    full_metadata
)


par(mfrow=c(1, 1), mar=rep(0.5, 4))
TernaryPlot(point="up", atip='DC = 0', btip='DC = 1', ctip='DC = 2 or 3', alab='% association with healthy trees', blab='% association with moderately healthy trees', clab='% association with declining trees')

data_points <- list(
R = c(1, 0, 0),
O = c(1, 0.5, 0),
Y = c(0.3, 0.2, 0.1),
G = c(0, 0.1, 0)
)

AddToTernary(points, data_points, pch=21, cex=2.8)
AddToTernary(text, data_points, names(data_points), cex=0.8, font=2)

#For ternary should use frequency in each class and not absolute number of trees in a class. The algorithm sums to 1, so if there is low occurence in trees in a class the proportion of uninhabited trees relative to the other classes will be taken into account

#three levels of tree cond (as numeric)

ASV_names = colnames(asv_tab.gt1K.ten_perc_ASVs)

ASV_freq_list = list()

for(i in 1:length(ASV_names)){
    temp_samples = filter(full_metadata.sorted, TreeCond == 0) %>% select(sample)
    temp_vec = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ASV_names[i]]
    TC0 = sum(temp_vec)/length(temp_vec)
    
    temp_samples = filter(full_metadata.sorted, TreeCond == 1) %>% select(sample)
    temp_vec = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ASV_names[i]]
    TC1 = sum(temp_vec)/length(temp_vec)
    
    temp_samples = filter(full_metadata.sorted, TreeCond == 2 | TreeCond == 3) %>% select(sample)
    temp_vec = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ASV_names[i]]
    TC23 = sum(temp_vec)/length(temp_vec)

    ASV_freq_list[[ASV_names[i]]] = c(TC0,TC1,TC23)

}

par(mfrow=c(1, 1), mar=rep(0.5, 4))
TernaryPlot(point="up", atip='DC = 0', btip='DC = 1', ctip='DC = 2 or 3', alab='% association with healthy trees', blab='% association with moderately healthy trees', clab='% association with declining trees')
AddToTernary(points, ASV_freq_list, pch=16, cex=1)

pdf("ISA_spp_tables/crown_dieback_ternary.pdf")
TernaryPlot(point="up", atip='DC = 0', btip='DC = 1', ctip='DC = 2 or 3', alab='% association with healthy trees', blab='% association with moderately healthy trees', clab='% association with declining trees')
AddToTernary(points, ASV_freq_list, pch=16, cex=1)
dev.off()



ASV_freq_list = list()

for(i in 1:length(ASV_names)){
    temp_samples = filter(full_metadata.sorted, TreeCond == 0) %>% select(sample)
    temp_vec = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ASV_names[i]]
    TC0 = sum(temp_vec)/length(temp_vec)
    
    temp_samples = filter(full_metadata.sorted, TreeCond == 1 | TreeCond == 2) %>% select(sample)
    temp_vec = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ASV_names[i]]
    TC12 = sum(temp_vec)/length(temp_vec)
    
    temp_samples = filter(full_metadata.sorted, TreeCond == 3) %>% select(sample)
    temp_vec = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ASV_names[i]]
    TC3 = sum(temp_vec)/length(temp_vec)
    
    ASV_freq_list[[ASV_names[i]]] = c(TC0,TC12,TC3)
    
}

TernaryPlot(point="up", atip='DC = 0', btip='DC = 1 or 2', ctip='DC = 3', alab='% association with healthy trees', blab='% association with moderately healthy trees', clab='% association with declining trees')
AddToTernary(points, ASV_freq_list, pch=16, cex=1)

pdf("ISA_spp_tables/crown_dieback_ternary_inermediate_class_groupings.pdf")
TernaryPlot(point="up", atip='DC = 0', btip='DC = 1 or 2', ctip='DC = 3', alab='% association with healthy trees', blab='% association with moderately healthy trees', clab='% association with declining trees')
AddToTernary(points, ASV_freq_list, pch=16, cex=1)
dev.off()


TernaryPlot(alab="Redder \u2192", blab="\u2190 Greener", clab="Bluer \u2192",
lab.col=c('red', 'darkgreen', 'blue'),
point='up', lab.cex=0.8, grid.minor.lines=0,
grid.lty='solid', col=rgb(0.9, 0.9, 0.9), grid.col='white',
axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6),
padding=0.08)
# Colour the background:
cols <- TernaryPointValues(rgb)
ColourTernary(cols, spectrum = NULL)
