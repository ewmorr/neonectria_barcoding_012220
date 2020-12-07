
require(tidyverse)
require(vegan)
require(indicspecies)

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

#three levels of tree cond (as numeric)
tree_cond.levels = c(0,1,2,3)
ASV_names = colnames(asv_tab.gt1K.ten_perc_ASVs)

tree_cond_freq_table = data.frame(
    treeCond = vector(mode = "numeric", length = 4*ncol(asv_tab.gt1K.ten_perc_ASVs)),
    ASV = vector(mode = "character", length = 4*ncol(asv_tab.gt1K.ten_perc_ASVs)),
    frequency = vector(mode = "numeric", length = 4*ncol(asv_tab.gt1K.ten_perc_ASVs)),
    stringsAsFactors = F
)

index = 0
for(i in 1:length(tree_cond.levels)){
    temp_samples = filter(full_metadata.sorted, TreeCond == tree_cond.levels[i]) %>% select(sample)
    temp_table = asv_tab.gt1K.ten_perc_ASVs[rownames(asv_tab.gt1K.ten_perc_ASVs) %in% temp_samples$sample, ]
    
    for(u in 1:length(ASV_names)){
        index = index+1
        tree_cond_freq_table$treeCond[index] = tree_cond.levels[i]
        tree_cond_freq_table$ASV[index] = ASV_names[u]
        tree_cond_freq_table$frequency[index] = sum(temp_table[,u])/nrow(temp_table)
    }
}

tree_cond_freq_table.no_absent = filter(tree_cond_freq_table, frequency != 0)

tree_cond_freq_table$frequency[tree_cond_freq_table$frequency == 0] = NA


tree_cond_freq_table.no_absent = full_join(
tree_cond_freq_table.no_absent,
tree_cond_freq_table.no_absent %>% group_by(ASV) %>% summarize(sum.treeCond = sum(treeCond))
)
tree_cond_freq_table.no_absent = full_join(
tree_cond_freq_table.no_absent,
tree_cond_freq_table.no_absent %>% group_by(ASV) %>% summarize(sum.freq = sum(frequency))
)

ggplot(tree_cond_freq_table.no_absent, aes(treeCond, reorder(ASV, treeCond), group = ASV, alpha = frequency)) +
geom_point(size = 1) +
geom_line(size = 1, na.rm = T) +
labs(x = "Crown dieback class", y = "ASVs") +
my_gg_theme +
theme(
axis.text.y = element_blank()
)

ggplot(tree_cond_freq_table.no_absent, aes(y = treeCond, x = reorder(ASV, sum.freq), group = ASV, alpha = frequency)) +
geom_point(size = 1) +
geom_line(size = 1, na.rm = T) +
labs(y = "Crown dieback class", x = "ASVs") +
my_gg_theme +
theme(
axis.text.x = element_blank()
)

p = ggplot(tree_cond_freq_table.no_absent, aes(y = treeCond, x = reorder(reorder(ASV, sum.freq), treeCond), group = ASV, alpha = frequency)) +
geom_point(size = 1) +
geom_line(size = 1, na.rm = T) +
labs(y = "Crown dieback class", x = "ASVs") +
my_gg_theme +
labs(alpha = "ASV frequncy") +
theme(
axis.text.x = element_blank(),
legend.title = element_text(size = 18)
)

pdf("ISA_spp_tables/ISA_crown_dieback_ASV_progression_all_ASVs_colored_by_frequency.pdf", width = 12, height = 6)
p
dev.off()
