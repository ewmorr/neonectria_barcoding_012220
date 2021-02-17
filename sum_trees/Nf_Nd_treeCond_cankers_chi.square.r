require(dplyr)
source("~/ggplot_theme.txt")

source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

lt_1K_samps = full_metadata %>%
filter(total_seqs < 1000) %>%
dplyr::select("sample")

#Species matrix on Nf and Nd occurence
Nf_v_Nd.bin.gt1K = Nf_v_Nd.bin %>%
filter(!sample %in% lt_1K_samps$sample)

full_metadata.sorted = left_join(
Nf_v_Nd.bin.gt1K,
full_metadata,
by = "sample"
)

full_metadata.sorted$occurence = factor(full_metadata.sorted$occurence, levels = c("none", "Nf", "Nd", "both"))

tree_cond_table.Nd = full_metadata.sorted %>% group_by(TreeCond, Nd) %>% summarize(n = n())
tree_cond_table.Nf = full_metadata.sorted %>% group_by(TreeCond, Nf) %>% summarize(n = n())
cankers_table.Nd = full_metadata.sorted %>% group_by(RaisedCanker, Nd) %>% summarize(n = n())
cankers_table.Nf = full_metadata.sorted %>% group_by(RaisedCanker, Nf) %>% summarize(n = n())

tree_cond_table.Nd.mat = matrix(tree_cond_table.Nd$n, 2, 4,
    dimnames = list(Nd = c("0", "1"), tree_cond = levels(as.factor(tree_cond_table.Nd$TreeCond)) )
)
tree_cond_table.Nf.mat = matrix(tree_cond_table.Nf$n, 2, 4,
dimnames = list(Nf = c("0", "1"), tree_cond = levels(as.factor(tree_cond_table.Nd$TreeCond)) )
)

tree_cond_table.Nd.test = chisq.test(tree_cond_table.Nd.mat)
tree_cond_table.Nf.test  = chisq.test(tree_cond_table.Nf.mat)

tree_cond_table.Nd.test = fisher.test(tree_cond_table.Nd.mat)
tree_cond_table.Nf.test  = fisher.test(tree_cond_table.Nf.mat)

#Neither is significant, however we already show sig positive association for Nd using HMSC and ISA. Probably if I split the crown dieback into 0-1 and 2-3 for Nd this would show positive (as it is for ISA).

#Nd cnakers 0-2 versus 3
chisq.test(matrix(c(25+14+17, 9+4+5, 14, 14), 2,2,
    dimnames = list(Nd = c("0", "1"), cankers = c("0-2", "3"))
)
)
#Nf cankers 0 vs 1-3
chisq.test(matrix(c(21, 13, 3+2+5, 15+20+23), 2,2,
dimnames = list(Nf = c("0", "1"), cankers = c("0-2", "3"))
)
)

#Nd crown 0-1 vs 2-3
chisq.test(matrix(c(15+32, 4+10, 18+5, 11+7), 2,2,
dimnames = list(Nd = c("0", "1"), cankers = c("0-1", "2-3"))
)
)
