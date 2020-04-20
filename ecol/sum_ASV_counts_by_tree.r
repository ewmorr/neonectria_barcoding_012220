require(tidyverse)

source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")

site_tree_label = full_metadata %>% filter(bench.control == "n" & locus == "ITS2") %>% select(c("Site", "Tree", "metadata.label"))

site_tree_label$Site.tree = interaction(site_tree_label$Site, site_tree_label$Tree)

#########################
#Make sure colnames match

b = colnames(asv_tab)
a = site_tree_label$metadata.label

not_in_tab = unique(a[! a %in% b])
site_tree_label = site_tree_label %>% filter(!metadata.label %in% not_in_tab)

####################################################
#new table indexed by tree instead of metadata.label
asv_tab.tree_sum = data.frame(matrix(
    data = 0,
    nrow = length(rownames(asv_tab)),
    ncol = site_tree_label$Site.tree %>% unique %>% length
))

rownames(asv_tab.tree_sum) = rownames(asv_tab)
colnames(asv_tab.tree_sum) = site_tree_label$Site.tree %>% unique

#add cols
for(i in 1:length(site_tree_label$Site.tree)){
    #tree.col = asv_tab.tree_sum[,colnames(asv_tab.tree_sum) == site_tree_label$Site.tree[i]]
    print(paste(site_tree_label$Site.tree[i], "plus", site_tree_label$metadata.label[i]))
    asv_tab.tree_sum[,colnames(asv_tab.tree_sum) == site_tree_label$Site.tree[i]] = asv_tab.tree_sum[,colnames(asv_tab.tree_sum) == site_tree_label$Site.tree[i]] + asv_tab[,colnames(asv_tab) == site_tree_label$metadata.label[i]]
}

#check that rowSums are equal
rowSums(asv_tab.tree_sum) - rowSums(asv_tab)


#################
#Testing routine#

#test_tab = data.frame(a = rep(1,5), b = seq(1, 5), c = rep(2,5), d = seq(5, 9))
#label_match = data.frame(to_col = c("x","x","z", "z"), from_col = c("a", "b", "c", "d"))

#test_tab.sum = data.frame(matrix(
#    data = 0,
#    nrow = length(rownames(test_tab)),
#    ncol = label_match$to_col %>% unique %>% length
#))
#colnames(test_tab.sum) = label_match$to_col %>% unique

#for(i in 1:length(label_match$to_col)){
#    test_tab.sum[,colnames(test_tab.sum) == label_match$to_col[i]] = test_tab.sum[,colnames(test_tab.sum) == label_match$to_col[i]] + test_tab[,colnames(test_tab) == label_match$from_col[i]]
#}
