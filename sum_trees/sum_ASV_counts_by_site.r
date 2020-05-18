require(tidyverse)

source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")

site_label = full_metadata %>% filter(bench.control == "n" & locus == "ITS2") %>% select(c("Site",  "metadata.label"))

#########################
#Make sure colnames match

b = colnames(asv_tab)
a = site_label$metadata.label

not_in_tab = unique(a[! a %in% b])
site_label = site_label %>% filter(!metadata.label %in% not_in_tab)

####################################################
#new table indexed by tree instead of metadata.label
asv_tab.site_sum = data.frame(matrix(
    data = 0,
    nrow = length(rownames(asv_tab)),
    ncol = site_label$Site %>% unique %>% length
))

rownames(asv_tab.site_sum) = rownames(asv_tab)
colnames(asv_tab.site_sum) = site_label$Site %>% unique

#add cols
for(i in 1:length(site_label$Site)){
    #tree.col = asv_tab.site_sum[,colnames(asv_tab.site_sum) == site_label$Site[i]]
    print(paste(site_label$Site[i], "plus", site_label$metadata.label[i]))
    asv_tab.site_sum[,colnames(asv_tab.site_sum) == site_label$Site[i]] = asv_tab.site_sum[,colnames(asv_tab.site_sum) == site_label$Site[i]] + asv_tab[,colnames(asv_tab) == site_label$metadata.label[i]]
}

#check that rowSums are equal
rowSums(asv_tab.site_sum) - rowSums(asv_tab)


