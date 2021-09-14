require(tidyverse)
require(Hmisc)
require(ppcor)

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

##############################
#Calculate Nf and Nd frequency as individual columns (i.e., no "both")

Nf_Nd_site_freq.spp = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
N.faginata = vector(mode = "numeric", length = length(site_info$Site)),
N.ditissima = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted %>% filter(Site == as.character(site_info$Site[i]))
    Nf_Nd_site_freq.spp$N.faginata[i] = ((filter(temp_tab, occurence == "Nf" | occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq.spp$N.ditissima[i] = ((filter(temp_tab, occurence == "Nd" | occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq.spp$Site[i] = as.character(site_info$Site[i])
}

Nf_Nd_site_freq.spp.metadata = full_join(Nf_Nd_site_freq.spp, site_info) %>% full_join(., site_climate.GDD)
pairwise_df = Nf_Nd_site_freq.spp.metadata %>% dplyr::select(c(N.faginata, duration_infection, HDD4.mean_nongrowing))

plot(N.faginata ~ duration_infection, data = pairwise_df)
plot(N.faginata ~ HDD4.mean_nongrowing, data = pairwise_df)

for(i in 1:ncol(pairwise_df)){
    print(colnames(pairwise_df)[i])
    print(shapiro.test(pairwise_df[,i]))
}
#Nf p = 0.03

#corr.matrix.pear = rcorr(as.matrix(pairwise_df), type = "pearson")
corr.matrix.spear = rcorr(as.matrix(pairwise_df), type = "spearman")

#cor.test(pairwise_df$N.faginata, pairwise_df$duration_infection, method = "spearman")
#cor.test(pairwise_df$N.faginata, pairwise_df$HDD4.mean_nongrowing, method = "spearman")

pairwise_partial_corrs = pcor(pairwise_df, method = "spearman")










