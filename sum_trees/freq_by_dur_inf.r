library(ggplot2)
library(tidyverse)
#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")


#filter by 1K min seqs
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
    temp_tab = full_metadata.sorted %>% filter(Site == site_info$Site[i])
    Nf_Nd_site_freq.spp$N.faginata[i] = ((filter(temp_tab, occurence == "Nf" | occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq.spp$N.ditissima[i] = ((filter(temp_tab, occurence == "Nd" | occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq.spp$Site[i] = as.character(site_info$Site[i])
}

###############################



Nf_Nd_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
neither = vector(mode = "numeric", length = length(site_info$Site)),
N.faginata = vector(mode = "numeric", length = length(site_info$Site)),
N.ditissima = vector(mode = "numeric", length = length(site_info$Site)),
both = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted %>% filter(Site == site_info$Site[i])
    Nf_Nd_site_freq$neither[i] = ((filter(temp_tab, occurence == "none")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$N.faginata[i] = ((filter(temp_tab, occurence == "Nf")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$N.ditissima[i] = ((filter(temp_tab, occurence == "Nd")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$both[i] = ((filter(temp_tab, occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(site_info$Site[i])
}

site_n = full_metadata.sorted %>% group_by(Site) %>% summarize(n = n())

Nf_Nd_site_freq.metadata = left_join(
Nf_Nd_site_freq,
site_info %>% select(Site, duration_infection)
)

plot(both ~ duration_infection, data = Nf_Nd_site_freq.metadata)
plot(both+N.ditissima ~ duration_infection, data = Nf_Nd_site_freq.metadata)
plot(both+N.faginata ~ duration_infection, data = Nf_Nd_site_freq.metadata)
summary(lm(both ~ duration_infection, data = Nf_Nd_site_freq.metadata))
