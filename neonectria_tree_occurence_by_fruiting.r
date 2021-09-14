#library(tidyverse)
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

plot(occurence ~ NeoFruiting, full_metadata.sorted)

full_metadata.sorted %>% select(c(occurence, NeoFruiting)) %>% group_by(occurence, NeoFruiting) %>% summarize (n())
full_metadata.sorted %>% select(c(occurence, NeoFruiting)) %>% group_by(NeoFruiting) %>% summarize (n())
