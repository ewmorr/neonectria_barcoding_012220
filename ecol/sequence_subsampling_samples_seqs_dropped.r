require(vegan)
require(tidyverse)
source("~/ggplot_theme.txt")
#read data
source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")
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


#################################################
#Filtering out taxa that occur in only one sample

asv_tab.no_singletons = asv_tab[rowSums(asv_tab > 0) > 1,]
#858 of 1173 rows remaining

sequence_counts = data.frame(rank = seq(1, length(colSums(asv_tab.no_singletons)), 1),
    sequences = sort(colSums(asv_tab.no_singletons), decreasing = T)
)

sequence_counts$rank %>% length
(sequence_counts %>% filter(sequences >= 1000))$rank %>% length
(sequence_counts %>% filter(sequences >= 5000))$rank %>% length

sequence_counts.gt1k = sequence_counts %>% filter(sequences >= 1000)
sequence_counts.gt5k = sequence_counts %>% filter(sequences >= 5000)

samples_dropped.1k = (sequence_counts$rank %>% length) - (sequence_counts %>% filter(sequences >= 1000))$rank %>% length

samples_dropped.5k = (sequence_counts$rank %>% length) - (sequence_counts %>% filter(sequences >= 5000))$rank %>% length

seqs_dropped.1k = sum(sequence_counts.gt1k$sequences - 1000)
seqs_dropped.5k = sum(sequence_counts.gt5k$sequences - 5000)
total_seqs = sum(sequence_counts$sequences)

asv_tab.1K = rrarefy(t(asv_tab.no_singletons[,colSums(asv_tab.no_singletons) >= 1000]), 1000)
asv_tab.5K = rrarefy(t(asv_tab.no_singletons[,colSums(asv_tab.no_singletons) >= 5000]), 5000)

taxa_remaining.1k = asv_tab.1K[,colSums(asv_tab.1K) > 0] %>% colnames %>% length
taxa_remaining.5k = asv_tab.1K[,colSums(asv_tab.5K) > 0] %>% colnames %>% length
total_taxa = rownames(asv_tab.no_singletons) %>% length

p = ggplot(sequence_counts, aes(rank, sequences)) +
geom_col() +
scale_y_log10(labels = fancy_scientific) +
geom_hline(yintercept = 1000) +
geom_hline(yintercept = 5000) +
my_gg_theme +
labs(title = paste("199 samples, ", total_seqs, "sequences, 858 taxa"), x = "sample rank" )+
annotate(x = 100, y = 75, geom = "text", size = 6, label = paste("rarefy 1K sequences: ", samples_dropped.1k, "samples dropped,\n", seqs_dropped.1k, "seqs dropped, ", total_taxa - taxa_remaining.1k, " taxa dropped")) +
annotate(x = 100, y = 10^5, geom = "text", size = 6, label = paste("rarefy 5K sequences: ", samples_dropped.5k, "samples dropped,\n", seqs_dropped.5k, "seqs dropped, ", total_taxa - taxa_remaining.5k, " taxa dropped"))

pdf("prelim_figs/sequence_subsampling.pdf", width = 8, height = 6)
print(p)
dev.off()

