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

#############################################################
#Avg global species richness at 1K and 5K subsampling depths#

spp_rich.1k = data.frame(sample_number = vector(mode = "numeric", length = 100), spp_rich = vector(mode = "numeric", length = 100))

spp_rich.5k = data.frame(sample_number = vector(mode = "numeric", length = 100), spp_rich = vector(mode = "numeric", length = 100))

for(i in 1:100){
    print(i)
    asv_tab.1K = rrarefy(t(asv_tab.no_singletons[,colSums(asv_tab.no_singletons) >= 1000]), 1000)
    asv_tab.5K = rrarefy(t(asv_tab.no_singletons[,colSums(asv_tab.no_singletons) >= 5000]), 5000)
    
    spp_rich.1k$sample_number[i] = i
    spp_rich.5k$sample_number[i] = i
    spp_rich.1k$spp_rich[i] = asv_tab.1K[,colSums(asv_tab.1K) > 0] %>% colnames %>% length
    spp_rich.5k$spp_rich[i] = asv_tab.5K[,colSums(asv_tab.5K) > 0] %>% colnames %>% length
}

avg_taxa.1k = sum(spp_rich.1k$spp_rich)/100
avg_taxa.5k = sum(spp_rich.5k$spp_rich)/100
sd_taxa.1k = sd(spp_rich.1k$spp_rich)
sd_taxa.5k = sd(spp_rich.5k$spp_rich)

total_taxa = rownames(asv_tab.no_singletons) %>% length
#############################################

#taxa_remaining.1k = asv_tab.1K[,colSums(asv_tab.1K) > 0] %>% colnames %>% length
#taxa_remaining.5k = asv_tab.5K[,colSums(asv_tab.5K) > 0] %>% colnames %>% length

######
#PLOT#
######

p = ggplot(sequence_counts, aes(rank, sequences)) +
geom_col() +
scale_y_log10(labels = fancy_scientific) +
geom_hline(yintercept = 1000) +
geom_hline(yintercept = 5000) +
my_gg_theme +
labs(title = paste("199 samples, ", total_seqs, "sequences, 858 taxa"), x = "sample rank" )+
annotate(x = 100, y = 75, geom = "text", size = 6, label = paste("rarefy 1K sequences: ", samples_dropped.1k, "samples dropped\n", seqs_dropped.1k, "seqs dropped, ", total_taxa - avg_taxa.1k, " taxa dropped\navg taxa", avg_taxa.1k, "+/-", signif(sd_taxa.1k,2), "stdev")) +
annotate(x = 100, y = 10^5, geom = "text", size = 6, label = paste("rarefy 5K sequences: ", samples_dropped.5k, "samples dropped\n", seqs_dropped.5k, "seqs dropped, ", total_taxa - avg_taxa.5k, " taxa dropped\navg taxa", avg_taxa.5k, "+/-", signif(sd_taxa.5k,2), "stdev"))

pdf("prelim_figs/sequence_subsampling.pdf", width = 8, height = 6)
print(p)
dev.off()

############################
#Rarefaction based analyses#

asv_tab.no_singletons.t = t(asv_tab.no_singletons)

asv_tab.no_singletons.rare_curve = rarecurve(asv_tab.no_singletons.t, step = 100, label = F)
names(asv_tab.no_singletons.rare_curve) = rownames(asv_tab.no_singletons.t)#paste("sample", 1:length(rownames(asv_tab.no_singletons.t)), sep = "")

protox <- mapply(FUN = function(x, y) {
    mydf <- as.data.frame(x)
    colnames(mydf) <- "value"
    mydf$sample <- y
    mydf$subsample <- attr(x, "Subsample")
    mydf
}, x = asv_tab.no_singletons.rare_curve, y = as.list(names(asv_tab.no_singletons.rare_curve)), SIMPLIFY = FALSE)

xy <- do.call(rbind, protox)
xy = rbind(protox)
rownames(xy) <- NULL  # pretty

p1 = ggplot(xy, aes(x = subsample, y = value, group = sample)) +
geom_vline(xintercept = 1000) +
geom_vline(xintercept = 5000) +
geom_line(alpha = 0.3) +
my_gg_theme +
#scale_x_continuous(limits= c(0,50000)) +
labs(y = "ASVs", x = "sequences subsampled")

p2 = ggplot(xy, aes(x = subsample, y = value, group = sample)) +
geom_vline(xintercept = 1000) +
geom_vline(xintercept = 5000) +
geom_line(alpha = 0.3) +
my_gg_theme +
scale_x_continuous(limits= c(0,50000)) +
labs(y = "ASVs", x = "sequences subsampled")

p3 = ggplot(xy, aes(x = subsample, y = value, group = sample)) +
geom_vline(xintercept = 1000) +
geom_vline(xintercept = 5000) +
geom_line(alpha = 0.3) +
my_gg_theme +
scale_x_continuous(limits= c(0,20000)) +
labs(y = "ASVs", x = "sequences subsampled")

p4 = ggplot(xy, aes(x = subsample, y = value, group = sample)) +
geom_vline(xintercept = 1000) +
geom_vline(xintercept = 5000) +
geom_line(alpha = 0.3) +
my_gg_theme +
scale_x_continuous(limits= c(0,10000)) +
labs(y = "ASVs", x = "sequences subsampled")

require(gridExtra)

pdf("prelim_figs/rarefaction_curves.pdf", width = 16, height = 10)
grid.arrange(p1,p2,p3,p4, ncol = 2)
dev.off()

#average richness per sample
spp_rich_rarefy.1k = rarefy(asv_tab.no_singletons.t[rowSums(asv_tab.no_singletons.t) >= 1000,], sample = 1000)

spp_rich_rarefy.5k = rarefy(asv_tab.no_singletons.t[rowSums(asv_tab.no_singletons.t) >= 5000,], sample = 5000)

