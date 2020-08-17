require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)
source("~/ggplot_theme.txt")

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

#site frequency table

Nf_Nd_site_freq = data.frame(
    Site = vector(mode = "character", length = length(site_info$Site)),
    N.faginata = vector(mode = "numeric", length = length(site_info$Site)),
    N.ditissima = vector(mode = "numeric", length = length(site_info$Site)),
#    total = vector(mode = "numeric", length = length(site_info$Site)),
    stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted %>% filter(Site == site_info$Site[i])
    Nf_Nd_site_freq$N.faginata[i] = ((filter(temp_tab, Nf > 0)) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$N.ditissima[i] = ((filter(temp_tab, Nd > 0)) %>% nrow)/nrow(temp_tab)
#    Nf_Nd_site_freq$total[i] = nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(site_info$Site[i])
}

Nf_Nd_site_freq.long = Nf_Nd_site_freq %>%
    pivot_longer(-Site, names_to = "spp", values_to = "frequency")

Nf_Nd_site_freq.metadata = left_join(
    Nf_Nd_site_freq.long,
    site_info %>% select(Site, duration_infection)
) %>% left_join(
., site_climate.GDD %>% select(Site, HDD4.mean_nongrowing)
)


p1 = ggplot(Nf_Nd_site_freq.metadata, aes(
    x = duration_infection,
    y = frequency,
    fill = HDD4.mean_nongrowing,
    shape = spp)
) +
geom_point(position = position_dodge(width = 3), size = 4) +
scale_fill_gradient2("Nongrowing\nseason GDD", low='black',mid='grey',high = 'white',midpoint=150) +
scale_shape_manual(values = c(21,22)) +
labs(shape = "", x = "Duration infection (yrs)", y = "Site-level species frequency\n(proportion trees detected)") +
my_gg_theme +
theme(
    legend.title = element_text(size = 20, hjust = 0)
)

pdf("prelim_figs/Nd_Nf_frequency_site.dur-inf-GDD.pdf", width = 8,height = 4)
p1
dev.off()
