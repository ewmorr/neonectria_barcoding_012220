require(tidyverse)
source("~/ggplot_theme.txt")
#read data
asv_tax = read.table("dada2_out/ASVs_taxonomy.tsv", header = T)
asv_tab = read.table("dada2_out/ASVs_counts.tsv", sep = "\t", header = T, row.names = 1)
track.long = read.csv("dada2_processing_tables_figs/read_processing_tracking.csv", row.names = 1)
id_bench_map = read.table("sample_data/sample_mapping.txt", header = T)
metadata_map = read.table("sample_data/metadata.txt", header = T)
survey_dat = read.table("sample_data/trees_site_survey_data.txt", header = T, sep = "\t")
neo_cov = read.table("sample_data/plug_neonectria_coverage.txt", header = T)
site_info = read.csv("sample_data/site_info.csv")
site_means = read.table("sample_data/BBD_survey_transect_data.site_mean.txt", header = T)
site_climate = read.table("sample_data/sites_climate.txt", header = T)

#read the mapping data
asv_tab.mapping = read.table("neo_map/all_seqs.derep.otu_tab.txt", header = T)
rownames(asv_tab.mapping) = asv_tab.mapping$OTUID
asv_tab.mapping$OTUID = NULL

#table joins
metadata_ordered = full_join(metadata_map, id_bench_map)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>%
    left_join(., site_info, by = "Site") %>%
    left_join(., site_means, by = "Site") %>%
    left_join(., site_climate, by = c("Site","lat","lon"))

full_metadata = full_join(metadata_ordered, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))
if("seq.rep" %in% colnames(full_metadata)){
    full_metadata = full_metadata %>% filter(seq.rep != "y")
}


#####################################
#Remove control samples from asv_tab#
asv_tab.negatives = semi_join(
data.frame(sample = rownames(t(asv_tab)), t(asv_tab)),
full_metadata %>% filter(bench.control != "n")
)
rownames(asv_tab.negatives) = asv_tab.negatives$sample
asv_tab.negatives$sample = NULL
asv_tab.negatives = t(asv_tab.negatives)

asv_tab = asv_tab[,!colnames(asv_tab) %in% colnames(asv_tab.negatives)]
asv_tab.mapping = asv_tab.mapping[,!colnames(asv_tab.mapping) %in% colnames(asv_tab.negatives)]

#############################
#Nf and Nd counts (sum ASVs)#

Nf_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Species == "s__faginata")$asv_name
Nd_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Species == "s__ditissima")$asv_name

Nf_counts = subset(asv_tab, rownames(asv_tab) %in% Nf_asvs) %>% colSums
Nd_counts = subset(asv_tab, rownames(asv_tab) %in% Nd_asvs) %>% colSums

Nf_v_Nd = full_join(
data.frame(Nf = Nf_counts, sample = names(Nf_counts)),
data.frame(Nd = Nd_counts, sample = names(Nd_counts)),
by = "sample")

Nf_v_Nd.long = gather(Nf_v_Nd, "spp", "Neonectria_count.asv", -sample)

#also do mapping file

Nf_counts.mapping = subset(asv_tab.mapping, rownames(asv_tab.mapping) %in% Nf_asvs) %>% colSums
Nd_counts.mapping = subset(asv_tab.mapping, rownames(asv_tab.mapping) %in% Nd_asvs) %>% colSums

Nf_v_Nd.mapping = full_join(
data.frame(Nf = Nf_counts.mapping, sample = names(Nf_counts.mapping)),
data.frame(Nd = Nd_counts.mapping, sample = names(Nd_counts.mapping)),
by = "sample")

Nf_v_Nd.mapping.long = gather(Nf_v_Nd.mapping, "spp", "Neonectria_count.mapping", -sample)

#join mapping and ASV counts
Nf_Nd.asv_mapping_comp = full_join(Nf_v_Nd.long, Nf_v_Nd.mapping.long, by = c("sample", "spp"))
#vserach does not write files with zero counts, so add zero for NAs
Nf_Nd.asv_mapping_comp[is.na(Nf_Nd.asv_mapping_comp)] = 0

asv_not_count = filter(Nf_Nd.asv_mapping_comp, Neonectria_count.asv == 0 & Neonectria_count.mapping > 0) %>%
    select(sample) %>% length
count_not_asv = filter(Nf_Nd.asv_mapping_comp, Neonectria_count.asv > 0 & Neonectria_count.mapping == 0) %>%
    select(sample) %>% length

p = ggplot(Nf_Nd.asv_mapping_comp, aes(Neonectria_count.asv +1, Neonectria_count.mapping +1)) +
geom_point(alpha = 0.5) +
scale_x_log10() +
scale_y_log10() +
my_gg_theme +
labs(x = "ASV counts + 1", y = "Mapping counts + 1")

pdf("neo_map/Neonectria_ASV_v_mapping_counts.pdf")
print(p)
dev.off()

