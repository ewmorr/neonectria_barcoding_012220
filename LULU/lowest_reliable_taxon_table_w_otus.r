require(tidyverse)

#read data
source("~/repo/neonectria_barcoding_012220/read_ASV_dat.r")
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

#and calculates neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata
#new table of ASV name by informative taxa

####################
#ASV count by taxon#
ASV_count_by_taxon = data.frame(ASV = row.names(asv_tax), asv_tax) %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarize(ASV_count = n())

##################################
#OTUs per taxon at dif sim levels#

#read OTU names files
OTUs.1 = read.table("seq_sim_ASV_cluster/ASVs.centroids.1.txt", header = F) %>% as.data.frame
colnames(OTUs.1) = "ASV"
OTUs.99 = read.table("seq_sim_ASV_cluster/ASVs.centroids.0.99.txt", header = F) %>% as.data.frame
colnames(OTUs.99) = "ASV"
OTUs.98 = read.table("seq_sim_ASV_cluster/ASVs.centroids.0.98.txt", header = F) %>% as.data.frame
colnames(OTUs.98) = "ASV"
OTUs.97 = read.table("seq_sim_ASV_cluster/ASVs.centroids.0.97.txt", header = F) %>% as.data.frame
colnames(OTUs.97) = "ASV"
OTUs.96 = read.table("seq_sim_ASV_cluster/ASVs.centroids.0.96.txt", header = F) %>% as.data.frame
colnames(OTUs.96) = "ASV"
OTUs.95 = read.table("seq_sim_ASV_cluster/ASVs.centroids.0.95.txt", header = F) %>% as.data.frame
colnames(OTUs.95) = "ASV"

#filter out plant and animal ASVs
OTUs.1 = OTUs.1 %>% filter(!ASV %in% plant_asvs & !ASV %in% animal_asvs)
OTUs.99 = OTUs.99 %>% filter(!ASV %in% plant_asvs & !ASV %in% animal_asvs)
OTUs.98 = OTUs.98 %>% filter(!ASV %in% plant_asvs & !ASV %in% animal_asvs)
OTUs.97 = OTUs.97 %>% filter(!ASV %in% plant_asvs & !ASV %in% animal_asvs)
OTUs.96 = OTUs.96 %>% filter(!ASV %in% plant_asvs & !ASV %in% animal_asvs)
OTUs.95 = OTUs.95 %>% filter(!ASV %in% plant_asvs & !ASV %in% animal_asvs)

#summarize OTU count by taxon
OTUs.1.by_taxon = left_join(
    OTUs.1,
    data.frame(ASV = row.names(asv_tax), asv_tax),
    by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(OTU_100_perc = n())

OTUs.99.by_taxon = left_join(
OTUs.99,
data.frame(ASV = row.names(asv_tax), asv_tax),
by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(OTU_99_perc = n())

OTUs.98.by_taxon = left_join(
OTUs.98,
data.frame(ASV = row.names(asv_tax), asv_tax),
by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(OTU_98_perc = n())

OTUs.97.by_taxon = left_join(
OTUs.97,
data.frame(ASV = row.names(asv_tax), asv_tax),
by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(OTU_97_perc = n())

OTUs.96.by_taxon = left_join(
OTUs.96,
data.frame(ASV = row.names(asv_tax), asv_tax),
by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(OTU_96_perc = n())

OTUs.95.by_taxon = left_join(
OTUs.95,
data.frame(ASV = row.names(asv_tax), asv_tax),
by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(OTU_95_perc = n())

#join table and print

ASV_OTU_count_by_taxon = full_join(ASV_count_by_taxon, OTUs.1.by_taxon) %>%
    full_join(., OTUs.99.by_taxon) %>%
    full_join(., OTUs.98.by_taxon) %>%
    full_join(., OTUs.97.by_taxon) %>%
    full_join(., OTUs.96.by_taxon) %>%
    full_join(., OTUs.95.by_taxon)

write.table(ASV_OTU_count_by_taxon, "tables/ASV_and_OTU_count_by_taxon.txt", row.names = F, quote = F, sep = "\t")

#####################################
#Sum sample and site occurence for each taxon#
samples_to_sum = full_metadata %>% filter(bench.control == "n" & locus == "ITS2") %>% select(sample, Site)

asv_tab.tax = inner_join(data.frame(ASV = row.names(asv_tax), asv_tax),
data.frame(ASV = rownames(asv_tab), asv_tab)
)

#summarize by taxon (sequence sums)
taxon_tab.tax = asv_tab.tax %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarize_if(is.numeric, sum)

#convert to binary and add taxon names back
taxon_tab = taxon_tab.tax[8:length(colnames(taxon_tab.tax))]
taxon_tab[taxon_tab > 0] = 1
#join Site data
taxon_tab_w_site.bin = inner_join(samples_to_sum, data.frame(sample = colnames(taxon_tab), t(taxon_tab)), by = "sample")

#Sum the number of sample that a taxon occurs in by site
asv_tab_w_site.sum_sample = taxon_tab_w_site.bin %>% group_by(Site) %>% summarize_if(is.numeric, sum) %>% as.data.frame
#add taxon names
rownames(asv_tab_w_site.sum_sample) = asv_tab_w_site.sum_sample$Site
asv_tab_w_site.sum_sample$Site = NULL
asv_tab_w_site.sum_sample.tax = data.frame(taxon_tab.tax[,1:7], t(asv_tab_w_site.sum_sample))

#binary by site
asv_tab_w_site.sum_sample.bin = t(asv_tab_w_site.sum_sample)
asv_tab_w_site.sum_sample.bin[asv_tab_w_site.sum_sample.bin>0] = 1

#make df

site_and_sample_counts_by_taxon = data.frame(
    asv_tab_w_site.sum_sample.tax[,1:7],
    sample_count = apply(asv_tab_w_site.sum_sample.tax[,8:length(colnames(asv_tab_w_site.sum_sample.tax))], 1, sum),
    site_count = apply(asv_tab_w_site.sum_sample.bin, 1, sum)
)

##############
#Join ASV and sample/site counts
taxon_all_counts = full_join(ASV_count_by_taxon, site_and_sample_counts_by_taxon)

write.table(taxon_all_counts, file = "tables/taxon_ASV_site_sample_counts.txt", sep = "\t", quote = F, row.names = F)

