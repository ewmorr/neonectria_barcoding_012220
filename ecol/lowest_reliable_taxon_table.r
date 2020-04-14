require(tidyverse)

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

#and calculates neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

system("mkdir tables")

#new table of ASV name by informative taxa

write.table(data.frame(ASV = names(asv_informative_taxa), taxon = unname(asv_informative_taxa)), file = "dada2_out/asv_lowest_informative_taxon.txt", row.names = F, col.names = F, quote = F, sep = "\t")

asv_informative_taxa.df = data.frame(ASV = names(asv_informative_taxa), taxon = unname(asv_informative_taxa))

####################
#ASV count by taxon#
ASV_count_by_taxon = data.frame(ASV = row.names(asv_tax), asv_tax) %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarize(ASV_count = n())

write.table(ASV_count_by_taxon, file = "tables/ASV_count_by_taxon.txt", sep = "\t", row.names = F, quote = F)

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

