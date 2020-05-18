require(tidyverse)
require(data.table)
get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

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

#Some gymnastics to get sample labels in the correct format. Could fix this upstream...
colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
track.long$sample <- unname(sapply(as.character(track.long$sample), get.sample.name))
#track.long$sample = paste0("X", track.long$sample)

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

#########################################
#Filter plants and animals out of tables#

plant_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Kingdom == "k__Viridiplantae")$asv_name
animal_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Kingdom == "k__Metazoa")$asv_name

asv_tab = subset(asv_tab, !rownames(asv_tab) %in% plant_asvs)
asv_tab = subset(asv_tab, !rownames(asv_tab) %in% animal_asvs)

asv_tax = subset(asv_tax, !rownames(asv_tax) %in% plant_asvs)
asv_tax = subset(asv_tax, !rownames(asv_tax) %in% animal_asvs)

################################################
#Negative & control samples table with taxonomy#

#make negatives only asv_tab (long format)
asv_tab.negatives = semi_join(
    data.frame(sample = rownames(t(asv_tab)), t(asv_tab)),
    full_metadata %>% filter(bench.control != "n")
)
rownames(asv_tab.negatives) = asv_tab.negatives$sample
asv_tab.negatives$sample = NULL
asv_tab.negatives = t(asv_tab.negatives)
#join with taxonomy
asv_tab.negatives.asvnames = data.frame(ASV = rownames(asv_tab.negatives), asv_tab.negatives)

asv_tab.negatives.long = melt(asv_tab.negatives.asvnames[rowSums(asv_tab.negatives) > 0,] %>% data.table,
id = "ASV", variable.name = "sample", value.name = "count") %>%
data.frame


#############################
#Nf and Nd counts (sum ASVs)#

#Also read mapping data in order to avg mapping counts and ASV counts for detection
#read the mapping data
asv_tab.mapping = read.table("neo_map/all_seqs.derep.otu_tab.txt", header = T)
rownames(asv_tab.mapping) = asv_tab.mapping$OTUID
asv_tab.mapping$OTUID = NULL
#remove negative controls from both tables
asv_tab = asv_tab[,!colnames(asv_tab) %in% colnames(asv_tab.negatives)]
asv_tab.mapping = asv_tab.mapping[,!colnames(asv_tab.mapping) %in% colnames(asv_tab.negatives)]
#BP225 did not make it trhough DADA2 so has to be removed manually
asv_tab.mapping$BP225 = NULL
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
#vsearch does not write files with zero counts, so add zero for NAs
Nf_Nd.asv_mapping_comp[is.na(Nf_Nd.asv_mapping_comp)] = 0
#add sample sequence totals from asv_tab
Nf_Nd.asv_mapping_comp.w_counts = left_join(Nf_Nd.asv_mapping_comp,
data.frame(sample = colnames(asv_tab), total_seqs = colSums(asv_tab))
)

#CAL AVG
Nf_Nd.asv_mapping_comp$Neo_avg = (Nf_Nd.asv_mapping_comp$Neonectria_count.asv + Nf_Nd.asv_mapping_comp$Neonectria_count.mapping) / 2

#remake count vectors
Nf_counts = as.vector((filter(Nf_Nd.asv_mapping_comp, spp == "Nf") %>% select(Neo_avg))[[1]] )
names(Nf_counts) = (filter(Nf_Nd.asv_mapping_comp, spp == "Nf") %>% select(sample))[[1]]

Nd_counts = as.vector((filter(Nf_Nd.asv_mapping_comp, spp == "Nd") %>% select(Neo_avg))[[1]] )
names(Nd_counts) = (filter(Nf_Nd.asv_mapping_comp, spp == "Nd") %>% select(sample))[[1]]

#############
#SUM by site#
#############

site_label = full_metadata %>% filter(bench.control == "n" & locus == "ITS2") %>% select(c("Site", "metadata.label"))

#########################
#Make sure colnames match

b = names(Nf_counts)
a = site_label$metadata.label

not_in_tab = unique(a[! a %in% b])
site_label = site_label %>% filter(!metadata.label %in% not_in_tab)

#New vectors indexed by tree instead of metadata.label
Nf_counts.site_sum = vector(mode = "numeric", length = site_label$Site %>% unique %>% length)
names(Nf_counts.site_sum) = site_label$Site %>% unique

Nd_counts.site_sum = vector(mode = "numeric", length = site_label$Site %>% unique %>% length)
names(Nd_counts.site_sum) = site_label$Site %>% unique

#add vals
for(i in 1:length(site_label$Site)){
    #print(paste(site_label$Site[i], "plus", site_label$metadata.label[i]))
    
    Nf_counts.site_sum[names(Nf_counts.site_sum) == site_label$Site[i]] = Nf_counts.site_sum[names(Nf_counts.site_sum) == site_label$Site[i]] + Nf_counts[names(Nf_counts) == site_label$metadata.label[i]]
    
    Nd_counts.site_sum[names(Nd_counts.site_sum) == site_label$Site[i]] = Nd_counts.site_sum[names(Nd_counts.site_sum) == site_label$Site[i]] + Nd_counts[names(Nd_counts) == site_label$metadata.label[i]]
}

#check that sums are equal
#sum(Nf_counts.site_sum) - sum(Nf_counts)
#sum(Nd_counts.site_sum) - sum(Nd_counts)

Nf_counts = Nf_counts.site_sum
Nd_counts = Nd_counts.site_sum

Nf_v_Nd = full_join(
data.frame(Nf = Nf_counts, Site = names(Nf_counts)),
data.frame(Nd = Nd_counts, Site = names(Nd_counts)),
by = "Site")

Nf_v_Nd.long = gather(Nf_v_Nd, "spp", "Neonectria_count.asv", -Site)

##############
#Add metadata#

####################################################################################
#First transform full_metadata to contain site level values and sample as Site#

full_metadata.site = full_metadata %>% filter(bench.control == "n" & locus == "ITS2" & !metadata.label %in% not_in_tab)
full_metadata.site$sample = as.character(full_metadata.site$Site)

#add number of plugs sampled per tree
plugs_per_sample = full_metadata.site %>% group_by(Site) %>% select(Site, Plug) %>% summarize(n_plugs = n())
full_metadata.site = full_metadata.site %>% select(Site, starts_with("Site_mean"), tmax, tmin, MAT, ppt, elev_m, lat, lon, duration_infection)
full_metadata.site = left_join(full_metadata.site %>% distinct, plugs_per_sample, by = "Site")

#then add metadata to dfs
Nf_v_Nd.long.metadata = left_join(Nf_v_Nd.long, full_metadata.site, by = "Site")

##########################################################
#Read in LULU filtered (93% minimum_simlarity) ASV table#

asv_tab = read.table("LULU/asv_tab.LULU_93.txt", header = T)
#remove negatives
asv_tab = asv_tab[,!colnames(asv_tab) %in% colnames(asv_tab.negatives)]
#remove zero sum asvs (though there actually are none)
asv_tab = asv_tab[rowSums(asv_tab) > 1,]

########################
#SUM ASV counts by site#
########################

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

asv_tab = asv_tab.site_sum

#filter out taxa that only occur in one site
asv_tab = asv_tab[rowSums(asv_tab > 0) > 1,]
#Also filter asv_tax based on retained ASVs
asv_tax = asv_tax[rownames(asv_tax) %in% rownames(asv_tab),]

############################################
#replace original matdata w/ transformed df#

full_metadata = full_metadata.site
total_seqs = data.frame(Site = colnames(asv_tab), total_seqs = colSums(asv_tab))
full_metadata = full_join(full_metadata, total_seqs, by = "Site")

