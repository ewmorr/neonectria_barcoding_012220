require(tidyverse)
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

#Some gymnastics to get sample labels in the correct format. Could fix this upstream...
colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
track.long$sample <- unname(sapply(as.character(track.long$sample), get.sample.name))
track.long$sample = paste0("X", track.long$sample)

#table joins
metadata_ordered = full_join(metadata_map, id_bench_map)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>% left_join(., site_info, by = "Site")

full_metadata = full_join(metadata_ordered, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))

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

Nf_v_Nd.long = gather(Nf_v_Nd, "spp", "Neonectria_count", -sample)

###########################
#Nf and Nd occurence (0/1)#

Nf_v_Nd.bin = full_join(
data.frame(Nf = as.numeric(as.matrix(Nf_counts) > 0), sample = names(Nf_counts)),
data.frame(Nd = as.numeric(as.matrix(Nd_counts) > 0), sample = names(Nd_counts)),
by = "sample")

for(i in 1:length(Nf_v_Nd.bin$sample)){
    if(Nf_v_Nd.bin$Nf[i] == 1 & Nf_v_Nd.bin$Nd[i] == 1){
        Nf_v_Nd.bin$occurence[i] = "both"
    }
    if(Nf_v_Nd.bin$Nf[i] == 1 & Nf_v_Nd.bin$Nd[i] == 0){
        Nf_v_Nd.bin$occurence[i] = "Nf"
    }
    if(Nf_v_Nd.bin$Nf[i] == 0 & Nf_v_Nd.bin$Nd[i] == 1){
        Nf_v_Nd.bin$occurence[i] = "Nd"
    }
    if(Nf_v_Nd.bin$Nf[i] == 0 & Nf_v_Nd.bin$Nd[i] == 0){
        Nf_v_Nd.bin$occurence[i] = "none"
    }
}

##############
#Add metadata#

Nf_v_Nd.long.metadata = left_join(Nf_v_Nd.long, full_metadata, by = "sample")

Nf_v_Nd.long.metadata = left_join(Nf_v_Nd.long.metadata,
data.frame(sample = track.long %>% filter(step == "nonchim") %>% select(sample),
total_seqs = (track.long %>% filter(step == "nonchim"))$count),
by = "sample"
)

Nf_v_Nd.bin.metadata = left_join(Nf_v_Nd.bin, full_metadata, by = "sample")

Nf_v_Nd.bin.metadata = left_join(Nf_v_Nd.bin.metadata,
data.frame(sample = track.long %>% filter(step == "nonchim") %>% select(sample),
total_seqs = (track.long %>% filter(step == "nonchim"))$count),
by = "sample"
)
