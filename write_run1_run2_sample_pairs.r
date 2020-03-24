require(tidyverse)
get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

#Read metadata
id_bench_map = read.table("sample_data/sample_mapping.txt", header = T)
metadata_map = read.table("sample_data/metadata.txt", header = T)
survey_dat = read.table("sample_data/trees_site_survey_data.txt", header = T, sep = "\t")
neo_cov = read.table("sample_data/plug_neonectria_coverage.txt", header = T)
site_info = read.csv("sample_data/site_info.csv")

colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
metadata_ordered = full_join(metadata_map, id_bench_map)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>% left_join(., site_info, by = "Site")

full_metadata = full_join(metadata_ordered, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))

sample_map.wide_run = full_metadata %>%
    filter(!is.na(metadata.label) & seq.rep == "n" & locus == "ITS2") %>%
    select(c(sample, run.seq, metadata.label))  %>%
    pivot_wider(., id_cols = metadata.label, names_from = run.seq, values_from = sample) %>%
    drop_na

write.table(sample_map.wide_run, file = "sample_data/run1_run2_sample_pairs.txt", sep = "\t", col.names = F, quote = F, row.names = F)
