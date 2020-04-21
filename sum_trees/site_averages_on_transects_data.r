require(tidyverse)

transect_dat = read.table("sample_data/BBD_survey_transect_data.txt", header = T)
transect_dat = filter(transect_dat, TreeSpecies == "Fg")
transect_dat$Tree = as.factor(transect_dat$Tree)
transect_dat.site_mean = transect_dat %>% group_by(Site) %>% summarize_if(is.numeric, mean)

colnames(transect_dat.site_mean)[2:length(colnames(transect_dat.site_mean))] = paste("Site_mean", colnames(transect_dat.site_mean)[2:length(colnames(transect_dat.site_mean))], sep = ".")

write.table(transect_dat.site_mean, file = "sample_data/BBD_survey_transect_data.site_mean.txt", quote = F, sep = "\t", row.names = F)
