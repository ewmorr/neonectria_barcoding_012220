
source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")
write.csv(asv_tab, "plug_and_tree_level_asv_and_metadata_tables/asv_tab.plug.csv", quote = F)
write.csv(Nf_v_Nd.bin.metadata, "plug_and_tree_level_asv_and_metadata_tables/metadata.plug.csv")
write.csv(asv_tax, "plug_and_tree_level_asv_and_metadata_tables/asv_tax.csv", quote = F)

source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
write.csv(asv_tab, "plug_and_tree_level_asv_and_metadata_tables/asv_tab.tree.csv", quote = F)
write.csv(Nf_v_Nd.bin.metadata, "plug_and_tree_level_asv_and_metadata_tables/metadata.tree.csv")

