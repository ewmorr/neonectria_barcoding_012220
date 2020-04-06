require(lulu)

asv_tab = read.table("dada2_out/ASVs_counts.tsv", header = T)
asv_matches = read.table("LULU/match_list.0.84min.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE) #note the flags here. algorithm does not work if these are not set (fails a couple of ways)

asv_tax = read.table("dada2_out/ASVs_taxonomy.tsv", header = T)
asv_tax.char = apply(asv_tax, 2, as.character)


lulu.asv_0.84 = lulu(
    otutable = asv_tab,
    matchlist = asv_matches,
    minimum_ratio_type = "min", #default = 'min', Min is recommended value. More conservative than `avg` in terms of designating an ASV as an error
    minimum_ratio = 1, #default = 1, i.e. the parent can never be less abundant than the error
    minimum_match = 84, #default = 84, consider increasing this above the default (e.g. 95%)
    minimum_relative_cooccurence = 0.95 #default = 0.95
)


lulu.asv_0.84 = lulu(
otutable = asv_tab,
matchlist = asv_matches,
minimum_match = 84
)

summary(lulu.asv_0.84)

#get all parents
parents = (lulu.asv_0.84$otu_map %>% filter(curated == "merged"))$parent_id %>% unique
#get all merged stats. need to do both manipulations to retain rownames
merged_asvs = data.frame(child = rownames(lulu.asv_0.84$otu_map),lulu.asv_0.84$otu_map)  %>% filter(parent_id %in% parents & curated != "parent")

merged_asvs_taxonomy_and_sim = data.frame(
    child_tax = vector(mode = "character", length = length(merged_asvs$child)),
    parent_tax = vector(mode = "character", length = length(merged_asvs$child)),
    seq_sim = vector(mode = "numeric", length = length(merged_asvs$child)),
    n_samples_child = vector(mode = "numeric", length = length(merged_asvs$child)),
    child_name = vector(mode = "character", length = length(merged_asvs$child)),
    parent_name = vector(mode = "character", length = length(merged_asvs$child)),
    stringsAsFactors=FALSE
)

for(i in 1:length(merged_asvs$child)){
    
    child_tax = as.matrix(data.frame(ASV = as.character(rownames(asv_tax)), asv_tax.char) %>% filter(ASV == merged_asvs$child[i]))
    parent_tax = as.matrix(data.frame(ASV = as.character(rownames(asv_tax)), asv_tax.char) %>% filter(ASV == merged_asvs$parent_id[i]))
    seq_sim = (asv_matches %>% filter(V1 == merged_asvs$child[i] & V2 == merged_asvs$parent_id[i]))$V3
 
    merged_asvs_taxonomy_and_sim[i,1] = paste0(child_tax[3:8], collapse = ";")
    merged_asvs_taxonomy_and_sim[i,2] = paste0(parent_tax[3:8], collapse = ";")
    merged_asvs_taxonomy_and_sim[i,3] = seq_sim
    merged_asvs_taxonomy_and_sim[i,4] = merged_asvs$spread[i]
    merged_asvs_taxonomy_and_sim[i,5] = merged_asvs$child[i] %>% as.character
    merged_asvs_taxonomy_and_sim[i,6] = merged_asvs$parent_id[i]
}

write.table(merged_asvs_taxonomy_and_sim, file = "LULU/merged_asvs_84.txt", row.names = F, sep = "\t", quote = F)
write.table(lulu.asv_0.84$curated_table, file = "LULU/asv_tab.LULU_84.txt", sep = "\t", quote = F)

################
#95 percent sim#
################

lulu.asv_0.95 = lulu(
otutable = asv_tab,
matchlist = asv_matches,
minimum_match = 95
)

summary(lulu.asv_0.95)


#get all parents
parents = (lulu.asv_0.95$otu_map %>% filter(curated == "merged"))$parent_id %>% unique
#get all merged stats. need to do both manipulations to retain rownames
merged_asvs = data.frame(child = rownames(lulu.asv_0.95$otu_map),lulu.asv_0.95$otu_map)  %>% filter(parent_id %in% parents & curated != "parent")

merged_asvs_taxonomy_and_sim = data.frame(
child_tax = vector(mode = "character", length = length(merged_asvs$child)),
parent_tax = vector(mode = "character", length = length(merged_asvs$child)),
seq_sim = vector(mode = "numeric", length = length(merged_asvs$child)),
n_samples_child = vector(mode = "numeric", length = length(merged_asvs$child)),
child_name = vector(mode = "character", length = length(merged_asvs$child)),
parent_name = vector(mode = "character", length = length(merged_asvs$child)),
stringsAsFactors=FALSE
)


for(i in 1:length(merged_asvs$child)){
    
    child_tax = as.matrix(data.frame(ASV = as.character(rownames(asv_tax)), asv_tax.char) %>% filter(ASV == merged_asvs$child[i]))
    parent_tax = as.matrix(data.frame(ASV = as.character(rownames(asv_tax)), asv_tax.char) %>% filter(ASV == merged_asvs$parent_id[i]))
    seq_sim = (asv_matches %>% filter(V1 == merged_asvs$child[i] & V2 == merged_asvs$parent_id[i]))$V3
    
    merged_asvs_taxonomy_and_sim[i,1] = paste0(child_tax[3:8], collapse = ";")
    merged_asvs_taxonomy_and_sim[i,2] = paste0(parent_tax[3:8], collapse = ";")
    merged_asvs_taxonomy_and_sim[i,3] = seq_sim
    merged_asvs_taxonomy_and_sim[i,4] = merged_asvs$spread[i]
    merged_asvs_taxonomy_and_sim[i,5] = merged_asvs$child[i] %>% as.character
    merged_asvs_taxonomy_and_sim[i,6] = merged_asvs$parent_id[i]
}

write.table(merged_asvs_taxonomy_and_sim, file = "LULU/merged_asvs_95.txt", row.names = F, sep = "\t", quote = F)
write.table(lulu.asv_0.95$curated_table, file = "LULU/asv_tab.LULU_95.txt", sep = "\t", quote = F)


#ASV count by taxon for dif filtering methods

####################
#ASV count by taxon#
ASV_count_by_taxon = data.frame(ASV = row.names(asv_tax), asv_tax) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(ASV_count = n())

ASVs_LULU_95 = left_join(
    data.frame(ASV = lulu.asv_0.95$curated_otus),
    data.frame(ASV = row.names(asv_tax), asv_tax),
    by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(LULU_95 = n())

ASVs_LULU_84 = left_join(
data.frame(ASV = lulu.asv_0.84$curated_otus),
data.frame(ASV = row.names(asv_tax), asv_tax),
by = "ASV"
) %>%
group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
summarize(LULU_84 = n())

ASV_OTU_count_by_taxon = full_join(ASV_count_by_taxon, ASVs_LULU_95) %>%
    full_join(., ASVs_LULU_84)

write.table(ASV_OTU_count_by_taxon, file = "LULU/ASV_count_by_taxon.txt", row.names = F, quote = F, sep = "\t")














