require(tidyverse)

asv_tax = read.table("ASVs_taxonomy.tsv", header = T)
asv_tab = read.table("ASVs_counts.tsv", sep = "\t", header = T) #need to fix this import bc column headers start with numbers
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
track.long = read.table("read_processing_tracking.tsv", header = T)
get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}-[ATCG]{8}", perl = T)[[1]][1]
sample.ids <- unname(sapply(colnames(asv_tab), get.sample.name))
colnames(asv_tab) = sample.ids
asv_tab = data.frame(asv_tab)
write.table(asv_tab, "ASVs_counts.new_headers.tsv", sep = "\t", quote = F)


id_bench_map = read.table("sample_labeling_bench_data.txt", header = T)
id_bench_map = filter(id_bench_map, !is.na(metadata.label))

metadata_map = read.table("metadata.txt", header = T)
#table manipulations

#this maintains ordering for microDecon, blank samples
metadata_ordered = full_join(metadata_map, id_bench_map)

track.long$sample = paste0("X", track.long$sample)
track.long.bench = full_join(track.long, id_bench_map, by = "sample")
track.long.bench.meta = full_join(track.long.bench, metadata_map, by = "metadata.label")

plot((count) ~ (QUBIT), filter(track.long.bench, QUBIT >= 0.4 & count >= 1000 & step == "input"))
plot((count) ~ QUBIT, filter(track.long.bench, QUBIT < 0.4 & count >= 1000 & step == "input"))
plot((count) ~ (QUBIT), filter(track.long.bench, QUBIT < 0.4 & count < 1000 & step == "input"))
plot((count) ~ QUBIT, filter(track.long.bench, QUBIT >= 1 & count < 500 & step == "input"))



pdf("sequence_counts_vs_dna_concentration.pdf", width = 8, height = 6)
ggplot(filter(track.long.bench, step == "input" & !is.na(metadata.label)), aes(QUBIT, count)) +
geom_point() +
scale_y_sqrt() +
scale_x_sqrt() +
facet_wrap(~plate.seq, ncol = 1) +
my_gg_theme
dev.off()


#
get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}-[ATCG]{8}", perl = T)[[1]][1]
sample.ids <- unname(sapply(colnames(asv_tab), get.sample.name))
colnames(asv_tab) = sample.ids


Nf_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Species == "s__faginata")$asv_name
Nd_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Species == "s__ditissima")$asv_name

Nf_counts = subset(asv_tab, rownames(asv_tab) %in% Nf_asvs) %>% colSums
Nd_counts = subset(asv_tab, rownames(asv_tab) %in% Nd_asvs) %>% colSums
colnames(asv_tab)[Nf_counts >= 1 & Nd_counts >= 1]

Nf_v_Nd = full_join(
    data.frame(Nf = Nf_counts, sample = names(Nf_counts)),
    data.frame(Nd = Nd_counts, sample = names(Nd_counts)),
    by = "sample")

Nf_v_Nd.long = gather(Nf_v_Nd, "spp", "count", -sample)


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

Nf_v_Nd.long.metadata = left_join(Nf_v_Nd.long, id_bench_map, by = "sample") %>%
left_join(., metadata_map, by = "metadata.label")

Nf_v_Nd.long.metadata = left_join(Nf_v_Nd.long.metadata,
    data.frame(sample = track.long %>% filter(step == "nonchim") %>% select(sample),
        total_seqs = (track.long %>% filter(step == "nonchim"))$count),
    by = "sample"
)

Nf_v_Nd.bin.metadata = left_join(Nf_v_Nd.bin, id_bench_map, by = "sample") %>%
left_join(., metadata_map, by = "metadata.label")

Nf_v_Nd.bin.metadata = left_join(Nf_v_Nd.bin.metadata,
data.frame(sample = track.long %>% filter(step == "nonchim") %>% select(sample),
total_seqs = (track.long %>% filter(step == "nonchim"))$count),
by = "sample"
)


#Neo occurence across sites, trees, plugs
pdf("Neo_occurence_min_1k_seqs_per_sample.pdf", width = 10)
ggplot(Nf_v_Nd.bin.metadata %>% filter(bench.control == "n" & seq.rep == "n" & total_seqs > 1000), aes(as.factor(Tree), ..count.., fill = occurence)) +
geom_histogram(stat = "count", width = 0.9, color = "black") +
facet_wrap(~Site, scales = "free_x", ncol = 5) +
scale_fill_manual(values = rev(cbPalette[1:4])) +
labs(x = "Trees", y = "count (plugs)") +
my_gg_theme +
theme(axis.text.x = element_blank())
dev.off()

ggplot(Nf_v_Nd.long.metadata %>% filter(bench.control == "n" & seq.rep == "n" & total_seqs > 1000), aes(as.factor(interaction(Tree, Plug)), count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~Site, scales = "free_x", ncol = 5) +
scale_fill_manual(values = rev(cbPalette[1:4])) +
labs(x = "Trees", y = "count (plugs)") +
my_gg_theme +
theme(axis.text.x = element_blank())

#ggplot(Nf_v_Nd.bin.metadata %>% filter(bench.control == "n" & resequenced == "y"), aes(sample, ..count.., fill = occurence)) +
#geom_histogram(stat = "count", width = 0.9, color = "black") +
#facet_wrap(~metadata.label, scales = "free_x", ncol = 5) +
#scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2])) +
#labs(x = "Trees", y = "count (plugs)") +
#my_gg_theme +
#theme(axis.text.x = element_blank())

#Sample resequencing
ggplot(Nf_v_Nd.long.metadata %>% filter(bench.control == "n" & resequenced == "y" & total_seqs > 1000), aes(sample, count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
#scale_y_log10() +
facet_wrap(~metadata.label, ncol = 4) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2])) +
labs(x = "Plate", y = "total sequences in sample") +
my_gg_theme


#ggplot(Nf_v_Nd.bin.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(bench.control, ..count.., fill = occurence)) +
#geom_histogram(stat = "count", width = 0.9, color = "black") +
##facet_wrap(~Site, scales = "free_x", ncol = 5) +
#scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
#labs(y = "count (plugs)") +
#my_gg_theme +
#theme(
#axis.text.x = element_text(angle = 35, hjust = 1))

#ggplot(Nf_v_Nd.bin.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(bench.control, ..count.., fill = occurence)) +
#geom_histogram(stat = "count", width = 0.9, color = "black") +
##facet_wrap(~Site, scales = "free_x", ncol = 5) +
#scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
#labs(y = "count (plugs)") +
#my_gg_theme +
#theme(
#axis.text.x = element_text(angle = 35, hjust = 1))


####
#Controls by seq depth
pdf("negative_controls.pdf", width = 9, height = 4)
ggplot(Nf_v_Nd.bin.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(sample, total_seqs, fill = occurence)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~bench.control, scales = "free_x", ncol = 3) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "total sequences in sample") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))+
theme(
axis.text.x = element_blank())
dev.off()


p1 = ggplot(Nf_v_Nd.long.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(sample, count, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~bench.control, scales = "free_x", ncol = 3) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "Neonectria counts") +
scale_y_log10() +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))+
theme(
axis.text.x = element_blank())

p2 = ggplot(Nf_v_Nd.long.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(sample, count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~bench.control, scales = "free_x", ncol = 3) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "proportion Neonectria") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))+
theme(
axis.text.x = element_blank())

pdf("bench_controls.pdf", width = 18, height = 8)
grid.arrange(p1, p2, ncol = 2)
dev.off()


#Control resequencing
ggplot(Nf_v_Nd.bin.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(plate.seq, total_seqs, fill = occurence)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_y_log10() +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "total sequences in sample") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))

p1 = ggplot(Nf_v_Nd.long.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(sample, count, fill = spp)) +
geom_col(width = 0.9, color = "black") +
scale_y_log10() +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "Neonectria counts") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 45, hjust = 1))

p2 = ggplot(Nf_v_Nd.long.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(sample, count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "proportion Neonectria") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 45, hjust = 1))

require(gridExtra)

pdf("bench_controls_resequencing.pdf", width = 15, height = 10)
grid.arrange(p1, p2, ncol = 1)
dev.off()


ggplot(Nf_v_Nd.long.metadata %>% filter(resequenced == "y"), aes(sample, count, fill = spp)) +
geom_col(width = 0.9, color = "black") +
scale_y_log10() +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "Neonectria counts") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(Nf_v_Nd.long.metadata %>% filter(resequenced == "y"), aes(interaction(plate.seq, sample), count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "proportion Neonectria") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 45, hjust = 1))

##################################
#COntaminant exploration

contam_ASVs = rownames(asv_tax[data.frame(asv_tab)$X232 > 1,])
colnames(asv_tab)[contam_ASVs >= 1]

asv_tab = read.table("ASVs_counts.new_headers.tsv", headers = T)
for(i in 1:length(contam_ASVs)){
    print(contam_ASVs[i])
    print(colnames(asv_tab)[as.matrix(asv_tab[contam_ASVs[i],]) > 0 ] )
}


asv_tab_ordered = asv_tab %>% select(one_of(as.character(metadata_ordered$sample)))
asv_tab_ordered = data.frame(OTUID = rownames(asv_tab_ordered), asv_tab_ordered)
#for true negatives first 5 columns are "blanks" and first 9 include planks plus plug negatives

#test microDecon on a site
asv_tab_ordered.ADN1 = data.frame(
    asv_tab_ordered[,c(1,6)],
    asv_tab_ordered %>% select(one_of(as.character((filter(metadata_ordered, Site == "ADN1" & bench.control == "n"))$sample)))
)
asv_tab_ordered.ADS1 = data.frame(
OTUID = asv_tab_ordered[,1], BlankMean = rowSums(asv_tab_ordered[,c(2,6)])/2,
asv_tab_ordered %>% select(one_of(as.character((filter(metadata_ordered, Site == "ADS1" & bench.control == "n"))$sample)))
)

deconatimanted.ADN1 = decon(data = asv_tab_ordered.ADN1, numb.blanks = 1, numb.ind = length(colnames(asv_tab_ordered.ADN1))-2, taxa = F)

#try on full set of samples

#INTIALIZE VARIABLES
asv_tab_ordered.decon = data.frame(
OTUID = asv_tab_ordered[,1], BlankMean = rowSums(asv_tab_ordered[,c(2,6)])/2
)
site_block_lens = vector("numeric", length = length(levels(as.factor(as.character((filter(metadata_ordered, bench.control != "extraction.negative"))$Site)))))
sites_list = levels(as.factor(as.character((filter(metadata_ordered, bench.control != "extraction.negative"))$Site)))

#Build table
for(i in 1:length(site_block_lens)){
 
    #filter search list using join
    samples = as.character((filter(metadata_ordered, Site == sites_list[i] & bench.control != "extraction.negative"))$sample)
    samples = inner_join(data.frame(sample = samples), data.frame(sample = colnames(asv_tab_ordered)))

    asv_tab_ordered.decon = data.frame(asv_tab_ordered.decon,
    asv_tab_ordered %>% select(one_of(samples$sample))
    )
    site_block_lens[i] = length(samples$sample)
    print(site_block_lens[i])
    print(c("number of cols ", ncol(asv_tab_ordered.decon)))
}

sum(site_block_lens)
ncol(asv_tab_ordered.decon)

deconatimanted.table = decon(data = asv_tab_ordered.decon, numb.blanks = 1, numb.ind = site_block_lens, taxa = F, thresh = 1, prop.thresh = 0)

decon.diff(asv_tab_ordered.decon, deconatimanted.table, numb.blanks = 1, numb.ind = site_block_lens)





decon_tab = deconatimanted.table$decon.table[,2:length(deconatimanted.table$decon.table)]
rownames(decon_tab) = deconatimanted.table$decon.table[,1]

Nf_counts.decon = subset(decon_tab, rownames(decon_tab) %in% Nf_asvs) %>% colSums
Nd_counts.decon = subset(decon_tab, rownames(decon_tab) %in% Nd_asvs) %>% colSums
colnames(decon_tab)[Nf_counts.decon >= 1 & Nd_counts.decon >= 1]

Nf_v_Nd.decon = full_join(
data.frame(Nf = Nf_counts.decon, sample = names(Nf_counts.decon)),
data.frame(Nd = Nd_counts.decon, sample = names(Nd_counts.decon)),
by = "sample")

Nf_v_Nd.decon.long = gather(Nf_v_Nd.decon, "spp", "count", -sample)


Nf_v_Nd.decon.bin = full_join(
data.frame(Nf = as.numeric(as.matrix(Nf_counts.decon) > 0), sample = names(Nf_counts.decon)),
data.frame(Nd = as.numeric(as.matrix(Nd_counts.decon) > 0), sample = names(Nd_counts.decon)),
by = "sample")

for(i in 1:length(Nf_v_Nd.decon.bin$sample)){
    if(Nf_v_Nd.decon.bin$Nf[i] == 1 & Nf_v_Nd.decon.bin$Nd[i] == 1){
        Nf_v_Nd.decon.bin$occurence[i] = "both"
    }
    if(Nf_v_Nd.decon.bin$Nf[i] == 1 & Nf_v_Nd.decon.bin$Nd[i] == 0){
        Nf_v_Nd.decon.bin$occurence[i] = "Nf"
    }
    if(Nf_v_Nd.decon.bin$Nf[i] == 0 & Nf_v_Nd.decon.bin$Nd[i] == 1){
        Nf_v_Nd.decon.bin$occurence[i] = "Nd"
    }
    if(Nf_v_Nd.decon.bin$Nf[i] == 0 & Nf_v_Nd.decon.bin$Nd[i] == 0){
        Nf_v_Nd.decon.bin$occurence[i] = "none"
    }
}

Nf_v_Nd.decon.long.metadata = left_join(Nf_v_Nd.decon.long, id_bench_map, by = "sample") %>%
left_join(., metadata_map, by = "metadata.label")

Nf_v_Nd.decon.long.metadata = left_join(Nf_v_Nd.decon.long.metadata,
data.frame(sample = track.long %>% filter(step == "nonchim") %>% select(sample),
total_seqs = (track.long %>% filter(step == "nonchim"))$count),
by = "sample"
)

Nf_v_Nd.decon.bin.metadata = left_join(Nf_v_Nd.decon.bin, id_bench_map, by = "sample") %>%
left_join(., metadata_map, by = "metadata.label")

Nf_v_Nd.decon.bin.metadata = left_join(Nf_v_Nd.decon.bin.metadata,
data.frame(sample = track.long %>% filter(step == "nonchim") %>% select(sample),
total_seqs = (track.long %>% filter(step == "nonchim"))$count),
by = "sample"
)


#Neo occurence across sites, trees, plugs
pdf("decon_Neo_occurence_min_1k_seqs_per_sample.pdf", width = 10)
ggplot(Nf_v_Nd.decon.bin.metadata %>% filter(bench.control == "n" & seq.rep == "n" & total_seqs > 1000), aes(as.factor(Tree), ..count.., fill = occurence)) +
geom_histogram(stat = "count", width = 0.9, color = "black") +
facet_wrap(~Site, scales = "free_x", ncol = 5) +
scale_fill_manual(values = rev(cbPalette[1:4])) +
labs(x = "Trees", y = "count (plugs)") +
my_gg_theme +
theme(axis.text.x = element_blank())
dev.off()

ggplot(Nf_v_Nd.decon.long.metadata %>% filter(bench.control == "n" & seq.rep == "n" & total_seqs > 1000), aes(as.factor(interaction(Tree, Plug)), count, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~Site, scales = "free_x", ncol = 5) +
scale_fill_manual(values = rev(cbPalette[1:4])) +
scale_y_log10() +
labs(x = "Trees", y = "count (plugs)") +
my_gg_theme +
theme(axis.text.x = element_blank())





####
#Controls by seq depth
pdf("negative_controls.decon.pdf", width = 9, height = 4)
ggplot(Nf_v_Nd.decon.bin.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(sample, total_seqs, fill = occurence)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~bench.control, scales = "free_x", ncol = 3) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "total sequences in sample") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))+
theme(
axis.text.x = element_blank())
dev.off()


p1 = ggplot(Nf_v_Nd.decon.long.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(sample, count, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~bench.control, scales = "free_x", ncol = 3) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "Neonectria counts") +
scale_y_log10() +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))+
theme(
axis.text.x = element_blank())

p2 = ggplot(Nf_v_Nd.decon.long.metadata %>% filter(bench.control != "n" & seq.rep == "n"), aes(sample, count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~bench.control, scales = "free_x", ncol = 3) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "proportion Neonectria") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))+
theme(
axis.text.x = element_blank())

pdf("bench_controls.decon.pdf", width = 18, height = 8)
grid.arrange(p1, p2, ncol = 2)
dev.off()


#Control resequencing
ggplot(Nf_v_Nd.decon.bin.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(plate.seq, total_seqs, fill = occurence)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_y_log10() +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "total sequences in sample") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 35, hjust = 1))

p1 = ggplot(Nf_v_Nd.decon.long.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(sample, count, fill = spp)) +
geom_col(width = 0.9, color = "black") +
scale_y_log10() +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "Neonectria counts") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 45, hjust = 1))

p2 = ggplot(Nf_v_Nd.decon.long.metadata %>% filter(bench.control != "n" & resequenced == "y"), aes(sample, count/total_seqs, fill = spp)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(metadata.label~bench.control, scales = "free_x", ncol = 5) +
scale_fill_manual(values = c("none" = cbPalette[1], "both" = cbPalette[4], "Nf" = cbPalette[2], "Nd" = cbPalette[3])) +
labs(y = "proportion Neonectria") +
my_gg_theme +
theme(
axis.text.x = element_text(angle = 45, hjust = 1))

require(gridExtra)

pdf("bench_controls_resequencing.decon.pdf", width = 15, height = 10)
grid.arrange(p1, p2, ncol = 1)
dev.off()


###########################################
#Build asv table with taxnomy for plotting#
###########################################

asv_tab.long = gather(data.frame(asv = rownames(asv_tab), asv_tab), "sample", "sequences", -asv)

asv_tab.long.tax = full_join(asv_tab.long, data.frame(asv = rownames(asv_tax), asv_tax))

asv_tab.long.tax = mutate(asv_tab.long.tax, fill_col = ifelse(is.na(Species), "white", ifelse(
    Species == "s__ditissima", cbPalette[3], ifelse(
    Species == "s__faginata", cbPalette[2], "white"
    )
)))

asv_cols = asv_tab.long.tax$fill_col
names(asv_cols) = asv_tab.long.tax$asv

asv_tab.long.tax.metadata = left_join(asv_tab.long.tax, id_bench_map, by = "sample") %>%
    left_join(., metadata_map, by = "metadata.label") %>%
    left_join(., data.frame(sample = asv_tab %>% colSums %>% names, total_seqs = asv_tab %>% colSums))

asv_tab.total_seqs = left_join(id_bench_map, data.frame(sample = asv_tab %>% colSums %>% names, total_seqs = asv_tab %>% colSums), by = "sample") %>%
left_join(., metadata_map, by = "metadata.label")

ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & resequenced == "y"),
    aes(sample, sequences, fill = asv)) +
geom_col(color = "black") +
facet_wrap(metadata.label~bench.control, scales = "free") +
scale_y_log10() +
scale_fill_manual(values = asv_cols, guide = F) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

####
#CONTROL SAMPLES COLORED BY NF AND ND
#pdf("control_samples_stacked_taxa.pdf", width = 12, height = 8)
p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n"),
aes(sample, sequences)) +
geom_col(color = "dark grey", aes(fill = asv)) +
geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs+1000, label = total_seqs), size = 5) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
scale_fill_manual(values = asv_cols, guide = F) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
#dev.off()

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n"),
aes(sample, sequences/total_seqs)) +
geom_col(color = "dark grey", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
scale_fill_manual(values = asv_cols, guide = F) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("control_samples_stacked_taxa_abs_and_rel_abd.pdf", width = 18, height = 14)
grid.arrange(p1,p2, ncol = 1)
dev.off()


p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n"),
aes(sample, sequences)) +
geom_col(color = "dark grey", aes(fill = asv)) +
geom_text(data = filter(asv_tab.total_seqs, bench.control == "n" & !is.na(total_seqs)), aes(y = total_seqs+1000, label = total_seqs), size = 5, angle = 90) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
scale_fill_manual(values = asv_cols, guide = F) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n"& !is.na(total_seqs)),
aes(sample, sequences/total_seqs)) +
geom_col(color = "dark grey", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
scale_fill_manual(values = asv_cols, guide = F) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

pdf("full_data_stacked_taxa_abs_and_rel_abd.pdf", width = 22, height = 16)
grid.arrange(p1,p2, ncol = 1)
dev.off()


#Plot only those asvs that appear in PBS control

clean.wash.asvs = asv_tab.long.tax.metadata %>% filter(bench.control == "clean.wash.negative" & sequences > 0) %>% select(asv) %>% unique


asv_tab.long.tax.metadata %>% filter(bench.control == "clean.wash.negative" & sequences > 0) %>% write.table("clean.wash.negative.taxa.txt", sep = "\t", row.names = F, quote = F)

p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & asv %in% clean.wash.asvs$asv),
aes(sample, sequences)) +
geom_col(color = "black", aes(fill = interaction(asv,Genus))) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
scale_fill_manual(values = cbPalette) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & asv %in% clean.wash.asvs$asv),
aes(sample, sequences/total_seqs)) +
geom_col(color = "black", aes(fill = interaction(asv,Genus))) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
scale_fill_manual(values = cbPalette) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("control_samples_stacked_taxa_abs_and_rel_abd_PBS_taxa_only.pdf", width = 18, height = 14)
grid.arrange(p1,p2, ncol = 1)
dev.off()

#Plot only those asvs that appear in extraction control

extraction.neg.asvs = asv_tab.long.tax.metadata %>% filter(bench.control == "extraction.negative" & sequences > 0) %>% select(asv) %>% unique

asv_tab.long.tax.metadata %>% filter(bench.control == "extraction.negative" & sequences > 0) %>% write.table("extraction.negative.taxa.txt", sep = "\t", row.names = F, quote = F)

p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & asv %in% extraction.neg.asvs$asv),
aes(sample, sequences)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & asv %in% extraction.neg.asvs$asv),
aes(sample, sequences/total_seqs)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("control_samples_stacked_taxa_abs_and_rel_abd_extraction_negative_taxa_only.pdf", width = 18, height = 14)
grid.arrange(p1,p2, ncol = 1)
dev.off()

#FOR SAMPLE DATA Plot only those asvs that appear in PBS control

clean.wash.asvs = asv_tab.long.tax.metadata %>% filter(bench.control == "clean.wash.negative" & sequences > 0) %>% select(asv) %>% unique

p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & asv %in% clean.wash.asvs$asv),
aes(sample, sequences)) +
geom_col(color = "black", aes(fill = interaction(asv,Genus))) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
scale_fill_manual(values = cbPalette) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & asv %in% clean.wash.asvs$asv),
aes(sample, sequences/total_seqs)) +
geom_col(color = "black", aes(fill = interaction(asv,Genus))) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
scale_fill_manual(values = cbPalette) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

pdf("full_data_stacked_taxa_abs_and_rel_abd_PBS_taxa_only.pdf", width = 22, height = 16)
grid.arrange(p1,p2, ncol = 1)
dev.off()

#FOR SAMPLE DATA Plot only those asvs that appear in PBS control

extraction.neg.asvs = asv_tab.long.tax.metadata %>% filter(bench.control == "extraction.negative" & sequences > 0) %>% select(asv) %>% unique

p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & asv %in% extraction.neg.asvs$asv),
aes(sample, sequences)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & asv %in% extraction.neg.asvs$asv),
aes(sample, sequences/total_seqs)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

pdf("full_data_stacked_taxa_abs_and_rel_abd_extraction_negative_taxa_only.pdf", width = 22, height = 16)
grid.arrange(p1,p2, ncol = 1)
dev.off()




#Plot only those asvs that appear in no.sample control

no.sample.asvs = asv_tab.long.tax.metadata %>% filter(bench.control == "no.sample" & sequences > 0) %>% select(asv) %>% unique
no.sample.asvs.gt10 = asv_tab.long.tax.metadata %>% filter(bench.control == "no.sample" & sequences > 10) %>% select(asv) %>% unique

asv_tab.long.tax.metadata %>% filter(bench.control == "no.sample" & sequences > 0) %>% write.table("no.sample.taxa.txt", sep = "\t", row.names = F, quote = F)
asv_tab.long.tax.metadata %>% filter(bench.control == "no.sample" & sequences > 10) %>% write.table("no.sample.taxa.gt10.txt", sep = "\t", row.names = F, quote = F)

p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & asv %in% no.sample.asvs.gt10$asv),
aes(sample, sequences)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control != "n" & asv %in% no.sample.asvs.gt10$asv),
aes(sample, sequences/total_seqs)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~bench.control, scales = "free_x", nrow = 2) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("control_samples_stacked_taxa_abs_and_rel_abd_no_sample_gt10_taxa_only.pdf", width = 18, height = 14)
grid.arrange(p1,p2, ncol = 1)
dev.off()

#SAMPLE DATA NO SAMPLE
p1 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & asv %in% no.sample.asvs.gt10$asv),
aes(sample, sequences)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

p2 = ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & asv %in% no.sample.asvs.gt10$asv),
aes(sample, sequences/total_seqs)) +
geom_col(color = "black", aes(fill = asv)) +
#geom_text(data = filter(asv_tab.total_seqs, bench.control != "n"), aes(y = total_seqs, label = total_seqs)) +
facet_wrap(~plate.seq, scales = "free_x", ncol = 1) +
#scale_y_log10() +
#scale_fill_manual(values = cbPalette) +
my_gg_theme +
labs(y = "Proportion sequences") +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank())

pdf("full_data_stacked_taxa_abs_and_rel_abd_no_sample_gt10_taxa_only.pdf", width = 22, height = 16)
grid.arrange(p1,p2, ncol = 1)
dev.off()

###########################################
pdf("control_samples_conc_vs_seqs.pdf", width = 8, height = 4)
ggplot(track.long.bench.meta %>% filter(bench.control != "n" & step == "input"), aes(QUBIT, count, color = bench.control)) +
geom_point() +
scale_color_manual(values = cbPalette) +
my_gg_theme
dev.off()

ggplot(asv_tab.long.tax.metadata %>% filter(bench.control == "n" & seq.rep == "n" & total_seqs >= 1000), aes(as.factor(interaction(Tree, Plug)), count/total_seqs, fill = asv)) +
geom_col(width = 0.9, color = "black") +
facet_wrap(~Site, scales = "free_x") +
scale_fill_manual(values = asv_cols, guide = F) +
#scale_y_log10() +
labs(x = "Trees", y = "count") +
my_gg_theme +
theme(axis.text.x = element_blank())

#Neonectria frequency by site

sites_nmds = asv_tab.gt5K.rare.mds.metadata.survey_dat.neo %>% filter(bench.control == "n") %>% select(Site) %>% unique

Nf_Nd_site_freq = data.frame(Site = vector(mode = "character", length = length(sites_nmds$Site)),
Nf = vector(mode = "numeric", length = length(sites_nmds$Site)),
Nd = vector(mode = "numeric", length = length(sites_nmds$Site)),
total = vector(mode = "numeric", length = length(sites_nmds$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(sites_nmds$Site)){
    temp_tab = asv_tab.gt5K.rare.mds.metadata.survey_dat.neo %>% filter(Site == sites_nmds$Site[i])
    Nf_Nd_site_freq$Nf[i] = (filter(temp_tab, Nf > 0)) %>% nrow
    Nf_Nd_site_freq$Nd[i] = (filter(temp_tab, Nd > 0)) %>% nrow
    Nf_Nd_site_freq$total[i] = nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(sites_nmds$Site[i])
}

Nf_Nd_site_freq.site_info = left_join(Nf_Nd_site_freq, site_info)


p1 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nd/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. ditissima frequency\n(proportion plugs)", x = "Latitude") +
scale_color_manual(values = cbPalette, guide = F) +
my_gg_theme

p2 = ggplot(Nf_Nd_site_freq.site_info , aes(lat, Nf/total, group = lat, color = state_prov)) +
geom_point(size = 3) +
labs(y = "N. faginata frequency\n(proportion plugs)", x = "Latitude") +
scale_color_manual(values = cbPalette) +
my_gg_theme

pdf("Neonectria_frequency_by_lat.pdf", width = 10, height = 4)
grid.arrange(p1,p2,ncol= 2, widths = c(0.45,0.55))
dev.off()




sample_richness = data.frame(sample = names(apply(asv_tab.gt5K.rare,1,function(x) sum(x > 0))), richness =  apply(asv_tab.gt5K.rare,1,function(x) sum(x > 0)))
sample_richness.metadata = left_join(sample_richness, id_bench_map) %>%
left_join(., metadata_map)


sample_richness.metadata.site_info = left_join(sample_richness.metadata, Nf_v_Nd.bin) %>% left_join(., site_info, by = "Site")

pdf("ASV_richness_by_lat.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(lat, richness, group = lat, color = state_prov)) +
geom_boxplot(width = .25) +
labs(y = "ASV richness\nper 5K sequences", x = "Latitude") +
scale_color_manual(values = cbPalette) +
my_gg_theme
dev.off()

pdf("ASV_richness_by_neo_occurence.pdf", width = 8, height = 4)
ggplot(sample_richness.metadata.site_info , aes(occurence, richness)) +
geom_boxplot(width = .25) +
labs(y = "ASV richness\nper 5K sequences", x = "Neonectria occurence") +
scale_color_manual(values = cbPalette, guide = F) +
scale_x_discrete(labels = c("both", "N. ditissima", "N. faginata", "none")) +
my_gg_theme
dev.off()


