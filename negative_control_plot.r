require(tidyverse)
require(data.table)
require(RColorBrewer)
get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

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
site_climate.GDD = read.table(file = "sample_data/site_info.GDD.txt", header = T)
#Some gymnastics to get sample labels in the correct format. Could fix this upstream...
colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
track.long$sample <- unname(sapply(as.character(track.long$sample), get.sample.name))
#track.long$sample = paste0("X", track.long$sample)

#table joins
metadata_ordered = full_join(metadata_map, id_bench_map)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>%
left_join(., site_info, by = "Site") %>%
left_join(., site_means, by = "Site") %>%
left_join(., site_climate, by = c("Site","lat","lon")) %>%
left_join(., site_climate.GDD, by = "Site")

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
#asv_tab.negatives = semi_join(
#data.frame(sample = rownames(t(asv_tab)), t(asv_tab)),
#full_metadata %>% filter(bench.control != "n")
#)
asv_tab.negatives = left_join(
full_metadata %>% filter(bench.control != "n") %>% select(sample),
data.frame(sample = rownames(t(asv_tab)), t(asv_tab))
)

asv_tab.negatives = asv_tab.negatives %>% filter(!is.na(sample))

rownames(asv_tab.negatives) = asv_tab.negatives$sample
asv_tab.negatives$sample = NULL
asv_tab.negatives = t(asv_tab.negatives)
asv_tab.negatives %>% head
asv_tab.negatives[is.na(asv_tab.negatives)] = 0
asv_tab.negatives = asv_tab.negatives[rowSums(asv_tab.negatives) > 0,]

total_sample_seq_counts = data.frame(sample = colnames(asv_tab.negatives), count = colSums(asv_tab.negatives))
#join with taxonomy

asv_tab.negatives.asvnames = data.frame(ASV = rownames(asv_tab.negatives), asv_tab.negatives)




#################################################
#Process asv_tax for lowest informative taxonomy#

asv_tax.char = apply(asv_tax, 2, as.character)
asv_tax.char[is.na(asv_tax.char)] = "unknown"
rownames(asv_tax.char) = rownames(asv_tax)

get_informative_taxa = function(x){
    found_info = 0
    for(i in length(x):2){
        if(x[i] != "unknown"){
            if(i == length(x)){
                return(paste(as.character(x[i-1]), sub("s__", "", as.character(x[i]) ), sep = " " ))
            }
            else{
                return(as.character(x[i]))
            }
            found_info = 1
            break
        }
    }
    if(found_info == 0){return(x[1])}
}

asv_informative_taxa = vector(mode = "character", length = length(rownames(asv_tax.char)))
asv_informative_taxa = apply(asv_tax.char, 1, function(x) get_informative_taxa(x))

#################
#Join with taxonomy and metadata

asv_tab.negatives.asvnames.taxa = left_join(
    asv_tab.negatives.asvnames,
    data.frame(ASV = names(asv_informative_taxa), taxonomy = asv_informative_taxa)

)

asv_tab.negatives.asvnames.taxa.sum = asv_tab.negatives.asvnames.taxa %>% group_by(taxonomy) %>% summarize(across(where(is.numeric), ~sum(.x)))

asv_tab.negatives.asvnames.taxa.sum.long = asv_tab.negatives.asvnames.taxa.sum %>% pivot_longer(cols = -taxonomy, names_to = "sample", values_to = "count")

plot_labels = read.csv("negatives_plot_labels.csv")

asv_tab.negatives.asvnames.taxa.sum.long.metadata = left_join(
asv_tab.negatives.asvnames.taxa.sum.long,
full_metadata %>% select(sample, bench.control)
) %>% left_join(.,
plot_labels
)




#To control y limits using https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel

facet_bounds <- read.table(header=TRUE,
text=
"bench.control ymin ymax
1     0     5000
2     0     10",
stringsAsFactors=FALSE)

ff <- with(facet_bounds,
data.frame(count=c(ymin,ymax),
bench.control=c(bench.control,bench.control)))

ggplot(asv_tab.negatives.asvnames.taxa.sum.long.metadata, aes(x = sample, y = count, fill = taxonomy)) +
geom_bar(stat = "identity") +
facet_wrap(~bench.control, scales = "free", ncol = 4) +
scale_fill_brewer(palette = "Paired") +
geom_point(data = ff, x = NA, fill = NA) +
#annotate(geom = "text", x = total_sample_seq_counts$sample, y = total_sample_seq_counts$count+2, label = total_sample_seq_counts$count) +
#my_gg_theme +
theme_bw() +
theme(
axis.text.x = element_blank()
)


facet_bounds <- read.table(header=TRUE,
text=
"panel ymin ymax
1     0     5000
2     0     15",
stringsAsFactors=FALSE)

ff <- with(facet_bounds,
data.frame(count=c(ymin,ymax),
panel=c(panel,panel)))

asv_tab.negatives.asvnames.taxa.sum.long.metadata$plot_label = factor(
    asv_tab.negatives.asvnames.taxa.sum.long.metadata$plot_label,
        levels = c(
            "N. ditissima positive 1",
            "N. ditissima positive 2",
            "N. faginata positive 1",
            "N. faginata positive 2",
            "plug negative 1",
            "plug negative 2",
            "plug negative 3",
            "clean wash negative",
            "used wash 1",
            "used wash 2",
            "used wash 3",
            "used wash 4",
            "used wash 5",
            "used wash 6",
            "used wash 7",
            "used wash 8",
            "used wash 9",
            "DNA extraction negative",
            "PCR negative 1",
            "PCR negative 2",
            "PCR negative 3"
    )
)

p1 = ggplot(asv_tab.negatives.asvnames.taxa.sum.long.metadata, aes(x = plot_label, y = count, fill = taxonomy)) +
geom_bar(stat = "identity") +
facet_wrap(~panel, scales = "free", ncol = 2) +
scale_fill_brewer(palette = "Paired") +
geom_point(data = ff, x = NA, fill = NA) +
labs(y = "Sequence count") +
theme_bw() +
theme(
axis.text.x = element_text(angle = 55, hjust = 1, color = "black"),
axis.text.y = element_text(color = "black"),
axis.title.x = element_blank(),
legend.title = element_blank(),
strip.text = element_blank(),
strip.background = element_rect(color = "white", fill = "white"),
panel.grid = element_blank()
)

#adjust grid widths
require(grid)
require(gtable)
gt = ggplot_gtable(ggplot_build(p1))
gtable_show_layout(gt)
gt
#these give the index in gt$widths that needs to be changed to adjust the widthof a panel
gt$layout$l[grep('panel-1-1', gt$layout$name)]
gt$layout$l[grep('panel-2-1', gt$layout$name)]


gt = ggplot_gtable(ggplot_build(p1))

#gt$widths[gt$layout$l[grep('panel-1-1', gt$layout$name)]] = gt$widths[gt$layout$l[grep('panel-1-1', gt$layout$name)]] * .15
gt$widths[gt$layout$l[grep('panel-2-1', gt$layout$name)]] = gt$widths[gt$layout$l[grep('panel-1-1', gt$layout$name)]] * 5.6

grid.draw(gt)

pdf("Fig_S1_negatives.pdf", width = 8, height = 4.5)
grid.draw(gt)
dev.off()


