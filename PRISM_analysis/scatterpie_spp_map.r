library(ggplot2)
library(tidyverse)
library(maps)
require(scatterpie)

source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

us_map = map_data(map='usa')

states <- map_data("state")
dim(states)

ggplot(data = states) +
geom_polygon(aes(x = long, y = lat, fill = region, group = group), color = "white") +
coord_fixed(1.3) +
guides(fill=FALSE)  # do this to leave off the color legend

states_include = c(
"maine",
"new hampshire",
"vermont",
"new york",
"massachusetts",
"connecticut",
"rhode island",
"pennsylvania",
"delaware",
"maryland",
"virginia",
"west virginia",
"north carolina",
"tennessee",
"kentucky",
"ohio",
"michigan",
"wisconsin",
"illinois",
"indiana"
)

states.study_area = subset(states, region %in% states_include)

ggplot(data = states.study_area) +
geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") +
coord_fixed(1.3) +
geom_point(data = site_info, aes(y = lat, x = lon)) +
my_gg_theme

#########################
#Get neonectria occurence

#filter by 1K min seqs
lt_1K_samps = full_metadata %>%
filter(total_seqs < 1000) %>%
dplyr::select("sample")

#Species matrix on Nf and Nd occurence
Nf_v_Nd.bin.gt1K = Nf_v_Nd.bin %>%
filter(!sample %in% lt_1K_samps$sample)

full_metadata.sorted = left_join(
Nf_v_Nd.bin.gt1K,
full_metadata,
by = "sample"
)

##############################
#Calculate Nf and Nd frequency as individual columns (i.e., no "both")

Nf_Nd_site_freq.spp = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
N.faginata = vector(mode = "numeric", length = length(site_info$Site)),
N.ditissima = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted %>% filter(Site == site_info$Site[i])
    Nf_Nd_site_freq.spp$N.faginata[i] = ((filter(temp_tab, occurence == "Nf" | occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq.spp$N.ditissima[i] = ((filter(temp_tab, occurence == "Nd" | occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq.spp$Site[i] = as.character(site_info$Site[i])
}

###############################



Nf_Nd_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
neither = vector(mode = "numeric", length = length(site_info$Site)),
N.faginata = vector(mode = "numeric", length = length(site_info$Site)),
N.ditissima = vector(mode = "numeric", length = length(site_info$Site)),
both = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted %>% filter(Site == site_info$Site[i])
    Nf_Nd_site_freq$neither[i] = ((filter(temp_tab, occurence == "none")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$N.faginata[i] = ((filter(temp_tab, occurence == "Nf")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$N.ditissima[i] = ((filter(temp_tab, occurence == "Nd")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$both[i] = ((filter(temp_tab, occurence == "both")) %>% nrow)/nrow(temp_tab)
    Nf_Nd_site_freq$Site[i] = as.character(site_info$Site[i])
}

Nf_Nd_site_freq.metadata = left_join(
    Nf_Nd_site_freq,
    site_info %>% select(Site, lat, lon)
)


######################
#Plot piecharts on map

p = ggplot() +
geom_polygon(
    data = states.study_area,
    aes(x = long, y = lat, group = group),
    fill = "light grey", color = "black"
) +
coord_fixed(1.0) +
geom_scatterpie(
    data = Nf_Nd_site_freq.metadata,
    aes(x = lon, y = lat),
    pie_scale = 2.5,
    cols = c("both", "N.faginata", "N.ditissima", "neither")
) +
scale_fill_manual(
    values = c("both" = "black", "N.faginata" = "grey42", "N.ditissima" = "grey72", "neither" = "white"),
    labels = c("both" = "both spp.", "N.faginata" = "N.faginata", "N.ditissima" = "N.ditissima", "neither" = "neither spp."),
    breaks = c("both", "N.faginata", "N.ditissima", "neither")
) +
labs(
x = NULL,
y = NULL,
fill = "Neonectria spp.\ndetected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/species_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()



