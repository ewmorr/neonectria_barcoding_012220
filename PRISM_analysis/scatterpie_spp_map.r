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
guides(fill=FALSE) # leave off the color legend
#scale_x_continuous(limits = c(-92, -67))
#scale_y_continuous(limits = c(35,80))

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

p1 = ggplot(data = states.study_area) +
geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") +
coord_fixed(1.3) +
#coord_quickmap()+
geom_point(data = site_info, aes(y = lat, x = lon)) +
my_gg_theme
p2 = ggplot(data = states.study_area) +
geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") +
#coord_fixed(1.3) +
coord_quickmap()+
geom_point(data = site_info, aes(y = lat, x = lon)) +
my_gg_theme
p3 = ggplot(data = states.study_area) +
geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") +
#coord_fixed(1.3) +
coord_map()+
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

site_n = full_metadata.sorted %>% group_by(Site) %>% summarize(n = n())

Nf_Nd_site_freq.metadata = left_join(
    Nf_Nd_site_freq,
    site_info %>% select(Site, lat, lon)
) %>% left_join(., site_n)


######################
#Plot piecharts on map

#b/w pies, green background
p = ggplot() +
geom_polygon(
data = states.study_area,
aes(x = long, y = lat*1.3, group = group),
fill = "#2ca25f", color = "black"
) +
coord_fixed(1) +
#coord_quickmap() +
geom_scatterpie(
data = Nf_Nd_site_freq.metadata,
    aes(x = lon, y = lat*1.3), pie_scale = 2.5,
cols = c("both", "N.faginata", "N.ditissima", "neither")
) +
#geom_scatterpie_legend(log10(c(5,10,20)), x=-70, y=48, n = 3, labeller = function(x){floor(10^x)}) +
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


pdf("PRISM_maps/species_occurence_in_sites.NEW.pdf", width = 8, height = 4)
p
dev.off()

#size scaling, b/w pies, green background
p = ggplot() +
geom_polygon(
    data = states.study_area,
    aes(x = long, y = lat*1.3, group = group),
    fill = "#2ca25f", color = "black"
) +
coord_fixed(1) +
#coord_quickmap() +
geom_scatterpie(
    data = Nf_Nd_site_freq.metadata,
    aes(x = lon, y = lat*1.3, r = log10(n)),
    pie_scale = 2.5,
    cols = c("both", "N.faginata", "N.ditissima", "neither")
) +
geom_scatterpie_legend(log10(c(5,10,20)), x=-70, y=48, n = 3, labeller = function(x){floor(10^x)}) +
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


pdf("PRISM_maps/species_occurence_in_sites.NEW.logSize.pdf", width = 8, height = 4)
p
dev.off()



####################

#size scaling, color pies, grey background
p = ggplot() +
geom_polygon(
data = states.study_area,
aes(x = long, y = lat*1.3, group = group),
fill = "light grey", color = "black"
) +
coord_fixed(1) +
#coord_quickmap() +
geom_scatterpie(
data = Nf_Nd_site_freq.metadata,
aes(x = lon, y = lat*1.3), pie_scale = 2.5,
cols = c("both", "N.faginata", "N.ditissima", "neither")
) +
#geom_scatterpie_legend(log10(c(5,10,20)), x=-70, y=48, n = 3, labeller = function(x){floor(10^x)}) +
scale_fill_manual(
values = c("both" = "#fc8d59", "N.faginata" = "#d7301f", "N.ditissima" = "#fdcc8a", "neither" = "#fef0d9"),
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


pdf("PRISM_maps/species_occurence_in_sites.NEW.colorPie.pdf", width = 8, height = 4)
p
dev.off()

# color pies, grey background
p = ggplot() +
geom_polygon(
data = states.study_area,
aes(x = long, y = lat*1.3, group = group),
fill = "light grey", color = "black"
) +
coord_fixed(1) +
#coord_quickmap() +
geom_scatterpie(
data = Nf_Nd_site_freq.metadata,
aes(x = lon, y = lat*1.3, r = log10(n)),
pie_scale = 2.5,
cols = c("both", "N.faginata", "N.ditissima", "neither")
) +
geom_scatterpie_legend(log10(c(5,10,20)), x=-70, y=48, n = 3, labeller = function(x){floor(10^x)}) +
scale_fill_manual(
values = c("both" = "#fc8d59", "N.faginata" = "#d7301f", "N.ditissima" = "#fdcc8a", "neither" = "#fef0d9"),
labels = c("both" = "both spp.", "N.faginata" = expression(italic("N. faginata")), "N.ditissima" = expression(italic("N. ditissima")), "neither" = "neither spp."),
breaks = c("both", "N.faginata", "N.ditissima", "neither")
) +
labs(
x = NULL,
y = NULL,
fill = expression(atop(paste(italic("Neonectria"), " spp.")~phantom(1), "(proportion trees)"))
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
legend.text = element_text(hjust = 0),
axis.text = element_blank(),
axis.ticks = element_blank()
)


pdf("PRISM_maps/species_occurence_in_sites.NEW.logSize.colorPie.pdf", width = 8, height = 4)
p
dev.off()



