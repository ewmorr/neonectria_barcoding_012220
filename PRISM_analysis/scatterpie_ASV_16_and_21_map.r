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
#Get ASV_16 occurence

#filter by 1K min seqs and set binary
lt_1K_samps = full_metadata %>%
filter(total_seqs < 1000) %>%
dplyr::select("sample")

asv_tab.gt1K.asv_16 = t(asv_tab["ASV_16",!colnames(asv_tab) %in% lt_1K_samps$sample])
asv_tab.gt1K.asv_21 = t(asv_tab["ASV_21",!colnames(asv_tab) %in% lt_1K_samps$sample])
asv_tab.gt1K.asv_16 = asv_tab.gt1K.asv_16 + asv_tab.gt1K.asv_21

asv_tab.gt1K.asv_16[asv_tab.gt1K.asv_16 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_16 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_16), asv_tab.gt1K.asv_16),
full_metadata,
by = "sample"
)


asv_16_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_16 %>% filter(Site == site_info$Site[i])
    asv_16_site_freq$absent[i] = ((filter(temp_tab, ASV_16 == 0)) %>% nrow)/nrow(temp_tab)
    asv_16_site_freq$present[i] = ((filter(temp_tab, ASV_16 == 1)) %>% nrow)/nrow(temp_tab)
    asv_16_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_16_site_freq.metadata = left_join(
    asv_16_site_freq,
    site_info %>% dplyr::select(Site, lat, lon)
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
    data = asv_16_site_freq.metadata,
    aes(x = lon, y = lat),
    pie_scale = 2.5,
    cols = c("present", "absent")
) +
scale_fill_manual(
    values = c("absent" = "white", "present" = "dark grey")
) +
labs(
x = NULL,
y = NULL,
fill = "Clonostachys rosea\ndetected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/Clonostachys_rosea_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()



