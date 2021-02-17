require(gridExtra)
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

states.study_area.zoom = subset(states, region %in% states_include.zoom)

states.study_area.zoom.wisc = states.study_area.zoom %>% filter(region == "wisconsin" & lat <= 46.325)
states.study_area.zoom.rest = states.study_area.zoom %>% filter(region != "wisconsin")

states.study_area.zoom = rbind(states.study_area.zoom.rest, states.study_area.zoom.wisc)


states.study_area.zoom[states.study_area.zoom$long < -90, "long"] = -90
states.study_area.zoom[states.study_area.zoom$lat < 34, "lat"] = 34

#filter by 1K min seqs and set binary
lt_1K_samps = full_metadata %>%
filter(total_seqs < 1000) %>%
dplyr::select("sample")

###############
#ASV_126
###############

#########################
#Get ASV_126 occurence

asv_tab.gt1K.ASV_126 = t(asv_tab["ASV_126",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_126[asv_tab.gt1K.ASV_126 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_126 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_126), asv_tab.gt1K.ASV_126),
full_metadata,
by = "sample"
)


ASV_126_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_126 %>% filter(Site == site_info$Site[i])
    ASV_126_site_freq$absent[i] = ((filter(temp_tab, ASV_126 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_126_site_freq$present[i] = ((filter(temp_tab, ASV_126 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_126_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_126_site_freq.metadata = left_join(
ASV_126_site_freq,
site_info %>% select(Site, lat, lon)
)

###############
#ASV_575
###############

#########################
#Get ASV_575 occurence

asv_tab.gt1K.ASV_575 = t(asv_tab["ASV_575",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_575[asv_tab.gt1K.ASV_575 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_575 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_575), asv_tab.gt1K.ASV_575),
full_metadata,
by = "sample"
)


ASV_575_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_575 %>% filter(Site == site_info$Site[i])
    ASV_575_site_freq$absent[i] = ((filter(temp_tab, ASV_575 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_575_site_freq$present[i] = ((filter(temp_tab, ASV_575 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_575_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_575_site_freq.metadata = left_join(
ASV_575_site_freq,
site_info %>% select(Site, lat, lon)
)

###############
#ASV_234
###############

#########################
#Get ASV_234 occurence

asv_tab.gt1K.ASV_234 = t(asv_tab["ASV_234",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_234[asv_tab.gt1K.ASV_234 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_234 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_234), asv_tab.gt1K.ASV_234),
full_metadata,
by = "sample"
)


ASV_234_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_234 %>% filter(Site == site_info$Site[i])
    ASV_234_site_freq$absent[i] = ((filter(temp_tab, ASV_234 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_234_site_freq$present[i] = ((filter(temp_tab, ASV_234 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_234_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_234_site_freq.metadata = left_join(
ASV_234_site_freq,
site_info %>% select(Site, lat, lon)
)

###############
#ASV_27
###############

#########################
#Get ASV_27 occurence

asv_tab.gt1K.ASV_27 = t(asv_tab["ASV_27",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_27[asv_tab.gt1K.ASV_27 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_27 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_27), asv_tab.gt1K.ASV_27),
full_metadata,
by = "sample"
)


ASV_27_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_27 %>% filter(Site == site_info$Site[i])
    ASV_27_site_freq$absent[i] = ((filter(temp_tab, ASV_27 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_27_site_freq$present[i] = ((filter(temp_tab, ASV_27 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_27_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_27_site_freq.metadata = left_join(
ASV_27_site_freq,
site_info %>% select(Site, lat, lon)
)

######################
######################
#Plot piecharts on map

p1 = ggplot() +
geom_polygon(
data = states.study_area.zoom,
aes(x = long, y = lat*1.3, group = group),
fill = "light grey", color = "black"
) +
coord_fixed(1.0) +
geom_scatterpie(
data = ASV_126_site_freq.metadata,
aes(x = lon, y = lat*1.3),
pie_scale = 2.5,
cols = c("present", "absent")
) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#016c59")
) +
labs(
x = NULL,
y = NULL,
fill = expression(atop(italic("A. alternatum")~phantom (10), "(proportion trees)")),
title = "A"
)+
my_gg_theme +
theme(
axis.text = element_blank(),
axis.ticks = element_blank(),
legend.background = element_rect(color = "black"),
legend.position = c(0.79,0.2),
legend.title = element_text(size = 16, hjust = 0),
legend.text = element_text(size = 16),
plot.title = element_text(margin = margin(b=-24,t=7.5), hjust = 0.01)
)


###############
p2 = ggplot() +
geom_polygon(
data = states.study_area.zoom,
aes(x = long, y = lat*1.3, group = group),
fill = "light grey", color = "black"
) +
coord_fixed(1.0) +
geom_scatterpie(
data = ASV_575_site_freq.metadata,
aes(x = lon, y = lat*1.3),
pie_scale = 2.5,
cols = c("present", "absent")
) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#1c9099")
) +
labs(
x = NULL,
y = NULL,
fill = expression(atop(italic("Microcera sp.")~phantom(10), "(proportion trees)")),
title = "B"
)+
my_gg_theme +
theme(
axis.text = element_blank(),
axis.ticks = element_blank(),
legend.background = element_rect(color = "black"),
legend.position = c(0.79,0.2),
legend.title = element_text(size = 16, hjust = 0),
legend.text = element_text(size = 16),
plot.title = element_text(margin = margin(b=-24,t=7.5), hjust = 0.01)
)


###############
p3 = ggplot() +
geom_polygon(
data = states.study_area.zoom,
aes(x = long, y = lat*1.3, group = group),
fill = "light grey", color = "black"
) +
coord_fixed(1.0) +
geom_scatterpie(
data = ASV_234_site_freq.metadata,
aes(x = lon, y = lat*1.3),
pie_scale = 2.5,
cols = c("present", "absent")
) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#67a9cf")
) +
labs(
x = NULL,
y = NULL,
fill = expression(atop(italic("H. surugaensis")~phantom(10), "(proportion trees)")),
title = "C"
)+
my_gg_theme +
theme(
axis.text = element_blank(),
axis.ticks = element_blank(),
legend.background = element_rect(color = "black"),
legend.position = c(0.79,0.2),
legend.title = element_text(size = 16, hjust = 0),
legend.text = element_text(size = 16),
plot.title = element_text(margin = margin(b=-24,t=7.5), hjust = 0.01)
)


###############
p4 = ggplot() +
geom_polygon(
data = states.study_area.zoom,
aes(x = long, y = lat*1.3, group = group),
fill = "light grey", color = "black"
) +
coord_fixed(1.0) +
geom_scatterpie(
data = ASV_27_site_freq.metadata,
aes(x = lon, y = lat*1.3),
pie_scale = 2.5,
cols = c("present", "absent")
) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#bdc9e1")
) +
labs(
x = NULL,
y = NULL,
fill = expression(atop(italic("V. taibaiensis")~phantom(100), "(proportion trees)")),
title = "D"
)+
my_gg_theme +
theme(
axis.text = element_blank(),
axis.ticks = element_blank(),
legend.background = element_rect(color = "black"),
legend.position = c(0.79,0.2),
legend.title = element_text(size = 16, hjust = 0),
legend.text = element_text(size = 16),
plot.title = element_text(margin = margin(b=-24,t=7.5), hjust = 0.01)
)


pdf("PRISM_maps/Mycoparasite_occurence_in_sites.pdf", width = 12, height = 8)
grid.arrange(p1,p2,p3,p4,ncol = 2)
dev.off()



