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

#filter by 1K min seqs and set binary
lt_1K_samps = full_metadata %>%
filter(total_seqs < 1000) %>%
dplyr::select("sample")

###############
#ASV_5
###############

#########################
#Get ASV_5 occurence


asv_tab.gt1K.asv_5 = t(asv_tab["ASV_5",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_5[asv_tab.gt1K.asv_5 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_5 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_5), asv_tab.gt1K.asv_5),
full_metadata,
by = "sample"
)


asv_5_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_5 %>% filter(Site == site_info$Site[i])
    asv_5_site_freq$absent[i] = ((filter(temp_tab, ASV_5 == 0)) %>% nrow)/nrow(temp_tab)
    asv_5_site_freq$present[i] = ((filter(temp_tab, ASV_5 == 1)) %>% nrow)/nrow(temp_tab)
    asv_5_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_5_site_freq.metadata = left_join(
    asv_5_site_freq,
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
    data = asv_5_site_freq.metadata,
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
fill = "ASV_5 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_5_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_27
###############

#########################
#Get ASV_27 occurence

asv_tab.gt1K.asv_27 = t(asv_tab["ASV_27",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_27[asv_tab.gt1K.asv_27 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_27 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_27), asv_tab.gt1K.asv_27),
full_metadata,
by = "sample"
)


asv_27_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_27 %>% filter(Site == site_info$Site[i])
    asv_27_site_freq$absent[i] = ((filter(temp_tab, ASV_27 == 0)) %>% nrow)/nrow(temp_tab)
    asv_27_site_freq$present[i] = ((filter(temp_tab, ASV_27 == 1)) %>% nrow)/nrow(temp_tab)
    asv_27_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_27_site_freq.metadata = left_join(
asv_27_site_freq,
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
data = asv_27_site_freq.metadata,
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
fill = "ASV_27 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_27_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_234
###############

#########################
#Get ASV_234 occurence

asv_tab.gt1K.asv_234 = t(asv_tab["ASV_234",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_234[asv_tab.gt1K.asv_234 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_234 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_234), asv_tab.gt1K.asv_234),
full_metadata,
by = "sample"
)


asv_234_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_234 %>% filter(Site == site_info$Site[i])
    asv_234_site_freq$absent[i] = ((filter(temp_tab, ASV_234 == 0)) %>% nrow)/nrow(temp_tab)
    asv_234_site_freq$present[i] = ((filter(temp_tab, ASV_234 == 1)) %>% nrow)/nrow(temp_tab)
    asv_234_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_234_site_freq.metadata = left_join(
asv_234_site_freq,
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
data = asv_234_site_freq.metadata,
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
fill = "ASV_234 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_234_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_575
###############

#########################
#Get ASV_575 occurence

asv_tab.gt1K.asv_575 = t(asv_tab["ASV_575",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_575[asv_tab.gt1K.asv_575 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_575 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_575), asv_tab.gt1K.asv_575),
full_metadata,
by = "sample"
)


asv_575_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_575 %>% filter(Site == site_info$Site[i])
    asv_575_site_freq$absent[i] = ((filter(temp_tab, ASV_575 == 0)) %>% nrow)/nrow(temp_tab)
    asv_575_site_freq$present[i] = ((filter(temp_tab, ASV_575 == 1)) %>% nrow)/nrow(temp_tab)
    asv_575_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_575_site_freq.metadata = left_join(
asv_575_site_freq,
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
data = asv_575_site_freq.metadata,
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
fill = "ASV_575 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_575_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()


###############
#ASV_94
###############

#########################
#Get ASV_94 occurence

asv_tab.gt1K.asv_94 = t(asv_tab["ASV_94",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_94[asv_tab.gt1K.asv_94 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_94 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_94), asv_tab.gt1K.asv_94),
full_metadata,
by = "sample"
)


asv_94_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_94 %>% filter(Site == site_info$Site[i])
    asv_94_site_freq$absent[i] = ((filter(temp_tab, ASV_94 == 0)) %>% nrow)/nrow(temp_tab)
    asv_94_site_freq$present[i] = ((filter(temp_tab, ASV_94 == 1)) %>% nrow)/nrow(temp_tab)
    asv_94_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_94_site_freq.metadata = left_join(
asv_94_site_freq,
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
data = asv_94_site_freq.metadata,
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
fill = "ASV_94 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_94_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()


###############
#ASV_135
###############

#########################
#Get ASV_135 occurence

asv_tab.gt1K.asv_135 = t(asv_tab["ASV_135",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_135[asv_tab.gt1K.asv_135 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_135 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_135), asv_tab.gt1K.asv_135),
full_metadata,
by = "sample"
)


asv_135_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_135 %>% filter(Site == site_info$Site[i])
    asv_135_site_freq$absent[i] = ((filter(temp_tab, ASV_135 == 0)) %>% nrow)/nrow(temp_tab)
    asv_135_site_freq$present[i] = ((filter(temp_tab, ASV_135 == 1)) %>% nrow)/nrow(temp_tab)
    asv_135_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_135_site_freq.metadata = left_join(
asv_135_site_freq,
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
data = asv_135_site_freq.metadata,
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
fill = "ASV_135 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_135_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_79
###############

#########################
#Get ASV_79 occurence

asv_tab.gt1K.asv_79 = t(asv_tab["ASV_79",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_79[asv_tab.gt1K.asv_79 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_79 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_79), asv_tab.gt1K.asv_79),
full_metadata,
by = "sample"
)


asv_79_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_79 %>% filter(Site == site_info$Site[i])
    asv_79_site_freq$absent[i] = ((filter(temp_tab, ASV_79 == 0)) %>% nrow)/nrow(temp_tab)
    asv_79_site_freq$present[i] = ((filter(temp_tab, ASV_79 == 1)) %>% nrow)/nrow(temp_tab)
    asv_79_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_79_site_freq.metadata = left_join(
asv_79_site_freq,
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
data = asv_79_site_freq.metadata,
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
fill = "ASV_79 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_79_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_10
###############

#########################
#Get ASV_10 occurence

asv_tab.gt1K.asv_10 = t(asv_tab["ASV_10",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_10[asv_tab.gt1K.asv_10 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_10 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_10), asv_tab.gt1K.asv_10),
full_metadata,
by = "sample"
)


asv_10_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_10 %>% filter(Site == site_info$Site[i])
    asv_10_site_freq$absent[i] = ((filter(temp_tab, ASV_10 == 0)) %>% nrow)/nrow(temp_tab)
    asv_10_site_freq$present[i] = ((filter(temp_tab, ASV_10 == 1)) %>% nrow)/nrow(temp_tab)
    asv_10_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_10_site_freq.metadata = left_join(
asv_10_site_freq,
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
data = asv_10_site_freq.metadata,
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
fill = "ASV_10 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_10_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()


###############
#ASV_12
###############

#########################
#Get ASV_12 occurence

asv_tab.gt1K.asv_12 = t(asv_tab["ASV_12",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_12[asv_tab.gt1K.asv_12 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_12 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.asv_12), asv_tab.gt1K.asv_12),
full_metadata,
by = "sample"
)


asv_12_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_12 %>% filter(Site == site_info$Site[i])
    asv_12_site_freq$absent[i] = ((filter(temp_tab, ASV_12 == 0)) %>% nrow)/nrow(temp_tab)
    asv_12_site_freq$present[i] = ((filter(temp_tab, ASV_12 == 1)) %>% nrow)/nrow(temp_tab)
    asv_12_site_freq$Site[i] = as.character(site_info$Site[i])
}

asv_12_site_freq.metadata = left_join(
asv_12_site_freq,
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
data = asv_12_site_freq.metadata,
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
fill = "ASV_12 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_12_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_807
###############

#########################
#Get ASV_807 occurence

asv_tab.gt1K.ASV_807 = t(asv_tab["ASV_807",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_807[asv_tab.gt1K.ASV_807 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_807 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_807), asv_tab.gt1K.ASV_807),
full_metadata,
by = "sample"
)


ASV_807_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_807 %>% filter(Site == site_info$Site[i])
    ASV_807_site_freq$absent[i] = ((filter(temp_tab, ASV_807 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_807_site_freq$present[i] = ((filter(temp_tab, ASV_807 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_807_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_807_site_freq.metadata = left_join(
ASV_807_site_freq,
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
data = ASV_807_site_freq.metadata,
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
fill = "ASV_807 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_807_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

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
data = ASV_126_site_freq.metadata,
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
fill = "ASV_126 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_126_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

###############
#ASV_69
###############

#########################
#Get ASV_69 occurence

asv_tab.gt1K.ASV_69 = t(asv_tab["ASV_69",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_69[asv_tab.gt1K.ASV_69 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_69 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_69), asv_tab.gt1K.ASV_69),
full_metadata,
by = "sample"
)


ASV_69_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_69 %>% filter(Site == site_info$Site[i])
    ASV_69_site_freq$absent[i] = ((filter(temp_tab, ASV_69 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_69_site_freq$present[i] = ((filter(temp_tab, ASV_69 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_69_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_69_site_freq.metadata = left_join(
ASV_69_site_freq,
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
data = ASV_69_site_freq.metadata,
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
fill = "ASV_69 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_69_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()


###############
#ASV_217
###############

#########################
#Get ASV_217 occurence

asv_tab.gt1K.ASV_217 = t(asv_tab["ASV_217",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_217[asv_tab.gt1K.ASV_217 > 0] = 1
#Species matrix on Nf and Nd occurence

full_metadata.sorted.ASV_217 = left_join(
data.frame(sample = rownames(asv_tab.gt1K.ASV_217), asv_tab.gt1K.ASV_217),
full_metadata,
by = "sample"
)


ASV_217_site_freq = data.frame(
Site = vector(mode = "character", length = length(site_info$Site)),
absent = vector(mode = "numeric", length = length(site_info$Site)),
present = vector(mode = "numeric", length = length(site_info$Site)),
stringsAsFactors = FALSE
)

for( i in 1:length(site_info$Site)){
    temp_tab = full_metadata.sorted.ASV_217 %>% filter(Site == site_info$Site[i])
    ASV_217_site_freq$absent[i] = ((filter(temp_tab, ASV_217 == 0)) %>% nrow)/nrow(temp_tab)
    ASV_217_site_freq$present[i] = ((filter(temp_tab, ASV_217 == 1)) %>% nrow)/nrow(temp_tab)
    ASV_217_site_freq$Site[i] = as.character(site_info$Site[i])
}

ASV_217_site_freq.metadata = left_join(
ASV_217_site_freq,
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
data = ASV_217_site_freq.metadata,
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
fill = "ASV_217 detected\n(proportion trees)"
)+
my_gg_theme +
theme(
legend.title = element_text(size = 20),
axis.text = element_blank(),
axis.ticks = element_blank()
)

pdf("PRISM_maps/ASV_217_occurence_in_sites.pdf", width = 8, height = 6)
p
dev.off()

