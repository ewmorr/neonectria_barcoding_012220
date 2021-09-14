require(tidyverse)
require(RColorBrewer)
require(gridExtra)
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
source("~/ggplot_theme.txt")

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


asv_tab.gt1K.asv_16 = t(asv_tab["ASV_16",!colnames(asv_tab) %in% lt_1K_samps$sample])
asv_tab.gt1K.asv_21 = t(asv_tab["ASV_21",!colnames(asv_tab) %in% lt_1K_samps$sample])
asv_tab.gt1K.asv_16 = asv_tab.gt1K.asv_16 + asv_tab.gt1K.asv_21

asv_tab.gt1K.asv_16[asv_tab.gt1K.asv_16 > 0] = 1


asv_tab.gt1K.ASV_807 = t(asv_tab["ASV_807",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.ASV_807[asv_tab.gt1K.ASV_807 > 0] = 1


asv_tab.gt1K.asv_19 = t(asv_tab["ASV_19",!colnames(asv_tab) %in% lt_1K_samps$sample])

asv_tab.gt1K.asv_19[asv_tab.gt1K.asv_19 > 0] = 1


all_spp_occurence = full_join(
    full_metadata.sorted %>% select(Nf, sample, Nd),
    data.frame(sample = rownames(asv_tab.gt1K.asv_16), asv_tab.gt1K.asv_16)
) %>% full_join(.,
    data.frame(sample = rownames(asv_tab.gt1K.ASV_807), asv_tab.gt1K.ASV_807)
) %>% full_join(.,
data.frame(sample = rownames(asv_tab.gt1K.asv_19), asv_tab.gt1K.asv_19)
)

colnames(all_spp_occurence) = c("N. faginata", "sample", "N. ditissima", "C. rosea", "N. ferrugineum", "F. babinda")

all_spp_occurence.long = pivot_longer(all_spp_occurence, -sample)

colnames(all_spp_occurence) = c("N.faginata", "sample", "N.ditissima", "C.rosea", "N.ferrugineum", "F.babinda")

all_spp_occurence.meta = left_join(all_spp_occurence, full_metadata.sorted %>% select(sample, Site))

#ggplots of point clouds
all_spp_occurence.meta$Site = factor(all_spp_occurence.meta$Site, levels = c("MEN1","MES1", "CW1", "ADN1", "ADS1","TSP1", "MI1", "WF1", "GK1j", "ASH2"))

###########
#Black points only with alpha

p1 = ggplot(all_spp_occurence.meta, aes(y = N.faginata, x = C.rosea)) +
geom_point(position = position_jitter(height = 0.15, width = 0.15), size = 3, alpha = 0.25, color = "black") +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
labs(y = expression(italic("N. faginata")), x = NULL) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6),
axis.text.x = element_blank()
)

p2 = ggplot(all_spp_occurence.meta, aes(y = N.faginata, x = N.ferrugineum)) +
geom_point(position = position_jitter(height = 0.15, width = 0.15), size = 3, alpha = 0.25, color = "black") +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
labs(y = "", x = NULL) +
my_gg_theme +
theme(
axis.text.x = element_blank(),
axis.text.y = element_blank()
)

p3 = ggplot(all_spp_occurence.meta, aes(y = N.faginata, x = F.babinda)) +
geom_point(position = position_jitter(height = 0.15, width = 0.15), size = 3, alpha = 0.25, color = "black") +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
labs(y = "", x = NULL) +
my_gg_theme +
guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) +
theme(
axis.text.x = element_blank(),
axis.text.y = element_blank()
)

p4 = ggplot(all_spp_occurence.meta, aes(y = N.ditissima, x = C.rosea)) +
geom_point(position = position_jitter(height = 0.15, width = 0.15), size = 3, alpha = 0.25, color = "black") +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
labs(y = expression(italic("N. ditissima")), x = expression(italic("C. rosea"))) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

p5 = ggplot(all_spp_occurence.meta, aes(N.ditissima, N.ferrugineum)) +
geom_point(position = position_jitter(height = 0.15, width = 0.15), size = 3, alpha = 0.25, color = "black") +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
labs(y = "", x = expression(italic("N. ferrugineum"))) +
my_gg_theme +
theme(
axis.text.y = element_blank()
)

p6 = ggplot(all_spp_occurence.meta, aes(N.ditissima, F.babinda)) +
geom_point(position = position_jitter(height = 0.15, width = 0.15), size = 3, alpha = 0.25, color = "black") +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
labs(y = "", x = expression(italic("F. babinda"))) +
my_gg_theme +
guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) +
theme(
axis.text.y = element_blank(),
legend.title = element_text(size = 15),
legend.text = element_text(size = 15)
)

pdf("neos_v_other_mapped_w_lines.black_points.pdf", width = 11, height = 6)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3, heights = c(0.45, 0.55))
dev.off()

#################
#Scatterpie plots
require(scatterpie)

C.rosea.Nf.Nd = all_spp_occurence.meta %>% group_by(N.faginata, N.ditissima, C.rosea) %>% summarize(n = n())
C.rosea.Nf.Nd.wide = C.rosea.Nf.Nd %>% pivot_wider(names_from = C.rosea, values_from = n)
colnames(C.rosea.Nf.Nd.wide) = c("N.faginata", "N.ditissima", "absent", "present")
C.rosea.Nf.Nd.wide$n = C.rosea.Nf.Nd.wide$absent + C.rosea.Nf.Nd.wide$present
C.rosea.Nf.Nd.wide$

N.ferrugineum.Nf.Nd = all_spp_occurence.meta %>% group_by(N.faginata, N.ditissima, N.ferrugineum) %>% summarize(n = n())
N.ferrugineum.Nf.Nd.wide = N.ferrugineum.Nf.Nd %>% pivot_wider(names_from = N.ferrugineum, values_from = n)
colnames(N.ferrugineum.Nf.Nd.wide) = c("N.faginata", "N.ditissima", "absent", "present")
N.ferrugineum.Nf.Nd.wide[is.na(N.ferrugineum.Nf.Nd.wide)] = 0
N.ferrugineum.Nf.Nd.wide$n = N.ferrugineum.Nf.Nd.wide$absent + N.ferrugineum.Nf.Nd.wide$present


F.babinda.Nf.Nd = all_spp_occurence.meta %>% group_by(N.faginata, N.ditissima, F.babinda) %>% summarize(n = n())
F.babinda.Nf.Nd.wide = F.babinda.Nf.Nd %>% pivot_wider(names_from = F.babinda, values_from = n)
colnames(F.babinda.Nf.Nd.wide) = c("N.faginata", "N.ditissima", "absent", "present")
F.babinda.Nf.Nd.wide[is.na(F.babinda.Nf.Nd.wide)] = 0
F.babinda.Nf.Nd.wide$n = F.babinda.Nf.Nd.wide$absent + F.babinda.Nf.Nd.wide$present


#scaled pies

p1 = ggplot(C.rosea.Nf.Nd.wide) +
geom_point(aes(x = N.faginata, y = N.ditissima)) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
coord_fixed(1.0) +
geom_scatterpie(
    data = C.rosea.Nf.Nd.wide,
    aes(x = N.faginata, y = N.ditissima,r = n/100),
    #pie_scale = 2.5,
    cols = c("present", "absent")
) +
geom_scatterpie_legend(c(5,15,25,45)/100, x=0, y=0, n = 4, labeller = function(x){x*100}) +
scale_fill_manual(
    values = c("absent" = "white", "present" = "#1c9099")
) +
labs(x = expression(italic("N. faginata")), y = expression(italic("N. ditissima")), fill = expression(italic("C. rosea"))) +
my_gg_theme +
theme(
legend.title = element_text(size = 18),
legend.text = element_text(size = 18)
)

ggplot(C.rosea.Nf.Nd.wide) +
geom_point(aes(x = N.faginata, y = N.ditissima)) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
coord_fixed(1.0) +
geom_scatterpie(
data = C.rosea.Nf.Nd.wide,
#aes(x = N.faginata, y = N.ditissima,r = log10(n)/5),
aes(x = N.faginata, y = N.ditissima,r = n/100),
#pie_scale = 2.5,
cols = c("present", "absent")
) +
#geom_scatterpie_legend(log10(c(5,15,25,45))/10, x=0, y=-0.3, n = 4, labeller = function(x){10^(x*10)}) +
geom_text(aes(x = N.faginata - 0.35 , y = N.ditissima + 0.45, label = paste("n =", n)), size = 5.5) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#1c9099")
) +
labs(x = expression(italic("N. faginata")), y = expression(italic("N. ditissima")), fill = expression(italic("C. rosea"))) +
my_gg_theme +
theme(
legend.title = element_text(size = 18),
legend.text = element_text(size = 18)
)


#not scaled pies
p1 = ggplot(C.rosea.Nf.Nd.wide) +
geom_point(aes(x = N.faginata, y = N.ditissima)) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
coord_fixed(1.0) +
geom_scatterpie(
    data = C.rosea.Nf.Nd.wide,
    aes(x = N.faginata, y = N.ditissima),
    pie_scale = 14,
    cols = c("present", "absent")
) +
scale_fill_manual(
    values = c("absent" = "white", "present" = "#1c9099")
) +
geom_text(aes(x = N.faginata , y = N.ditissima + 0.4, label = paste("n =", n)), size = 5.5) +
labs(x = expression(italic("N. faginata")), y = expression(italic("N. ditissima")), fill = expression(italic("C. rosea"))) +
my_gg_theme +
theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 17)
)



#########################
#Scaled pies with text labels for sample size

p1 =
ggplot(C.rosea.Nf.Nd.wide) +
geom_point(aes(x = N.faginata, y = N.ditissima)) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
coord_fixed(1.0) +
geom_scatterpie(
    data = C.rosea.Nf.Nd.wide,
    aes(x = N.faginata, y = N.ditissima,r = n/110),
    cols = c("present", "absent")
) +
#geom_text(aes(x = N.faginata - 0.1, y = N.ditissima + 0.45, label = paste(paste("n =", n), paste(", present =", present), sep = "")), size = 5) +
geom_text(aes(x = N.faginata, y = N.ditissima + 0.45, label = paste(paste("n =", n), paste(", present =", present), sep = "")), size = 5) +
scale_fill_manual(
    values = c("absent" = "white", "present" = "#1c9099"), breaks = c("absent", "present")
) +
labs(x = "", y = expression(italic("N. ditissima")), fill = expression(italic("C. rosea"))) +
my_gg_theme +
theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    #legend.position = c(x = 0.3575, y = 0.61),
    legend.position = c(x = 0.138, y = 0.615),
    #legend.position = "top",
#    legend.background = element_rect(color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.6)
)


p2 =
ggplot(N.ferrugineum.Nf.Nd.wide) +
geom_point(aes(x = N.faginata, y = N.ditissima)) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
coord_fixed(1.0) +
geom_scatterpie(
data = N.ferrugineum.Nf.Nd.wide,
aes(x = N.faginata, y = N.ditissima,r = n/110),
cols = c("present", "absent")
) +
#geom_text(aes(x = N.faginata, y = N.ditissima + 0.45, label = paste(paste("n =", n), paste(", present =", present), sep = "")), size = 5) +
#geom_text(aes(x = N.faginata - 0.25, y = N.ditissima + 0.45, label = paste("present =", present)), size = 5) +
geom_text(aes(x = N.faginata, y = N.ditissima + 0.45, label = paste("present =", present)), size = 5) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#67a9cf"), breaks = c("absent", "present")
) +
labs(x = expression(italic("N. faginata")), y = "", fill = expression(italic("N. ferrugineum"))) +
my_gg_theme +
theme(
legend.title = element_text(size = 13),
legend.text = element_text(size = 14),
#legend.position = c(x = 0.3575, y = 0.61),
legend.position = c(x = 0.165, y = 0.61375),
#legend.position = "top",
#legend.background = element_rect(color = "black"),
axis.text.y = element_text(angle = 90, hjust = 0.6)
)


p3 =
ggplot(F.babinda.Nf.Nd.wide) +
geom_point(aes(x = N.faginata, y = N.ditissima)) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5), expand = c(0,0)) +
coord_fixed(1.0) +
geom_scatterpie(
data = F.babinda.Nf.Nd.wide,
aes(x = N.faginata, y = N.ditissima,r = n/110),
cols = c("present", "absent")
) +
#geom_text(aes(x = N.faginata, y = N.ditissima + 0.45, label = paste(paste("n =", n), paste(", present =", present), sep = "")), size = 5) +
#geom_text(aes(x = N.faginata - 0.23, y = N.ditissima + 0.45, label = paste("present =", present)), size = 5) +
geom_text(aes(x = N.faginata, y = N.ditissima + 0.45, label = paste("present =", present)), size = 5) +
scale_fill_manual(
values = c("absent" = "white", "present" = "#016c59"), breaks = c("absent", "present")
) +
labs(x = "", y = "", fill = expression(italic("F. babinda"))) +
my_gg_theme +
theme(
legend.title = element_text(size = 14),
legend.text = element_text(size = 14),
#legend.position = c(x = 0.3575, y = 0.61),
legend.position = c(x = 0.138, y = 0.615),
#legend.position = "top",
#legend.background = element_rect(color = "black"),
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

require(gridExtra)

pdf("mapped_spp_occurence_w_Nf_Nd.pdf", width = 15, height = 6)
grid.arrange(p1,p2,p3, ncol = 3)
dev.off()

#############
#Older plots#
#############
#With site colors
p1 = ggplot(all_spp_occurence.meta, aes(N.faginata, C.rosea, color = Site, shape = as.factor(N.ditissima))) +
geom_point(position = position_jitter(height = 0.1, width = 0.1), size = 2.5) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present"), guide = F) +
scale_color_brewer(palette = "Paired",
    labels = c(
    "ADN1" = "New York N",
    "ADS1" = "New York S",
    "ASH2" = "North Carolina",
    "CW1" = "New Hampshire",
    "GK1j" = "West Virginia",
    "MEN1" = "Maine N",
    "MES1" = "Maine S",
    "MI1" = "Michigan",
    "TSP1" = "Pennsylvania",
    "WF1" = "Wisconsin"
    ), guide = F
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.1,1.1)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
labs(x = expression(italic("N. faginata")), y = expression(italic("N. ferrugineum"))) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

p2 = ggplot(all_spp_occurence.meta, aes(N.faginata, N.ferrugineum, color = Site, shape = as.factor(N.ditissima))) +
geom_point(position = position_jitter(height = 0.1, width = 0.1), size = 2.5) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present"), guide = F) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
), guide = F
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.1,1.1)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
labs(x = expression(italic("N. faginata")), y = expression(italic("N. ferrugineum"))) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

p3 = ggplot(all_spp_occurence.meta, aes(N.faginata, F.babinda, color = Site, shape = as.factor(N.ditissima))) +
geom_point(position = position_jitter(height = 0.1, width = 0.1), size = 2.5) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present")) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
)
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.1,1.1)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
labs(color = NULL, shape = expression(italic("N. ditissima")), x = expression(italic("N. faginata")), y = expression(italic("F. babinda"))) +
my_gg_theme +
guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) +
theme(
    axis.text.y = element_text(angle = 90, hjust = 0.6),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
)

require(gridExtra)

pdf("neos_v_other_mapped.pdf", width = 16, height = 4.5)
grid.arrange(p1,p2,p3, ncol = 3, widths = c(0.29, 0.29, 0.42))
dev.off()


p4 = ggplot(all_spp_occurence.meta, aes(N.faginata, C.rosea, color = Site, shape = as.factor(N.ditissima))) +
geom_point(position = position_jitter(height = 0.45, width = 0.45), size = 3.2) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present"), guide = F) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
), guide = F
) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
labs(x = "", y = expression(italic("C. rosea"))) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

p5 = ggplot(all_spp_occurence.meta, aes(N.faginata, N.ferrugineum, color = Site, shape = as.factor(N.ditissima))) +
geom_point(position = position_jitter(height = 0.45, width = 0.45), size = 3.2) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present"), guide = F) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
), guide = F
) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
labs(x = expression(italic("N. faginata")), y = expression(italic("N. ferrugineum"))) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

p6 = ggplot(all_spp_occurence.meta, aes(N.faginata, F.babinda, color = Site, shape = as.factor(N.ditissima))) +
geom_point(position = position_jitter(height = 0.45, width = 0.45), size = 3.2) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present")) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
)
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present"), limits = c(-0.5,1.5)) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
geom_hline(yintercept = 0.5, color = "black") +
geom_vline(xintercept = 0.5, color = "black") +
annotate(geom= "rect", xmin = 0.5, xmax = 0.5, ymin = -0.5, ymax = 1.5) +
labs(color = NULL, shape = expression(italic("N. ditissima")), x = "", y = expression(italic("F. babinda"))) +
my_gg_theme +
guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6),
legend.title = element_text(size = 15),
legend.text = element_text(size = 15)
)

require(gridExtra)

pdf("neos_v_other_mapped_w_lines.pdf", width = 16, height = 4.5)
grid.arrange(p4,p5,p6, ncol = 3, widths = c(0.29, 0.29, 0.42))
dev.off()



ggplot(all_spp_occurence.meta, aes(N.ditissima, C.rosea, color = Site, shape = as.factor(N.faginata))) +
geom_point(position = position_jitter(height = 0.05, width = 0.05), alpha = 0.75) +
scale_shape_manual(values = c("0" = 1, "1" = 19), labels = c("0" = "absent", "1" = "present"), guide = F) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
), guide = F
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present")) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

ggplot(all_spp_occurence.meta, aes(N.ditissima, N.ferrugineum, color = Site)) +
geom_point(position = position_jitter(height = 0.05, width = 0.05), alpha = 0.75) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
)
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present")) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)

ggplot(all_spp_occurence.meta, aes(N.ditissima, F.babinda, color = Site)) +
geom_point(position = position_jitter(height = 0.05, width = 0.05), alpha = 0.75) +
scale_color_brewer(palette = "Paired",
labels = c(
"ADN1" = "New York N",
"ADS1" = "New York S",
"ASH2" = "North Carolina",
"CW1" = "New Hampshire",
"GK1j" = "West Virginia",
"MEN1" = "Maine N",
"MES1" = "Maine S",
"MI1" = "Michigan",
"TSP1" = "Pennsylvania",
"WF1" = "Wisconsin"
)
) +
scale_y_continuous(breaks = c(0,1), labels = c("absent", "present")) +
scale_x_continuous(breaks = c(0,1), labels = c("absent", "present")) +
my_gg_theme +
theme(
axis.text.y = element_text(angle = 90, hjust = 0.6)
)


#Chi squared

nf_cr = all_spp_occurence.meta %>% group_by(N.faginata, C.rosea) %>% summarize(n = n())
nf_nefe = all_spp_occurence.meta %>% group_by(N.faginata, N.ferrugineum) %>% summarize(n = n())
nf_fb = all_spp_occurence.meta %>% group_by(N.faginata, F.babinda) %>% summarize(n = n())

nd_cr = all_spp_occurence.meta %>% group_by(N.ditissima, C.rosea) %>% summarize(n = n())
nd_nefe = all_spp_occurence.meta %>% group_by(N.ditissima, N.ferrugineum) %>% summarize(n = n())
nd_fb = all_spp_occurence.meta %>% group_by(N.ditissima, F.babinda) %>% summarize(n = n())


chisq.test(matrix(c(nf_cr$n[1], nf_cr$n[3], nf_cr$n[2], nf_cr$n[4]), 2,2,
dimnames = list(Nf = c("0", "1"), Cr = c("0", "1"))
)
)

fisher.test(matrix(c(nf_cr$n[1], nf_cr$n[3], nf_cr$n[2], nf_cr$n[4]), 2,2,
dimnames = list(Nf = c("0", "1"), Cr = c("0", "1"))
), alternative = "two.sided"
)

fisher.test(matrix(c(nf_nefe$n[1], nf_nefe$n[3], nf_nefe$n[2], 0), 2,2,
dimnames = list(Nf = c("0", "1"), nefe = c("0", "1"))
), alternative = "two.sided"
)

fisher.test(matrix(c(nf_fb$n[1], nf_fb$n[3], nf_fb$n[2], nf_fb$n[4]), 2,2,
dimnames = list(Nf = c("0", "1"), Cr = c("0", "1"))
), alternative = "two.sided"
)


fisher.test(matrix(c(nd_cr$n[1], nd_cr$n[3], nd_cr$n[2], nd_cr$n[4]), 2,2,
dimnames = list(nd = c("0", "1"), Cr = c("0", "1"))
), alternative = "two.sided"
)

fisher.test(matrix(c(nd_nefe$n[1], nd_nefe$n[3], nd_nefe$n[2], 0), 2,2,
dimnames = list(nd = c("0", "1"), nefe = c("0", "1"))
), alternative = "two.sided"
)

fisher.test(matrix(c(nd_fb$n[1], nd_fb$n[3], nd_fb$n[2], nd_fb$n[4]), 2,2,
dimnames = list(nd = c("0", "1"), Cr = c("0", "1"))
), alternative = "two.sided"
)





