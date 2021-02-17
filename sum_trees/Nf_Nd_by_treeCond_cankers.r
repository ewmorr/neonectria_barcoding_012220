library(tidyverse)
require(gridExtra)
source("~/ggplot_theme.txt")

source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

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

full_metadata.sorted$occurence = factor(full_metadata.sorted$occurence, levels = c("none", "Nf", "Nd", "both"))


tree_cond_table = full_metadata.sorted %>% group_by(TreeCond, occurence) %>% summarize(n = n())
cankers_table = full_metadata.sorted %>% group_by(RaisedCanker, occurence) %>% summarize(n = n())

summary(aov(n ~ RaisedCanker, cankers_table %>% filter(occurence == "none")))

summary(aov(n ~ RaisedCanker, cankers_table %>% filter(occurence == "both")))
summary(aov(n ~ TreeCond, tree_cond_table %>% filter(occurence == "both")))


#######################################
#Offset bars (Jeff's suggested layout)#
#######################################
full_metadata.sorted$NfNd.sum = full_metadata.sorted[,"Nf"] + full_metadata.sorted[,"Nd"]

for(i in 1:nrow(full_metadata.sorted)){
    if(full_metadata.sorted$NfNd.sum[i] == 1 & full_metadata.sorted$Nd[i] == 1){
        full_metadata.sorted$NfNd.sum[i] = 0.5
    }
}


tree_cond.presAbs = full_metadata.sorted %>%
select(c("Nf", "Nd", "Site", "sample", "TreeCond", "NfNd.sum")) %>%
pivot_longer(c(-Site, -sample, -TreeCond, -NfNd.sum))


ggplot(tree_cond.presAbs %>% filter(value != 0 & Site == "WF1"), aes(x = as.factor(sample), y = name, group = name, color = name)) +
geom_segment(size = 5) +
scale_y_discrete(expand = c(0, 1), breaks = NULL) +
#xlim(
coord_polar() +
my_gg_theme

ggplot(tree_cond.presAbs %>% filter(Site == "WF1"), aes(x = reorder(as.factor(sample), NfNd.sum), y = name, group = name, fill = as.factor(value))) +
geom_tile(color = "black") +
scale_y_discrete(expand = c(0, 1), breaks = NULL) +
#xlim(
coord_polar() +
scale_fill_manual(values = c("white", "black"), labels = c("absent", "present")) +
my_gg_theme +
theme(
axis.text.x = element_blank()
)

ggplot(ml[value != 0], aes(x = fct_month, y = variable, group = variable, color = variable)) +
geom_line(size = 5) +
scale_y_discrete(expand = c(0, 1), breaks = NULL) +
my_gg_theme

ggplot(ml[value != 0], aes(x = fct_month, y = variable,
group = variable, colour = variable)) +
geom_line(size = 5) +
scale_y_discrete(expand = c(0, 1), breaks = NULL) +
#xlim(month.abb) +
coord_polar() +
theme_bw() + xlab(NULL) + ylab(NULL)

p1 = ggplot(full_metadata.sorted, aes(x = TreeCond, fill = occurence)) +
#geom_bar(position = position_dodge(width = 0.9), color = "black") +
geom_bar(color = "black") +
scale_fill_manual(
    values = c("both" = "black", "Nf" = "grey42", "Nd" = "grey72", "none" = "white"),
    labels = c("both" = "both spp.", "Nf" = "N.faginata", "Nd" = "N.ditissima", "none" = "neither spp."),
    breaks = c("both", "Nf", "Nd", "none")
) +
my_gg_theme +
labs(title = "A", x = "Crown dieback", y = "Trees (count)") +
theme(
plot.title = element_text(hjust = -0.175)
)


p4 = ggplot(full_metadata.sorted, aes(x = RaisedCanker, fill = occurence)) +
#geom_bar(position = position_dodge(width = 0.9), color = "black") +
geom_bar(color = "black") +
scale_fill_manual(
values = c("both" = "black", "Nf" = "grey42", "Nd" = "grey72", "none" = "white"),
labels = c("both" = "both spp.", "Nf" = "N.faginata", "Nd" = "N.ditissima", "none" = "neither spp."),
breaks = c("both", "Nf", "Nd", "none")
) +
my_gg_theme +
labs(title = "D", x = "Cankers", y = "Trees (count)")+
theme(
plot.title = element_text(hjust = -0.175)
)


p2 = ggplot(full_metadata.sorted, aes(y = Nf, x = TreeCond)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
my_gg_theme +
scale_y_continuous(breaks = c(0,1)) +
labs(title = "B", x = "Crown dieback", y = "N. faginata\noccurrence")+
theme(
plot.title = element_text(hjust = -0.325)
)


p3 = ggplot(full_metadata.sorted, aes(y = Nd, x = TreeCond)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
scale_y_continuous(breaks = c(0,1)) +
my_gg_theme +
labs(title = "C", x = "Crown dieback", y = "N. ditissima\noccurrence")+
theme(
plot.title = element_text(hjust = -0.325)
)


p5 = ggplot(full_metadata.sorted, aes(y = Nf, x = RaisedCanker)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
scale_y_continuous(breaks = c(0,1)) +
my_gg_theme +
labs(title = "E", x = "Cankers", y = "N. faginata\noccurrence")+
theme(
plot.title = element_text(hjust = -0.325)
)


p6 = ggplot(full_metadata.sorted, aes(y = Nd, x = RaisedCanker)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
scale_y_continuous(breaks = c(0,1)) +
my_gg_theme +
labs(title = "F", x = "Crown dieback", y = "N. ditissima\noccurrence")+
theme(
plot.title = element_text(hjust = -0.325)
)


pdf("prelim_figs/Nf_Nd_by_treeCond_cankers.pdf", width = 14, height = 6)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 3, widths = c(0.5,0.25,0.25))
dev.off()


#vertical
treeCond_n = (full_metadata.sorted %>% group_by(TreeCond) %>% summarize(n = n()))$n
canker_n = (full_metadata.sorted %>% group_by(RaisedCanker) %>% summarize(n = n()))$n

p1 = ggplot(full_metadata.sorted, aes(x = TreeCond, fill = occurence)) +
#geom_bar(position = position_dodge(width = 0.9), color = "black") +
geom_bar(position = "fill", color = "black", width = 0.9) +
annotate(geom = "text", x = c(0,1,2,3), y = rep(1.05, 4), label = treeCond_n, size = 5) +
scale_fill_manual(
values = c("both" = "black", "Nf" = "grey42", "Nd" = "grey72", "none" = "white"),
labels = c("both" = "both spp.", "Nf" = "N.faginata", "Nd" = "N.ditissima", "none" = "neither spp."),
breaks = c("both", "Nf", "Nd", "none"),
guide = F
) +
scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
my_gg_theme +
labs(title = "A", x = "Crown dieback", y = "Trees (proportion)") +
theme(
plot.title = element_text(hjust = -0.275),
#axis.title.x = element_blank(),
#legend.position = c(0.775, 0.76),
legend.background = element_rect(color = "black")
)


p2 = ggplot(full_metadata.sorted, aes(x = RaisedCanker, fill = occurence)) +
#geom_bar(position = position_dodge(width = 0.9), color = "black") +
geom_bar(position = "fill", color = "black", width = 0.9) +
annotate(geom = "text", x = c(0,1,2,3), y = rep(1.05, 4), label = canker_n, size = 5) +
scale_fill_manual(
values = c("both" = "black", "Nf" = "grey42", "Nd" = "grey72", "none" = "white"),
labels = c("both" = "both spp.", "Nf" = "N.faginata", "Nd" = "N.ditissima", "none" = "neither spp."),
breaks = c("none", "Nf", "Nd", "both")#,
#guide = F
) +
scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
my_gg_theme +
labs(title = "B", x = "Cankers", y = "Trees (count)")+
theme(
#plot.title = element_text(hjust = -0.175),
#axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text.y = element_blank()
)


p3 = ggplot(full_metadata.sorted, aes(y = Nf, x = TreeCond)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
my_gg_theme +
scale_y_continuous(breaks = c(0,1)) +
labs(title = "C", x = "", y = "N. faginata (0/1)")+
theme(
plot.title = element_text(hjust = -0.15),
axis.title.x = element_blank()
)


p5 = ggplot(full_metadata.sorted, aes(y = Nd, x = TreeCond)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
scale_y_continuous(breaks = c(0,1)) +
my_gg_theme +
labs(title = "E", x = "Crown dieback", y = "N. ditissima (0/1)")+
theme(
plot.title = element_text(hjust = -0.15)
)


p4 = ggplot(full_metadata.sorted, aes(y = Nf, x = RaisedCanker)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
scale_y_continuous(breaks = c(0,1)) +
my_gg_theme +
labs(title = "D", x = "", y = "N. faginata (0/1)")+
theme(
plot.title = element_text(hjust = -0.05),
axis.title.y = element_blank(),
axis.title.x = element_blank()
)


p6 = ggplot(full_metadata.sorted, aes(y = Nd, x = RaisedCanker)) +
geom_point(position = position_jitter(width = 0.1, height = 0)) +
stat_smooth(method = "glm", method.args = list(family = "binomial"), color = "black", linetype = 2, se = F) +
#coord_flip() +
scale_y_continuous(breaks = c(0,1)) +
my_gg_theme +
labs(title = "F", x = "Cankers", y = "N. ditissima (0/1))")+
theme(
plot.title = element_text(hjust = -0.05),
axis.title.y = element_blank()
)


#pdf("prelim_figs/Nf_Nd_by_treeCond_cankers.vertical.pdf", width = 10, height = 10)
#grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 2)
#dev.off()

pdf("prelim_figs/Nf_Nd_by_treeCond_cankers.proportion.pdf", width = 10, height = 4)
grid.arrange(p1,p2,ncol = 2, widths = c (0.44,0.56))
dev.off()


pdf("prelim_figs/Nf_Nd_by_treeCond_cankers.glm_logistic.vertical.pdf", width = 10, height = 8)
grid.arrange(p3,p4,p5,p6,p6ncol = 2)
dev.off()

