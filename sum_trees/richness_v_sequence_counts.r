require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)
source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
#this pulls in objects:
#asv_tab
#asv_tax
#id_bench_map

#joins metadata files to get metadata object:
#full_metadata

#creates negatives and controls only asv_tab (long format):
#asv_tab.negatives.long

#creates new object with lowest informative taxonomic level and a character asv_tax with unknown instead of NA
#asv_informative_taxa
#asv_tax.char

#and calcuduration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

asv_tab.gt1k = asv_tab[,colSums(asv_tab) >= 1000]

richness = data.frame(richness = colSums(asv_tab.gt1k > 0), total_seqs = colSums(asv_tab.gt1k))

mod1 = lm(richness ~ (total_seqs), data = richness)
mod2 = lm(richness ~ log10(total_seqs), data = richness)
mod3 = lm(richness ~ sqrt(total_seqs), data = richness)
mod4 = lm(richness ~ poly(total_seqs, 2), data = richness)

anova(mod1, mod2)
anova(mod1, mod3)
anova(mod1, mod4)

qqnorm(residuals(mod1))
plot(residuals(mod1))

qqnorm(residuals(mod4))
plot(residuals(mod4))

ggplot(richness, aes(total_seqs, richness)) +
geom_point()+
geom_smooth(method = "lm", formula = y ~ poly(x,2)) +
my_gg_theme

ggplot(richness, aes(total_seqs, richness)) +
geom_point()+
geom_smooth(method = "lm") +
my_gg_theme

ggplot(richness, aes(log10(total_seqs), richness)) +
geom_point()+
geom_smooth(method = "lm") +
my_gg_theme

ggplot(richness, aes(sqrt(total_seqs), richness)) +
geom_point()+
geom_smooth(method = "lm") +
my_gg_theme

p.val = summary(mod3)[[4]][2,4]
r2 = summary(mod3)[[8]]

p1 = ggplot(richness, aes(total_seqs, richness)) +
geom_point()+
geom_smooth(method = "lm") +
scale_x_continuous(
    trans = "sqrt",
    labels = fancy_scientific, breaks = c(10^3, 10^4, 5*10^4, 10^5, 2*10^5, 3*10^5)) +
labs(title = paste("R2 =", signif(r2,3), ", P =",signif(p.val,2)),
x = "total sequences (sqrt)") +
my_gg_theme

pdf("prelim_figs/richness_vs_seqs.tree_sum.pdf", width = 8, height = 6)
p1
dev.off()

