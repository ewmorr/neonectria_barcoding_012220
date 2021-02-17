require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)
set.seed(6)
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

#and calc duration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

#######
#PLOTS#
######
#NMDS#

##########################
#1000 seqs per sample min#
##########################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
asv_tab.gt1K.rare = asv_tab.gt1K.rare[,colSums(asv_tab.gt1K.rare > 0) > 0]

full_metadata.sorted = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare)),
full_metadata
) %>% left_join(., Nf_v_Nd.bin)

#score ASVs by frequency
#score sampes by presence absence of species occurence
#run adonis by presence absence and extract R2

asv_tab.gt1K.rare.incidence = colSums(asv_tab.gt1K.rare > 0)/nrow(asv_tab.gt1K.rare)
plot(sort(asv_tab.gt1K.rare.incidence))

asv_tab.gt1K.rare.geTenPercInc = asv_tab.gt1K.rare[,asv_tab.gt1K.rare.incidence >= 0.1]

asv_adonis_coefs = data.frame(
    R2 = vector(mode = "numeric", length = ncol(asv_tab.gt1K.rare.geTenPercInc)),
    p.val = vector(mode = "numeric", length = ncol(asv_tab.gt1K.rare.geTenPercInc)),
    ASV = vector(mode = "character", length = ncol(asv_tab.gt1K.rare.geTenPercInc)),
    stringsAsFactors = F
)

foo.adonis = adonis(asv_tab.gt1K.rare ~ full_metadata.sorted$Nf)

for(i in 1:nrow(asv_adonis_coefs)){
    print(i)
    asv_nam = colnames(asv_tab.gt1K.rare.geTenPercInc)[i]
    temp.adonis = adonis(log10(asv_tab.gt1K.rare[,colnames(asv_tab.gt1K.rare) != asv_nam]+1) ~ (asv_tab.gt1K.rare.geTenPercInc[,i] > 0) )
    asv_adonis_coefs$R2[i] = temp.adonis$aov.tab[1,5]
    asv_adonis_coefs$p.val[i] = temp.adonis$aov.tab[1,6]
    asv_adonis_coefs$ASV[i] = asv_nam
}

asv_adonis_coefs.incidence = left_join(
    asv_adonis_coefs,
    data.frame(ASV = names(asv_tab.gt1K.rare.incidence), freq = asv_tab.gt1K.rare.incidence)

)

asv_adonis_coefs.incidence$color = vector(mode = "character", length = nrow(asv_adonis_coefs.incidence))

for(i in 1:nrow(asv_adonis_coefs.incidence)){
    if(asv_adonis_coefs.incidence$ASV[i] %in% Nf_asvs){
        asv_adonis_coefs.incidence$color[i] = "N.faginata"
    }else if(asv_adonis_coefs.incidence$ASV[i] %in% Nd_asvs){
        asv_adonis_coefs.incidence$color[i] = "N.ditissima"
    }else{
        asv_adonis_coefs.incidence$color[i] = "other spp."
    }
}

#lm of R2 ~ freq
summary(lm(R2 ~ freq, data = asv_adonis_coefs.incidence))

#plot
p2 = ggplot(asv_adonis_coefs.incidence, aes(freq, R2, fill = color, size = color)) +
geom_point(alpha = 0.75, shape = 21) +
geom_smooth(method = "lm", se = F, color = "black", linetype = 2) +
annotate(geom = "text", label = expression(paste("R"^2, "= 0.165, P = 0.001")), x = 0.35, y = 0.12, size = 6) +
labs(x = "sample incidence", y = expression(paste("adnois R"^2)), title = "B") +
scale_fill_manual(
    values = c("other spp." = "black", "N.ditissima" = "light grey", "N.faginata" = "grey54"),
    labels = c("other spp." = "other spp.", "N.ditissima" = expression(italic("N.ditissima")), "N.faginata" = expression(italic("N.faginata")))
) +
scale_size_manual(values = c("other spp." = 2, "N.ditissima" = 4, "N.faginata" = 4), guide = F) +
scale_y_continuous(breaks = c(0,0.05,0.1)) +
guides(fill = guide_legend(override.aes = list(shape = 22, size = 7.5, linetype = 0))) +
my_gg_theme +
theme(
axis.text.x = element_blank(),
axis.title.y = element_blank(),
plot.title = element_text(hjust = -0.12)
)


p1 = ggplot(asv_adonis_coefs.incidence, aes(reorder(ASV, -R2), R2, fill = color)) +
geom_bar(stat = "identity",color = "black") +
labs(x = "rank", y = expression(paste("adnois R"^2)), title = "A") +
scale_fill_manual(values = c("other spp." = "grey24", "N.ditissima" = "light grey", "N.faginata" = "grey54"), guide = F) +
scale_y_continuous(breaks = c(0,0.05,0.1)) +
my_gg_theme +
theme(
axis.text.x = element_blank(),
plot.title = element_text(hjust = -0.175)
)

require(gridExtra)

pdf("prelim_figs/adonis_rank.pdf", width = 12, height = 4)
grid.arrange(p1,p2,ncol = 2, widths = c(0.5,0.5))
dev.off()

asv_adonis_coefs.incidence[order(asv_adonis_coefs.incidence$R2),]

top_asvs = c("ASV_12", "ASV_10", "ASV_27", "ASV_5", "ASV_1")

asv_tax[top_asvs,]

write.table(asv_adonis_coefs.incidence, "adonis_rank.txt", row.names = F, col.names = T, sep = "\t", quote = F)

