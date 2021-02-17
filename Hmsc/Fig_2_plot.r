
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
require(RColorBrewer)
source("~/ggplot_theme.txt")
VP.vals.support = read.table(file = "HMSC/Nf_Nd_variance_partitioning.bin.spatial.VP_table.txt", header = T, sep = "\t")

VP.vals.support$ASV = factor(VP.vals.support$ASV, levels = c("N. faginata", "N. ditissima"))

p1 = ggplot(VP.vals.support %>% filter(variable != "intercept"), aes(ASV, variable, size = R2, fill = betaMean, color = P.val)) +
geom_point(shape = 21, stroke = 1.5) +
scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0) +
#scale_fill_gradientn(colours = c("#2c7bb6", "white", "#d7191c"),
#values = scales::rescale(c(-0.25, -.15,-0.1,-0.05,-0.025, 0 ,0.025, 0.05, 0.1,0.25)))  +
scale_color_manual(values = c("P<0.05" = "black", "n.s." = "white"), labels = c("P<0.05" = "P>0.95", "n.s." = "n.s.")) +
guides(color = guide_legend(override.aes = list(fill = "grey"))) +
scale_size_continuous(range= c(5,25), breaks = c(0,0.05,0.15), limits = c(0,0.2))+
my_gg_theme +
guides(
    color = guide_legend(override.aes = list(size = 5, shape = 21, fill = "dark grey")),
    size = guide_legend(override.aes = list(shape = 16, color = "dark grey"))
) +
labs(
x = "",
y = "",
size = expression(paste("R"^2)),
color = "Support",
fill = "Slope"
#title = "HMSC variance partitioning"
) +
theme(
legend.title = element_text(size = 20),
axis.text = element_text(size = 18),
axis.text.x = element_text(face = "italic"),
legend.key = element_rect(fill = "white")
)

pdf("HMSC/Fig_2.pdf", width = 8.5, height = 7.5)
p1
dev.off()


