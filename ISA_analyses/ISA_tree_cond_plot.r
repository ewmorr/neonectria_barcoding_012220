
require(tidyverse)
require(vegan)
require(indicspecies)

source("~/ggplot_theme.txt")

tree_cond_table = read.table("ISA_spp_tables/TreeCond.txt", header = T)
tree_cond_table$ASV = rownames(tree_cond_table)
tree_cond_table.sig = filter(tree_cond_table, p.value < 0.05)

tree_cond_table.sig[
    with(tree_cond_table.sig, order("s.0", "s.1", "s.2", "s.3"))
]

tree_cond_table.sig.long = pivot_longer(tree_cond_table.sig,
    cols = starts_with("s."),
    names_to = c("prepend", "cat"),
    names_sep = "\\.",
    values_to = "indicator"
)

tree_cond_table.sig.long = filter(tree_cond_table.sig.long, indicator > 0)
tree_cond_table.sig.long$cat = tree_cond_table.sig.long$cat %>% as.numeric
tree_cond_table.sig.long = full_join(
tree_cond_table.sig.long,
tree_cond_table.sig.long %>% group_by(ASV) %>% summarize(sum.cat = sum(cat))
)

p = ggplot(tree_cond_table.sig.long, aes(cat, reorder(ASV, cat), group = ASV)) +
geom_point() +
geom_line(size = 2) +
labs(x = "Crown dieback class", y = "Indicator ASVs") +
my_gg_theme +
theme(
axis.text.y = element_blank()
)

pdf("ISA_spp_tables/ISA_crown_dieback_ASV_progression.pdf", width = 6, height = 6)
p
dev.off()
