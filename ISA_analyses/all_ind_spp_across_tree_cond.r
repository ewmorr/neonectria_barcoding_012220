
require(tidyverse)
require(vegan)
require(indicspecies)

source("~/ggplot_theme.txt")

tree_cond_table = read.table("ISA_spp_tables/TreeCond.txt", header = T)
tree_cond_table$ASV = rownames(tree_cond_table)
tree_cond_table.is_na = filter(tree_cond_table, !is.na(p.value))

tree_cond_table.is_na[
    with(tree_cond_table.is_na, order("s.0", "s.1", "s.2", "s.3"))
]

for(i in 1:nrow(tree_cond_table.is_na)){
    if(tree_cond_table.is_na$p.value[i] < 0.05){
        tree_cond_table.is_na$sig[i] = "P<0.05"
    }else{
        tree_cond_table.is_na$sig[i] = "n.s."
    }
}
tree_cond_table.is_na.long = pivot_longer(tree_cond_table.is_na,
    cols = starts_with("s."),
    names_to = c("prepend", "cat"),
    names_sep = "\\.",
    values_to = "indicator"
)

tree_cond_table.is_na.long = filter(tree_cond_table.is_na.long, indicator > 0)
tree_cond_table.is_na.long$cat = tree_cond_table.is_na.long$cat %>% as.numeric

tree_cond_table.is_na.long = full_join(
tree_cond_table.is_na.long,
tree_cond_table.is_na.long %>% group_by(ASV) %>% summarize(sum.cat = sum(cat))
)

p = ggplot(tree_cond_table.is_na.long, aes(cat, reorder(ASV, cat), group = ASV, color = sig, alpha = stat)) +
geom_point(size = 1) +
geom_line(size = 1) +
labs(x = "Crown dieback class", y = "Indicator ASVs") +
scale_color_manual(values = c("black", "red")) +
my_gg_theme +
theme(
axis.text.y = element_blank()
)

pdf("ISA_spp_tables/ISA_crown_dieback_ASV_progression_all_ASVs_colored_by_sig.pdf", width = 6, height = 6)
p
dev.off()
