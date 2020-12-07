require(dplyr)
require(ggplot2)
source("~/ggplot_theme.txt")
merged_93 = read.table("LULU/merged_ASVs_93.txt", header = T)

mean(merged_93$seq_sim)
median(merged_93$seq_sim)
sd(merged_93$seq_sim)
range(merged_93$seq_sim)
quantile(merged_93$seq_sim)

mean(merged_93$n_samples_child)
median(merged_93$n_samples_child)
sd(merged_93$n_samples_child)
range(merged_93$n_samples_child)
quantile(merged_93$n_samples_child)

ggplot(merged_93, aes(x = seq_sim)) +
 geom_histogram() +
my_gg_theme

ggplot(merged_93, aes(x = n_samples_child)) +
geom_histogram() +
my_gg_theme
