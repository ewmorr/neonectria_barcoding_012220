require(tidyverse)
require(GGally)
source("~/ggplot_theme.txt")

source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

full_dat = left_join(
data.frame(sample = Nf_v_Nd.bin$sample),
full_metadata,
by = "sample"
)


pairwise_df = data.frame(
    Site = full_dat$Site,
    dbh = full_dat$dbh,
    neo_fruit = full_dat$NeoFruiting,
    cankers = full_dat$RaisedCanker,
#    tarry_spot = full_dat$Site_mean.TarrySpots,
    wax = full_dat$Wax,
#    xyloc = full_dat$Site_mean.Xylococcus,
    tree_cond = full_dat$TreeCond
)

ggplot(pairwise_df, aes(x = neo_fruit)) +
geom_histogram() +
facet_wrap(~Site, ncol = 5) +
my_gg_theme


pairwise_df %>% group_by(Site) %>% summarize(median(neo_fruit), mean(neo_fruit))


ggplot(pairwise_df, aes(x = cankers)) +
geom_histogram() +
facet_wrap(~Site) +
my_gg_theme

pairwise_df %>% group_by(Site) %>% summarize(median(cankers), mean(cankers))


ggplot(pairwise_df, aes(x = wax)) +
geom_histogram() +
facet_wrap(~Site) +
my_gg_theme

pairwise_df %>% group_by(Site) %>% summarize(median(wax), mean(wax))


ggplot(pairwise_df, aes(x = tree_cond)) +
geom_histogram() +
facet_wrap(~Site) +
my_gg_theme

pairwise_df %>% group_by(Site) %>% summarize(median(tree_cond), mean(tree_cond))

