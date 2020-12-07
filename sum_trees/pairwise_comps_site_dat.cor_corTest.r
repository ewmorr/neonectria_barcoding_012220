require(tidyverse)
require(Hmisc)
#source("~/ggplot_theme.txt")

site_info = read.csv("sample_data/site_info.csv")
site_means = read.table("sample_data/BBD_survey_transect_data.site_mean.txt", header = T)
site_climate = read.table("sample_data/site_info.GDD.txt", header = T)
site_climate.normals = read.table("sample_data/sites_climate.txt", header = T)

full_dat = full_join(site_info, site_climate, by = "Site") %>%
    full_join(., site_means, by = "Site") %>%
    full_join(., site_climate.normals, by = c("Site", "lat", "lon")
)


pairwise_df = data.frame(
lat = full_dat$lat,
lon = full_dat$lon,
elev = full_dat$elev_m,

GDD4.summer = full_dat$HDD4.mean_growing,
GDD4.winter = full_dat$HDD4.mean_nongrowing,
freeze_thaw.summer = full_dat$freezeThaw.mean_growing,
freeze_thaw.winter = full_dat$freezeThaw.mean_nongrowing,
precip.summer = full_dat$ppt.mean_growing,
precip.winter = full_dat$ppt.mean_nongrowing,
MAT = full_dat$MAT,
tmin = full_dat$tmin,
tmax = full_dat$tmax
)

pairwise_df = data.frame(
    lat = full_dat$lat,
    lon = full_dat$lon,
    elev = full_dat$elev_m,

    GDD4.summer = full_dat$HDD4.mean_growing,
    GDD4.winter = full_dat$HDD4.mean_nongrowing,
    freeze_thaw.summer = full_dat$freezeThaw.mean_growing,
    freeze_thaw.winter = full_dat$freezeThaw.mean_nongrowing,
    precip.summer = full_dat$ppt.mean_growing,
    precip.winter = full_dat$ppt.mean_nongrowing,
    infec_dur = full_dat$duration_infection,
    dbh = full_dat$Site_mean.dbh,
    neo_fruit = full_dat$Site_mean.NeoFruiting,
    cankers = full_dat$Site_mean.RaisedCanker,
#    tarry_spot = full_dat$Site_mean.TarrySpots,
    wax = full_dat$Site_mean.Wax,
#    xyloc = full_dat$Site_mean.Xylococcus,
    tree_cond = full_dat$Site_mean.TreeCond
)

for(i in 1:ncol(pairwise_df)){
    print(colnames(pairwise_df)[i])
    print(shapiro.test(pairwise_df[,i]))
}

corr.matrix = rcorr(as.matrix(pairwise_df), type = "pearson")
corr.matrix.spear = rcorr(as.matrix(pairwise_df), type = "spearman")

write.table(data.frame(row.var = rownames(corr.matrix$r), corr.matrix$r), "tables/site_level_cor.r.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(data.frame(row.var = rownames(corr.matrix$P), corr.matrix$P), "tables/site_level_cor.P.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(data.frame(row.var = rownames(corr.matrix.spear$r), corr.matrix.spear$r), "tables/site_level_cor_spear.r.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(data.frame(row.var = rownames(corr.matrix.spear$P), corr.matrix.spear$P), "tables/site_level_cor_spear.P.txt", sep = "\t", quote = F, row.names = F, col.names = T)



cor.test(pairwise_df$lat, pairwise_df$elev, method = "pearson")
cor.test(pairwise_df$lat, pairwise_df$GDD4.summer, method = "pearson")
cor.test(pairwise_df$lat, pairwise_df$GDD4.winter, method = "pearson")
cor.test(pairwise_df$lat, pairwise_df$precip.summer, method = "pearson")
cor.test(pairwise_df$lon, pairwise_df$elev, method = "pearson")

#P.vals and r are the same from cor.test versus rcorr



