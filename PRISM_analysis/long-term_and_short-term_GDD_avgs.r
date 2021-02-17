require(tidyverse)

site_info = read.csv("sample_Data/site_info.csv")
sites_climate.seasonal_summaries = read.table(file = "sample_data/summarized_12_year_climate_dailys.txt", sep = "\t", header = T)

#For all sites except TSP1 use previous years data for short-term. TSP1 was sampled in mid-Dec so can use the data for the year of sampling, but all other sites were sampled in either January or spring (apr may jun) so should use previous year.
#OR ACTUALLY going to use YEAR OF sampling data for all sites except CW1 and ASH2 (these two were sampled in mid-January so doesn't make sense to use year of)

site_info$year_sampled[site_info$Site == "CW1"] = 2018
site_info$year_sampled[site_info$Site == "ASH2"] = 2018
colnames(sites_climate.seasonal_summaries)[1] = "year_sampled"

site_info.seasonal_summaries = left_join(
    site_info %>% select(Site, year_sampled),
    sites_climate.seasonal_summaries %>% select(-doy.grow_start, -doy.grow_end, -year_len),
    by = c("Site", "year_sampled")
)

#Just use average from 2009-2018 foe everything for long-term

site_means = sites_climate.seasonal_summaries %>%
    filter(year_sampled %in% seq(2009, 2018)) %>%
    group_by(Site, season) %>%
    summarize(
        HDD4.mean = mean(HDD4.sum),
        freezeThaw.mean = mean(freezeThaw.sum),
        HDD4.sd = sd(HDD4.sum),
        freezeThaw.sd = sd(freezeThaw.sum),
        ppt.mean = mean(ppt.sum),
        ppt.sd = sd(ppt.sum),
        length.mean = mean(length),
        length.sd = sd(length)
    )

site_means.doy_start_end = sites_climate.seasonal_summaries %>%
filter(year_sampled %in% seq(2009, 2018) & season == "growing") %>%
group_by(Site) %>%
summarize(
doy.grow_start.mean = mean(doy.grow_start),
doy.grow_end.mean = mean(doy.grow_end)
)


site_info.seasonal_summaries = left_join(site_info.seasonal_summaries, site_means, by = c("Site", "season"))
#pivot_wider for one row per site
site_info.seasonal_summaries.wide = site_info.seasonal_summaries %>%
    pivot_wider(names_from = season, values_from = c(-Site, -year_sampled, -season)) %>%
    select(-year_sampled)

write.table(site_info.seasonal_summaries.wide, file = "sample_data/site_info.GDD.txt", sep = "\t", quote = F, row.names = F)
