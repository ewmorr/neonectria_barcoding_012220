
source("~/repo/neonectria_barcoding_012220/ecol/read_ASV_dat.LULU_tab.r")


BioSample.metadata = full_metadata %>% select(c("sample", "Site_tree_plug", "dateSamp", "state_prov", "county", "lat", "lon", "elev_m", "dbh", "NeoFruiting", "RaisedCanker", "Wax", "TreeCond", "Site", "Tree", "Plug"))

#geo_loc_name edits
BioSample.metadata$state_prov = paste("USA:", BioSample.metadata$state_prov, sep = " ")
BioSample.metadata$state_prov = paste(BioSample.metadata$state_prov, ",", sep = "")
BioSample.metadata$county = paste(BioSample.metadata$county, "county,", sep = " ")
BioSample.metadata$geo_loc_name = paste(BioSample.metadata$state_prov, BioSample.metadata$county)

#description column edits
BioSample.metadata$Site = paste("Site:", BioSample.metadata$Site)
BioSample.metadata$Tree = paste("Tree:", BioSample.metadata$Tree)
BioSample.metadata$Plug = paste("Plug:", BioSample.metadata$Plug)
BioSample.metadata$description = paste(BioSample.metadata$Site, BioSample.metadata$Tree)
BioSample.metadata$description = paste(BioSample.metadata$description, BioSample.metadata$Plug)

BioSample.metadata$geo_loc_name_site_tree_plug = paste(BioSample.metadata$geo_loc_name, BioSample.metadata$description)

#lat lon edits
BioSample.metadata$lat = paste(round(BioSample.metadata$lat, 4), "N")
BioSample.metadata$lon = paste(sub("-", "", round(BioSample.metadata$lon, 4)), "W")
BioSample.metadata$lat_lon = paste(BioSample.metadata$lat, BioSample.metadata$lon)


#Edit misc. parameters -- DBH, Wax, Neonectria, TreeCond, Cankers
BioSample.metadata$dbh = paste("DBH (cm):", BioSample.metadata$dbh)
BioSample.metadata$Wax = paste("Beech scale density:", BioSample.metadata$Wax)
BioSample.metadata$NeoFruiting = paste("Neonectria perithecia density:", BioSample.metadata$NeoFruiting)
BioSample.metadata$RaisedCanker = paste("Canker density:", BioSample.metadata$RaisedCanker)
BioSample.metadata$TreeCond = paste("Crown dieback:", BioSample.metadata$TreeCond)


BioSample.metadata.columns = BioSample.metadata %>% select(c(-state_prov, -county, -lat, -lon))

write.table(BioSample.metadata.columns, "BioSample_metadata.txt", row.names = F, sep = "\t",
quote = F)


#Trying with site as biosample

BioSample.site.metadata = full_metadata %>% select(c("Site", "dateSamp", "state_prov", "county", "lat", "lon", "elev_m"))
BioSample.site.metadata.uniqs = BioSample.site.metadata %>% unique %>% filter(!is.na(dateSamp))

BioSample.site.metadata.uniqs$state_prov = paste("USA:", BioSample.site.metadata.uniqs$state_prov, sep = " ")
BioSample.site.metadata.uniqs$state_prov = paste(BioSample.site.metadata.uniqs$state_prov, ",", sep = "")
BioSample.site.metadata.uniqs$county = paste(BioSample.site.metadata.uniqs$county, "county", sep = " ")
BioSample.site.metadata.uniqs$geo_loc_name = paste(BioSample.site.metadata.uniqs$state_prov, BioSample.site.metadata.uniqs$county)

BioSample.site.metadata.uniqs$lat = paste(round(BioSample.site.metadata.uniqs$lat, 4), "N")
BioSample.site.metadata.uniqs$lon = paste(sub("-", "", round(BioSample.site.metadata.uniqs$lon, 4)), "W")
BioSample.site.metadata.uniqs$lat_lon = paste(BioSample.site.metadata.uniqs$lat, BioSample.site.metadata.uniqs$lon)

write.table(BioSample.site.metadata.uniqs, "BioSample_metadata.site.txt", row.names = F, sep = "\t",
quote = F)



