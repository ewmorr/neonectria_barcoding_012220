require(devtools)
require(reshape2)
require(dplyr)
require(raster)

require(hutilscpp)
#install_github(repo = "ropensci/prism")
require(prism)

#prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_tmin")
#get_prism_dailys(type = "tmin", minDate="2007-01-01", maxDate="2018-12-31", keepZip = TRUE)
#prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_tmax")
#get_prism_dailys(type = "tmax", minDate="2007-01-01", maxDate="2018-12-31", keepZip = TRUE)
#prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_ppt")
#get_prism_dailys(type = "ppt", minDate="2007-01-01", maxDate="2018-12-31", keepZip = TRUE)

#options(prism.path = "/Users/ericmorrison/PRISM_data")

site_coords = read.table("sample_data/site_coords.txt", header = T)
sites_climate = site_coords


#####
#this processes the first file of daily data, i.e. the first day#
new_file<-1#this number corresponds to the row of the file of interest
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- reshape2::melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns
prism_coords = m.df[,1:2]

################
#nearest coords#

nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
m.df[nrst_coords$pos,]

sites_climate$tmax = m.df[nrst_coords$pos,colnames(m.df) == "value"]

coord_day = regmatches(names(RS),regexpr("[0-9]{8,}",names(RS), perl = T)) %>% as.numeric
sites_climate$day = coord_day

###########################
#Start of loop to loop the files for daily values

prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_tmax")

num_files = ls_prism_data()$files %>% length


sites_climate_list = list()

for(i in 1:num_files){
    
    print(i)
    
    prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_tmax")

    new_file<-i#this number corresponds to the row of the file of interest
    RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    
    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    prism_coords = m.df[,1:2]
    
    #nearest coords#
    
    nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
    m.df[nrst_coords$pos,]
    
    #Get tmax and date
    sites_climate$tmax = m.df[nrst_coords$pos,colnames(m.df) == "value"]
    coord_day = regmatches(names(RS),regexpr("[0-9]{8,}",names(RS), perl = T)) %>% as.numeric
    sites_climate$day = coord_day
    
    #reset to get tmin
    prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_tmin")
    
    new_file<-i#this number corresponds to the row of the file of interest
    RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    
    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    prism_coords = m.df[,1:2]
    
    #nearest coords#
    
    nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
    m.df[nrst_coords$pos,]
    
    #Get tmax and date
    sites_climate$tmin = m.df[nrst_coords$pos,colnames(m.df) == "value"]

    #Save to list indexed by file number
    sites_climate_list[[i]] = sites_climate

}

#coerce list to df
sites_climate_list.df = do.call(rbind.data.frame, sites_climate_list)

saveRDS(sites_climate_list, "intermediate_RDS/sites_daily_tmin_tmax_20072018.list.rds")
saveRDS(sites_climate_list.df, "intermediate_RDS/sites_daily_tmin_tmax_20072018.df.rds")


####################
#PPT data done separartely due to hd storage constraints

site_coords = read.table("sample_data/site_coords.txt", header = T)
sites_climate = site_coords

prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_ppt")

num_files = ls_prism_data()$files %>% length


sites_climate_list = list()

for(i in 1:num_files){
    
    print(i)
    
    #prism_set_dl_dir("/Users/ericmorrison/PRISM_dailys_neonectria_sampling_ppt")
    
    new_file<-i#this number corresponds to the row of the file of interest
    RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    
    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    prism_coords = m.df[,1:2]
    
    #nearest coords#
    
    nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
    m.df[nrst_coords$pos,]
    
    #Get tmax and date
    sites_climate$ppt = m.df[nrst_coords$pos,colnames(m.df) == "value"]
    coord_day = regmatches(names(RS),regexpr("[0-9]{8,}",names(RS), perl = T)) %>% as.numeric
    sites_climate$day = coord_day
    #Save to list indexed by file number
    sites_climate_list[[i]] = sites_climate
    
}

#coerce list to df
sites_climate_list.df = do.call(rbind.data.frame, sites_climate_list)

saveRDS(sites_climate_list, "intermediate_RDS/sites_daily_ppt_20072018.list.rds")
saveRDS(sites_climate_list.df, "intermediate_RDS/sites_daily_ppt_20072018.df.rds")


