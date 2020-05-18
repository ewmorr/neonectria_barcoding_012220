require(devtools)
require(reshape2)
require(dplyr)
require(raster)

require(hutilscpp)
#install_github(repo = "ropensci/prism")
require(prism)

prism_set_dl_dir("/Users/ericmorrison/PRISM_data")
#options(prism.path = "/Users/ericmorrison/PRISM_data")
get_prism_normals(type = 'ppt', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmean', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmax', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmin', resolution = '4km', annual = T, keepZip = TRUE)


site_coords = read.table("sample_data/site_coords.txt", header = T)
sites_climate = site_coords

#####
#ppt#
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

sites_climate$ppt = m.df[nrst_coords$pos,colnames(m.df) == "value"]

#####
#tmax#

new_file<-2#this number corresponds to the row of the file of interest
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- reshape2::melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns

################
#nearest coords#

nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
m.df[nrst_coords$pos,]

sites_climate$tmax = m.df[nrst_coords$pos,colnames(m.df) == "value"]

#####
#MAT#

new_file<-3#this number corresponds to the row of the file of interest
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- reshape2::melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns

################
#nearest coords#

nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
m.df[nrst_coords$pos,]

sites_climate$MAT = m.df[nrst_coords$pos,colnames(m.df) == "value"]

#####
#tmin#

new_file<-4#this number corresponds to the row of the file of interest
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- reshape2::melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns

################
#nearest coords#

nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
m.df[nrst_coords$pos,]

sites_climate$tmin = m.df[nrst_coords$pos,colnames(m.df) == "value"]

#######
#write#

write.table(sites_climate, file = "sample_data/sites_climate.txt", quote = F, sep = "\t", row.names = F)
