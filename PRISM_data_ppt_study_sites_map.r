require(devtools)
require(reshape2)
require(dplyr)
require(raster)

install_github(repo = "ropensci/prism")
require(prism)

prism_set_dl_dir("/Users/ericmorrison/PRISM_data")
#options(prism.path = "/Users/ericmorrison/PRISM_data")
get_prism_normals(type = 'ppt', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmean', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmax', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmin', resolution = '4km', annual = T, keepZip = TRUE)


#Convert raster to point data

new_file<-1#this number corresponds to the row of the file of interest
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") ##assign projection info

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns


site_coords = read.table("sample_data/site_coords.txt", header = T)

minLat=min(site_coords$lat)-2
maxLat=max(site_coords$lat)+2
minLon=min(site_coords$lon)-2
maxLon=max(site_coords$lon)+2


require(ggplot2)
require(ggmap)

#ggplot()+
#geom_raster(data=m.df, aes(x=lon, y=lat, fill=value))+
#my_gg_theme +

m.df.study_area<-m.df%>%filter(minLat < lat, lat < maxLat, minLon < lon, lon <maxLon)%>%
mutate(ppt = value)%>%
dplyr::select(-value)

#ggplot()+
#geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=ppt))+
#my_gg_theme +


p1 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=ppt))+
geom_point(data=site_coords, aes(x=lon, y = lat), color = "black") +
my_gg_theme +
theme(
legend.title = element_text(size = 20)
) +
scale_fill_gradient2("precip.\n(mm)", high='darkslateblue',mid='lightblue',low = 'red',midpoint = 1500)

plot_height = (maxLat-minLat)/2
plot_width = (maxLon-minLon)/2

pdf("site_ppt_map.pdf", height = plot_height, width = plot_width)
p1
dev.off()
