require(devtools)
require(reshape2)
require(dplyr)
require(raster)

#install_github(repo = "ropensci/prism")
require(prism)

prism_set_dl_dir("/Users/ericmorrison/PRISM_data")
#options(prism.path = "/Users/ericmorrison/PRISM_data")
get_prism_normals(type = 'ppt', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmax', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmean', resolution = '4km', annual = T, keepZip = TRUE)
get_prism_normals(type = 'tmin', resolution = '4km', annual = T, keepZip = TRUE)


#Convert raster to point data

new_file<-3#this number corresponds to the row of the file of interest
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

m.df.study_area<-m.df%>%filter(minLat < lat, lat < maxLat, minLon < lon, lon <maxLon)%>%
mutate(ppt = value)%>%
dplyr::select(-value)

p1 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=ppt))+
geom_point(data=site_coords, aes(x=lon, y = lat), color = "black") +
my_gg_theme +
theme(
legend.title = element_text(size = 20)
) +
scale_fill_gradient2("MAT", low='darkslateblue',mid='lightblue',high = 'red',midpoint=10)

plot_height = (maxLat-minLat)/2
plot_width = (maxLon-minLon)/2

pdf("PRISM_maps/site_MAT_map.pdf", height = plot_height, width = plot_width)
print(p1)
dev.off()

#big text

p1 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=ppt))+
geom_point(data=site_coords, aes(x=lon, y = lat), color = "black", size = 5) +
my_gg_theme +
labs(y = NULL, x = NULL) +
theme(
legend.title = element_text(size = 35),
legend.text = element_text(size = 30),
axis.text = element_blank(),
axis.ticks = element_blank(),
legend.position = c(0.825,0.225)
) +
scale_fill_gradient2(expression("MAT ("*degree*C*")"), low='darkslateblue',mid='lightblue',high = 'red',midpoint=10)

plot_height = (maxLat-minLat)/2
plot_width = (maxLon-minLon)/2

pdf("PRISM_maps/site_MAT_map.big_text.pdf", height = plot_height, width = plot_width)
print(p1)
dev.off()

