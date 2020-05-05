require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)

source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

#1000 seqs per sample table
#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
asv_tab.gt1K.rare = asv_tab.gt1K.rare[,colSums(asv_tab.gt1K.rare) > 0]

##########################
#within site comm. cmomp.#

site_levels = full_metadata$Site %>% unique

######
#NMDS#
nmds_list = vector("list", 10)

for(i in 1:10){
    row_names = full_metadata %>% filter(Site == site_levels[i]) %>% select(sample)
    temp_tab = asv_tab.gt1K.rare[rownames(asv_tab.gt1K.rare) %in% row_names$sample, ]
    temp_tab = temp_tab[,colSums(temp_tab) > 0]
    temp_nmds = metaMDS(log10(temp_tab +1), autotransform = F, k = 2, zerodist = "add")
    nmds_list[[i]] = temp_nmds
}

nmds_list.points = vector("list", 10)

for(i in 1:10){
    nmds_list.points[[i]] = nmds_list[[i]]$points
}

nmds_list.points.df = do.call("rbind", nmds_list.points)

nmds_list.points.df.metadata = left_join(
data.frame(sample = rownames(nmds_list.points.df), nmds_list.points.df),
full_metadata,
by = "sample"
)

######
#PCOA#

#nmds_list = vector("list", 10)
#
#for(i in 1:10){
#    row_names = full_metadata %>% filter(Site == site_levels[i]) %>% select(sample)
#    temp_tab = asv_tab.gt1K.rare[rownames(asv_tab.gt1K.rare) %in% row_names$sample, ]
#    temp_tab = temp_tab[,colSums(temp_tab) > 0]
#    temp_nmds = capscale(log10(temp_tab +1) ~ 1, method = "bray")
#    nmds_list[[i]] = temp_nmds
#}
#
#nmds_list.points = vector("list", 10)
#
#for(i in 1:10){
#    nmds_list.points[[i]] = nmds_list[[i]]$CA$u[,1:2]
#}
#
#nmds_list.points.df = do.call("rbind", nmds_list.points)
#
#nmds_list.points.df.metadata = left_join(
#data.frame(sample = rownames(nmds_list.points.df), nmds_list.points.df),
#full_metadata,
#by = "sample"
#)

#Run GAMs

surf_vars = nmds_list.points.df.metadata %>%
select("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "Xylococcus", "TreeCond", "sample", "Site")

surf_list = vector("list")
surf_list.anova = vector("list")

surf_names = c("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "Xylococcus", "TreeCond")

for(i in 1:10){
    print(i)
    surf_vars.temp = surf_vars %>% filter(Site == site_levels[i])

    for(u in 1:7){
        k = (surf_vars.temp[,u] %>% unique %>% length)
        if(k < 3){
            next
        }
        surf_list[[site_levels[i]]][[surf_names[u]]] = ordisurf(nmds_list[[i]] ~ surf_vars.temp[,u], k = k)
        surf_list.anova[[site_levels[i]]][[surf_names[u]]] = anova(surf_list[[site_levels[i]]][[surf_names[u]]])
    }
}

plot_list = vector("list")

#pull Site names from list
site_levels = names(surf_list)

#make plots
for(i in 1:length(site_levels)){
    #get NMDS points and metadata for site
    surf_names = names(surf_list[[ site_levels[[i]] ]])
    temp_points = nmds_list.points.df.metadata %>% filter(Site == site_levels[i])
    
    for(u in 1:length(surf_list[[i]])){
        #extract plottable from surf
        surf.grid = expand.grid(x = surf_list[[site_levels[i]]][[surf_names[u]]]$grid$x, y = surf_list[[site_levels[i]]][[surf_names[u]]]$grid$y)
        surf.grid$z = as.vector(surf_list[[site_levels[i]]][[surf_names[u]]]$grid$z)
        surf.grid = data.frame(na.omit(surf.grid))
        
        p = ggplot() +
        geom_point(data = temp_points, aes(MDS1, MDS2, fill = !!parse_expr(surf_names[u]) ), shape = 22, size = 2.5) +
        stat_contour(data = surf.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
        scale_color_gradient(name ="GAM", low = "light grey", high = "black") +
        scale_fill_gradient(name = "obs", low = "light grey", high = "black") +
        my_gg_theme +
        theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
        labs(
            title = paste(
                site_levels[i], surf_names[u], "\nR2 =",
                signif( as.numeric(surf_list.anova[[site_levels[i]]][[surf_names[u]]]$r.sq), 2),
                ", P =",
                signif( as.numeric(surf_list.anova[[site_levels[i]]][[surf_names[u]]]$s.pv), 2)
                ),
            x = "PCoA 1",
            y = "PCoA 2")

        plot_list[[site_levels[i]]][[surf_names[u]]] = p
    }
}

#print plots
for(i in 1:length(plot_list)){
    name_site = names(plot_list)[i]
    n_col = floor(length(names(plot_list[[i]]))/2+0.5)
    w = n_col*6
    h = 10
    pdf(file = paste("GAM_fits/", name_site, "w_n_site.pdf", sep = ""), width = w, height = h)
    grid.arrange(grobs = plot_list[[i]], ncol = n_col)
    dev.off()
}

#coefs table

tot_elements = lapply(surf_list, function(x) length(names(x))) %>% unlist %>% sum

coefs_table = data.frame(
    site = vector(mode = "character", length = tot_elements),
    var = vector(mode = "character", length = tot_elements),
    R2 = vector(mode = "numeric", length = tot_elements),
    p.val = vector(mode = "numeric", length = tot_elements),
    sig = vector(mode = "character", length = tot_elements)
    stringsAsFactors = FALSE
)
index_count = 0

for(i in 1:length(site_levels)){
    #get NMDS points and metadata for site
    surf_names = names(surf_list[[ site_levels[[i]] ]])

    for(u in 1:length(surf_list[[i]])){
        index_count = index_count + 1
        
        coefs_table$site[index_count] = site_levels[i] %>% as.character
        coefs_table$var[index_count] = surf_names[u] %>% as.character
        coefs_table$R2[index_count] = signif( as.numeric(surf_list.anova[[site_levels[i]]][[surf_names[u]]]$r.sq), 2)
        coefs_table$p.val[index_count] = signif( as.numeric(surf_list.anova[[site_levels[i]]][[surf_names[u]]]$s.pv), 2)
        if(surf_list.anova[[site_levels[i]]][[surf_names[u]]]$s.pv < 0.05){
            coefs_table$sig[index_count] = "*"
        }else{
            coefs_table$sig[index_count] = "ns"
        }
    }
}

write.table(coefs_table, file = "GAM_fits/within_site_coefs.txt", sep = "\t", quote = F, row.names = F)
