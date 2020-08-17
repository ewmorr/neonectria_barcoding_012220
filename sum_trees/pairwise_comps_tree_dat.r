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
    dbh = full_dat$dbh,
    neo_fruit = full_dat$NeoFruiting,
    cankers = full_dat$RaisedCanker,
#    tarry_spot = full_dat$Site_mean.TarrySpots,
    wax = full_dat$Wax,
#    xyloc = full_dat$Site_mean.Xylococcus,
    tree_cond = full_dat$TreeCond
)

##################
#pairwise lms

#This script contains some additional functions that will be useful for performing pairwise LMs and plotting
source("~/repo/neonectria_barcoding_012220/ecol/general_paired_simpleLM.r")

#To get a first look at the significant relationships in the matrix
pairwise_df.lm_coefs = pairwise_simpleLM(pairwise_df) %>% filter(beta.pv < 0.05, LHS != RHS) %>% select(LHS, RHS, beta.pv, R2)

pairwise_df.lm_coefs %>% filter(R2 >= 0.5)

#######################################################

#################################
#Custon GGally plotting function#


GGscatterPlot <- function(data, mapping, ...,
method = "pearson") {
    
    #Get correlation coefficient
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    
    cor <- cor(x, y, method = method)
    
    fit <- lm(y ~ x)
    fs <- summary(fit)$fstatistic
    pValue <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)

    #print(r.squared)
    #print(pValue)
   
    
    #define corColor
    corColors <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
    
#    print(corColors)
    print(cor)
    
    if (cor <= -0.8) {
        corCol <- corColors[9]
        print(corCol)
    } else if (cor <= -0.6) {
        corCol <- corColors[8]
        print(corCol)
    } else if (cor <= -0.4) {
        corCol <- corColors[7]
        print(corCol)
    } else if (cor <= -0.2) {
        corCol <- corColors[6]
        print(corCol)
    } else if (cor < 0.2) {
        corCol <- corColors[5]
        print(corCol)
    } else if (cor < 0.4) {
        corCol <- corColors[4]
        print(corCol)
    } else if (cor < 0.6) {
        corCol <- corColors[3]
        print(corCol)
    } else if (cor < 0.8) {
        corCol <- corColors[2]
        print(corCol)
    } else {
        corCol <- corColors[1]
        print(corCol)
    }

    #Define linetype
    if(pValue < 0.05){
        my_line_type = "solid"
    }else{
        my_line_type = "dashed"
    }
    #Assemble data frame
     
    df <- data.frame(x = x, y = y)
    # PCA
    nonNull <- x!=0 & y!=0
    dfpc <- prcomp(~x+y, df[nonNull,])
    df$cols <- predict(dfpc, df)[,1]
    # Define the direction of color range based on PC1 orientation:
    dfsum <- x+y
    colDirection <- ifelse(dfsum[which.max(df$cols)] <
    dfsum[which.min(df$cols)],
    1,
    -1)
    #Get 2D density for alpha
    dens2D <- MASS::kde2d(df$x, df$y)
    df$density <- fields::interp.surface(dens2D ,
    df[,c("x", "y")])
    
    if (any(df$density==0)) {
        mini2D = min(df$density[df$density!=0]) #smallest non zero value
        df$density[df$density==0] <- mini2D
    }
    #Prepare plot
    pp <- ggplot(df, aes(x=x, y=y)) +
    #ggplot2::geom_label(
    #data = data.frame(
    #xlabel = 0,#max(x, na.rm = TRUE),
    #ylabel = min(y, na.rm = TRUE),
    #lab = round(cor, digits = 3)),
    #mapping = ggplot2::aes(x = xlabel,
    #y = ylabel,
    #label = lab),
    #hjust = 0, vjust = 0.5,
    #size = 3, fontface = "bold",
    #inherit.aes = FALSE # do not inherit anything from the ...
    #) +
    geom_smooth(method = "lm", color = "dark grey") +
    ggplot2::geom_point(shape=16, show.legend = FALSE) +
    theme(
        panel.background = element_rect(color = "black", linetype = my_line_type, fill = corCol, size = 1),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")
    )
    
    return(pp)
}

##################################
#Function for only lower triangle#

gpairs_lower <- function(g){
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow -1
    
    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1
    
    g
}

#plot
p1 = ggpairs(pairwise_df,
    lower = list(continuous = GGscatterPlot),
    upper = "blank",
    diag = "blank",
    switch = "both",
    columnLabels = c("lat", "lon", "elev" , "GDD4\nsummer", "GDD4\nwinter", "freeze_thaw\nsummer", "freeze_thaw\nwinter", "precip\nsummer", "precip\nwinter", "infec_dur","dbh","neo_fruit","cankers", "wax","tree_cond")) +
theme(
strip.text = element_text(size = 11),
plot.title = element_text(hjust = 0.85, vjust = -17, size = 25)
)

g1 = gpairs_lower(p1)

pdf("prelim_figs/pairwise_comps_site_dat.GDD_based.tree_level_vars.pdf", width = 16, height = 16)
print(g1)
dev.off()



