library(ggplot2)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
#source("~/ggplot_theme.txt")

#source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

#lt_1K_samps = full_metadata %>%
#filter(total_seqs < 1000) %>%
#dplyr::select("sample")

#Species matrix on Nf and Nd occurence
#Nf_v_Nd.bin.gt1K = Nf_v_Nd.bin %>%
#filter(!sample %in% lt_1K_samps$sample)

#full_metadata.sorted = left_join(
#Nf_v_Nd.bin.gt1K,
#full_metadata,
#by = "sample"
#)

#full_metadata.sorted$occurence = factor(full_metadata.sorted$occurence, levels = c("none", "Nf", "Nd", "both"))

#tree_cond_table = full_metadata.sorted %>% group_by(TreeCond, occurence) %>% summarize(n = n())
#cankers_table = full_metadata.sorted %>% group_by(RaisedCanker, occurence) %>% summarize(n = n())

#dput(tree_cond_table)
#dput(cankers_table)

#############
#GGPLOT THEME
my_gg_theme = theme(panel.background = element_rect(fill='white', colour='black'),
panel.grid.major=element_blank(),
panel.grid.minor= element_blank(),
text=element_text(family="sans"),
axis.text=element_text(size=15, color="black"),
axis.ticks = element_line(color = "black"),
plot.title = element_text(hjust=0, size=20),
axis.title = element_text(size=17),
legend.title = element_blank(),
legend.text = element_text(size = 19),
strip.text = element_text(size = 15),
axis.title.x = element_text(margin = margin(t= 10)),
axis.title.y = element_text(margin = margin(r=10))    )

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    #reformat zeros
    l <- gsub("0e\\+00","0",l)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
}
###############################

#####################
#Data
tree_cond_table = structure(list(TreeCond = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 2L,
2L, 2L, 2L, 3L, 3L, 3L, 3L), occurence = structure(c(1L, 2L,
3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), .Label = c("none",
"Nf", "Nd", "both"), class = "factor"), n = c(3L, 12L, 1L, 3L,
11L, 21L, 1L, 9L, 9L, 9L, 2L, 9L, 3L, 2L, 1L, 6L)), row.names = c(NA,
-16L), class = c("grouped_df", "tbl_df", "tbl", "data.frame"), groups = structure(list(
TreeCond = 0:3, .rows = list(1:4, 5:8, 9:12, 13:16)), row.names = c(NA,
-4L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE))

cankers_table = structure(list(RaisedCanker = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L,
2L, 2L, 2L, 3L, 3L, 3L, 3L), occurence = structure(c(1L, 2L,
3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 4L, 1L, 2L, 3L, 4L), .Label = c("none",
"Nf", "Nd", "both"), class = "factor"), n = c(18L, 7L, 3L, 6L,
2L, 12L, 1L, 3L, 2L, 15L, 5L, 4L, 10L, 1L, 13L)), row.names = c(NA,
-15L), class = c("grouped_df", "tbl_df", "tbl", "data.frame"), groups = structure(list(
RaisedCanker = 0:3, .rows = list(1:4, 5:8, 9:11, 12:15)), row.names = c(NA,
-4L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE))
#####################
tree_cond_table %>% group_by(TreeCond) %>% summarize(sum(n))
cankers_table %>% group_by(RaisedCanker) %>% summarize(sum(n))

#Data formatting

#tree cond
tree_cond_table.start_end = data.frame(
    start = vector(mode = "numeric", length = 8),
    end = vector(mode = "numeric", length = 8),
    cond = c(0,0,1,1,2,2,3,3),
    spp = rep(c("Nf", "Nd"), 4),
    stringsAsFactors = F
)

conds = unique(tree_cond_table$TreeCond)

##########
#Start for
for(i in 1:length(conds)){
    temp.tab = tree_cond_table[tree_cond_table$TreeCond == conds[i],]
    
    ###
    #There are no cases in this data where Nd is greater than Nf, but if so these lines would be necessary
    spp = vector(mode = "character", length = 2)
    if(nrow(temp.tab[temp.tab$occurence == "Nd",]) == 0){
        spp = c("Nf", "Nd")
    }else
    if(temp.tab[temp.tab$occurence == "Nf", "n"] > temp.tab[temp.tab$occurence == "Nd", "n"]){
        spp = c("Nf", "Nd")
    }else{
        spp = c("Nd", "Nf")
    }
    ###

    #Calcs for start and end of spp bars
    tree_cond_table.start_end[
        tree_cond_table.start_end$cond == conds[i] & tree_cond_table.start_end$spp == spp[1], "start"
        ] = 0
    
    tree_cond_table.start_end[
        tree_cond_table.start_end$cond == conds[i] & tree_cond_table.start_end$spp == spp[1], "end"
        ] =
        (temp.tab[temp.tab$occurence == spp[1],"n"] + temp.tab[temp.tab$occurence == "both","n"])/sum(temp.tab$n)
    
    tree_cond_table.start_end[
        tree_cond_table.start_end$cond == conds[i] & tree_cond_table.start_end$spp == spp[2], "start"
        ] =
        temp.tab[temp.tab$occurence == spp[1],"n"]/sum(temp.tab$n)
    
    tree_cond_table.start_end[
        tree_cond_table.start_end$cond == conds[i] & tree_cond_table.start_end$spp == spp[2], "end"
        ] =
        (temp.tab[temp.tab$occurence == spp[1],"n"] + temp.tab[temp.tab$occurence == "both","n"] + temp.tab[temp.tab$occurence == spp[2],"n"])/sum(temp.tab$n)
}
##########

treeCond_n = tree_cond_table %>% group_by(TreeCond) %>% summarize(n = sum(n))
colnames(treeCond_n) = c("cond", "n")
tree_cond_table.start_end = full_join(tree_cond_table.start_end, treeCond_n, by = "cond")

tree_cond_table.start_end$spp = factor(tree_cond_table.start_end$spp, levels = c("Nf", "Nd"))


#tree cond
canker_table.start_end = data.frame(
start = vector(mode = "numeric", length = 8),
end = vector(mode = "numeric", length = 8),
cond = c(0,0,1,1,2,2,3,3),
spp = rep(c("Nf", "Nd"), 4),
stringsAsFactors = F
)

conds = unique(cankers_table$RaisedCanker)

##########
#Start for
for(i in 1:length(conds)){
    temp.tab = cankers_table[cankers_table$RaisedCanker == conds[i],]
    
    ###
    #There are no cases in this data where Nd is greater than Nf, but if so these lines would be necessary
    spp = vector(mode = "character", length = 2)
    if(nrow(temp.tab[temp.tab$occurence == "Nd",]) == 0){
        spp = c("Nf", "Nd")
    }else
    if(temp.tab[temp.tab$occurence == "Nf", "n"] > temp.tab[temp.tab$occurence == "Nd", "n"]){
        spp = c("Nf", "Nd")
    }else{
        spp = c("Nd", "Nf")
    }
    ###
    
    #Calcs for start and end of spp bars
    canker_table.start_end[
    canker_table.start_end$cond == conds[i] & canker_table.start_end$spp == spp[1], "start"
    ] = 0
    
    canker_table.start_end[
    canker_table.start_end$cond == conds[i] & canker_table.start_end$spp == spp[1], "end"
    ] =
    (temp.tab[temp.tab$occurence == spp[1],"n"] + temp.tab[temp.tab$occurence == "both","n"])/sum(temp.tab$n)
    
    canker_table.start_end[
    canker_table.start_end$cond == conds[i] & canker_table.start_end$spp == spp[2], "start"
    ] =
    temp.tab[temp.tab$occurence == spp[1],"n"]/sum(temp.tab$n)
    
    if(nrow(temp.tab[temp.tab$occurence == "Nd",]) == 0){
        canker_table.start_end[
        canker_table.start_end$cond == conds[i] & canker_table.start_end$spp == spp[2], "end"
        ] =
        (temp.tab[temp.tab$occurence == spp[1],"n"] + temp.tab[temp.tab$occurence == "both","n"])/sum(temp.tab$n)
    }else{
    canker_table.start_end[
    canker_table.start_end$cond == conds[i] & canker_table.start_end$spp == spp[2], "end"
    ] =
    (temp.tab[temp.tab$occurence == spp[1],"n"] + temp.tab[temp.tab$occurence == "both","n"] + temp.tab[temp.tab$occurence == spp[2],"n"])/sum(temp.tab$n)
    }
}
##########

cankers_n = cankers_table %>% group_by(RaisedCanker) %>% summarize(n = sum(n))
colnames(cankers_n) = c("cond", "n")
canker_table.start_end = full_join(canker_table.start_end, cankers_n, by = "cond")


canker_table.start_end$spp = factor(canker_table.start_end$spp, levels = c("Nf", "Nd"))


####################
#PLOTS

p1 = ggplot(tree_cond_table.start_end, aes(ymin = start, ymax = end, x = as.factor(cond), color = spp, label = n)) +
    geom_linerange(position = position_dodge(width = 0.5), size = 6) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), expand = expansion(mult = c(0.003,0.1))) +
    scale_color_manual(
        values = c("#d7301f", "#fdcc8a"),
        labels = c("Nf" = "N. faginata", "Nd" = "N. ditissima"),
        guide = F
    ) +
    geom_text(y = 1.05, size = 6, color = "black") +
    annotate(geom="segment", x = 0.4, xend = 4.6, y = 1, yend = 1) +
#    annotate(geom = "text", x = c(1,2,3,4), y = rep(1.05, 4), label = treeCond_n$n, size = 6) +
#    coord_cartesian(ylim = c(0,1), clip = "off") +
    labs(x = "Tree condition", y = "Trees (proportion)", title = "A") +
    my_gg_theme +
    theme(
        plot.title = element_text(hjust = -0.25, margin = margin(b = -20))
    )


p2 = ggplot(canker_table.start_end, aes(ymin = start, ymax = end, x = as.factor(cond), color = spp, label = n)) +
    geom_linerange(position = position_dodge(width = 0.5), size = 6) +
scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), expand = expansion(mult = c(0.003,0.1))) +
    scale_color_manual(values = c("#d7301f", "#fdcc8a"), labels = c("Nf" = expression(italic("N. faginata")), "Nd" = expression(italic("N. ditissima"))), guide = F ) +
geom_text(y = 1.05, size = 6, color = "black") +
annotate(geom="segment", x = 0.4, xend = 4.6, y = 1, yend = 1) +
#    annotate(geom = "text", x = c(1,2,3,4), y = rep(1.05, 4), label = treeCond_n$n, size = 6) +
#    coord_cartesian(ylim = c(0,1), clip = "off") +
    labs(x = "Cankers", y = "Trees (proportion)", title = "B") +
    my_gg_theme +
    theme(
        plot.margin = margin(l = 25, t = 5, ,r = 5,b = 5),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.11, margin = margin(b = -20)),
        axis.text.y = element_blank(),
        legend.text = element_text(hjust = 0)
    )


pdf("offset_bars_tree_cond_cankers.pdf", width = 8, height = 4)
grid.arrange(p1, p2, widths = c(0.525, 0.475) )
dev.off()









####################
#Example plots

pdf("offset_bars_tree_condition.pdf", height = 6, width = 8)
ggplot(tree_cond_table.start_end, aes(ymin = start, ymax = end, x = as.factor(cond), color = spp)) +
geom_linerange(position = position_dodge(width = 0.5), size = 6) +
scale_y_continuous(limits = c(0,1)) +
scale_color_brewer(palette = "Dark2", labels = c("Nf" = "N. faginata", "Nd" = "N.ditissima")) +
annotate(geom = "text", x = c(1,2,3,4), y = rep(0.9, 4), label = treeCond_n$n, size = 6) +
labs(x = "Tree condition", y = "Trees (proportion)") +
my_gg_theme
dev.off()

tree_cond_table.start_end$spp = factor(tree_cond_table.start_end$spp, levels = c("Nd", "Nf"))

pdf("offset_circles.pdf")
ggplot(tree_cond_table.start_end %>% filter(cond == 0), aes(xmin = start, xmax = end, y = spp, color = spp)) +
geom_linerange(position = position_dodge(width = 2.5), size = 15) +
scale_x_continuous(limits = c(0,1)) +
scale_color_brewer(palette = "Dark2", labels = c("Nf" = "N. faginata", "Nd" = "N.ditissima")) +
#annotate(geom = "text", x = c(1,2,3,4), y = rep(0.9, 4), label = treeCond_n$n, size = 5) +
coord_polar() +
my_gg_theme +
theme(
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text = element_blank(),
axis.ticks = element_blank()
)
dev.off()

















