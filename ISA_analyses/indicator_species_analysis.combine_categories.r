#Before running this script need to set memory limit for environment#in terminal run
#cd ~
#touch .Renviron
#vim .Renviron
#then insert and save the following line to .Renviron
#R_MAX_VSIZE=100Gb
#Memory size can be adjusted

require(tidyverse)
require(vegan)
require(indicspecies)

rm(list = ls())
#Sys.setenv('R_MAX_VSIZE'=128000000000)
Sys.getenv('R_MAX_VSIZE')

source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
#this pulls in objects:
#asv_tab
#asv_tax
#id_bench_map

#joins metadata files to get metadata object:
#full_metadata

#creates negatives and controls only asv_tab (long format):
#asv_tab.negatives.long

#creates new object with lowest informative taxonomic level and a character asv_tax with unknown instead of NA
#asv_informative_taxa
#asv_tax.char

#and calc duration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

##########################
#1000 seqs per sample min#
##########################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
asv_tab.gt1K.rare.bin = asv_tab.gt1K.rare
asv_tab.gt1K.rare.bin[asv_tab.gt1K.rare.bin > 0] = 1

full_metadata.sorted = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare)),
full_metadata
) %>% left_join(., Nf_v_Nd.bin)

full_metadata.sorted$Wax.lowHigh = vector(mode = "character", length = nrow(full_metadata.sorted))
full_metadata.sorted$TreeCond.lowHigh = vector(mode = "character", length = nrow(full_metadata.sorted))
full_metadata.sorted$NeoFruiting.lowHigh = vector(mode = "character", length = nrow(full_metadata.sorted))
full_metadata.sorted$RaisedCanker.lowHigh = vector(mode = "character", length = nrow(full_metadata.sorted))

for(i in 1:nrow(full_metadata.sorted)){
    if(full_metadata.sorted$Wax[i] <= 1){
        full_metadata.sorted$Wax.lowHigh[i] = "lowWax"
    }else{
        full_metadata.sorted$Wax.lowHigh[i] = "highWax"
    }
    if(full_metadata.sorted$TreeCond[i] <= 1){
        full_metadata.sorted$TreeCond.lowHigh[i] = "lowCrownDie"
    }else{
        full_metadata.sorted$TreeCond.lowHigh[i] = "highCrownDie"
    }
    if(full_metadata.sorted$NeoFruiting[i] <= 1){
        full_metadata.sorted$NeoFruiting.lowHigh[i] = "lowNeoFruit"
    }else{
        full_metadata.sorted$NeoFruiting.lowHigh[i] = "highNeoFruit"
    }
    if(full_metadata.sorted$RaisedCanker[i] <= 1){
        full_metadata.sorted$RaisedCanker.lowHigh[i] = "lowCanker"
    }else{
        full_metadata.sorted$RaisedCanker.lowHigh[i] = "highCanker"
    }
}

full_metadata.sorted$diseaseSeverityCombined = interaction(
    full_metadata.sorted$Wax.lowHigh,
    full_metadata.sorted$TreeCond.lowHigh,
    full_metadata.sorted$NeoFruiting.lowHigh,
    full_metadata.sorted$RaisedCanker.lowHigh
)

full_metadata.sorted %>% group_by(diseaseSeverityCombined) %>% summarize(n())

#for ISA only consider single factor level (not combinations) for combined disease severity vars
duleg = T


#vars to use from full_metadata.sorted
#duration_infection NeoFruiting RaisedCanker  Wax Xylococcus TreeCond quartile occurence

combined_disease_severity_ISA = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$diseaseSeverityCombined, duleg = T)
summary(combined_disease_severity_ISA)
summary(combined_disease_severity_ISA, indvalcom = T)


full_metadata.sorted$diseaseSeverityCombined = interaction(
full_metadata.sorted$Wax.lowHigh,
full_metadata.sorted$TreeCond.lowHigh,
#full_metadata.sorted$NeoFruiting.lowHigh,
full_metadata.sorted$RaisedCanker.lowHigh
)

full_metadata.sorted %>% group_by(diseaseSeverityCombined) %>% summarize(n())

combined_disease_severity_ISA = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$diseaseSeverityCombined, duleg = T)
summary(combined_disease_severity_ISA)


full_metadata.sorted$diseaseSeverityCombined = interaction(
#full_metadata.sorted$Wax.lowHigh,
full_metadata.sorted$TreeCond.lowHigh,
full_metadata.sorted$NeoFruiting.lowHigh,
full_metadata.sorted$RaisedCanker.lowHigh
)

full_metadata.sorted %>% group_by(diseaseSeverityCombined) %>% summarize(n())

combined_disease_severity_ISA = multipatt(asv_tab.gt1K.rare.bin, cluster = full_metadata.sorted$diseaseSeverityCombined, duleg = T)
summary(combined_disease_severity_ISA)
