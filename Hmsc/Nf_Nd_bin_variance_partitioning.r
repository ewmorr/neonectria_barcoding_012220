require(Hmsc)
require(vegan)
library(parallel)
library(corrplot)
require(MASS)

set.seed(1)
#vignette_3_multivariate_high
#Getting started with HMSC-R: high-dimensional
#multivariate models (source, pdf)
#vignette_2_multivariate_low
#Getting started with HMSC-R: low-dimensional
#multivariate models (source, pdf)
#vignette_4_spatial      Getting started with HMSC-R: spatial models
#(source, pdf)
#vignette_1_univariate   Getting started with HMSC-R: univariate models
#(source, pdf)
#vignette_5_performance
#Testing the performance of Hmsc with simulated
#data (source, pdf)


source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

#working with only asv_1 and asv_4 which are the dominant Nf and Nd asvs respectively

asv_tab.neonectria = asv_tab[rownames(asv_tab) %in% c("ASV_1", "ASV_4"),] %>% t()
asv_tab.neonectria[asv_tab.neonectria > 0] = 1

full_metadata.sorted = left_join(
    data.frame(sample = rownames(asv_tab.neonectria)),
    full_metadata,
    by = "sample"
)

#############################
#Get spatial distance matrix#

convert_lat_lon_to_cartesian = function(lat, lon){
    #https://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates
    lat_rad = lat*pi/180
    lon_rad = lon*pi/180
    
    R = 6371
    x = R * cos(lat_rad) * cos(lon_rad)
    y = R * cos(lat_rad) * sin(lon_rad)
    z = R *sin(lat_rad)
    return(c(x,y,z))
}

site_coords = data.frame(lat = full_metadata.sorted$lat, lon = full_metadata.sorted$lon)

site_cart_coords = data.frame(
    x = vector("numeric", length(site_coords$lat)),
    y = vector("numeric", length(site_coords$lat)),
    z = vector("numeric", length(site_coords$lat))
)

for(i in 1:length(site_coords$lat)){
    site_cart_coords[i,1:3] = convert_lat_lon_to_cartesian(site_coords[i,1], site_coords[i,2])
}

rownames(site_cart_coords) = full_metadata.sorted$sample
spatial_dist_mat = vegdist(site_cart_coords, method = "euclidean")
#Set spatial matrix as random effect
rL.spatial = HmscRandomLevel(distMat = spatial_dist_mat)
#set latent variables to one (This wold not be ideal in a community with more than 2 spp
#rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1)

#Or use lat lon
xycoords = data.frame(lat = full_metadata.sorted$lat, lon = full_metadata.sorted$lon)
xycoords$lon = xycoords$lon - min(xycoords$lon)
#rownames(xycoords) = full_metadata.sorted$sample
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1)

####################
#Also set a random level of total_seqs
rL.seqs = HmscRandomLevel(xData = full_metadata.sorted$total_seqs)

############################
#Begin of model specification

#Species matrix
Y = as.matrix(asv_tab.neonectria)
#rownames(Y) = NULL

#covariate data fixed effects
XData = full_metadata.sorted %>% dplyr::select("HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "RaisedCanker", "Wax", "total_seqs")
#total_seqs should be supplied as part of the xData as described here https://github.com/hmsc-r/HMSC/issues/10 may also need to adjust priors

#Supply sample designation to study design
studyDesign = data.frame(
    spatial = full_metadata.sorted$Site
)
#It's very importnat to set above to factor and not character, otherwise the preds are not recoverable

#Would include this as a random effect if not setting spatial mat
#rL = HmscRandomLevel(units = studyDesign$sample)
#rL$nfMax = 15


#model specification
m.spatial = Hmsc(Y=Y, XData=XData, XFormula=~.,
studyDesign=studyDesign, ranLevels=list(spatial = rL.spatial),distr="probit")

#########################################################
#This block works through MCMC
rL = HmscRandomLevel(units = full_metadata.sorted$sample)

studyDesign = data.frame(
sample = full_metadata.sorted$sample
)

m.spatial = Hmsc(Y=Y, XData=XData, XFormula=~.,
studyDesign=studyDesign, ranLevels=list(sample = rL),distr="probit")
#########################################################










#########################################################
#This block does not through MCMC
rL = HmscRandomLevel(units = full_metadata.sorted$sample)

xycoords = data.frame(lat = full_metadata.sorted$lat, lon = full_metadata.sorted$lon)
xycoords$lon = xycoords$lon - min(xycoords$lon)
rownames(xycoords) = sprintf('spatial_%.2d',1:117)
#rownames(xycoords) = full_metadata.sorted$sample
rL.spatial = HmscRandomLevel(sData = xycoords)
#rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1)

rL.site = HmscRandomLevel(units = full_metadata.sorted$Site)

studyDesign = data.frame(
sample = full_metadata.sorted$sample %>% as.factor,
site = full_metadata.sorted$Site %>% as.factor,
spatial = sprintf('spatial_%.2d',1:117) %>% as.factor
)

m.spatial = Hmsc(Y=Y, XData=XData, XFormula=~.,
studyDesign=studyDesign, ranLevels=list(sample = rL, spatial = rL.spatial, site = rL.site),distr="probit")
#########################################################


###########################
#MCMC sampling
#if running fast set test.run = T, but these results will not be reliable
nChains = 2
test.run = T
if (test.run){
    # with this option, the vignette runs fast but results are not reliable
    thin = 1
    samples = 10
    transient = 5
    verbose = 0
} else {
    # with this option, the vignette evaluates slow but it reproduces the results of
    # the .pdf version
    thin = 10
    samples = 1000
    transient = 1000
    verbose = 0
}

m.spatial = sampleMcmc(m.spatial, thin = thin, samples = samples, transient = transient,
nChains = nChains, verbose = verbose,updater=list(GammaEta=FALSE))

######################
#MCMC convergence

#ess == effective sample size
#psrf == potential scale reduction factor (closeto 1 indicates good MCMC convergence)
#beta == species niches (as realted to covarites)
#gamma == influence of traits on species niches
#omega == residual species association
#rho == phylogenetic signal

mpost = convertToCodaObject(m.spatial)
par(mfrow=c(3,2))
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)
ns = 50
sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
    tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

##########################
#Model fit and partioning

preds = computePredictedValues(m.spatial)
MF = evaluateModelFit(hM=m.spatial, predY=preds)
hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$TjurR2),2)))

#######################
#Variance partioning
#the group vairable assigns X covariates to different covariate groups
#so first look at the design matrix

head(m.spatial$X)
#For our real data we fit an intercept and three coninuous variables, so they can each be assigned separate groups
#If instead we had a categorical variable the levels could be assigned to a single group along with the intercept

VP = computeVariancePartitioning(m.spatial, group = c(1,2,3,4,5,6,7,8), groupnames = c("intercept","HDD4.nongrowing", "freezeThaw.nongrowing", "ppt.nongrowing", "duration_infection", "RaisedCanker", "Wax", "total_seqs"))
plotVariancePartitioning(m.spatial, VP = VP)

####################
#Plot variance partioning
#If included a square of continuous variable, negative response would indicate intermediate niche optimum (i.e., abundance goes up initially but then goes down)
postBeta = getPostEstimate(m.spatial, parName = "Beta")
plotBeta(m.spatial, post = postBeta, param = "Support",
plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
#This can also be mapped on a tree with plotTree = T, but then tree must be included in model

############################################
#transform VP and postBeta object for ggplot
#VP cpntains R2 vals and postBeta contains support (i.e. alpha)

VP.vals = data.frame(VP$vals)
colnames(VP.vals) = c("Nf", "Nd")
VP.vals = VP.vals[1:8,] #this removes the random effect

VP.vals$variable = c("intercept","HDD4.nongrowing", "freezeThaw.nongrowing", "ppt.nongrowing", "duration_infection", "cankers", "wax", "total_seqs")
VP.vals.long = VP.vals %>%
    pivot_longer(-variable, names_to = "ASV", values_to = "R2")


#Transform R2 based on positive or negative response
postBeta.mean = data.frame(postBeta$mean)
colnames(postBeta.mean) = c("Nf", "Nd")
postBeta.mean$variable = c("intercept","HDD4.nongrowing", "freezeThaw.nongrowing", "ppt.nongrowing", "duration_infection", "cankers", "wax", "total_seqs")
postBeta.mean.long = postBeta.mean %>%
pivot_longer(-variable, names_to = "ASV", values_to = "mean")

for(i in 1:length(postBeta.mean.long$mean)){
    if(postBeta.mean.long$mean[i] < 0){
        VP.vals.long$R2[i] = VP.vals.long$R2[i] * -1
    }
}


postBeta.support = data.frame(postBeta$support)
colnames(postBeta.support) = c("Nf", "Nd")
postBeta.support$variable = c("intercept","HDD4.nongrowing", "freezeThaw.nongrowing", "ppt.nongrowing", "duration_infection", "cankers", "wax", "total_seqs")
postBeta.support.long = postBeta.support %>%
pivot_longer(-variable, names_to = "ASV", values_to = "support")


postBeta.supportNeg = data.frame(postBeta$supportNeg)
colnames(postBeta.supportNeg) = c("Nf", "Nd")
postBeta.supportNeg$variable = c("intercept","HDD4.nongrowing", "freezeThaw.nongrowing", "ppt.nongrowing", "duration_infection", "cankers", "wax", "total_seqs")
postBeta.supportNeg.long = postBeta.supportNeg %>%
pivot_longer(-variable, names_to = "ASV", values_to = "supportNeg")


VP.vals.support = full_join(VP.vals.long, postBeta.support.long) %>%
    full_join(., postBeta.supportNeg.long)
VP.vals.support = data.frame(VP.vals.support)

VP.vals.support$P.val = vector(mode = "character", length = length(VP.vals.support$support))

for(i in 1:length(VP.vals.support$P.val)){
    if(VP.vals.support$supportNeg[i] > 0.95 || VP.vals.support$support[i] > 0.95){
        VP.vals.support$P.val[i] = "P<0.05"
    }else{
        VP.vals.support$P.val[i] = "n.s."
    }
}

require(RColorBrewer)
source("~/ggplot_theme.txt")

VP.vals.support$variable = factor(VP.vals.support$variable, levels = c("intercept","HDD4.nongrowing", "freezeThaw.nongrowing", "ppt.nongrowing", "duration_infection", "cankers", "wax", "total_seqs"))

p1 = ggplot(VP.vals.support, aes(ASV, variable, fill = R2, color = P.val)) +
geom_tile(size = 1, height = 0.975, width = 0.975) +
scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0) +
scale_color_manual(values = c("P<0.05" = "black", "n.s." = "white"), ) +
my_gg_theme +
labs(
    x = "",
    y = "",
    fill = "R2",
    color = "alpha",
    title = "HMSC variance partitioning\nNeonectria presence-absence"
) +
theme(legend.title = element_text(size = 20))

pdf("HMSC/Nf_Nd_variance_partitioning.bin.pdf", width = 8, height = 8)
p1
dev.off()

#scale_fill_continuous(limits = c(-1,1))

####################
####################
#Estimated residual var between spp

OmegaCor = computeAssociations(m.spatial)
supportLevel = 0.55 #this is alpha = 0.05
toPlot = ((OmegaCor[[1]]$support>supportLevel)
+ (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
col=colorRampPalette(c("blue","white","red"))(500),
tl.cex=.6, tl.col="black",
title=paste("random effect level:", m.spatial$rLNames[1]), mar=c(0,0,1,0))


#################
#Refit model with additional aprams


#covariate data fixed effects
XData = full_metadata.sorted %>% dplyr::select("HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "RaisedCanker", "Wax", "dbh", "TreeCond", "total_seqs")



#########################################################
#This block works through MCMC
rL = HmscRandomLevel(units = full_metadata.sorted$sample)

studyDesign = data.frame(
sample = full_metadata.sorted$sample
)

m.spatial = Hmsc(Y=Y, XData=XData, XFormula=~.,
studyDesign=studyDesign, ranLevels=list(sample = rL),distr="probit")
#########################################################



m.spatial = sampleMcmc(m.spatial, thin = thin, samples = samples, transient = transient,
nChains = nChains, verbose = verbose,updater=list(GammaEta=FALSE))

######################
#MCMC convergence

#ess == effective sample size
#psrf == potential scale reduction factor (closeto 1 indicates good MCMC convergence)
#beta == species niches (as realted to covarites)
#gamma == influence of traits on species niches
#omega == residual species association
#rho == phylogenetic signal

mpost = convertToCodaObject(m.spatial)
par(mfrow=c(3,2))
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)
ns = 50
sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
    tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

##########################
#Model fit and partioning

preds = computePredictedValues(m.spatial)
MF = evaluateModelFit(hM=m.spatial, predY=preds)
hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$TjurR2),2)))

#######################
#Variance partioning
#the group vairable assigns X covariates to different covariate groups
#so first look at the design matrix

head(m.spatial$X)
#For our real data we fit an intercept and three coninuous variables, so they can each be assigned separate groups
#If instead we had a categorical variable the levels could be assigned to a single group along with the intercept

VP = computeVariancePartitioning(m.spatial, group = c(1,2,3,4,5,6,7,8,9,10), groupnames = c("intercept","HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "RaisedCanker", "Wax", "dbh", "TreeCond", "total_seqs"))
plotVariancePartitioning(m.spatial, VP = VP)

####################
#Plot variance partioning
#If included a square of continuous variable, negative response would indicate intermediate niche optimum (i.e., abundance goes up initially but then goes down)
postBeta = getPostEstimate(m.spatial, parName = "Beta")
plotBeta(m.spatial, post = postBeta, param = "Support",
plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
#This can also be mapped on a tree with plotTree = T, but then tree must be included in model


############################################
#transform VP and postBeta object for ggplot
#VP cpntains R2 vals and postBeta contains support (i.e. alpha)

VP.vals = data.frame(VP$vals)
colnames(VP.vals) = c("Nf", "Nd")
VP.vals = VP.vals[1:10,] #this removes the random effect

VP.vals$variable = c("intercept","HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "cankers", "wax", "dbh", "treeCond", "total_seqs")
VP.vals.long = VP.vals %>%
pivot_longer(-variable, names_to = "ASV", values_to = "R2")


#Transform R2 based on positive or negative response
postBeta.mean = data.frame(postBeta$mean)
colnames(postBeta.mean) = c("Nf", "Nd")
postBeta.mean$variable = c("intercept","HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "cankers", "wax", "dbh", "treeCond", "total_seqs")
postBeta.mean.long = postBeta.mean %>%
pivot_longer(-variable, names_to = "ASV", values_to = "mean")

for(i in 1:length(postBeta.mean.long$mean)){
    if(postBeta.mean.long$mean[i] < 0){
        VP.vals.long$R2[i] = VP.vals.long$R2[i] * -1
    }
}


postBeta.support = data.frame(postBeta$support)
colnames(postBeta.support) = c("Nf", "Nd")
postBeta.support$variable = c("intercept","HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "cankers", "wax", "dbh", "treeCond", "total_seqs")
postBeta.support.long = postBeta.support %>%
pivot_longer(-variable, names_to = "ASV", values_to = "support")


postBeta.supportNeg = data.frame(postBeta$supportNeg)
colnames(postBeta.supportNeg) = c("Nf", "Nd")
postBeta.supportNeg$variable = c("intercept","HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "cankers", "wax", "dbh", "treeCond", "total_seqs")
postBeta.supportNeg.long = postBeta.supportNeg %>%
pivot_longer(-variable, names_to = "ASV", values_to = "supportNeg")


VP.vals.support = full_join(VP.vals.long, postBeta.support.long) %>%
full_join(., postBeta.supportNeg.long)
VP.vals.support = data.frame(VP.vals.support)

VP.vals.support$P.val = vector(mode = "character", length = length(VP.vals.support$support))

for(i in 1:length(VP.vals.support$P.val)){
    if(VP.vals.support$supportNeg[i] > 0.95 || VP.vals.support$support[i] > 0.95){
        VP.vals.support$P.val[i] = "P<0.05"
    }else{
        VP.vals.support$P.val[i] = "n.s."
    }
}

require(RColorBrewer)
source("~/ggplot_theme.txt")

VP.vals.support$variable = factor(VP.vals.support$variable, levels = c("intercept","HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "duration_infection", "cankers", "wax", "dbh", "treeCond", "total_seqs"))

p1 = ggplot(VP.vals.support, aes(ASV, variable, fill = R2, color = P.val)) +
geom_tile(size = 1, height = 0.975, width = 0.975) +
scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0) +
scale_color_manual(values = c("P<0.05" = "black", "n.s." = "white"), ) +
my_gg_theme +
labs(
x = "",
y = "",
fill = "R2",
color = "alpha",
title = "HMSC variance partitioning\nNeonectria presence-absence"
) +
theme(legend.title = element_text(size = 20))

pdf("HMSC/Nf_Nd_variance_partitioning.bin.plus_vars.pdf", width = 8, height = 8)
p1
dev.off()
