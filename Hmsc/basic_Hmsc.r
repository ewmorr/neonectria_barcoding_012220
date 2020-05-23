require(Hmsc)
require(vegan)
library(parallel)
library(corrplot)
require(MASS)
require(ape)
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


#This is workflow in vignette 3

source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

#asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")

#Filter table to include only those sample with greater than 1000 seqs
asv_tab.gt3freq = asv_tab[rowSums(asv_tab > 0) >= 3,]
asv_tab.gt1k = asv_tab.gt3freq[,colSums(asv_tab.gt3freq) > 1000]
#minSamps = floor(length(colSums(asv_tab.gt1k))/10)
asv_tab.gt1k.rare = rrarefy(t(asv_tab.gt1k), 1000)

asv_tab.gt1k.rare = asv_tab.gt1k.rare[,colSums(asv_tab.gt1k.rare > 0) >= 2]

top_50 = (colSums(asv_tab.gt1k.rare) %>% sort(.,decreasing = T))[1:50] %>% names

asv_tab.gt1k.rare.top_50.log10 = log10(asv_tab.gt1k.rare[,colnames(asv_tab.gt1k.rare) %in% top_50] + 1)

asv_tab.gt1k.rare.top_50 = asv_tab.gt1k.rare[,colnames(asv_tab.gt1k.rare) %in% top_50]

full_metadata.sorted = left_join(
    data.frame(sample = rownames(asv_tab.gt1k.rare)),
    full_metadata,
    by = "sample"
)

############################
#Begin of model specification
Y = as.matrix(asv_tab.gt1k.rare.top_50.log10)
rownames(Y) = NULL
#covariate data for model
XData = full_metadata.sorted %>% dplyr::select("NeoFruiting", "ppt", "duration_infection")
#XData = full_metadata.sorted %>% select("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "TreeCond", "ppt", "tmin", "tmax", "duration_infection")
#other needed
studyDesign = full_metadata.sorted %>% dplyr::select(sample)
#It's very importnat to set this to factor and not character, otherwise the preds are not recoverable
studyDesign$sample = as.factor(studyDesign$sample)

rL = HmscRandomLevel(units = studyDesign$sample)
rL$nfMax = 15

#model specification
m = Hmsc(Y = Y, XData = XData, XFormula = ~.,
studyDesign = studyDesign, ranLevels = list(sample = rL))


######################
#Or use simulated data

ns = 50
phy = ape::rcoal(n=ns, tip.label = sprintf('species_%.3d',1:ns), br = "coalescent")
plot(phy, show.tip.label = FALSE, no.margin = TRUE)
C = vcv(phy, model = "Brownian", corr = TRUE)
spnames = colnames(C)
traits = matrix(NA,ncol =2,nrow = ns)
for (i in 1:2){
    traits[,i] = matrix(mvrnorm(n = 1, mu = rep(0,ns), Sigma=C))
}
rownames(traits) = spnames
colnames(traits) = c("habitat.use","thermal.optimum")
traits = as.data.frame(traits)
par(fig = c(0,0.6,0,0.8), mar=c(6,0,2,0))
plot(phy, show.tip.label = FALSE)
par(fig = c(0.6,0.9,0.025,0.775), mar=c(6,0,2,0), new=T)
plot.new()
image.plot(t(traits),axes=FALSE,legend.width = 3,legend.shrink=1,
col = colorRampPalette(c("blue","white","red"))(200))
text(x=1.1, y=0.72, srt = 90, "H", cex=0.9, pos = 4)
text(x=1.4, y=0.72, srt = 90, "T", cex=0.9, pos = 4)
n = 101
habitat = factor(sample(x = c("forest","open"), size = n, replace=TRUE))
climate = rnorm(n)
nc = 4
mu = matrix(0,nrow=nc,ncol=ns)
#expected niche of each species related to the "covariate" intercept
mu[1, ] = -traits$thermal.optimum^2/4-traits$habitat.use
#expected niche of each species related to the covariate forest
#(open area as reference level, so included in intercept)
mu[2, ] = 2*traits$habitat.use
#expected niche of each species related to the covariate climate
mu[3, ] = traits$thermal.optimum/2
#expected niche of each species related to the covariate climate*climate
mu[4, ] = -1/4
beta = mu + 0.25*matrix(rnorm(n = ns*nc), ncol=ns)
X = cbind(rep(1,ns), as.numeric(habitat=="forest"), climate, climate*climate)
L = X%*%beta
Y = L + mvrnorm(n=n, mu=rep(0,ns), Sigma=diag(ns))
colnames(Y) = spnames

#model specification
XData = data.frame(climate = climate, habitat = habitat)
XFormula = ~habitat + poly(climate,degree = 2,raw = TRUE)
TrFormula = ~habitat.use + thermal.optimum
studyDesign = data.frame(sample = sprintf('sample_%.3d',1:n))
rL = HmscRandomLevel(units = studyDesign$sample)
rL$nfMax = 15

m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
#TrData = traits, TrFormula = TrFormula,
#phyloTree = phy,
studyDesign = studyDesign, ranLevels = list(sample = rL))

###########################
#MCMC sampling
#if running fast set test.run = T, but these results will not be reliable

nChains = 2
test.run = TRUE
if (test.run){
    #with this option, the vignette evaluates in ca. 10 minutes in a laptop
    thin = 1
    samples = 100
    transient = 50
} else {
    #with this option, the vignette evaluates in ca. 2 hrs in a laptop
    thin = 10
    samples = 1000
    transient = 500
}
verbose = 0

m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
nChains = nChains, nParallel = nChains, verbose = verbose)

######################
#MCMC convergence

#ess == effective sample size
#psrf == potential scale reduction factor (closeto 1 indicates good MCMC convergence)
#beta == species niches (as realted to covarites)
#gamma == influence of traits on species niches
#omega == residual species association
#rho == phylogenetic signal

mpost = convertToCodaObject(m)
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

#These will be NULL without a phylo tree fit (rho)
print("ess.rho:")
effectiveSize(mpost$Rho)
print("psrf.rho:")
gelman.diag(mpost$Rho)$psrf

##########################
#Model fit and partioning

preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))

#######################
#Variance partioning
#the group vairable assigns X covariates to different covariate groups
#so first look at the design matrix

head(m$X)
#For our real data we fit an intercept and three coninuous variables, so they can each be assigned separate groups
#If instead we had a categorical variable the levels could be assigned to a single group along with the intercept

VP = computeVariancePartitioning(m, group = c(1,2,3,4), groupnames = c("intercept","neoFruit", "ppt", "durInf"))
plotVariancePartitioning(m, VP = VP)


########################
#Traits
#Let us then ask how much the traits explain out of the variation among the species in their responses to environmental covariates

#We did not have any explicit trait data, so this will be NA

#kable(VP$R2T$Beta)
VP$R2T$Beta

#We next ask how the influence of traits on species niches propagates into the influence of traits on species abundances.
VP$R2T$Y

####################
#Plot variance partioning
#If included a square of continuous variable, negative response would indicate intermediate niche optimum (i.e., abundance goes up initially but then goes down)
postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support",
plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
#This can also be mapped on a tree with plotTree = T, but then tree must be included in model

#Gamma models how species traits influence niche (so this is NULL)
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post=postGamma, param="Support", supportLevel = 0.95)


####################
####################
#Estimated residual var between spp

OmegaCor = computeAssociations(m)
supportLevel = 0.95 #this is alpha = 0.05
toPlot = ((OmegaCor[[1]]$support>supportLevel)
+ (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
col=colorRampPalette(c("blue","white","red"))(500),
tl.cex=.6, tl.col="black",
title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))


#####################
#Strength of phylo signal

summary(mpost$Rho)


#######################
#Plotting variation of gradients

#this is breaking, but may not be that useful
Gradient = constructGradient(m,focalVariable = "duration_infection",
non.focalVariables = list("NeoFruiting", "ppt")
