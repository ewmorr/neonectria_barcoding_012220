
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")

library(Hmsc)
library(parallel)
library(corrplot)
library(MASS)
set.seed(6)

###############################
#Y matrix

#remove samples with less than 1000 seqs
lt_1K_samps = full_metadata %>%
    filter(total_seqs < 1000) %>%
    dplyr::select("sample")

#Species matrix on Nf and Nd occurence
Ytable = Nf_v_Nd.bin %>%
    filter(!sample %in% lt_1K_samps$sample) %>%
    dplyr::select(Nf, Nd)
rownames(Ytable) = ( Nf_v_Nd.bin %>% filter(!sample %in% lt_1K_samps$sample) )$sample
yprob = as.matrix(Ytable)

#species matrix if using asv_tab
#working with only asv_1 and asv_4 which are the dominant Nf and Nd asvs respectively
#asv_tab.neonectria = asv_tab[rownames(asv_tab) %in% c("ASV_1", "ASV_4"),] %>% t()
#asv_tab.neonectria[asv_tab.neonectria > 0] = 1
#yprob = asv_tab.neonectria

full_metadata.sorted = left_join(
data.frame(sample = rownames(Ytable)),
full_metadata,
by = "sample"
)


#####################
#Define random levels
xycoords = matrix(c(site_info$lat, site_info$lon),ncol=2)

colnames(xycoords) = c("x-coordinate","y-coordinate")
rownames(xycoords) = site_info$Site

studyDesign = data.frame(sample = as.factor(full_metadata.sorted$sample), site_spatial = as.factor(full_metadata.sorted$Site))

rL.sample = HmscRandomLevel(units = as.factor(full_metadata.sorted$sample))
rL.site_spatial = HmscRandomLevel(sData = xycoords)

#################
#covariate matrix

#first running full set of vars (including wax and excluding winter precip [they have 0.5 R^2], including infection duration and excluding cankers [cankers is correalted with inf_dur, tree_cond, wax, winter freeze-thaw)
XData = full_metadata.sorted %>% dplyr::select("HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "duration_infection", "dbh", "Wax", "TreeCond", "total_seqs")
XData[,colnames(XData) == "total_seqs"] = log(XData[,colnames(XData) == "total_seqs"])

XData = apply(XData,2,scale)
XData = data.frame(XData)

m.spatial.full_covars = Hmsc(Y=yprob, XData=XData, XFormula=~.,
studyDesign=studyDesign, ranLevels=list("sample"=rL.sample, "site_spatial" = rL.site_spatial),distr="probit")

#set priors for total_seqs so that it is used as an offset for sampling effort
#This is failing with error saying V0 is not positive definite...
#K = ncol(XData)
#m.spatial <- setPriors(m.spatial,
#V0 = diag(c(rep(1, K-1), 0.00001)),
#f0 = K+1,
#mGamma = c(rep(0, K-1), 1),
#UGamma = diag(c(rep(1, K-1), 1e-8))
#)

#Set MCMC sampling paramteres and run
nChains = 2
test.run = F
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


m.spatial.full_covars = sampleMcmc(m.spatial.full_covars, thin = thin, samples = samples, transient = transient,
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

preds = computePredictedValues(m.spatial.full_covars)
MF = evaluateModelFit(hM=m.spatial.full_covars, predY=preds)
hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$TjurR2),2)))

#######################
#Variance partioning
#the group variable assigns X covariates to different covariate groups
#so first look at the design matrix

head(m.spatial.full_covars$X)
#For our real data we fit an intercept and three continuous variables, so they can each be assigned separate groups
#If instead we had a categorical variable the levels could be assigned to a single group along with the intercept

VP = computeVariancePartitioning(m.spatial.full_covars, group = c(1,2,3,4,5,6,7,8), groupnames = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences"))
plotVariancePartitioning(m.spatial, VP = VP)

VP.2 = computeVariancePartitioning(m.spatial.full_covars, group = c(1,2,2,3,4,5,6,7), groupnames = c("intercept","climate", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences"))

VP.3 = computeVariancePartitioning(m.spatial.full_covars, group = c(1,2,2,3,4,5,5,6), groupnames = c("intercept","climate", "duration infection", "DBH", "disease", "total sequences"))

####################
#Plot variance partioning
#If included a square of continuous variable, negative response would indicate intermediate niche optimum (i.e., abundance goes up initially but then goes down)
postBeta = getPostEstimate(m.spatial.full_covars, parName = "Beta")
plotBeta(m.spatial.full_covars, post = postBeta, param = "Support",
plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
#This can also be mapped on a tree with plotTree = T, but then tree must be included in model

############################################
#transform VP and postBeta object for ggplot
#VP cpntains R2 vals and postBeta contains support (i.e. alpha)

VP.vals = data.frame(VP$vals)
colnames(VP.vals) = c("Nf", "Nd")
VP.vals = VP.vals[1:8,] #this removes the random effect

##################
#Multiply proportion variance by Turj R2
VP.vals$Nf = VP.vals$Nf * MF$TjurR2[1]
VP.vals$Nd = VP.vals$Nd * MF$TjurR2[2]


VP.vals$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences")
VP.vals.long = VP.vals %>%
    pivot_longer(-variable, names_to = "ASV", values_to = "R2")



#Transform R2 based on positive or negative response
postBeta.mean = data.frame(postBeta$mean)
colnames(postBeta.mean) = c("Nf", "Nd")
postBeta.mean$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences")
postBeta.mean.long = postBeta.mean %>%
pivot_longer(-variable, names_to = "ASV", values_to = "betaMean")

#for(i in 1:length(postBeta.mean.long$mean)){
#    if(postBeta.mean.long$mean[i] < 0){
#        VP.vals.long$R2[i] = VP.vals.long$direction[i] = -1
#    }
#}


postBeta.support = data.frame(postBeta$support)
colnames(postBeta.support) = c("Nf", "Nd")
postBeta.support$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences")
postBeta.support.long = postBeta.support %>%
pivot_longer(-variable, names_to = "ASV", values_to = "support")


postBeta.supportNeg = data.frame(postBeta$supportNeg)
colnames(postBeta.supportNeg) = c("Nf", "Nd")
postBeta.supportNeg$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences")
postBeta.supportNeg.long = postBeta.supportNeg %>%
pivot_longer(-variable, names_to = "ASV", values_to = "supportNeg")


VP.vals.support = full_join(VP.vals.long, postBeta.support.long) %>%
    full_join(., postBeta.supportNeg.long) %>%
    full_join(., postBeta.mean.long)
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

VP.vals.support$variable = factor(VP.vals.support$variable, levels = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "duration infection", "DBH", "beech scale", "crown dieback", "total sequences"))

VP.vals.support[VP.vals.support["ASV"] == "Nf", "ASV"] = "N. faginata"
VP.vals.support[VP.vals.support["ASV"] == "Nd", "ASV"] = "N. ditissima"

#p1 = ggplot(VP.vals.support %>% filter(variable != "intercept"), aes(ASV, variable, fill = R2, color = P.val)) +
#geom_tile(size = 1, height = 0.975, width = 0.975) +
#scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0) +
#scale_color_manual(values = c("P<0.05" = "black", "n.s." = "white")) +
#guides(color = guide_legend(override.aes = list(fill = "grey"))) +
#my_gg_theme +
#labs(
#    x = "",
#    y = "",
#    fill = expression(paste("R"^2)),
#    color = "support",
#    title = "HMSC variance partitioning"
#) +
#theme(
#    legend.title = element_text(size = 20),
#    axis.text = element_text(size = 18)
#)


p1 = ggplot(VP.vals.support %>% filter(variable != "intercept"), aes(ASV, variable, size = R2, fill = betaMean, color = P.val)) +
geom_point(shape = 21) +
scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0) +
#scale_fill_gradientn(colours = c("#2c7bb6", "white", "#d7191c"),
#values = scales::rescale(c(-0.25, -.15,-0.1,-0.05,-0.025, 0 ,0.025, 0.05, 0.1,0.25)))  +
scale_color_manual(values = c("P<0.05" = "black", "n.s." = "white")) +
guides(color = guide_legend(override.aes = list(fill = "grey"))) +
scale_size_continuous(range= c(5,25), breaks = c(0,0.05,0.15), limits = c(0,0.2))+
my_gg_theme +
guides(
    color = guide_legend(override.aes = list(size = 5, shape = 21, fill = "dark grey")),
    size = guide_legend(override.aes = list(shape = 16, color = "dark grey"))
) +
labs(
x = "",
y = "",
size = expression(paste("R"^2)),
color = "Support",
fill = "Slope"
#title = "HMSC variance partitioning"
) +
theme(
legend.title = element_text(size = 20),
axis.text = element_text(size = 18),
legend.key = element_rect(fill = "white")
)

pdf("HMSC/Nf_Nd_variance_partitioning.bin.spatial.mod.pdf", width = 8.5, height = 8)
p1
dev.off()

write.table(VP.vals.support, "Nf_Nd_variance_partitioning.bin.spatial.VP_table.txt", row.names = F, sep = "\t", quote = F)
#scale_fill_continuous(limits = c(-1,1))

####################
####################
#Estimated residual var between spp

OmegaCor = computeAssociations(m.spatial.full_covars)
supportLevel = 0.5 #this is alpha = 0.05
toPlot = ((OmegaCor[[1]]$support>supportLevel)
+ (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
col=colorRampPalette(c("blue","white","red"))(500),
tl.cex=.6, tl.col="black",
title=paste("random effect level:", m.spatial$rLNames[1]), mar=c(0,0,1,0))


