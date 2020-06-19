
library(Hmsc)
library(MASS)
set.seed(6)

n = 100
ns = 2

plot_coords = matrix(runif(2*10), ncol = 2)
xycoords = matrix(rep(t(plot_coords),10),ncol=2, byrow = T)

colnames(xycoords) = c("x-coordinate","y-coordinate")
rownames(xycoords) = paste0("sample_", seq(1:n))
colnames(plot_coords) = c("x-coordinate","y-coordinate")
rownames(plot_coords) = paste0("site_", letters[seq(1:10)])

beta1 = c(-1,1)
alpha = rep(0,ns)
beta = cbind(alpha,beta1)
x = cbind(rep(1,n),rnorm(n))
Lf = x%*%t(beta)
sigma.spatial = c(2)
alpha.spatial = c(0.35)
Sigma = sigma.spatial^2*exp(-as.matrix(dist(xycoords))/alpha.spatial)
eta1 = mvrnorm(mu=rep(0,n), Sigma=Sigma)
lambda1 = c(1,2)
Lr = eta1%*%t(lambda1)
L = Lf + Lr
#y = as.matrix(L + matrix(rnorm(n*ns),ncol=ns))
yprob = 1*((L +matrix(rnorm(n*ns),ncol=ns))>0)
XData = data.frame(x1=x[,2])


#Study desgin, random effect, and model specification

studyDesign = data.frame(
    sample = as.factor(paste0("sample_", seq(1:n))),
    site = as.factor(paste0("site_", letters[seq(1:10)])),
    spatial = as.factor(paste0("sample_", seq(1:n)))
)

rL.sample = HmscRandomLevel(units = as.factor(studyDesign$sample))
rL.site = HmscRandomLevel(units = as.factor(studyDesign$Site))
rL.spatial = HmscRandomLevel(sData = xycoords)


m.spatial = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
studyDesign=studyDesign, ranLevels=list("sample"=rL.sample, "site" = rL.site, "spatial" = rL.spatial),distr="probit")

m.spatial = sampleMcmc(m.spatial, thin = 1, samples = 10, transient = 5,
nChains = 2, verbose = 0,updater=list(GammaEta=FALSE))

#spatial effect defined at level of site
rL.spatial = HmscRandomLevel(sData = plot_coords)


studyDesign = data.frame(
sample = as.factor(paste0("sample_", seq(1:n))),
site = as.factor(paste0("site_", letters[seq(1:10)])),
spatial = as.factor(paste0("site_", letters[seq(1:10)]))
)


m.spatial = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
studyDesign=studyDesign, ranLevels=list("sample"=rL.sample, "site" = rL.site, "spatial" = rL.spatial),distr="probit")

m.spatial = sampleMcmc(m.spatial, thin = 1, samples = 10, transient = 5,
nChains = 2, verbose = 0,updater=list(GammaEta=FALSE))


#Updated from #47

plot_coords = matrix(runif(2*10), ncol = 2)
rownames(plot_coords) = paste0("site_", letters[seq(1:10)])

studyDesign = data.frame(
    sample = as.factor(paste0("sample_", seq(1:n))),
    site_spatial = as.factor(
        rep(paste0("site_", letters[seq(1:10)]), 10)
    )
)

rL.sample = HmscRandomLevel(units = paste0("sample_", seq(1:n)))
rL.site_spatial = HmscRandomLevel(sData = plot_coords)

m.spatial = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
studyDesign=studyDesign, ranLevels=list("sample"=rL.sample, "site_spatial" = rL.site_spatial),distr="probit")

m.spatial = sampleMcmc(m.spatial, thin = 1, samples = 10, transient = 5,
nChains = 2, verbose = 0,updater=list(GammaEta=FALSE))


