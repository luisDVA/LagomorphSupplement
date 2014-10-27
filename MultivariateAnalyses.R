## Multivariate Analyses

# Principal Components Analysis of environmental variables
## Note: set the working directory path to a folder with all the files provided in the supplementary data

# load table of environmental data  (some variables are already transformed for normality)  
fullClimT <- read.csv ("environmentalData.csv", row.names=1)

# visualize transformed climate data
par (mfrow=c(4,4))
lapply (names(fullClimT), function(xn) hist(fullClimT[[xn]], main = paste("Histogram of", xn))) 
dev.off ()

# PCA
library (FactoMineR)
lagPCA <- PCA (fullClimT)
pcaCors <- round (lagPCA$var$cor, 4) # loadings
pcaContr <- round (lagPCA$var$contrib, 4) # contributions
# explore PCA results
pcaCors
pcaContr
summary (lagPCA)
# first two component axes already merged with the complete dataset

## Fitting phylogenetic generalized linear mixed models 
## aka  Bayesian Phylogenetic Mixed Models (BPMM) following Botero et al. (2013) 

require (MCMCglmm)
require (plyr)
require (geiger)
require (picante)
## DATA PREPARATION
#  read full dataset

fullDataset <- read.csv ("fullDataset.csv")
row.names (fullDataset) <-fullDataset$inPhylogeny 

# drop data deficient species from dataset for extinction risk model
fullDataBPMM  <- fullDataset[!is.na(fullDataset$ERisk), ]
# read phylogeny
lphylo  <- read.nexus ("lphyloST.nex")
# drop data deficient species from phylogeny
lphyloR <- treedata (lphylo,fullDataBPMM)$phy

# create phylogenetic covariance matrix
Ainv <- inverseA(lphyloR)$Ainv

#  define priors
#  B prior follows suggested prior for ordinal regression by Gelman (2008), modified by J. Hadfield see function documentation
prior <- list (B=list(mu=rep(0,6), V=gelman.prior(~ Dim.1+Dim.2+MinHPD15+pcntConverted+bsizes, data=fullDataBPMM, scale=1+pi^2/3)),
               R = list(V = 1, fix= 1), 
               G = list(G1= list(V=1e-8, nu=-1)))                            
# run model 
fullmodERisk <- MCMCglmm(ERisk ~ Dim.1+Dim.2+MinHPD15+pcntConverted+bsizes, data=fullDataBPMM, random=~inPhylogeny, ginverse=list(inPhylogeny=Ainv), family = "ordinal", prior= prior, nitt=2500000, burnin=500000, thin=2000)

# evaluate convergence
heidel.diag (fullmodERisk$Sol) # Heidelberg and Welch (1983) diagnostic test
autocorr (fullmodERisk$Sol) # autocorrelation of succesive samples (should be near 0)
plot (fullmodERisk) # visual inspection of mixing properties of the MCMC chain

# explore results
summary (fullmodERisk)


## Correlates of evolutionry distinctiveness
# Calculate ED for the given phylogeny
evoldist <- evol.distinct (lphylo, type="fair.proportion")
evoldist$ED <- evoldist$w
evoldist$w <- NULL
fullDataset <- merge (fullDataset, evoldist, by.x="inPhylogeny", by.y="Species")

# phylogenetic covariance matrix
AinvED <- inverseA (lphylo)$Ainv
# prior specification
priorED = list (B = list(mu = rep(0,5), V = diag(5)), R = list(V= (var(fullDataset$ED)/2),nu=1), G = list( G1 = list(V=1, nu=1, alpha.mu=0, alpha.V=1000 )))              
# run model
modED <- MCMCglmm (ED ~ MammRich+Dim.1+Dim.2+log(EOO), data=fullDataset, random=~inPhylogeny, ginverse=list(inPhylogeny=AinvED), prior= priorED, nitt=2500000, burnin=500000, thin=2000)


## Predicting threat status with model object
#  using code and the package postMCMCglmm by Joshua Wiley 
#  Obtains Predicted probabilities from MCMCglmm ordinal probit models
require (devtools)
install_github ("postMCMCglmm", "JWiley")
require (postMCMCglmm)

# Predict probabilities
yhat <- predict2 (fullmodERisk, type="response", Z=NULL)
# summary of predicted probabilities
sumyhat <- summary (yhat)
# combine
longsum <- do.call (rbind.data.frame, sumyhat)
#round 
longsum <- round (longsum, 3)
# create a level indicator
longsum$Level <- factor (rep(1:5, each = nrow(sumyhat[[1]])))
spnams<-as.character(fullDataBPMM$inPhylogeny)
longsum$sp <- rep (spnams, 5)
# Arrange data into Data Frame
splif <- split (longsum, longsum$sp)
# Which category has the highest probability
predictedRL <- lapply (splif, function(splif)splif$Level[which.max(splif$M)])
predictedRL <- ldply (predictedRL)              
predictedRL
