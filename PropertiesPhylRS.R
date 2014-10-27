## The function for grafting missing species to a phylogeny (expanded tree approach)
## by François Michonneau is provided at:

https://stat.ethz.ch/pipermail/r-sig-phylo/2012-January/001826.html
# Clade membership from taxonomic data:
##  Nesolagus_timminsi: sister to Nesolagus_netscheri
##  Ochotona_hoffmanni: sister to Ochotona_hyperborea
##  Pronolagus_saundersiae: sister to Pronolagus_rupestris
##  Sylvilagus_obscurus: sister to Sylvilagus_floridanus
##  Sylvilagus_varynaensis: sister to Sylvilagus_brasiliensis


## Note: set the working directory path to a folder with all the files provided in the supplementary data
## Properties of the lagomorph phylogeny

## Comparing the number of monotypic species with the number of genera 
## for all mammalian Orders

library (plyr)

# count species per genus from the IUCN mammal taxonomy
Mammals_Higher_Taxonomy <- read.csv ("Mammals_Higher_Taxonomy.csv")
freqTable <- count (Mammals_Higher_Taxonomy, c ("Order","Genus"))
# count no. of genera per Order
ngen <- ddply (freqTable, .(Order), summarize, noGen=length(unique(Genus)))
# identify monotypic genera
mono1 <- ddply (freqTable, .(Order), summarize, monot1=length(freq[freq<2]))
# plotting
plot (ngen$noGen, mono1$monot1, xlab="genera in order", ylab="number of monotypic species",xlim=c(0,200),ylim=c(0,80),cex=0.7)
title(main="a)",adj=0,cex.main=1, font.main = 1 )
# Values for Lagomorphs
points (12,7,pch=24, bg="black")


## Comparing clade size distributions

library (apTreeshape)
library (sm)

# Generate ERM trees (Yule model, equal probabilities)
scspec <- lapply (rtreeshape(100, tip.number=87, model="yule"),
                 FUN=smaller.clade.spectrum)
scspecDF <- ldply (scspec)
vecSim <- scspecDF[, 1]

# Read phylo tree using "read.nexus" from the "ape" package
# both files provided in supplement
## "lphyloST" is the lagomorph clade from the Rolland et al. (2014) redated tree with five species grafted 
## "lphyloPASTIS" is the Ge et al. (2013) mitochondrial phylogeny completed using PASTIS package

lphylo  <- read.nexus("lphyloST.nex") 

# Count clade sizes in Lagomorph phylogeny
csizeST <- smaller.clade.spectrum (as.treeshape(lphylo))
CSphST <- csizeST[,1]

# Real and simulated values
allCsizes <- c(CSphST,vecSim)
labelsVec <-c (rep("tree", 86), rep("simulated", 8600))
labelsVecF <- factor (labelsVec)

# Comparing densities
sm.options (col.band="gray")
sm.density.compare (allCsizes, labelsVecF, model="equal", xlab="clade sizes", col=c('black', 'black'), lty=c(1, 3), col.band="gray")
title(main="b)",adj=0,cex.main=1, font.main = 1 )
# title (main="Density distribution of clade sizes")
# add legend to plot
legend (x ='topright', lty=c(1, 3), legend=c("simulated", "tree"))

# Testing tree imbalance

colless.test (as.treeshape(lphylo), model="yule", alternative="greater", n.mc=1000)

# Phylogenetic signal in extinction risk 

library (caper)

# read full dataset

fullDataset <- read.csv ("fullDataset.csv")
row.names (fullDataset) <-fullDataset$inPhylogeny 

# Create comparative.data object
lagomsDFc <- comparative.data (lphylo, fullDataset, names.col= "inPhylogeny", na.omit=TRUE)

# D test 
DstatLagoms <- phylo.d (lagomsDFc, binvar=Threatened)
DstatLagoms

# Repeat test assuming Data Deficient species are either threatened or not threatened

DDspecies <- which (is.na (fullDataset$Threatened))

# Assume DD species are threatened
fullDataDDth <- fullDataset
fullDataDDth$Threatened[DDspecies] <- 1
lagomsDFcDDth <- comparative.data (lphylo, fullDataDDth, names.col= "inPhylogeny", na.omit=TRUE)
DstatLagomsDDth <- phylo.d (lagomsDFcDDth, binvar=Threatened)
DstatLagomsDDth

# Assume DD species are not threatened
fullDataDDnt <- fullDataset
fullDataDDnt$Threatened[DDspecies] <- 0
lagomsDFcDDnt <- comparative.data (lphylo, fullDataDDnt, names.col= "inPhylogeny", na.omit=TRUE)
DstatLagomsDDnt <- phylo.d (lagomsDFcDDnt, binvar=Threatened)
DstatLagomsDDnt

