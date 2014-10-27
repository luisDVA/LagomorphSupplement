# Diversification analysis

library(geiger)
library(zoo)
library(plyr)
require(devtools)

install_github ("PANDA", "hmorlon")

# list of diversification models to be fitted, from Morlon et al. (2011)
# B constant
# BD constant
# B variable exp
# B variable lin
# B variable exp, D constant 
# B variable lin, D constant
# B constant, D variable exp
# B constant, D variable lin
# BD variable exp
# BD variable lin



# sylvilagus
# sylvilagus species list from full dataset
sylv <- subset(fullDataset,Genus == "Sylvilagus")
row.names(sylv) <- sylv$inPhylogeny
# extract subtree
sylvPhyl <- treedata(lagPS,sylv)$phy
# relative total time for use in all the following diversification analyses 
tot_time<-max(node.age(sylvPhyl)$ages)

# for all remaining species repeat the same code but using a tree that excludes Sylvilagus
# all lagomorphs except sylvilagus
# get species list from full dataset
noSyl <- subset(fullDataset,Genus != "Sylvilagus")
row.names(noSyl) <- noSyl$inPhylogeny
noSylPhyl <- treedata(lagPS,noSyl)$phy


# 1 # B constant
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.09)
mu_par<-c()
Bconstant <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                    cst.lamb = TRUE, fix.mu=TRUE)


# 2 # BD constant
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.09)
mu_par<-c(0.01)
BDconstant<- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                    cst.lamb = TRUE, cst.mu=T,fix.mu=FALSE)


# 3 # B variable exponential
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
BvariableE <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                     expo.lamb = TRUE, fix.mu=T)

# 4 # B variable lin
f.lamb <-function(t,y){y[1] + y[2] * t}
f.mu<-function(t,y){0}
lamb_par<-c(0.09, 0.001)
mu_par<-c()
BvariableL<- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par, 
                    fix.mu=T)


# 5 # B variable exp, D constant 
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.05, 0.01)
mu_par<-c(0.01)
BvariableEDconstant <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                              expo.lamb = TRUE,cst.mu=T, fix.mu=F)


# 6 # B variable lin, D constant
f.lamb <-function(t,y){y[1] + y[2] * t}
f.mu<-function(t,y){0}
lamb_par<-c(0.09, 0.001)
mu_par<-c(0.01)
BvariableLDconstant <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                              ,cst.mu=T,fix.mu=F)

# 7 # B constant, D variable exp                 
f.lamb <- function(t,y){y[1]}
f.mu <- function(t,y){y[1] * exp(y[2] * t)}
lamb_par<-c(0.09)
mu_par<-c(0.05, 0.01)        
BconstantDvariableE <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                              cst.lamb=T,expo.mu=T,fix.mu=F)

# 8 # B constant, D variable lin
f.lamb <- function(t,y){y[1]}
f.mu <- function(t,y){y[1] + y[2] * t}
lamb_par<-c(0.09)
mu_par<-c(0.05, 0.01)        
BconstantDvariableL <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                              cst.lamb=T,fix.mu=F)

# 9 # BD variable exp
f.lamb <- function(t,y){y[1] * exp(y[2] * t)}
f.mu <- function(t,y){y[1] * exp(y[2] * t)}
lamb_par<-c(0.05, 0.01) 
mu_par<-c(0.05, 0.01)        
BDvariableE <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                      expo.lamb=T,expo.mu=T,fix.mu=F)

# 10 # BD variable lin
f.lamb <- function(t,y){y[1] + y[2] * t}
f.mu <- function(t,y){y[1] + y[2] * t}
lamb_par<-c(0.09, 0.001)
mu_par<-c(0.05, 0.01)   
BDvariableLin <- fit_bd(sylvPhyl,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                        fix.mu=F)

# Find the best model
modelAiccVals <- c(Bconstant$aicc, BDconstant$aicc, BvariableE$aicc,
                   BvariableL$aicc, BvariableEDconstant$aicc, BvariableLDconstant$aicc,
                   BconstantDvariableE$aicc,BconstantDvariableL$aicc,BDvariableE$aicc,
                   BDvariableLin$aicc)

# figure 3 was created using the outputs from a modification to the diversity through time function

# This function stores the inferred number of species at each time unit 

forplotdtt <- function (fit.bd, tot_time, N0 = 1) 
{
  t <- seq(0, tot_time, length.out = 100)
  r <- function(t) {
    -fit.bd$f.lamb(t) + fit.bd$f.mu(t)
  }
  R <- function(s) {
    .Integrate(r, 0, s)
  }
  N <- N0 * exp(Vectorize(R)(t))
  plottingargs <- list(t,N)
  return (plottingargs)
  
}

environment(forplotdtt) <- asNamespace('PANDA')

# diversity through time curve for Sylvilagus

syltN <- forplotdtt(sylvilagusFinalBvarE,18,17)


# remaining lagomorphs

lagremtN <- forplotdtt(lagRemBDvarL,50,70)


# manipulate data for plotting
syldt <- cbind.data.frame(syltN)
names(syldt) <- c("t","Nsyl")
lagremdt <- cbind.data.frame(lagremtN)
names(lagremdt) <- c("t","Nlagrem")
divertt <- merge(lagremdt,syldt,by="t",all.x=T,all.y=T)
divertt$paraacc <- divertt$Nsyl
divertt$paraacc[133:196] <- 0
divertt$paraacc <- na.approx(divertt$paraacc)
divertt$accsum <- rowSums(divertt[,c(2,4)])
diverttplot <- divertt
diverttplot$paraacc <- NULL
melted <- melt(diverttplot,id.vars = "t" )

# plot inferred diversity 

ggplot(data=melted,aes(x=t,y=value,color=variable,linetype=variable))+
  geom_smooth()+scale_x_reverse()+ylab("number of species")+
  xlab("time (Mya)")+theme_classic()+scale_linetype_manual(values=c(1,1,5))+
  scale_color_manual(values=c("#fdb863","#b2abd2","#5e3c99"),labels=c("all other lagomorphs","Sylvilagus","cumulative diversity"))+
  theme(legend.title=element_blank())+guides(linetype=FALSE)
