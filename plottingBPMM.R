# Code to plot significant results from extinction risk model
require (ggplot2)
require (postMCMCglmm)
require(devtools)

# using code and the package "postMCMCglmm" by Joshua Wiley 
# Obtains Predicted probabilities from MCMCglmm ordinal probit models
install_github ("postMCMCglmm", "JWiley")

# Extinction risk model object must be in the workspace
# new data for prediction
# vectors for each variable

Dim1 = seq (min(fullDataBPMM$Dim.1), max(fullDataBPMM$Dim.1), length.out=500)
Dim2 = seq (min(fullDataBPMM$Dim.2), max(fullDataBPMM$Dim.2), length.out=500)
minhpd = seq (min(fullDataBPMM$MinHPD15), max(fullDataBPMM$MinHPD15), length.out=500)
converted = seq (min(fullDataBPMM$pcntConverted), max(fullDataBPMM$pcntConverted),length.out=500)
bsizes = seq (min(fullDataBPMM$bsizes), max(fullDataBPMM$bsizes),length.out=500)

# Plotting effects of % points in converted habitat
# other varaibles held constant

datlagPcT <- as.matrix (data.frame("(Intercept)"= 1,
                               "Dim.1"= mean (Dim1),
                               "Dim.2" = mean (Dim2), 
                               "MinHPD15" = mean (minhpd),
                               "pcntConverted" = converted,
                               "bsizes" = mean (bsizes)))                               

## calculate predicted values (fixed effects only)
predlagPcT <- predict2(fullmodERisk, X=datlagPcT, Z = NULL, use = "all", type = "response", varepsilon = 1)

## summarize values
spredlagPcT <- summary (predlagPcT)

## combine predicted probs + HPD intervals with prediction data
preddatlagPcT <- as.data.frame (cbind(do.call(rbind, rep(list(datlagPcT), 5)), do.call (rbind, spredlagPcT)))
## indicator for which level of the outcome
preddatlagPcT$outcome <- factor (rep(c("LC","NT","VU","EN","CR"), each = nrow(datlagPcT)), 
                                 levels=c("LC", "NT", "VU", "EN", "CR"))

## plot it

pctConvplot <- ggplot (preddatlagPcT, aes(x = pcntConverted, y = M, colour = outcome)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = outcome), alpha = .25) +
  geom_line(size=2)+ylab("probability")+xlab("percent converted habitat")+
  scale_color_manual(values=c("#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"),name="Red List Status")+
  scale_fill_manual(values=c("#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), guide=FALSE)+
  theme_classic()+theme(legend.position = "none")

############################
# plot effect of Human Population Density
# other variables held constant

datlagHPD <- as.matrix(data.frame("(Intercept)"= 1,
                                  "Dim.1"= mean(Dim1),
                                  "Dim.2" = mean(Dim2), 
                                  "MinHPD15" = minhpd,
                                  "pcntConverted" = mean(converted),
                                  "bsizes" = mean (bsizes)))




## get the predicted values (fixed effects only)
predlagHPD <- predict2 (fullmodERisk, X=datlagHPD, Z = NULL, use = "all", type = "response", varepsilon = 1)

## summarize them
spredlagHPD <- summary (predlagHPD)

## combine predicted probs + HPD intervals with prediction data
preddatlagHPD <- as.data.frame (cbind(do.call(rbind, rep(list(datlagHPD), 5)), do.call(rbind, spredlagHPD)))
## indicator for which level of the outcome
preddatlagHPD$outcome <- factor (rep(c("LC", "NT", "VU", "EN", "CR"), 
                                     each = nrow(datlagHPD)), levels=c("LC", "NT", "VU", "EN","CR"))


## plot it

hpdplot <- ggplot (preddatlagHPD, aes(x = MinHPD15, y = M, colour = outcome)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = outcome), alpha = .25) +
  geom_line(size=2)+ylab("probability")+xlab("minimum human population density")+
  scale_color_manual(values=c("#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"),name="Red List Status")+
  scale_fill_manual(values=c("#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), guide=FALSE)+
  theme_classic()+theme(legend.position = "none")

# side by side plots with the multiplot function from the R Graphics Cookbook by Winston Chang

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# side by side plots (legend can be un-suppressed or drawn manually)
multiplot(hpdplot,pctConvplot,cols=2)
  