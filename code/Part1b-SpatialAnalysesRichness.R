
load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")

library(rgeos)
library(maptools)

#### Spatially plot the positives, negatives, and neutrals ####

#these are sites with no changes in diversity through time
spatialNeutrals<- latlongs[match(colnames(richness)[nonSig], latlongs$Handle),]
coordinates(spatialNeutrals)<- ~Longitude+Latitude
neutralCentroid<- gCentroid(spatialNeutrals)
plot(spatialNeutrals, col="lightgray")
points(neutralCentroid, col="darkgray", pch=16, cex=1.5)

#these are sites with decreases in diversity through time
spatialNegatives<- latlongs[match(colnames(richness)[sigNeg], latlongs$Handle),]
coordinates(spatialNegatives)<- ~Longitude+Latitude
negativeCentroid<- gCentroid(spatialNegatives)
plot(spatialNegatives, col="red", add=T)
points(negativeCentroid, col="darkred", pch=16, cex=1.5)

#these are sites with increases in diversity through time
spatialPositives<- latlongs[match(colnames(richness)[sigPos], latlongs$Handle),]
coordinates(spatialPositives)<- ~Longitude+Latitude
positiveCentroid<- gCentroid(spatialPositives)
plot(spatialPositives, col="blue", add=T)
points(positiveCentroid, col="darkblue", pch=16, cex=1.5)


warmPalette<- colorRampPalette(c("red", "orange"))(43)
negCols<- apply(richness[,sigNeg], 2, function(x) max(which(!is.na(x))))-1
#reds= time series starts younger, orange=time series starts older

coolPalette<- colorRampPalette(c("green", "blue"))(43)
posCols<- apply(richness[,sigPos], 2, function(x) max(which(!is.na(x))))-1
#greens= time series starts younger, blues=time series starts older


data(wrld_simpl)
pdf(file="figures/Map-RichnessThruTime.pdf", height=9, width=10)
par(mfrow=c(1,1))
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)

plot(spatialNeutrals, pch=3, cex=0.75, col="lightgray", add=T)
plot(spatialNegatives, pch=15, col=warmPalette[negCols], add=T)
plot(spatialPositives, pch=17, col=coolPalette[posCols], add=T)

legend("bottom", pch=c(rep(15, 8), 3, rep(17, 8)), col=c(warmPalette[seq(43,1, -6)], "gray", (coolPalette[seq(1,43, 6)])), 
       legend=timeKYR[c(seq(43,1, -6)], 0, seq(1,43, 6))], cex=0.75, horiz=T)


dev.off()


pdf(file="figures/Map-RichnessThruTime-All.pdf", height=9, width=10)
par(mfrow=c(2,2), mar=c(4,4,4,4)+0.1)

plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (i in 1:nrow(richness)){
  lines(richness[i,]~allTimes, col="gray") 
}
lines(richnessMeans~allTimes, col="black", lwd=2)

plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[nonSig[k],]~allTimes, col="gray")
}
lines(nonSigRichnessMeans~allTimes, col="black", lwd=2)

coolPalette<- colorRampPalette(c("green", "blue"))(21)
plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigPos)){
  lines(richness[sigPos[k],]~allTimes, col=coolPalette[max(which(!is.na(richness[sigPos[k],])))-1])
}
lines(posRichnessMeans~allTimes, col="black", lwd=2)

warmPalette<- colorRampPalette(c("red", "orange"))(21)
plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigNeg)){
  lines(richness[sigNeg[k],]~allTimes, col=warmPalette[max(which(!is.na(richness[sigNeg[k],])))-1])
}
lines(negRichnessMeans~allTimes, col="black", lwd=2)


dev.off()




#### plot spatial patterns of richness change ####

par(mfrow=c(5,5))
for (j in 1:nrow(siteRichChanges)){
  specificLocs<- siteLocs[match(colnames(siteRichChanges), sites)]
  plot(siteRichChanges[j,]~specificLocs@coords[,2], pch=16) #plot richness change as a function of latitude
  summary(lm(siteRichChanges[j,]~specificLocs@coords[,2]))  
}


#### New analyses ####
# superimpose present-day biome map to delimit regions, then examine richness changes within them
# determine whether same set of species are involved in majority of colonizations or extirpations

