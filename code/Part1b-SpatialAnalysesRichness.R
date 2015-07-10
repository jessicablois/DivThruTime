load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")

#spatially plot the positives, negatives, and neutrals
library(sp)
library(rgeos)

#siteLocs<- read.delim("~/Dropbox/Research/Community Paleomodels/projects/pollen/data/All.site.data-withagemodel-finalv2.txt", sep="\t", header=T)

# sigPos #these are sites that show DECREASING RICHNESS through time
# sigNeg #these are sites that show INCREASING RICHNESS through time

#these are sites with no changes in diversity through time
spatialNeutrals<- latlongs[match(rownames(richness)[nonSig], latlongs$Handle),]
coordinates(spatialNeutrals)<- ~Longitude+Latitude
neutralCentroid<- gCentroid(spatialNeutrals)
plot(spatialNeutrals, col="lightgray")
#points(neutralCentroid, col="darkgray", pch=16, cex=1.5)

#these are sites with increases in diversity through time
spatialNegatives<- latlongs[match(rownames(richness)[sigNeg], latlongs$Handle),]
coordinates(spatialNegatives)<- ~Longitude+Latitude
negativeCentroid<- gCentroid(spatialNegatives)
plot(spatialNegatives, col="blue", add=T)
#points(negativeCentroid, col="darkblue", pch=16, cex=1.5)

#these are sites with decreases in diversity through time
spatialPositives<- latlongs[match(rownames(richness)[sigPos], latlongs$Handle),]
coordinates(spatialPositives)<- ~Longitude+Latitude
positiveCentroid<- gCentroid(spatialPositives)
plot(spatialPositives, col="red", add=T)
#points(positiveCentroid, col="red", pch=16, cex=1.5)


warmPalette<- colorRampPalette(c("red", "orange"))(21)
negCols<- apply(richness[sigNeg,], 1, function(x) max(which(!is.na(x))))-1
#reds= time series starts younger, orange=time series starts older

coolPalette<- colorRampPalette(c("green", "blue"))(21)
posCols<- apply(richness[sigPos,], 1, function(x) max(which(!is.na(x))))-1
#greens= time series starts younger, blues=time series starts older

library(maptools)
data(wrld_simpl)
pdf(file="figures/Map-RichnessThruTime.pdf", height=9, width=10)
par(mfrow=c(1,1))
plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)

plot(spatialNeutrals, pch=3, cex=0.75, col="lightgray", add=T)
plot(spatialNegatives, pch=15, col=warmPalette[negCols], add=T)
plot(spatialPositives, pch=17, col=coolPalette[posCols], add=T)

legend("bottom", pch=c(rep(15, 7), 3, rep(17, 7)), col=c(warmPalette[seq(21,1, -3)], "gray", (coolPalette[seq(3,21, 3)])), 
       legend=c(seq(21,1, -3), 0, seq(3,21, 3)), cex=0.75, horiz=T)


# legend("bottomright", pch=16, cex=0.75,
#        col=c("red", "blue", "gray"), 
#        c("Richness increases",
#          "Richness decreases", 
#          "No richness change"))

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


