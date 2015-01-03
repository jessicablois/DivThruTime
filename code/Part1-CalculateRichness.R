#### Part 1: BASIC PATTERNS ####
source("code/Part0-GlobalParam.R")
source("code/DivThruTimeFunctions.R")

# Calculate richness at each site
richness<- matrix(ncol=length(allTimes), nrow=length(sites))
colnames(richness)<- allTimes
rownames(richness)<- sites

for (i in 1:length(sites)){
  sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  datTimes<- dat[,1] # pull out time periods
  
  if (length(datTimes)<=1){
    richness[i,]<- NA     
  }else{
    matchedTimes<- allTimes[na.omit(match(datTimes, allTimes))]
    
    minTime<- min(matchedTimes, na.rm=T)
    maxTime<- max(matchedTimes, na.rm=T)
    
    richness[i,]<- calcSiteRichness(dat, minTime, maxTime, pollenThreshold, interval) 
  }
}

richnessMeans<- colMeans(richness, na.rm=T)

# Calculate sample size
sampSize<- apply(richness, 2, sampleSize)

#### at which sites is diversity increasing or decreasing? ####
timesKYR<- allTimes/1000
linearSlopes<- vector(length=nrow(richness))
pvals<- vector(length=nrow(richness))
for (i in 1:nrow(richness)){
  tempSpp<- richness[i,]
  if (length(which(!is.na(richness[i,])))>2){ #if there are more than two points, fit a model and store the slope
    lin.mod <- lm(tempSpp~timesKYR)
    plot(tempSpp~timesKYR, pch=16, ylim=c(0, max(richness, na.rm=T)))
    abline(lin.mod)
    linearSlopes[i]<- lin.mod$coefficients[2]
    pvals[i]<- summary(lin.mod)$coefficients[8]
  }else{
    linearSlopes[i]<- NaN
    pvals[i]<- NaN
  }
}

significants<- which(pvals<=0.05)
positives<- which(linearSlopes>0) #these are sites that show DECREASING RICHNESS through time
negatives<- which(linearSlopes< 0) #these are sites that show INCREASING RICHNESS through time
neutrals<- which(linearSlopes==0)

sigPos<- intersect(significants, positives) #these are sites that show DECREASING RICHNESS through time
sigNeg<- intersect(significants, negatives) #these are sites that show INCREASING RICHNESS through time
nonSig<- seq(1:length(linearSlopes))[-union(sigPos, sigNeg)]

posRichnessMeans<- colMeans(richness[sigPos,], na.rm=T)
negRichnessMeans<- colMeans(richness[sigNeg,], na.rm=T)
nonSigRichnessMeans<- colMeans(richness[nonSig,], na.rm=T)


#### plot spatial patterns of richness change ####
siteRichChanges<- richness[,2:22]-richness[,1:21]
for (j in 1:ncol(siteRichChanges)){
  specificLocs<- siteLocs[match(rownames(siteRichChanges), sites)]
  plot(siteRichChanges[,j]~specificLocs@coords[,2], pch=16)
  summary(lm(siteRichChanges[,j]~specificLocs@coords[,2]))  
}


save(list=c("richness", "richnessMeans", "sampSize", "sigPos", "sigNeg", "nonSig", "siteRichChanges"), file="workspaces/richness.RData")


#### Plotting ####

pdf(file="figures/Correlation-Richness-NoSites.pdf", height=4, width=6)
plot(richnessMeans ~ sampSize, pch=16, xlab="Number of sites", ylab="Mean richness")
test1<- cor.test(richnessMeans, sampSize) #is richness correlated with sample size?
legend("bottomright", bty="n", paste("cor=",round(test1$estimate, 4), "; p=", round(test1$p.value, 4), sep=""), cex=0.75)
dev.off()

# Plot the mean genus richness across all sites and overlay sample size
pdf(file="figures/RichnessSampSizeThruTime-all.pdf", height=4, width=6)
par(mar=c(4,4,4,4)+0.1)
plot(richnessMeans~allTimes, type="l", xlim=c(21000,0), xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
par(new=T)
plot(sampSize~allTimes, pch=16, col="red", ylim=c(0, max(sampSize)), xlim=c(21000,0), axes=F, xlab="", ylab="")
axis(4, at=seq(0, max(sampSize), by=100), col.axis="red")
mtext("Number of Sites", 4, line=3, col="red")
dev.off()

#Plot each individual line
pdf(file="figures/RichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (i in 1:nrow(richness)){
  lines(richness[i,]~allTimes, col=rainbow(nrow(richness))[i]) 
}
lines(richnessMeans~allTimes, col="black", lwd=2)
dev.off()

#Plot each individual line- sig pos vs sig neg
pdf(file="figures/RichnessThruTime-all-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[nonSig[k],]~allTimes, col="gray")
}
for (k in 1:length(sigPos)){
  lines(richness[sigPos[k],]~allTimes, col="blue")
}
for (k in 1:length(sigNeg)){
  lines(richness[sigNeg[k],]~allTimes, col="red")
}
lines(richnessMeans~allTimes, col="black", lwd=2)
dev.off()

pdf(file="figures/RichnessThruTime-threePanels-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,3))
plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[nonSig[k],]~allTimes, col="gray")
}
lines(nonSigRichnessMeans~allTimes, col="black", lwd=2)

plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigPos)){
  lines(richness[sigPos[k],]~allTimes, col="blue")
}
lines(posRichnessMeans~allTimes, col="black", lwd=2)


plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigNeg)){
  lines(richness[sigNeg[k],]~allTimes, col="red")
}
lines(negRichnessMeans~allTimes, col="black", lwd=2)

dev.off()


#Plot each individual line

pdf(file="figures/RichnessThruTime-upper-vs-lower-withLines.pdf", height=4, width=6)
plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
for (i in 1:nrow(upperRich)){
  lines(upperRich[i,]~allTimes, col="red", lty=2) 
}
lines(upperRichMean~allTimes, col="red", lwd=3)

for (i in 1:nrow(lowerRich)){
  lines(lowerRich[i,]~allTimes, col="blue", lty=2) 
}
lines(lowerRichMean~allTimes, col="blue", lwd=3)
lines(richnessMeans~allTimes, col="black", lwd=2)
dev.off()

colorRampPalette(c("#FF0000FF", "#00FF00FF"))(11)
colorRampPalette(c("#00FFFFFF", "#FF00FFFF"))(11)

#### Plot the increasing vs decreasing sites ####
pdf(file="figures/RichnessThruTime-DECREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigPos)){
  plot(richness[sigPos[k],]~allTimes, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)),
       xlab="", ylab="Site Richness")
  mtext(rownames(richness)[sigPos[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-INCREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigNeg)){
  plot(richness[sigNeg[k],]~allTimes, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)),
       xlab="", ylab="Site Richness")
  mtext(rownames(richness)[sigNeg[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-nonSig.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(nonSig)){
  plot(richness[nonSig[k],]~allTimes, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)),
       xlab="", ylab="Site Richness")
  mtext(rownames(richness)[nonSig[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()



#spatially plot the positives, negatives, and neutrals
library(sp)
library(rgeos)

siteLocs<- read.delim("~/Dropbox/Research/Community Paleomodels/projects/pollen/data/All.site.data-withagemodel-finalv2.txt", sep="\t", header=T)

#these are sites with no changes in diversity through time
spatialNeutrals<- siteLocs[match(rownames(richness)[nonSig], siteLocs$Handle),]
coordinates(spatialNeutrals)<- ~Longitude+Latitude
neutralCentroid<- gCentroid(spatialNeutrals)
plot(spatialNeutrals, col="lightgray", add=T)
#points(neutralCentroid, col="darkgray", pch=16, cex=1.5)

#these are sites with increases in diversity through time
spatialNegatives<- siteLocs[match(rownames(richness)[sigNeg], siteLocs$Handle),]
coordinates(spatialNegatives)<- ~Longitude+Latitude
negativeCentroid<- gCentroid(spatialNegatives)
plot(spatialNegatives, col="blue", add=T)
#points(negativeCentroid, col="darkblue", pch=16, cex=1.5)

#these are sites with decreases in diversity through time
spatialPositives<- siteLocs[match(rownames(richness)[sigPos], siteLocs$Handle),]
coordinates(spatialPositives)<- ~Longitude+Latitude
positiveCentroid<- gCentroid(spatialPositives)
plot(spatialPositives, col="red", add=T)
#points(positiveCentroid, col="red", pch=16, cex=1.5)



libary(maptools)
data(wrld_simpl)
pdf(file="figures/Map-RichnessThruTime.pdf", height=9, width=10)

plot(wrld_simpl[grep("United States", wrld_simpl@data$NAME),], xlim=c(-160, -60), ylim=c(15, 80), axes=FALSE, col='light yellow')
plot(wrld_simpl[grep("Canada", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)
plot(wrld_simpl[grep("Mexico", wrld_simpl@data$NAME),], axes=FALSE, col='light yellow', add=T)

plot(spatialNeutrals, pch=16, cex=0.75, col="lightgray", add=T)
plot(spatialNegatives, pch=16, col="blue", add=T)
plot(spatialPositives, pch=16, col="red", add=T)

legend("bottomright", pch=16, cex=0.75,
       col=c("red", "blue", "gray"), 
       c("Richness increases",
         "Richness decreases", 
         "No richness change"))

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

plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigPos)){
  lines(richness[sigPos[k],]~allTimes, col="blue")
}
lines(posRichnessMeans~allTimes, col="black", lwd=2)


plot(richnessMeans ~ allTimes, 
     xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigNeg)){
  lines(richness[sigNeg[k],]~allTimes, col="red")
}
lines(negRichnessMeans~allTimes, col="black", lwd=2)

dev.off()



