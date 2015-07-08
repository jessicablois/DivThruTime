#### Part 1: BASIC PATTERNS ####
source("code/Part0-GlobalParam.R")
source("code/DivThruTimeFunctions.R")

# Do you want to restrict analysis to sites with X number of samples?
restrict <- "yes"
if (restrict =="yes"){
  minSamp <- 6
}

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

richness<- richness[,ncol(richness):1]  # switch richness dataframe to run from past to present
richness<- t(richness)

#### restrict analysis to sites with at least 6 samples ####
if (restrict == "yes"){
  pa <- richness
  pa[which(pa > 0)] <- 1
  noSamp <- colSums(pa, na.rm=T)
  richness<- richness[,which(noSamp >= minSamp)]
}

richnessMeans<- rowMeans(richness, na.rm=T)

# Calculate sample size
sampSize<- apply(richness, 1, sampleSize)

#### at which sites is diversity increasing or decreasing? ####
timeKYR<- rev(allTimes/1000)
richMat<- cbind(timeKYR, richness)

linearSlopes<- vector(length=ncol(richness)) # each site has a slope
pvals<- vector(length=ncol(richness))  # each site has a p value

for (i in 2:ncol(richMat)){
  if (length(which(!is.na(richMat[,i])))>2){ #if there are more than two points, fit a model and store the slope
    lin.mod <- lm(richMat[,i]~richMat[,1])
    plot(richMat[,i]~richMat[,1], xlim=c(21,0), pch=16, ylim=c(0, max(richMat[,-1], na.rm=T)))
    abline(lin.mod)
    linearSlopes[i-1]<- -(lin.mod$coefficients[2])  #need to make this the opposite of original to have time running forward
    pvals[i-1]<- summary(lin.mod)$coefficients[8]
  }else{
    linearSlopes[i-1]<- NaN
    pvals[i-1]<- NaN
  }
}

significants<- which(pvals<=0.05)
positives<- which(linearSlopes>0) #these are sites that show INCREASING RICHNESS through time
negatives<- which(linearSlopes< 0) #these are sites that show DECREASING RICHNESS through time
neutrals<- which(linearSlopes==0)

sigPos<- intersect(significants, positives) #these are sites that show INCREASING RICHNESS through time
sigNeg<- intersect(significants, negatives) #these are sites that show DECREASING RICHNESS through time
nonSig<- seq(1:length(linearSlopes))[-union(sigPos, sigNeg)]

posRichnessMeans<- rowMeans(richness[,sigPos], na.rm=T)
negRichnessMeans<- rowMeans(richness[,sigNeg], na.rm=T)
nonSigRichnessMeans<- rowMeans(richness[,nonSig], na.rm=T)

richZeros <- richness
richZeros[which(is.na(richZeros))]<- 0
siteRichChanges<- -(richZeros[1:(nrow(richZeros)-1),]-richZeros[2:nrow(richZeros),])
propRichChanges<- siteRichChanges/richness[1:(nrow(richZeros)-1), ]

propRichMeans<- rowMeans(propRichChanges, na.rm=T)

save(list=c("richness", "richnessMeans", "sampSize", "sigPos", "sigNeg", "nonSig", "siteRichChanges", "propRichChanges", "propRichMeans"), file="workspaces/richness.RData")



#### Plotting ####

pdf(file="figures/Correlation-Richness-NoSites.pdf", height=4, width=6)
plot(richnessMeans ~ sampSize, pch=16, xlab="Number of sites", ylab="Mean richness")
test1<- cor.test(richnessMeans, sampSize) #is richness correlated with sample size?
legend("bottomright", bty="n", paste("cor=",round(test1$estimate, 4), "; p=", round(test1$p.value, 4), sep=""), cex=0.75)
dev.off()

# Plot the mean genus richness across all sites and overlay sample size
pdf(file="figures/RichnessSampSizeThruTime-all.pdf", height=4, width=6)
par(mar=c(4,4,4,4)+0.1)
plot(richnessMeans~timeKYR, type="l", xlim=c(21,0), xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
par(new=T)
plot(sampSize~timeKYR, pch=16, col="darkgray", ylim=c(0, max(sampSize)), xlim=c(21,0), axes=F, xlab="", ylab="", cex=0.75)
axis(4, at=seq(0, max(sampSize), by=100), col.axis="darkgray")
mtext("Number of Sites", 4, line=3, col="darkgray")
dev.off()

#Plot each individual line
pdf(file="figures/RichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
for (i in 1:ncol(richness)){
  lines(richness[,i]~timeKYR, col=rainbow(nrow(richness))[i]) 
}
lines(richnessMeans~timeKYR, col="black", lwd=2)
dev.off()

#Plot each individual line- sig pos vs sig neg
pdf(file="figures/RichnessThruTime-all-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeKYR, col="gray")
}
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~timeKYR, col="blue")
}
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~timeKYR, col="red")
}
lines(richnessMeans~timeKYR, col="black", lwd=2)
dev.off()

pdf(file="figures/RichnessThruTime-threePanels-withPosNegLines.pdf", height=4, width=6)
par(mfrow=c(1,3))
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(nonSig)){
  lines(richness[,nonSig[k]]~timeKYR, col="gray")
}
lines(nonSigRichnessMeans~timeKYR, col="black", lwd=2)

coolPalette<- colorRampPalette(c("green", "blue"))(43)
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigPos)){
  lines(richness[,sigPos[k]]~timeKYR, col=coolPalette[max(which(!is.na(richness[,sigPos[k]])))-1])
}
lines(posRichnessMeans~timeKYR, col="black", lwd=2)

warmPalette<- colorRampPalette(c("red", "orange"))(43)
plot(richnessMeans ~ timeKYR, 
     xlim=c(21,0), ylim=c(0, max(richness, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Genus Richness")
for (k in 1:length(sigNeg)){
  lines(richness[,sigNeg[k]]~timeKYR, col=warmPalette[max(which(!is.na(richness[,sigNeg[k]])))-1])
}
lines(negRichnessMeans~timeKYR, col="black", lwd=2)

dev.off()



#### Plot the increasing vs decreasing sites ####
pdf(file="figures/RichnessThruTime-DECREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigPos)){
  plot(richness[,sigPos[k]]~timeKYR, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)), xlim=c(21,0),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[sigPos[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-INCREASING.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(sigNeg)){
  plot(richness[,sigNeg[k]]~timeKYR, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)), xlim=c(21,0),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[sigNeg[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

pdf(file="figures/RichnessThruTime-nonSig.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
for (k in 1:length(nonSig)){
  plot(richness[,nonSig[k]]~timeKYR, 
       type="l", 
       ylim=c(0, max(richness, na.rm=T)),
       xlab="", ylab="Site Richness")
  mtext(colnames(richness)[nonSig[k]], side=3, line=-1.5, adj=1, cex=0.75)
}
dev.off()

#### Question: Do sites show most significant proportional change at the same times? ####
#### And do these changes correspond to times with greatest climate change?

timeStep<- as.numeric(rownames(siteRichChanges))

#Plot each individual line
pdf(file="figures/PropRichnessThruTime-all-withLines.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(propRichMeans ~ timeStep, 
     xlim=c(21000,0), ylim=c(min(propRichChanges, na.rm=T), max(propRichChanges, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Mean Proportional Richness Change")
for (i in 1:ncol(propRichChanges)){
  lines(propRichChanges[,i]~timeStep, col=rainbow(nrow(propRichChanges))[i]) 
}
lines(propRichMeans~timeStep, col="black", lwd=2)
dev.off()



# Determine number of sites showing >= 25% change in richness
noOverallSites <- vector(length=length(timeStep))
noPosSites <- vector(length=length(timeStep))
noNegSites <- vector(length=length(timeStep))
propPosSites <- vector(length=length(timeStep))
propNegSites <- vector(length=length(timeStep))

for (j in 1:nrow(propRichChanges)){
  noOverallSites[j] <- length(which(!is.na(propRichChanges[j, ])))
  noPosSites[j] <- length(which(propRichChanges[j, ] >= .25))
  propPosSites[j] <- noPosSites[j] / noOverallSites[j]
  noNegSites[j] <- length(which(propRichChanges[j, ] <= -.25))
  propNegSites[j] <- noNegSites[j] / noOverallSites[j]
}

par(mfrow=c(2,1))
plot(noPosSites~timeStep,
     xlim=c(21000,0), ylim=c(0, max(noOverallSites)), 
     pch=16, col="red", 
     xlab="Time step starting (yr BP)", ylab="Number of sites with >25% change")
points(noNegSites~timeStep,
       xlim=c(21000,0), ylim=c(0, max(noPosSites, noNegSites)), 
       pch=16, col="blue")
points(noOverallSites~timeStep, pch=16, col="gray")

plot(propPosSites~timeStep,
     xlim=c(21000,0), ylim=c(0, max(propPosSites, propNegSites)), 
     type="l", col="red", 
     xlab="Time step starting (yr BP)", ylab="Proportion of sites with >25% change")
points(propNegSites~timeStep,
       xlim=c(21000,0), ylim=c(0, max(propPosSites, propNegSites)), 
       type="l", col="blue")

# Is there a relationship between times with lots of positive and negative change?
cor.test(propNegSites, propPosSites)
plot(propNegSites~propPosSites, pch=16)
abline(lm(propNegSites~propPosSites))

### NEXT STEPS: plot the density of different rates of change with time periods ####
# need a big two column dataframe: column 1 is the time step, column two is the change.  
# maybe three columns: time step, positive proportional changes, negative proportional changes.


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




