load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")
source("code/DivThruTimeFunctions.R")
library(reshape)
library(vegan)

## Question: despite no changes in richness, is there still compositional change through time at the sites? ####

#positives, negatives, neutrals
#sigPos, sigNeg, nonSig

sigVector<- vector(length=ncol(richness))
sigVector[sigPos]<- "sigPos"
sigVector[sigNeg]<- "sigNeg"
sigVector[nonSig]<- "nonSig"

# calculate dissimilarity

siteChangesPair<- matrix(ncol=ncol(richness), nrow=nrow(richness))
siteChangesPair<- as.data.frame(siteChangesPair)
colnames(siteChangesPair)<- colnames(richness)
rownames(siteChangesPair)<- rownames(richness)
siteChangesSeq<- siteChangesPair

for (i in 1:ncol(richness)){
  timeSeries <- richness[,i]
  timeSeries <- na.omit(timeSeries)
  ages <- as.numeric(names(timeSeries))
  d<- dist(timeSeries)
  d<- as.matrix(d)
  d2 <- melt(d)[melt(lower.tri(d))$value,]
  names(d2) <- c("t1", "t2", "distance")
  
  pairwisematches <- vector(length=length(ages))
  sequentialmatches <- vector(length=length(ages))
  for (f in 2:length(ages)){
    pairwisematches[f-1] <- intersect(which(d2[,1]==ages[f]), which(d2[,2]==ages[f-1]))
    sequentialmatches[f-1] <- intersect(which(d2[,1]==ages[f]), which(d2[,2]==ages[1]))
  }
  
  pair<- d2[pairwisematches,]
  seq<- d2[sequentialmatches,]
  
  siteChangesPair[match(pair$t1, rownames(siteChangesPair)),i] <- pair$distance
  siteChangesSeq[match(seq$t1, rownames(siteChangesPair)),i] <- seq$distance
}

pairMeans <- apply(siteChangesPair, 1, mean, na.rm=T)
seqMeans <- apply(siteChangesSeq, 1, mean, na.rm=T)

# Plot the changes ####
pdf(file="figures/RichnessChangeThruTime-all-withLines.pdf", height=6, width=10)
par(mfrow=c(1,2))
plot(pairMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(siteChangesPair, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Euclidean Distance of Richness Change - pairwise")
for (i in 1:ncol(siteChangesPair)){
  lines(siteChangesPair[,i]~timeNeg, col="gray") 
}
lines(pairMeans~timeNeg, col="black", lwd=2)

plot(seqMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(siteChangesSeq, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Euclidean Distance of Richness Change - from youngest")
for (i in 1:ncol(siteChangesSeq)){
  lines(siteChangesSeq[,i]~timeNeg, col="gray") 
}
lines(seqMeans~timeNeg, col="black", lwd=2)
dev.off()


### What about keeping track of positive vs negative changes? ####
siteChangesPair<- matrix(ncol=ncol(richness), nrow=nrow(richness))
siteChangesPair<- as.data.frame(siteChangesPair)
colnames(siteChangesPair)<- colnames(richness)
rownames(siteChangesPair)<- rownames(richness)
siteChangesSeq<- siteChangesPair

for (i in 1:ncol(richness)){
  siteName<- colnames(richness)[i]
  sitePath<- files[match(siteName, gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  dat<- na.omit(dat)
  
  y<- calcSiteRichness(dat, minTime, maxTime, pollenThreshold, interval) 

  ages<- dat[,1] # pull out time periods
  datPA<- y$datPA
  rownames(datPA)<- ages
  d<- vegdist(datPA, method="jaccard")
  d<- as.matrix(d)
  d2 <- melt(d)[melt(lower.tri(d))$value,]
  names(d2) <- c("t1", "t2", "jaccard")
  
  pairwisematches <- vector(length=length(ages))
  sequentialmatches <- vector(length=length(ages))
  for (f in 2:length(ages)){
    pairwisematches[f-1] <- intersect(which(d2[,1]==ages[f]), which(d2[,2]==ages[f-1]))
    sequentialmatches[f-1] <- intersect(which(d2[,1]==ages[f]), which(d2[,2]==ages[1]))
  }
  
  pair<- d2[pairwisematches,]
  seq<- d2[sequentialmatches,]
  
  siteChangesPair[match(pair$t1, rownames(siteChangesPair)),i] <- pair$jaccard
  siteChangesSeq[match(seq$t1, rownames(siteChangesPair)),i] <- seq$jaccard
  
}

pairMeans <- apply(siteChangesPair, 1, mean, na.rm=T)
seqMeans <- apply(siteChangesSeq, 1, mean, na.rm=T)

# Plot the changes ####
pdf(file="figures/JaccardThruTime-all-withLines.pdf", height=6, width=10)
par(mfrow=c(1,2))
plot(pairMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(siteChangesPair, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Jaccard Distance - pairwise")
for (i in 1:ncol(siteChangesPair)){
  lines(siteChangesPair[,i]~timeNeg, col="gray") 
}
lines(pairMeans~timeNeg, col="black", lwd=2)

plot(seqMeans ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(siteChangesSeq, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Jaccard Distance - from youngest")
for (i in 1:ncol(siteChangesSeq)){
  lines(siteChangesSeq[,i]~timeNeg, col="gray") 
}
lines(seqMeans~timeNeg, col="black", lwd=2)
dev.off()


### OK, back to looking at overall compositional change, with abundance, not just richness ####
siteChangesPairAbund<- matrix(ncol=ncol(richness), nrow=nrow(richness))
siteChangesPairAbund<- as.data.frame(siteChangesPairAbund)
colnames(siteChangesPairAbund)<- colnames(richness)
rownames(siteChangesPairAbund)<- rownames(richness)
siteChangesSeqAbund<- siteChangesPairAbund

for (i in 1:ncol(richness)){
  siteName<- colnames(richness)[i]
  sitePath<- files[match(siteName, gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  dat<- na.omit(dat)
  
  ages<- dat[,1] # pull out time periods
  rownames(dat)<- ages
  dat<- dat[,-1]
  
  d<- vegdist(dat, method="jaccard")
  d<- as.matrix(d)
  d2 <- melt(d)[melt(lower.tri(d))$value,]
  names(d2) <- c("t1", "t2", "jaccard")
  
  pairwisematches <- vector(length=length(ages))
  sequentialmatches <- vector(length=length(ages))
  for (f in 2:length(ages)){
    pairwisematches[f-1] <- intersect(which(d2[,1]==ages[f]), which(d2[,2]==ages[f-1]))
    sequentialmatches[f-1] <- intersect(which(d2[,1]==ages[f]), which(d2[,2]==ages[1]))
  }
  
  pair<- d2[pairwisematches,]
  seq<- d2[sequentialmatches,]
  
  siteChangesPairAbund[match(pair$t1, rownames(siteChangesPairAbund)),i] <- pair$jaccard
  siteChangesSeqAbund[match(seq$t1, rownames(siteChangesPairAbund)),i] <- seq$jaccard
  
}

pairMeansAbund <- apply(siteChangesPairAbund, 1, mean, na.rm=T)
seqMeansAbund <- apply(siteChangesSeqAbund, 1, mean, na.rm=T)

# Plot the changes ####
pdf(file="figures/JaccardThruTime-abund-all-withLines.pdf", height=6, width=10)
par(mfrow=c(1,2))
plot(pairMeansAbund ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(siteChangesPairAbund, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Jaccard Distance - pairwise")
for (i in 1:ncol(siteChangesPairAbund)){
  lines(siteChangesPairAbund[,i]~timeNeg, col="gray") 
}
lines(pairMeansAbund~timeNeg, col="black", lwd=2)

plot(seqMeansAbund ~ timeNeg, 
     xlim=c(-21000,0), ylim=c(0, max(siteChangesSeqAbund, na.rm=T)), 
     type="n", 
     xlab="Time slice (yr BP)", ylab="Jaccard Distance - from youngest")
for (i in 1:ncol(siteChangesSeqAbund)){
  lines(siteChangesSeqAbund[,i]~timeNeg, col="gray") 
}
lines(seqMeansAbund~timeNeg, col="black", lwd=2)
dev.off()

# Does compositional change correlate to climate change? either rates or magnitudes? ####
# load and process climate data ####
library(raster)
library(MASS)

# Do sites have higher richness in warm vs cool, dry vs mesic times? ####
climateDir<- "/Volumes/bloisgroup/bloislab/Data/Climate/Paleo/CCSM3_500/With_PaleoShorelines/"
tVar<- "tmax_year_ave"
pVar<- "prcp_year_ave"

tStack<- stack(paste(climateDir, allTimes[1], "BP/", tVar, ".tif", sep=""))
pStack<- stack(paste(climateDir, allTimes[1], "BP/", pVar, ".tif", sep=""))

for (i in 2:length(allTimes)){
  tStack <- addLayer(tStack, paste(climateDir, allTimes[i], "BP/", tVar, ".tif", sep=""))
  pStack <- addLayer(pStack, paste(climateDir, allTimes[i], "BP/", pVar, ".tif", sep=""))
}

names(tStack)<- paste("time", allTimes, sep="_")
names(pStack)<- paste("time", allTimes, sep="_")

siteLatLongs <- latlongs[match(colnames(richness), latlongs$Handle), c("Longitude", "Latitude")]
siteT <- extract(tStack, siteLatLongs)
siteP <- extract(pStack, siteLatLongs)

rownames(siteT) <- rownames(siteP) <- colnames(richness)
colnames(siteT) <- colnames(siteP) <- allTimes

siteT <- t(siteT)
siteP <- t(siteP)

## match climate velocity to site changes ####
# siteChangesPair, siteChangesSeq...21000 = change between 21000 and 20500, 500 = change between 500 and 0

siteTempChangesPair<- siteTempVelocPair<- matrix(nrow=nrow(siteChangesPair), ncol=ncol(siteChangesPair))
sitePrecipChangesPair<- sitePrecipVelocPair<- matrix(nrow=nrow(siteChangesPair), ncol=ncol(siteChangesPair))

rownames(siteTempChangesPair) <- rownames(siteTempVelocPair) <- rownames(sitePrecipChangesPair) <- rownames(sitePrecipVelocPair) <- rownames(siteChangesPair)
colnames(siteTempChangesPair) <- colnames(siteTempVelocPair) <- colnames(sitePrecipChangesPair) <- colnames(sitePrecipVelocPair) <- colnames(siteChangesPair)

# total climate change
siteTempChangesPair[2:43,]<- siteT[1:42,] - siteT[2:43,]
sitePrecipChangesPair[2:43,]<- siteP[1:42,] - siteP[2:43,]

for (k in 2:nrow(siteChangesPair)){
  t2<- as.numeric(rownames(siteChangesPair)[k])
  t1<- t2-500
  
  #read in climate data and extract velocity values from shared sites
  tempVeloc<- stack(paste(climateDir, "Climate Velocity/", tVar, "-", t2, "-", t1, ".tif", sep=""))
  names(tempVeloc)<- c("temporalGrad", "spatialGrad", "Velocity", "BRNG")
  
  precipVeloc<- stack(paste(climateDir, "Climate Velocity prcp/", pVar, "-", t2, "-", t1, ".tif", sep=""))
  names(precipVeloc)<- c("temporalGrad", "spatialGrad", "Velocity", "BRNG")

  siteTempVelocPair[k,]<- extract(tempVeloc$Velocity, siteLocs[match(colnames(siteChangesPair), sites)])
  sitePrecipVelocPair[k,]<- extract(precipVeloc$Velocity, siteLocs[match(colnames(siteChangesPair), sites)])

}


# Plot:
# Amount of temperature change, amount of precipitation change, on average
# Rates of temperature change, rates of precipitation change, on average

par(mfrow=c(2,2))
  plot(rowMeans(siteChangesPairAbund, na.rm=T) ~ timeNeg, 
       type="l", ylim=c(0.1, 0.5), xlab="Time (years BP)", ylab="Mean Compositional Change")
  plot(rowMeans(siteChangesPairAbund, na.rm=T) ~ timeNeg, 
       type="l", ylim=c(0.1, 0.5),  xlab="Time (years BP)", ylab="Mean Compositional Change")

  plot(rowMeans(siteTempChangesPair, na.rm=T) ~ timeNeg, 
       type="l", xlab="Time (years BP)", ylab="Mean Max Temperature Change")
  segments(-21000, 0, 0, 0)
  plot(rowMeans(sitePrecipChangesPair, na.rm=T) ~ timeNeg, 
       type="l", xlab="Time (years BP)", ylab="Mean Precipitation Change")
  segments(-21000, 0, 0, 0)
  
#   plot(rowMeans(siteTempVelocPair, na.rm=T) ~ timeNeg,
#        type="l", xlab="Time (years BP)", ylab="Mean Max Temperature Velocity")
#   segments(-21000, 0, 0, 0)
#   plot(rowMeans(sitePrecipVelocPair, na.rm=T) ~ timeNeg, 
#        type="l", xlab="Time (years BP)", ylab="Mean Precipitation Velocity")
#   segments(-21000, 0, 0, 0)
#   

## Plot mean compositional change against climate change
pdf(file="figures/Jaccard-climatechange.pdf", width=10, height=6)
  par(mfrow=c(2,2))
  plot(rowMeans(siteChangesPair, na.rm=T) ~ rowMeans(siteTempChangesPair, na.rm=T),
       pch=16, ylab="Mean Compositional Change - PA", xlab="Mean Max Temperature Change")
  abline(lm(rowMeans(siteChangesPair, na.rm=T) ~ rowMeans(siteTempChangesPair, na.rm=T)))
  summary(lm(rowMeans(siteChangesPair, na.rm=T) ~ rowMeans(siteTempChangesPair, na.rm=T)))
  legend("bottomright", legend="NS", bty="n")
  
  plot(rowMeans(siteChangesPair, na.rm=T) ~ rowMeans(sitePrecipChangesPair, na.rm=T),
       pch=16, ylab="Mean Compositional Change - PA", xlab="Mean Precipitation Change")
  abline(lm(rowMeans(siteChangesPair, na.rm=T) ~ rowMeans(sitePrecipChangesPair, na.rm=T)))
  summary(lm(rowMeans(siteChangesPair, na.rm=T) ~ rowMeans(sitePrecipChangesPair, na.rm=T)))
  legend("bottomright", legend="NS", bty="n")
  
  plot(rowMeans(siteChangesPairAbund, na.rm=T) ~ rowMeans(siteTempChangesPair, na.rm=T),
       pch=16, ylab="Mean Compositional Change - Abund", xlab="Mean Max Temperature Change")
  abline(lm(rowMeans(siteChangesPairAbund, na.rm=T) ~ rowMeans(siteTempChangesPair, na.rm=T)))
  summary(lm(rowMeans(siteChangesPairAbund, na.rm=T) ~ rowMeans(siteTempChangesPair, na.rm=T)))
  legend("bottomright", legend="R2adj= 0.14, p=.009", bty="n")
  
  
  plot(rowMeans(siteChangesPairAbund, na.rm=T) ~ rowMeans(sitePrecipChangesPair, na.rm=T),
       pch=16, ylab="Mean Compositional Change - Jaccard", xlab="Mean Precipitation Change")
  abline(lm(rowMeans(siteChangesPairAbund, na.rm=T) ~ rowMeans(sitePrecipChangesPair, na.rm=T)))
  summary(lm(rowMeans(siteChangesPairAbund, na.rm=T) ~ rowMeans(sitePrecipChangesPair, na.rm=T)))
  legend("bottomright", legend="NS", bty="n")
  
dev.off()  
  
  # No significant associations with richness (siteChangesPair) ###
  # Sig associations with temperature with abundance (siteChangesPairAbund) ###
  
  
  ## What about individual sites, not means? ####
  convertDataFrame <- function(x){ y <- as.matrix(x); z <- as.vector(y); return(z) }

  compchangevector <- convertDataFrame(siteChangesPairAbund)
  tempchangevector <- convertDataFrame(siteTempChangesPair)
  precipchangevector <- convertDataFrame(sitePrecipChangesPair)
  par(mfrow=c(1,2))
  plot(compchangevector ~ tempchangevector, pch=16, cex=0.5)
  abline(lm(compchangevector ~ tempchangevector))
  summary(lm(compchangevector ~ tempchangevector))
  plot(compchangevector ~ precipchangevector, pch=16, cex=0.5)
  abline(lm(compchangevector ~ precipchangevector))
  summary(lm(compchangevector ~ precipchangevector))
  # yes, significant, but very little variation explained
  
  ## What about just the sig pos or sig neg sites?  e.g., the sites that show richness change? 
  compchangevectorSigPos <- convertDataFrame(siteChangesPairAbund[,sigPos])
  tempchangevectorSigPos <- convertDataFrame(siteTempChangesPair[,sigPos])
  precipchangevectorSigPos <- convertDataFrame(sitePrecipChangesPair[,sigPos])
  par(mfrow=c(1,2))
  plot(compchangevectorSigPos ~ tempchangevectorSigPos, pch=16, cex=0.5)
  abline(lm(compchangevectorSigPos ~ tempchangevectorSigPos))
  summary(lm(compchangevectorSigPos ~ tempchangevectorSigPos))
  plot(compchangevectorSigPos ~ precipchangevectorSigPos, pch=16, cex=0.5)
  abline(lm(compchangevectorSigPos ~ precipchangevectorSigPos))
  summary(lm(compchangevectorSigPos ~ precipchangevectorSigPos))
  
  compchangevectorsigNeg <- convertDataFrame(siteChangesPairAbund[,sigNeg])
  tempchangevectorsigNeg <- convertDataFrame(siteTempChangesPair[,sigNeg])
  precipchangevectorsigNeg <- convertDataFrame(sitePrecipChangesPair[,sigNeg])
  par(mfrow=c(1,2))
  plot(compchangevectorsigNeg ~ tempchangevectorsigNeg, pch=16, cex=0.5)
  abline(lm(compchangevectorsigNeg ~ tempchangevectorsigNeg))
  summary(lm(compchangevectorsigNeg ~ tempchangevectorsigNeg))
  plot(compchangevectorsigNeg ~ precipchangevectorsigNeg, pch=16, cex=0.5)
  abline(lm(compchangevectorsigNeg ~ precipchangevectorsigNeg))
  summary(lm(compchangevectorsigNeg ~ precipchangevectorsigNeg))
  # yes, some are significant, but very little variation explained
  
  ## what about dividing change between Pleistocene and Holocene
  siteChangesPairPsigPos <- siteChangesPairAbund[25:nrow(siteChangesPairAbund),sigPos]
  siteTempChangesPairPsigPos <- siteTempChangesPair[25:nrow(siteChangesPair),sigPos]
  sitePrecipChangesPairPsigPos <- sitePrecipChangesPair[25:nrow(siteChangesPair),sigPos]
  
  siteChangesPairHsigPos <- siteChangesPairAbund[1:24,sigPos]
  siteTempChangesPairHsigPos <- siteTempChangesPair[1:24,sigPos]
  sitePrecipChangesPairHsigPos <- sitePrecipChangesPair[1:24,sigPos]
  
  siteChangesPairPsigNeg <- siteChangesPairAbund[25:nrow(siteChangesPairAbund),sigNeg]
  siteTempChangesPairPsigNeg <- siteTempChangesPair[25:nrow(siteChangesPair),sigNeg]
  sitePrecipChangesPairPsigNeg <- sitePrecipChangesPair[25:nrow(siteChangesPair),sigNeg]
  
  siteChangesPairHsigNeg <- siteChangesPairAbund[1:24, sigNeg]
  siteTempChangesPairHsigNeg <- siteTempChangesPair[1:24, sigNeg]
  sitePrecipChangesPairHsigNeg <- sitePrecipChangesPair[1:24, sigNeg]

  # sigPos
  par(mfrow=c(2,2))
  plot(convertDataFrame(siteChangesPairPsigPos) ~ convertDataFrame(siteTempChangesPairPsigPos), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairPsigPos) ~ convertDataFrame(siteTempChangesPairPsigPos)))
  summary(lm(convertDataFrame(siteChangesPairPsigPos) ~ convertDataFrame(siteTempChangesPairPsigPos)))
  plot(convertDataFrame(siteChangesPairPsigPos) ~ convertDataFrame(sitePrecipChangesPairPsigPos), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairPsigPos) ~ convertDataFrame(sitePrecipChangesPairPsigPos)))
  summary(lm(convertDataFrame(siteChangesPairPsigPos) ~ convertDataFrame(sitePrecipChangesPairPsigPos)))
  
  plot(convertDataFrame(siteChangesPairHsigPos) ~ convertDataFrame(siteTempChangesPairHsigPos), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairHsigPos) ~ convertDataFrame(siteTempChangesPairHsigPos)))
  summary(lm(convertDataFrame(siteChangesPairHsigPos) ~ convertDataFrame(siteTempChangesPairHsigPos)))
  plot(convertDataFrame(siteChangesPairHsigPos) ~ convertDataFrame(sitePrecipChangesPairHsigPos), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairHsigPos) ~ convertDataFrame(sitePrecipChangesPairHsigPos)))
  summary(lm(convertDataFrame(siteChangesPairHsigPos) ~ convertDataFrame(sitePrecipChangesPairHsigPos)))
  
  # sigNeg
  par(mfrow=c(2,2))
  plot(convertDataFrame(siteChangesPairPsigNeg) ~ convertDataFrame(siteTempChangesPairPsigNeg), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairPsigNeg) ~ convertDataFrame(siteTempChangesPairPsigNeg)))
  summary(lm(convertDataFrame(siteChangesPairPsigNeg) ~ convertDataFrame(siteTempChangesPairPsigNeg)))
  
  plot(convertDataFrame(siteChangesPairPsigNeg) ~ convertDataFrame(sitePrecipChangesPairPsigNeg), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairPsigNeg) ~ convertDataFrame(sitePrecipChangesPairPsigNeg)))
  summary(lm(convertDataFrame(siteChangesPairPsigNeg) ~ convertDataFrame(sitePrecipChangesPairPsigNeg)))
  
  plot(convertDataFrame(siteChangesPairHsigNeg) ~ convertDataFrame(siteTempChangesPairHsigNeg), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairHsigNeg) ~ convertDataFrame(siteTempChangesPairHsigNeg)))
  summary(lm(convertDataFrame(siteChangesPairHsigNeg) ~ convertDataFrame(siteTempChangesPairHsigNeg)))
  
  plot(convertDataFrame(siteChangesPairHsigNeg) ~ convertDataFrame(sitePrecipChangesPairHsigNeg), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairHsigNeg) ~ convertDataFrame(sitePrecipChangesPairHsigNeg)))
  summary(lm(convertDataFrame(siteChangesPairHsigNeg) ~ convertDataFrame(sitePrecipChangesPairHsigNeg)))
  
  #targeted periods: end of bolling allerod, end of younger dryas
  # warm: 15 - 14.5 , 12 - 11.5
  # cold: 14.5 - 14, 13-12.5
  
  siteChangesPairTarget <- siteChangesPair[c(which(rownames(siteChangesPair)=="13000"), which(rownames(siteChangesPair)=="14500")), ]
  siteTempChangesPairTarget <- siteTempChangesPair[c(which(rownames(siteTempChangesPair)=="13000"), which(rownames(siteTempChangesPair)=="14500")), ]
  sitePrecipChangesPairTarget <- sitePrecipChangesPair[c(which(rownames(sitePrecipChangesPair)=="13000"), which(rownames(sitePrecipChangesPair)=="14500")), ]
  
  par(mfrow=c(1,2))
  plot(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(siteTempChangesPairTarget), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(siteTempChangesPairTarget)))
  summary(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(siteTempChangesPairTarget)))
  
  plot(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(sitePrecipChangesPairTarget), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(sitePrecipChangesPairTarget)))
  summary(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(sitePrecipChangesPairTarget)))
  
  
  
 
### Plot richness and abundance dissim, pairwise and sequential, with lines ####
  # Plot the changes ####
  pdf(file="figures/JaccardThruTime-RichAndAbund-withLines.pdf", height=8, width=10)
  par(mfrow=c(2,2))
  plot(pairMeans ~ timeNeg, 
       xlim=c(-21000,0), ylim=c(0, max(siteChangesPair, na.rm=T)), 
       type="n", 
       xlab="Time slice (yr BP)", ylab="Jaccard Distance (richness) - pairwise")
  for (i in 1:ncol(siteChangesPair)){
    lines(siteChangesPair[,i]~timeNeg, col="gray") 
  }
  lines(pairMeans~timeNeg, col="black", lwd=2)
  
  plot(seqMeans ~ timeNeg, 
       xlim=c(-21000,0), ylim=c(0, max(siteChangesSeq, na.rm=T)), 
       type="n", 
       xlab="Time slice (yr BP)", ylab="Jaccard Distance (richness) - from youngest")
  for (i in 1:ncol(siteChangesSeq)){
    lines(siteChangesSeq[,i]~timeNeg, col="gray") 
  }
  lines(seqMeans~timeNeg, col="black", lwd=2)
  
  plot(pairMeansAbund ~ timeNeg, 
       xlim=c(-21000,0), ylim=c(0, max(siteChangesPairAbund, na.rm=T)), 
       type="n", 
       xlab="Time slice (yr BP)", ylab="Jaccard Distance (abundance) - pairwise")
  for (i in 1:ncol(siteChangesPairAbund)){
    lines(siteChangesPairAbund[,i]~timeNeg, col="gray") 
  }
  lines(pairMeansAbund~timeNeg, col="black", lwd=2)
  
  plot(seqMeansAbund ~ timeNeg, 
       xlim=c(-21000,0), ylim=c(0, max(siteChangesSeqAbund, na.rm=T)), 
       type="n", 
       xlab="Time slice (yr BP)", ylab="Jaccard Distance (abundance) - from youngest")
  for (i in 1:ncol(siteChangesSeqAbund)){
    lines(siteChangesSeqAbund[,i]~timeNeg, col="gray") 
  }
  lines(seqMeansAbund~timeNeg, col="black", lwd=2)
  
  dev.off()
  
  #targeted periods: end of bolling allerod, end of younger dryas ####
  # warm: 15 - 14.5 , 12 - 11.5
  # cold: 14.5 - 14, 13-12.5
  
  siteChangesPairTarget <- siteChangesPair[c(which(rownames(siteChangesPair)=="12000"), which(rownames(siteChangesPair)=="15000")), ]
  siteTempChangesPairTarget <- siteTempChangesPair[c(which(rownames(siteTempChangesPair)=="12000"), which(rownames(siteTempChangesPair)=="15000")), ]
  sitePrecipChangesPairTarget <- sitePrecipChangesPair[c(which(rownames(sitePrecipChangesPair)=="12000"), which(rownames(sitePrecipChangesPair)=="15000")), ]
  
  par(mfrow=c(1,2))
  plot(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(siteTempChangesPairTarget), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(siteTempChangesPairTarget)))
  summary(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(siteTempChangesPairTarget)))
  
  plot(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(sitePrecipChangesPairTarget), pch=16, cex=0.5)
  abline(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(sitePrecipChangesPairTarget)))
  summary(lm(convertDataFrame(siteChangesPairTarget) ~ convertDataFrame(sitePrecipChangesPairTarget)))
  
  
  