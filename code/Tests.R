load("workspaces/richness.RData")
source("code/Part0-GlobalParam.R")

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

