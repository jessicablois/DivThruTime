#### Preliminary Code ####
source("code/DivThruTimeFunctions.R")

# Note also that this is calling the original data
# Do I want to use the pollen thresholds established by Kaitlin and Diego?

# pollenDir<- "/Volumes/bloisgroup/bloislab/Data/Biological" #fix this later to connect to server instead of local computer
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"
files<- list.files(path=pollenDir, pattern="gdm.data")
sites<- gsub(".gdm.data.csv", "", files)
times<- seq(0, 21000, by=500)

richness<- matrix(ncol=length(times), nrow=length(sites))
colnames(richness)<- times
rownames(richness)<- sites

for (i in 1:length(sites)){
  richness[i,]<- calcSiteRichness(sites[i], 0, 21000)
}

# Calculate sample size
sampSize<- apply(richness, 2, sampleSize)

# Plot the mean genus richness across all sites and overlay sample size
pdf(file="figures/RichnessSampSizeThruTime-all.pdf", height=4, width=6)
  par(mar=c(4,4,4,4)+0.1)
  plot(colMeans(richness, na.rm=T)~times, type="l", xlim=c(0, 21000), xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
  par(new=T)
  plot(sampSize~times, pch=16, col="red", ylim=c(0, max(sampSize)), xlim=c(0, 21000), axes=F, xlab="", ylab="")
  axis(4, at=seq(0, max(sampSize), by=100))
  mtext("Number of Sites", 4, line=3, col="red")
dev.off()

#Plot each individual line
richnessMeans<- colMeans(richness, na.rm=T)
pdf(file="figures/RichnessThruTime-all-withLines.pdf", height=4, width=6)
  plot(colMeans(richness, na.rm=T)~times, 
       xlim=c(0, 21000), ylim=c(0, max(richness, na.rm=T)), 
       type="n", 
       xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
    for (i in 1:nrow(richness)){
      lines(richness[i,]~times, col=rainbow(nrow(richness))[i]) 
    }
    lines(richnessMeans~times, col="black", lwd=2)
dev.off()

#Are there breakpoints in the lines?
#And if so, are they generally in the same place?
linearSlopes<- vector(length=nrow(richness))
for (i in 1:nrow(richness)){
  tempSpp<- richness[i,]
  if (length(which(!is.na(richness[i,])))>2){ #if there are more than two points, fit a model and store the slope
    lin.mod <- lm(tempSpp~times)
    plot(tempSpp~times, pch=16, ylim=c(0, max(richness, na.rm=T)))
    abline(lin.mod)
    linearSlopes[i]<- lin.mod$coefficients[2]
  }else{
    linearSlopes[i]<- NaN
  }
}

positives<- which(linearSlopes>=0.0005)
negatives<- which(linearSlopes<= -0.0005)
neutrals<- intersect(which(linearSlopes<0.0005),which(linearSlopes> -0.0005))


pdf(file="figures/RichnessThruTime-positives.pdf", height=9, width=10)
  par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
  for (k in 1:length(positives)){
    plot(richness[positives[k],]~times, 
         type="l", 
         ylim=c(0, max(richness, na.rm=T)),
         xlab="", ylab="Site Richness")
    mtext(rownames(richness)[positives[k]], side=3, line=-1.5, adj=1, cex=0.75)
  }
dev.off()

pdf(file="figures/RichnessThruTime-negatives.pdf", height=9, width=10)
  par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
  for (k in 1:length(negatives)){
    plot(richness[negatives[k],]~times, 
         type="l", 
         ylim=c(0, max(richness, na.rm=T)),
         xlab="", ylab="Site Richness")
    mtext(rownames(richness)[negatives[k]], side=3, line=-1.5, adj=1, cex=0.75)
  }
dev.off()

pdf(file="figures/RichnessThruTime-neutrals.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
  for (k in 1:length(neutrals)){
    plot(richness[neutrals[k],]~times, 
         type="l", 
         ylim=c(0, max(richness, na.rm=T)),
         xlab="", ylab="Site Richness")
    mtext(rownames(richness)[neutrals[k]], side=3, line=-1.5, adj=1, cex=0.75)
  }
dev.off()

#spatially plot the positives, negatives, and neutrals
#these are sites with no changes in diversity through time
spatialNeutrals<- siteLocs[match(rownames(richness)[neutrals], siteLocs$Handle),]
coordinates(spatialNeutrals)<- ~lon_alb_km+lat_alb_km
neutralCentroid<- gCentroid(spatialNeutrals)
plot(spatialNeutrals, col="lightgray")
points(neutralCentroid, col="darkgray", pch=16, cex=1.5)

#these are sites with increases in diversity through time
spatialNegatives<- siteLocs[match(rownames(richness)[negatives], siteLocs$Handle),]
coordinates(spatialNegatives)<- ~lon_alb_km+lat_alb_km
negativeCentroid<- gCentroid(spatialNegatives)
plot(spatialNegatives, col="blue", add=T)
points(negativeCentroid, col="darkblue", pch=16, cex=1.5)

#these are sites with decreases in diversity through time
spatialPositives<- siteLocs[match(rownames(richness)[positives], siteLocs$Handle),]
coordinates(spatialPositives)<- ~lon_alb_km+lat_alb_km
positiveCentroid<- gCentroid(spatialPositives)
plot(spatialPositives, col="red", add=T)
points(positiveCentroid, col="red", pch=16, cex=1.5)


segmented.mod <- segmented(lin.mod, seg.Z = ~times, psi=10500)
plot(segmented.mod, add=T)


