#### Arrange pollen data so that presence-absence corresponds to the 5% threshold ####

dat0<- read.csv("~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by time/PollenAbund_Q_0bp.csv", header=T, row.names=1)
pollenMax<- apply(dat0, 2, max)
pollenThreshold<- pollenMax*0.05


#### Part 1: BASIC PATTERNS ####
# Calculate richness at each site

source("code/DivThruTimeFunctions.R")

# pollenDir<- "/Volumes/bloisgroup/bloislab/Data/Biological" #fix this later to connect to server instead of local computer
pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"
files<- list.files(path=pollenDir, pattern="gdm.data")

sites<- gsub(".gdm.data.csv", "", files)
allTimes<- seq(0, 21000, by=500)

richness<- matrix(ncol=length(allTimes), nrow=length(sites))
colnames(richness)<- allTimes
rownames(richness)<- sites

for (i in 1:length(sites)){
  sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  datTimes<- dat[,1] # pull out time periods
  
  minTime<- min(datTimes)
  maxTime<- max(datTimes)

  richness[i,]<- calcSiteRichness(dat, minTime, maxTime, pollenThreshold)
   
}


richnessMeans<- colMeans(richness, na.rm=T)

# Calculate sample size
sampSize<- apply(richness, 2, sampleSize)

pdf(file="figures/Correlation-Richness-NoSites.pdf", height=4, width=6)
  plot(richnessMeans ~ sampSize, pch=16, xlab="Number of sites", ylab="Mean richness")
  test1<- cor.test(richnessMeans, sampSize) #is richness correlated with sample size?
  legend("bottomright", bty="n", paste("cor=",round(test1$estimate, 4), "; p=", round(test1$p.value, 4), sep=""), cex=0.75)
dev.off()

# Plot the mean genus richness across all sites and overlay sample size
pdf(file="figures/RichnessSampSizeThruTime-all.pdf", height=4, width=6)
  par(mar=c(4,4,4,4)+0.1)
  plot(richnessMeans~allTimes, type="l", xlim=c(21000,0), xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
  par(new=T)
  plot(sampSize~allTimes, pch=16, col="red", ylim=c(0, max(sampSize)), xlim=c(21000,0), axes=F, xlab="", ylab="")
  axis(4, at=seq(0, max(sampSize), by=100), col.axis="red")
  mtext("Number of Sites", 4, line=3, col="red")
dev.off()

#Plot each individual line
pdf(file="figures/RichnessThruTime-all-withLines.pdf", height=4, width=6)
  plot(richnessMeans ~ allTimes, 
       xlim=c(21000,0), ylim=c(0, max(richness, na.rm=T)), 
       type="n", 
       xlab="Time slice (kyr BP)", ylab="Mean Genus Richness")
    for (i in 1:nrow(richness)){
      lines(richness[i,]~allTimes, col=rainbow(nrow(richness))[i]) 
    }
    lines(richnessMeans~allTimes, col="black", lwd=2)
dev.off()

#### Part 2: COMPOSITIONAL CHANGE VERSUS CLIMATE CHANGE ####
# Did more compositional change occurred between periods with more rapid climate change?
library(raster)
albersCRS<- CRS("+proj=aea +lat_1=20 +lat2=60 +lat0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# extract amount of compositional change
richChange<- richnessMeans[43:2]-richnessMeans[42:1]

#read in climate velocity stacks.
# the temporalGrad layer has the magnitudes of climate change

climateDir<- "/Volumes/bloisgroup/bloislab/Data/Climate/Paleo/CCSM3_500/With_PaleoShorelines"

var<- "tmax_year_ave"

climateChange<- vector(length=length(richChange))
sites<- read.delim("~/Dropbox/Research/Community Paleomodels/projects/pollen/data/All.site.data-withagemodel-finalv2.txt", sep="\t", header=T)
siteLocs<- SpatialPoints(sites[,c('Longitude','Latitude')], proj4string = CRS("+proj=longlat +datum=WGS84"))
siteLocs<- spTransform(siteLocs, albersCRS)

for (i in length(allTimes):2){
  t1<- allTimes[i]
  t2<- allTimes[i-1]
  
  climVeloc<- stack(paste(climateDir, "/Climate Velocity/", var, "-", t1, "-", t2, ".tif", sep=""))
  names(climVeloc)<- c("temporalGrad", "spatialGrad", "Velocity", "BRNG")
  
  climChange<- extract(climVeloc$temporalGrad, siteLocs)
  climateChange[i-1]<- cellStats(climVeloc$temporalGrad, stat="mean", na.rm=T)
}
climateChange<- climateChange[42:1]

plot(richChange~climateChange, pch=16)
summary(lm(richChange~climateChange))
# no relationship.  But, this ic climate change averaged across the whole US.  
# What about at the actual locations of the sites?
# OK, getting closer, but it's still calling all points for all times, not the points used in each slice











#Are there breakpoints in the lines?
#And if so, are they generally in the same place?
linearSlopes<- vector(length=nrow(richness))
for (i in 1:nrow(richness)){
  tempSpp<- richness[i,]
  if (length(which(!is.na(richness[i,])))>2){ #if there are more than two points, fit a model and store the slope
    lin.mod <- lm(tempSpp~allTimes)
    plot(tempSpp~allTimes, pch=16, ylim=c(0, max(richness, na.rm=T)))
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
    plot(richness[positives[k],]~allTimes, 
         type="l", 
         ylim=c(0, max(richness, na.rm=T)),
         xlab="", ylab="Site Richness")
    mtext(rownames(richness)[positives[k]], side=3, line=-1.5, adj=1, cex=0.75)
  }
dev.off()

pdf(file="figures/RichnessThruTime-negatives.pdf", height=9, width=10)
  par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
  for (k in 1:length(negatives)){
    plot(richness[negatives[k],]~allTimes, 
         type="l", 
         ylim=c(0, max(richness, na.rm=T)),
         xlab="", ylab="Site Richness")
    mtext(rownames(richness)[negatives[k]], side=3, line=-1.5, adj=1, cex=0.75)
  }
dev.off()

pdf(file="figures/RichnessThruTime-neutrals.pdf", height=9, width=10)
par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)
  for (k in 1:length(neutrals)){
    plot(richness[neutrals[k],]~allTimes, 
         type="l", 
         ylim=c(0, max(richness, na.rm=T)),
         xlab="", ylab="Site Richness")
    mtext(rownames(richness)[neutrals[k]], side=3, line=-1.5, adj=1, cex=0.75)
  }
dev.off()

#spatially plot the positives, negatives, and neutrals
library(sp)
library(rgeos)

siteLocs<- read.delim("~/Dropbox/Research/Community Paleomodels/projects/pollen/data/All.site.data-withagemodel-finalv2.txt", sep="\t", header=T)

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


#### Calc additions and subtractions ####

## Note: this works in principle, but need to make 'changes' matrix the size of all taxa, and match taxa at sites to colnames

pollenDir<- "~/Dropbox/Research/Community Paleomodels/projects/pollen/output/all data by site/"
files<- list.files(path=pollenDir, pattern="gdm.data")
sites<- gsub(".gdm.data.csv", "", files)
times<- seq(0, 21000, by=500)

times<- seq(minTime, maxTime, by=500)
nTax<- 26
taxa<- c( "Abies",           "Betula",          "Larix",           "Ostrya.Carpinus", "Picea",          
          "Quercus",         "Alnus",           "Ambrosia.type",   "Fraxinus",        "Populus",        
          "Acer",            "Carya",           "Fagus",           "Juglans",         "Salix",          
          "Tsuga",           "Ulmus",           "Artemisia",       "Galium",         "Ephedra",        
          "Platanus",        "Tilia",           "Corylus",         "Iva",             "Thalictrum",     
          "Pinus")

changes<- as.data.frame(matrix(ncol=nTax+2, nrow=length(times)-1))
changes[,1]<- times[length(times):2]
changes[,2]<- times[(length(times)-1):1]
colnames(changes)[1:2]<- c("Btime", "Etime")
colnames(changes)[3:ncol(additions)]<- taxa
changesAll<- list(changes)

for (i in 1:length(sites)){
  sitePath<- files[match(sites[i], gsub(".gdm.data.csv", "", files))]
  dat<- read.csv(paste(pollenDir, sitePath, sep="")) #read data
  datTimes<- dat[,1] # pull out time periods
  
  datRev<- as.matrix(dat[,-1]) # remove times and convert to matrix
  datRev[which(datRev>0)]<- 1 #replace abundance with presence-absence   
  
  changes<- as.data.frame(matrix(ncol=nTax+2, nrow=length(times)-1))
  changes[,1]<- times[length(times):2]
  changes[,2]<- times[(length(times)-1):1]
  colnames(changes)[1:2]<- c("Btime", "Etime")
  colnames(changes)[3:ncol(additions)]<- taxa
  
  for (j in length(datTimes):2){
    changesTemp<- datRev[j-1,]-datRev[j,]  #0 means no change, -1 means species loss, +1 means species gain
    changes[which(changes[,1]==datTimes[j]), 3:ncol(changes)]<- changesTemp
  }
  
  changesAll[[i]]<-changes 
}